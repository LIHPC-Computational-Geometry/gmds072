/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux (2015)
 *
 * franck.ledoux@cea.fr
 *
 * The FRAME software is a computer program whose purpose is to provide a set
 * of algorithms to build 2D and 3D meshes using frame field concept. The main
 * focus of these algorithms is quadrilateral and hexahedral meshing.
 *
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL-C
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
/*---------------------------------------------------------------------------*/
// GMDS File Headers
#include <GMDS/Utils/Timer.h>
#include <GMDS/Utils/Log.h>
#include <GMDS/IG/IGMesh.h>
#include <GMDS/IO/MeditReader.h>
#include <GMDS/Algo/BoundaryOperator.h>
#include <GMDS/Math/Plane.h>
#include <GMDS/IO/VTKWriter.h>
#include <GMDS/Math/Numerics.h>
/*---------------------------------------------------------------------------*/
// FHeDo File Headers
#include "FHeDo.h"
#include "FieldGenerator.h"
#include "OpenNLFieldSolverStrategy.h"
#include "EigenFieldSolverStrategy.h"
#include "HDMeshGenerator.h"
#include "FrameFieldSmoother.h"
#include "SingularityRegionMesher.h"
#include <StreamlineBuilder.h>

/*---------------------------------------------------------------------------*/
using namespace fhedo;
using namespace gmds;
/*---------------------------------------------------------------------------*/
FHeDo::FHeDo()
:m_mesh(MeshModel(DIM3 | R | F | E | N | R2N | F2N | E2N |
                  R2F | F2R | F2E | E2F | R2E | N2R | N2F | N2E))
{
    initParams();
    m_pg=NULL;
}
/*---------------------------------------------------------------------------*/
FHeDo::~FHeDo()
{
    cleanMeshMarks();
    if(m_pg!=NULL){
        delete m_pg;
    }
}
/*---------------------------------------------------------------------------*/
void FHeDo::execute()
{
    //==================================================================
    // Timers definition
    //==================================================================
    TimePoint start_mesh_reading,
    start_mesh_repair,
    start_boundary_extraction,
    start_FF_generation,
    start_PGP,
    start_HDom, end_algo;

    //==================================================================
    // MESH READING
    start_mesh_reading.update();
    readFile();

    //==================================================================
    // MESH PREPARATION
    start_mesh_repair.update();
    prepareMesh();

    //Boolean marks used in the algorithm must be initialized
    initMeshMarks();

    //==================================================================
    // REMOVE BOUNDARY SLIVERS
    // such slivers create special configuration that lead to specific
    // cases that are not handled by next algorithm (point generation,
    // hex generation)
    // WARNING THE IMPLEMENTATION IS WRONG - MUST BE UPDATED
    // removeBoundarySlivers();


    //==================================================================
    // BOUNDARY DATA BUILDING
    start_boundary_extraction.update();
    buildBoundaryData();

    computeTargetSize();

    //==================================================================
    // FRAME FIELD GENERATION
    start_FF_generation.update();

    computeCotangentWeights();
    generateFrameField();

    if(m_param_ff.generate_streamlines)
    {
        IGMesh streamline_mesh(m_mesh.getModel());
        StreamlineBuilder builder(m_mesh,m_param_gl,streamline_mesh,m_bm);
        builder.execute();
    }

    if(m_param_gl.stop_at==ParamsGlobal::FF_GEN){
        end_algo.update();

        TimeInterval time_mesh_reading     = start_mesh_repair-start_mesh_reading;
        TimeInterval time_mesh_repair      = start_boundary_extraction-start_mesh_repair;
        TimeInterval time_mesh_bnd_extract = start_FF_generation-start_boundary_extraction;
        TimeInterval time_ff               = start_PGP-start_FF_generation;
        TimeInterval time_full             = end_algo-start_mesh_reading;

        Log::mng()<<"-------------------------------------------------\n";
        Log::mng()<<"finished computation at " << end_algo<<"s\n";
        Log::mng()<< "elapsed time: " << time_full << "s\n";
        Log::mng()<< "\t mesh reading: " << time_mesh_reading<< "s\n";
        Log::mng()<< "\t mesh repair : " << time_mesh_repair<< "s\n";
        Log::mng()<< "\t bnd extract.: " << time_mesh_bnd_extract<< "\n";
        Log::mng()<< "\t fr. field g.: " << time_ff<< "s\n";
        Log::mng()<<"-------------------------------------------------\n";

        //Boolean marks cleaning
        cleanMeshMarks();
        return;
    }
    //==================================================================
    // POINT FIELD GENERATION
    start_PGP.update();
    generatePointField();
    if(m_param_gl.stop_at==ParamsGlobal::PF_GEN){
        end_algo.update();

        TimeInterval time_mesh_reading     = start_mesh_repair-start_mesh_reading;
        TimeInterval time_mesh_repair      = start_boundary_extraction-start_mesh_repair;
        TimeInterval time_mesh_bnd_extract = start_FF_generation-start_boundary_extraction;
        TimeInterval time_ff               = start_PGP-start_FF_generation;
        TimeInterval time_pgp              = start_HDom-start_PGP;
        TimeInterval time_full             = end_algo-start_mesh_reading;

        Log::mng()<<"-------------------------------------------------\n";
        Log::mng()<<"finished computation at " << end_algo<<"s\n";
        Log::mng()<< "elapsed time: " << time_full << "s\n";
        Log::mng()<< "\t mesh reading: " << time_mesh_reading<< "s\n";
        Log::mng()<< "\t mesh repair : " << time_mesh_repair<< "s\n";
        Log::mng()<< "\t bnd extract.: " << time_mesh_bnd_extract<< "\n";
        Log::mng()<< "\t fr. field g.: " << time_ff<< "s\n";
        Log::mng()<< "\t points g.   : " << time_pgp<< "s\n";
        Log::mng()<<"-------------------------------------------------\n";

        //Boolean marks cleaning
        cleanMeshMarks();
        return;

    }
    //==================================================================
    // HEX DOM GENERATION
    start_HDom.update();
    generateHexDomMesh();

    end_algo.update();
    TimeInterval time_mesh_reading     = start_mesh_repair-start_mesh_reading;
    TimeInterval time_mesh_repair      = start_boundary_extraction-start_mesh_repair;
    TimeInterval time_mesh_bnd_extract = start_FF_generation-start_boundary_extraction;
    TimeInterval time_ff               = start_PGP-start_FF_generation;
    TimeInterval time_pgp              = start_HDom-start_PGP;
    TimeInterval time_hdom             = end_algo-start_HDom;
    TimeInterval time_full             = end_algo-start_mesh_reading;

    Log::mng()<<"-------------------------------------------------\n";
    Log::mng()<<"finished computation at " << end_algo<<"s\n";
    Log::mng()<< "elapsed time: " << time_full << "s\n";
    Log::mng()<< "\t mesh reading: " << time_mesh_reading<< "s\n";
    Log::mng()<< "\t mesh repair : " << time_mesh_repair<< "s\n";
    Log::mng()<< "\t bnd extract.: " << time_mesh_bnd_extract<< "\n";
    Log::mng()<< "\t fr. field g.: " << time_ff<< "s\n";
    Log::mng()<< "\t points g.   : " << time_pgp<< "s\n";
    Log::mng()<< "\t hex-dom mesh: " << time_hdom<< "s\n";
    Log::mng()<<"-------------------------------------------------\n";

    //Boolean marks cleaning
    cleanMeshMarks();
}

/*---------------------------------------------------------------------------*/
void FHeDo::readFile()
{
    Log::mng()<< "Reading...\n ";
    MeditReader<IGMesh> reader(m_mesh);

    reader.read(m_param_gl.mesh_file, DIM3|R|N);

    Log::mng()<< "(R,F,E,N) = ("
    << m_mesh.getNbRegions()
    << ", " << m_mesh.getNbFaces()
    << ", " << m_mesh.getNbEdges()
    << ", " << m_mesh.getNbNodes() <<")\n";
    Log::mng().flush();
}

/*---------------------------------------------------------------------------*/
void FHeDo::prepareMesh()
{
    Log::mng()<< "Build topology connectivity...\n ";
    IGMeshDoctor doctor(&m_mesh);

    doctor.buildFacesAndR2F();
    doctor.buildEdgesAndX2E();
    doctor.updateUpwardConnectivity();

    Log::mng()<<"\t done ";
    Log::mng()<< " - (R,F,E,N) = ("<< m_mesh.getNbRegions()
    << ", " << m_mesh.getNbFaces()<< ", " << m_mesh.getNbEdges()
    << ", " << m_mesh.getNbNodes() <<")\n";
    Log::mng().flush();

}
/*---------------------------------------------------------------------------*/
void FHeDo::initMeshMarks()
{
    m_bm.mark_node_on_surf  = m_mesh.getNewMark<Node>();
    m_bm.mark_node_on_curv  = m_mesh.getNewMark<Node>();
    m_bm.mark_node_on_pnt   = m_mesh.getNewMark<Node>();
    m_bm.mark_node_isolated = m_mesh.getNewMark<Node>();

    m_bm.mark_edge_on_surf  = m_mesh.getNewMark<Edge>();
    m_bm.mark_edge_on_curv  = m_mesh.getNewMark<Edge>();

    m_bm.mark_face_on_surf  = m_mesh.getNewMark<Face>();
}

/*---------------------------------------------------------------------------*/
void FHeDo::cleanMeshMarks()
{
    //==================================================================
    //clean marks
    m_mesh.unmarkAll<Node>(m_bm.mark_node_on_surf);
    m_mesh.unmarkAll<Node>(m_bm.mark_node_on_curv);
    m_mesh.unmarkAll<Node>(m_bm.mark_node_on_pnt);
    m_mesh.unmarkAll<Node>(m_bm.mark_node_isolated);

    m_mesh.unmarkAll<Edge>(m_bm.mark_edge_on_surf);
    m_mesh.unmarkAll<Edge>(m_bm.mark_edge_on_curv);

    m_mesh.unmarkAll<Face>(m_bm.mark_face_on_surf);

    //==================================================================
    //free marks
    m_mesh.freeMark<Node>(m_bm.mark_node_on_surf);
    m_mesh.freeMark<Node>(m_bm.mark_node_on_curv);
    m_mesh.freeMark<Node>(m_bm.mark_node_on_pnt);
    m_mesh.freeMark<Node>(m_bm.mark_node_isolated);

    m_mesh.freeMark<Edge>(m_bm.mark_edge_on_surf);
    m_mesh.freeMark<Edge>(m_bm.mark_edge_on_curv);

    m_mesh.freeMark<Face>(m_bm.mark_face_on_surf);
}
/*---------------------------------------------------------------------------*/
void FHeDo::computeTargetSize()
{
    //==================================================================
    //We compute the target spacing size if the user didn't provide it
    double spacing = m_param_pf.spacing;
    if(m_param_pf.with_user_spacing){
        Log::mng()<<"\t User-defined spacing: ";
    }
    else
    {
        Log::mng()<<"\t Computed spacing (PGP Ray's value): ";
        double mesh_volume = 0;
        for(IGMesh::region_iterator itr = m_mesh.regions_begin();
            !itr.isDone();
            itr.next()) {
            mesh_volume+= std::abs(itr.value().volume());
        }
        spacing = 2*pow(mesh_volume/m_mesh.getNbNodes(), 1.0/3.0);
    }
    Log::mng()<<spacing<<"\n";
    Log::mng().flush();
    m_param_pf.spacing=spacing;

}
/*---------------------------------------------------------------------------*/
void FHeDo::buildBoundaryData()
{
    Log::mng()<<"Boundary data building (marks, face normals)\n";
    BoundaryOperator boundaryOp(&m_mesh);


    if (!boundaryOp.isValid()) {
        std::cout << "Invalid model for boundary operations" << std::endl;
        throw GMDSException("Invalid model for boundary operations");
    }

    //==================================================================
    // Mark boundary cells
    boundaryOp.markCellOnGeometry(m_bm.mark_face_on_surf,
                                  m_bm.mark_edge_on_surf,
                                  m_bm.mark_node_on_surf,
                                  m_bm.mark_edge_on_curv,
                                  m_bm.mark_node_on_curv,
                                  m_bm.mark_node_on_pnt,
                                  m_bm.mark_node_isolated);


    //==================================================================
    // COLOR SURFACE AND CURVE NODES AND ASSIGN BND NORMALS

    //Color names used to name the mesh variables are defined in
    //the boundary operator class
    boundaryOp.colorEdges(m_bm.mark_edge_on_curv, m_bm.mark_node_on_pnt);

    boundaryOp.colorNodes(m_bm.mark_node_on_pnt);

    for(IGMesh::node_iterator itn = m_mesh.nodes_begin(); !itn.isDone();
        itn.next()){
        Node n = itn.value();
        if(m_mesh.isMarked(n,m_bm.mark_node_on_surf)){
            math::Vector nv= boundaryOp.getOutputNormalOfABoundaryNode(n);
            m_bnd_normals[n.getID()]=math::Vector3d(nv.X(),nv.Y(),nv.Z());
        }
    }

    //now we color nodes on curves and surfaces
    Variable<int>* color_f = m_mesh.getVariable<int>(GMDS_FACE,"BND_SURFACE_COLOR");
    Variable<int>* color_c = m_mesh.getVariable<int>(GMDS_EDGE,"BND_CURVE_COLOR");


    Variable<int>* color_nf= m_mesh.newVariable<int>(GMDS_NODE,"BND_SURFACE_COLOR");
    Variable<int>* color_nc= m_mesh.newVariable<int>(GMDS_NODE,"BND_CURVE_COLOR");

    for(IGMesh::face_iterator itf = m_mesh.faces_begin(); !itf.isDone();
        itf.next()){
        Face f = itf.value();

        //only faces on surface are of interest
        if(!m_mesh.isMarked(f, m_bm.mark_face_on_surf))
            continue;

        std::vector<Node> f_nodes=f.get<Node>();
        for(auto ni:f_nodes){
            if( m_mesh.isMarked(ni,m_bm.mark_node_on_surf) &&
               !m_mesh.isMarked(ni,m_bm.mark_node_on_curv) &&
               !m_mesh.isMarked(ni,m_bm.mark_node_on_pnt)){
                (*color_nf)[ni.getID()]=(*color_f)[f.getID()];
            }
        }
    }
    for(IGMesh::edge_iterator ite = m_mesh.edges_begin(); !ite.isDone();
        ite.next()){
        Edge e = ite.value();

        //only edges on surface are of interest
        if(!m_mesh.isMarked(e, m_bm.mark_edge_on_curv))
            continue;

        std::vector<Node> e_nodes=e.get<Node>();
        for(auto ni:e_nodes){
            if( m_mesh.isMarked(ni,m_bm.mark_node_on_curv) &&
               !m_mesh.isMarked(ni,m_bm.mark_node_on_pnt)){
                (*color_nc)[ni.getID()]=(*color_c)[e.getID()];
            }
        }
    }
}

/*---------------------------------------------------------------------------*/
bool FHeDo::removeBoundarySlivers()
{
    Log::mng()<<"Removal of boundary slivers\n";

    std::set<Region> slivers;
    for(IGMesh::face_iterator itf = m_mesh.faces_begin();
        !itf.isDone(); itf.next()){
        Face f = itf.value();
        std::vector<TCellID> reg_f = f.getIDs<Region>();
        if(reg_f.size()!=1)
            continue;

        //so f is a boundary face.
        std::vector<Node> n = f.get<Node>();
        //We get is only adjacent tetrahedral element
        Region r = f.get<Region>()[0];

        Node n_opp;
        std::vector<Node> r_nodes = r.get<Node>();
        auto found_opp = false;
        for(int i=0; i<4 && !found_opp; i++){
            if(r_nodes[i].getID()!=n[0].getID() &&
               r_nodes[i].getID()!=n[1].getID() &&
               r_nodes[i].getID()!=n[2].getID()){
                n_opp = r_nodes[i];
                found_opp=true;
            }
        }
        if(!found_opp)
            throw GMDSException("Topological error in bnd sliver deletion");

        math::Point p =  n_opp.getPoint();
        math::Plane pl(n[0].getPoint(), n[1].getPoint(), n[2].getPoint());
        math::Point pr = pl.project(p);
        math::Vector3d v01(n[0].getPoint(), n[1].getPoint());
        math::Vector3d v02(n[0].getPoint(), n[2].getPoint());
        math::Vector3d v12(n[1].getPoint(), n[2].getPoint());

        double min_dist = math::min3(v01.norm(),v02.norm(),v12.norm());

        if(p.distance(pr)<0.1*min_dist)
            slivers.insert(r);
    }
    Log::mng()<<"\t boundary slivers: "<<slivers.size()<<"\n";
    std::set<TCellID> faces_to_remove;
    std::set<TCellID> edges_to_remove;
    for(auto s:slivers){

        std::vector<Node> n = s.get<Node>();
        for(auto ni:n){
            ni.remove(s);
        }

        std::vector<Face> f = s.get<Face>();
        for(auto fi:f){
            std::vector<TCellID> reg_fi = fi.getIDs<Region>();
            if(reg_fi.size()<2){
                faces_to_remove.insert(fi.getID());
                std::vector<Edge> e = fi.get<Edge>();
                for(auto ei:e){
                    edges_to_remove.insert(ei.getID());
                }
            }

            fi.remove(s);
        }

        m_mesh.deleteRegion(s);
    }

    for(auto fi:faces_to_remove){
        Face fdi = m_mesh.get<Face>(fi);
        std::vector<Node> n = fdi.get<Node>();
        for(auto ni:n){
            ni.remove(fdi);
        }

        std::vector<Edge> e = fdi.get<Edge>();
        for(auto ei:e){
            ei.remove(fdi);
        }

        m_mesh.deleteFace(fi);

    }

    for(auto ei:edges_to_remove){
        Edge edi = m_mesh.get<Edge>(ei);
        std::vector<Face> f = edi.get<Face>();
        if(f.empty()){
            std::vector<Node> n = edi.get<Node>();
            for(auto ni:n){
                ni.remove(edi);
            }
            m_mesh.deleteEdge(ei);
        }

    }
    return (!slivers.empty());
}
/*---------------------------------------------------------------------------*/
void FHeDo::generateFrameField()
{

    Log::mng()<<"> FRAME FIELD GENERATION (";
    //==================================================================
    // STEP 1 - We generate the frame field on the initial mesh
    //==================================================================
    FieldSolverStrategyItf* solver=0;

    if(m_param_ff.solver_type==ParamsFrameField::OPENNL){
        Log::mng()<<"OPENNL";
        solver = new OpenNLFieldSolverStrategy();
    }
    else if(m_param_ff.solver_type==ParamsFrameField::EIGEN){
        Log::mng()<<"EIGEN";
        Log::mng().flush();
        throw GMDSException("Not yet implemented");
    }
    else{
        Log::mng().flush();
        throw GMDSException("Wrong solver choice");
    }
    Log::mng()<<")\n";
    Log::mng().flush();
    solver = new OpenNLFieldSolverStrategy();
    FieldGenerator ffg(solver,
                       &m_mesh,
                       m_param_gl,
                       m_param_ff,
                       m_bm);
    ffg.execute();
    //==================================================================
    // STEP 2 - Optional - If we want the premesh the singularity lines
    // we have to extract them and provide a new domain to be meshed,
    // i.e. initial one minus singularity regions, and compute the
    // frame field onto it.
    //==================================================================
    gmds::Variable<gmds::math::Chart>* chart_field = ffg.chartField();
    //==================================================================
    // BUILDING OF A AxisAngleRotation Field
    gmds::Variable<math::AxisAngleRotation>* rot_field=NULL;

    rot_field = m_mesh.newVariable<math::AxisAngleRotation>(GMDS_NODE,
                                                            "rotation_field");


    for(IGMesh::node_iterator itn = m_mesh.nodes_begin(); !itn.isDone();
        itn.next()){
        Node n = itn.value();
        math::Chart ci = (*chart_field)[n.getID()];
        (*rot_field)[n.getID()] = math::AxisAngleRotation(ci);
    }

    if(m_param_ff.premeshing_sing_lines){
        Log::mng()<<"Pre-meshing of singularity regions";
        Log::mng().flush();

        SingularityRegionMesher srm(&m_mesh,
                                    m_param_gl,
                                    m_param_pf.spacing,
                                    m_bm);
        srm.execute();
        // m_mesh is modified so.
        //we compute the frame field again on the new mesh
      //  ffg.execute();

    }



    delete solver;

}
/*---------------------------------------------------------------------------*/
void FHeDo::generatePointField()
{
    Log::mng()<<"> POINT FIELD GENERATION \n";
    Log::mng().flush();


    //==================================================================
    //We call the point generation algorithm
    m_pg = new PointGenerator(&m_mesh,
                              m_param_gl,
                              m_bnd_normals,
                              m_bm,
                              m_param_pf.spacing,
                              m_param_pf.curl);
    m_pg->execute();

    if(m_param_gl.with_debug_files){
        writePoints(m_pg->points());
    }
}
/*---------------------------------------------------------------------------*/
void FHeDo::generateHexDomMesh()
{

    Log::mng()<<"> HEX-DOMINANT MESH GENERATION\n";
    Log::mng().flush();

    HDMeshGenerator hex_dom(&m_mesh,
                            m_param_gl,
                            m_bm,
                            m_param_hd,
                            m_pg->points(),
                            m_pg->charts(),
                            m_pg->pointMeshData(),
                            m_pg->pointTypes(),
                            m_pg->pointClassification(),
                            m_pg->pointCurveNumbering(),
                            m_pg->pointSurfaceNumbering(),
                            m_pg->pointSurfaceNormal(),
                            m_param_pf.spacing);
    hex_dom.execute();

}
/*---------------------------------------------------------------------------*/
void FHeDo::computeCotangentWeights()
{
    Variable<double>* cot_w = m_mesh.newVariable<double>(GMDS_EDGE,"cot_weight");

    if(m_param_ff.with_cotangent_weights){

        //=================================================================
        // First, we initialize all to zero
        for(IGMesh::edge_iterator ite = m_mesh.edges_begin(); !ite.isDone();
            ite.next()){
            Edge e = ite.value();
            (*cot_w)[e.getID()] = 0;
        }

        //=================================================================
        // Then we get the contribution from each adjacent tet
        for(IGMesh::region_iterator itr = m_mesh.regions_begin();
            !itr.isDone(); itr.next()){
            Region r = itr.value();
            std::vector<Edge> r_edges = r.get<Edge>();
            std::vector<Node> r_nodes = r.get<Node>();
            for(auto e: r_edges){
                std::vector<Node> e_nodes = e.get<Node>();
                std::vector<Node> others;
                others.reserve(2);
                for(auto n:r_nodes){
                    if(n.getID()!=e_nodes[0].getID() &&
                       n.getID()!=e_nodes[1].getID()){
                        others.push_back(n);
                    }
                }

                (*cot_w)[e.getID()]+= math::cotangentWeight(e_nodes[0].getPoint(),
                                                            e_nodes[1].getPoint(),
                                                            others[0].getPoint(),
                                                            others[1].getPoint());
            }
        }
    }
    else{
        //All weights equal to 1
        for(IGMesh::edge_iterator ite = m_mesh.edges_begin(); !ite.isDone();
            ite.next()){
            Edge e = ite.value();
            (*cot_w)[e.getID()] = 1;
        }

    }
}
/*---------------------------------------------------------------------------*/
bool FHeDo::setParameters(const std::string& AFileName)
{
    m_param.parseIni(AFileName);
    //=================================================================
    // GLOBAL PARAMETERS
    //=================================================================
    std::string val_string;
    int         val_int;
    double      val_double;
    bool        val_bool;

    if(m_param.get("Global","stop_at", val_string)){
        if(val_string=="FF_GEN"){
            m_param_gl.stop_at = ParamsGlobal::FF_GEN;
        }
        else if(val_string=="FF_SMOOTH")
            m_param_gl.stop_at = ParamsGlobal::FF_GEN;
        else if(val_string=="PF_GEN")
            m_param_gl.stop_at = ParamsGlobal::PF_GEN;
        else if(val_string=="HEX_DOM")
            m_param_gl.stop_at = ParamsGlobal::HEX_DOM;
        else //DEFAULT
            m_param_gl.stop_at = ParamsGlobal::HEX_CAV;
    }
    if(m_param.get("Global","with_debug_files", val_bool)){
        m_param_gl.with_debug_files=val_bool;
    }

    if(m_param.get("Global","input_file", val_string)){
        m_param_gl.mesh_file=val_string;
    }
    if(m_param.get("Global","output_dir", val_string)){
        m_param_gl.output_dir=val_string;
    }

    //=================================================================
    // FRAME FIELD PARAMETERS
    //=================================================================
    if(m_param.get("FrameField","solver", val_string)){
        if(val_string=="EIGEN")
            m_param_ff.solver_type = ParamsFrameField::EIGEN;
        else //DEFAULT
            m_param_ff.solver_type = ParamsFrameField::OPENNL;
    }
    if(m_param.get("FrameField","with_cotangent_weights", val_bool)){
        m_param_ff.with_cotangent_weights=val_bool;
    }
    if(m_param.get("FrameField","with_smoothing", val_bool)){
        m_param_ff.with_smoothing=val_bool;
    }
    if(m_param.get("FrameField","smoothing_algo", val_string)){
        if(val_string=="LIU")
            m_param_ff.smoothing_algo = ParamsFrameField::LIU;
        else //DEFAULT
            m_param_ff.smoothing_algo = ParamsFrameField::RAY;
    }
    if(m_param.get("FrameField","smoothing_nb_iterations", val_int)){
        m_param_ff.smoothing_nb_iter=val_int;
    }
    if(m_param.get("FrameField","smoothing_epsilon", val_double)){
        m_param_ff.smoothing_epsilon=val_double;
    }
    if(m_param.get("FrameField","with_mesh_adaptation", val_bool)){
        m_param_ff.with_mesh_adaptation=val_bool;
    }
    if(m_param.get("FrameField","premesh_sing_lines", val_bool)){
        m_param_ff.premeshing_sing_lines=val_bool;
    }
    if(m_param.get("FrameField","generate_streamlines", val_bool)){
        m_param_ff.generate_streamlines=val_bool;
    }




    //=================================================================
    // POINT FIELD PARAMETERS
    //=================================================================
    if(m_param.get("PointField","curl", val_double)){
        m_param_pf.curl=val_double;
    }
    if(m_param.get("PointField","with_user_spacing", val_bool)){
        m_param_pf.with_user_spacing=val_bool;
    }
    if(m_param.get("PointField","spacing", val_double)){
        m_param_pf.spacing=val_double;
    }
    //=================================================================
    // HEX DOM PARAMETERS
    //=================================================================
    if(m_param.get("HexDom","with_edge_interpolation", val_bool)){
        m_param_hd.with_edge_interpolation=val_bool;
    }
    if(m_param.get("HexDom","edge_cone_tolerance", val_double)){
        m_param_hd.edge_cone_tolerance=val_double;
    }
    if(m_param.get("HexDom","with_quad_surface", val_bool)){
        m_param_hd.with_quad_surface=val_bool;
    }
    if(m_param.get("HexDom","with_whisker_weaving", val_bool)){
        m_param_hd.with_whisker_weaving=val_bool;
    }
    if(m_param.get("HexDom","with_pyramids", val_bool)){
        m_param_hd.with_pyramids=val_bool;
    }

    return true;

}
/*---------------------------------------------------------------------------*/
void FHeDo::initParams()
{
    m_param.add("Global","stop_at"                      , Parameters::STRING_P);
    m_param.add("Global","with_debug_files"             , Parameters::BOOL_P);
    m_param.add("Global","input_file"                   , Parameters::STRING_P);
    m_param.add("Global","output_dir"                   , Parameters::STRING_P);

    m_param.add("FrameField","solver"                   , Parameters::STRING_P);
    m_param.add("FrameField","with_cotangent_weights"    , Parameters::BOOL_P);
    m_param.add("FrameField","with_smoothing"            , Parameters::BOOL_P);
    m_param.add("FrameField","smoothing_algo"           , Parameters::STRING_P);
    m_param.add("FrameField","smoothing_nb_iterations"  , Parameters::INT_P);
    m_param.add("FrameField","smoothing_epsilon"        , Parameters::DOUBLE_P);
    m_param.add("FrameField","with_mesh_adaptation"      , Parameters::BOOL_P);
    m_param.add("FrameField","premesh_sing_lines"       , Parameters::BOOL_P);
    m_param.add("FrameField","generate_streamlines"     , Parameters::BOOL_P);

    m_param.add("PointField","curl"                     , Parameters::DOUBLE_P);
    m_param.add("PointField","with_user_spacing"        , Parameters::BOOL_P);
    m_param.add("PointField","spacing"                  , Parameters::DOUBLE_P);

    m_param.add("HexDom"    ,"with_edge_interpolation"  , Parameters::BOOL_P);
    m_param.add("HexDom"    ,"edge_cone_tolerance"      , Parameters::DOUBLE_P);
    m_param.add("HexDom"    ,"with_quad_surface"        , Parameters::BOOL_P);
    m_param.add("HexDom"    ,"with_whisker_weaving"     , Parameters::BOOL_P);
    m_param.add("HexDom"    ,"with_pyramids"            , Parameters::BOOL_P);

    m_param_gl.stop_at           = ParamsGlobal::HEX_CAV;
    m_param_gl.with_debug_files  = false;
    m_param_gl.mesh_file         = "";
    m_param_gl.output_dir        = "";

    m_param_ff.solver_type            = ParamsFrameField::OPENNL;
    m_param_ff.smoothing_algo         = ParamsFrameField::RAY;
    m_param_ff.with_cotangent_weights = true;
    m_param_ff.with_smoothing         = false;
    m_param_ff.smoothing_nb_iter      = 100;
    m_param_ff.smoothing_epsilon      = 1e-4;
    m_param_ff.with_mesh_adaptation   = false;
    m_param_ff.premeshing_sing_lines  = false;
    m_param_ff.generate_streamlines   = false;


    m_param_pf.curl             = 0.35;
    m_param_pf.with_user_spacing= false;
    m_param_pf.spacing          =1;

    m_param_hd.with_edge_interpolation= false;
    m_param_hd.edge_cone_tolerance    = 25.0;//tolerance in degree
    m_param_hd.with_quad_surface      = false;
    m_param_hd.with_whisker_weaving   = false;
    m_param_hd.with_pyramids          = false;

}

/*---------------------------------------------------------------------------*/
void FHeDo::writePoints(const std::vector<math::Point>& APnts)
{
    MeshModel model(DIM3 | R  | N | R2N );
    IGMesh mesh(model);
    for(unsigned int i=0; i<APnts.size();i++){
        math::Point pi = APnts[i];
        Node ni = mesh.newNode(pi);
        mesh.newTet(ni,ni,ni,ni);
    }
    VTKWriter<IGMesh> writer(mesh);
    std::string file_name = m_param_gl.output_dir+"PG";
    writer.write(file_name, DIM3 | R | N);
}

/*---------------------------------------------------------------------------*/
