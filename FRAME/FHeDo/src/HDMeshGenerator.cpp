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
#include <GMDS/IG/IGMesh.h>
#include <GMDS/IG/IGMeshDoctor.h>
#include <GMDS/IO/VTKWriter.h>
#include <GMDS/Math/Triangle.h>
#include <GMDS/Math/Plane.h>
#include <GMDS/Utils/OrientedGraph.h>
#include <GMDS/Utils/Log.h>
#include <GMDS/Algo/BoundaryOperator.h>
#include <GMDS/Utils/Timer.h>
/*---------------------------------------------------------------------------*/
// STL File Headers
#include <set>
#include <algorithm>
/*---------------------------------------------------------------------------*/
// FHeDo File Headers
#include "Tools.h"
#include "TetMeshManipulator.h"
#include "HDMeshGenerator.h"
#include "CavitySurfacePaver.h"
#include "TriangularSurfaceManipulator.h"
#include "CMLSurfEvalImpl.h"
#include "WhiskerWeaving.h"
/*---------------------------------------------------------------------------*/
// TETGEN File Headers
#include "tetgen.h"
/*---------------------------------------------------------------------------*/
// CAMAL File Headers
#include "CMLPaver.hpp"
/*---------------------------------------------------------------------------*/
using namespace gmds;
using namespace fhedo;
/*---------------------------------------------------------------------------*/
void computeCombinationsRec(std::vector<int>& cmb, int n, int p , int i, int k,
                            std::vector<std::vector<int> >& solutions)
{
    if (k == p) {
        solutions.push_back(cmb);
    }
    else if (i < n) {
        computeCombinationsRec(cmb,n,p,i+1,k,solutions);
        cmb[k] = i;
        computeCombinationsRec(cmb,n,p,i+1,k+1,solutions);
    }
}
/*---------------------------------------------------------------------------*/
// Returns in ASol the Cnp combinations of AP integer into [0;AN]
void computeCombinations(int AN, int AP,
                         std::vector<std::vector<int> >& ASol)
{
    std::vector<int> witness;
    witness.resize(3);
    computeCombinationsRec(witness, AN, AP,0,0, ASol);
}

/*---------------------------------------------------------------------------*/
HDMeshGenerator::
HDMeshGenerator(IGMesh*                             AMesh,
                const ParamsGlobal&                 AParamGL,
                const ParamsMark&                   AMarks,
                const ParamsHexDom&                 AParamHD,
                const std::vector<math::Point>&     APnts,
                const std::vector<math::Chart>&     ACharts,
                const std::vector<Cell::Data>&      AData,
                const std::vector<int>&             ATypes,
                const std::vector<int>&             AClass,
                const std::vector<int>&             ACurv,
                const std::vector<int>&             ASurf,
                const std::map<int, math::Vector3d>&ANormal,
                const double                        ASpacing):
m_mesh(AMesh),
m_param_gl(AParamGL),
m_bm(AMarks),
m_param_hexdom(AParamHD),
m_pnt(APnts),
m_chart(ACharts),
m_mesh_data(AData),
m_type(ATypes),
m_classification(AClass),
m_curve(ACurv),
m_surface(ASurf),
m_normal(ANormal),
m_spacing(ASpacing),
m_hexdom(IGMesh(MeshModel(DIM3 | R | F | N | R2N | R2F | F2N | F2R | N2F))),
m_dot_tolerance(0.75)
{
    m_rot_field = m_mesh->getVariable<math::AxisAngleRotation>(GMDS_NODE,
                                                               "rotation_field");
    
    
    if(m_pnt.size()!=m_chart.size() ||
       m_pnt.size()!=m_type.size() ||
       m_pnt.size()!=m_classification.size()||
       m_pnt.size()!=m_surface.size()) {
        throw GMDSException("Incompatible size vector as input!");
    }
    
    m_used.resize(m_pnt.size(), 0);
}
/*---------------------------------------------------------------------------*/
void HDMeshGenerator::execute()
{
    if(m_param_gl.with_debug_files){
        writeInput();
    }
    
    //======================================================================
    // STEP 1 - PREFILTER ON THE POINT DISTANCE
    // SHOULD BE IMPROVED
    //======================================================================
    TimePoint t_start;
    createDistanceFilter();
    
    TimePoint t_prefilter;
    Log::mng()<<t_prefilter-t_start<<"\n";
    Log::mng().flush();

    //======================================================================
    // STEP 2 - Compute the exact relation between points and cells of the
    //          tetrahedral background mesh. Indeed, this relation is loose
    //          at the end of the previous algorithm (point generation). We
    //          only have a close tet for each point (the tet in which the
    //          point was generated)
    //======================================================================
    Log::mng()<<"Mesh association - ";
    TimePoint t_ma_start;
    
    computeMeshAssociation();
    
    TimePoint t_ma_end;
    Log::mng()<<t_ma_end-t_ma_start<<"s \n";
    Log::mng().flush();

    
    if(m_param_gl.with_debug_files){
        writeInput();
    }

    //======================================================================
    // STEP 3 - For each point, build best-fit oriented edges to "adjacent"
    //          points
    //======================================================================
    Log::mng()<<"Build local edges - ";
    TimePoint t_le_start;
    
    std::vector<std::vector<OrientedEdge> > oriented_edges_init;
    buildOrientedEdges(oriented_edges_init);

    TimePoint t_le_end;
    Log::mng()<<t_le_end-t_le_start<<"s \n";
    Log::mng().flush();

    
    if(m_param_gl.with_debug_files)
        writeEdges(oriented_edges_init, m_param_gl.output_dir+"/EDGES_LOCAL");

    //======================================================================
    // STEP 4 - Correction of the oriented-edges to create edges
    //======================================================================
    Log::mng()<<"Build global edges - ";
    TimePoint t_ge_start;
    std::vector<std::vector<OrientedEdge> > oriented_edges;
    buildEdges(oriented_edges_init,oriented_edges);
    
    TimePoint t_ge_end;
    Log::mng()<<t_ge_end-t_ge_start<<"s \n";
    Log::mng().flush();

    if(m_param_gl.with_debug_files)
        writeEdges(oriented_edges, m_param_gl.output_dir+"/EDGES_GLOBAL");

    //======================================================================
    // STEP 5 - For each stable point compute and store its hex-corners
    //======================================================================
    Log::mng()<<"Build hex corners - ";
    TimePoint t_hc_start;
    
    buildHexCorners(oriented_edges);
    
    TimePoint t_hc_end;
    Log::mng()<<t_hc_end-t_hc_start<<"s \n";
    Log::mng().flush();

    //======================================================================
    // STEP 6 - Build stable hexahedral elts
    //======================================================================
    Log::mng()<<"Build stable hexahedra - ";
    TimePoint t_h_start;

    buildHexahedral();
    TimePoint t_h_end;
    Log::mng()<<t_h_end-t_h_start<<"s \n";
    Log::mng().flush();

    //======================================================================
    // STEP 7 - Build 3-sided base prism elts
    //======================================================================
    //    buildPrism3(oriented_edges);
    //    std::cout<<"Prism3 built"<<std::endl;
    if(m_param_gl.with_debug_files)
        writeInput();
    
    Log::mng()<<"Fill in cavities - ";
    TimePoint t_c_start;

    m_hexdom.newVariable<int>(GMDS_FACE, "Cavity");
    if(m_param_gl.stop_at==ParamsGlobal::HEX_DOM){
        if(m_param_gl.with_debug_files)
            writeOutput();
        return;
    }
    
    //======================================================================
    // STEP 8- Cavities filling
    //======================================================================

    buildCavities();
    
    TimePoint t_c_end;
    Log::mng()<<t_c_end-t_c_start<<"s \n";
    Log::mng().flush();
    
    if(m_param_gl.with_debug_files)
        writeOutput();
    
}


/*---------------------------------------------------------------------------*/
void HDMeshGenerator::createDistanceFilter()
{
    double max_distance = 2*m_spacing;

    Log::mng()<<"Global distance: "<<max_distance<<"\n";
    Log::mng().flush();
    
    for(auto i=0; i<m_pnt.size(); i++){
        math::Point pi = m_pnt[i];
        if(m_type[i]!=FRAME_SING){
            for(unsigned int j=i+1; j<m_pnt.size(); j++){
                math::Point pj = m_pnt[j];
                if(m_type[j]!=FRAME_SING){
                    if(pi.distance(pj)<max_distance){
                        m_filter[i].push_back(j);
                        m_filter[j].push_back(i);
                    }
                }
            }
            
        }
    }
}

/*---------------------------------------------------------------------------*/
void HDMeshGenerator::computeMeshAssociation()
{
    for(auto i=0; i<m_pnt.size(); i++){
        computeMeshAssociation(i);
    }
}
/*---------------------------------------------------------------------------*/
void HDMeshGenerator::computeMeshAssociation(const int AID)
{
    int geom_dim = m_classification[AID];
    TCellID start_tet_id = m_mesh_data[AID].id;
    double distance = m_spacing;
    
    bool global_surface = false;
    if(geom_dim==3){
        m_mesh_data[AID] = getRegionContaining(m_pnt[AID],
                                               start_tet_id);
    }
    else if (geom_dim==2){
        int surf_id = m_surface[AID];
        if(global_surface){
            m_mesh_data[AID]= Cell::Data(2,closestFace(m_pnt[AID],surf_id).getID());
        }
        else{
            m_mesh_data[AID] = getBoundaryFaceContaining(m_pnt[AID],
                                                         start_tet_id,
                                                         distance,
                                                         surf_id);
        }
        
    }
    else if (geom_dim==1){
        int curv_id = m_curve[AID];
        if(global_surface){
            m_mesh_data[AID]= Cell::Data(1,closestEdge(m_pnt[AID],curv_id).getID());
        }
        else{
            m_mesh_data[AID] = getBoundaryEdgeContaining(m_pnt[AID],
                                                         start_tet_id,
                                                         distance,
                                                         curv_id);
        }
        
    }
    // Otherwise we do nothing since it is useless for our algorithm
    // since oriented edges starting from points classified on curves and
    // vertices are not built on frame interpolation
}
/*---------------------------------------------------------------------------*/
Face HDMeshGenerator::closestFace(math::Point& AP, const int ASurfID)
{
    Variable<int>* color = m_mesh->getVariable<int>(GMDS_FACE,"BND_SURFACE_COLOR");
    Face closest_face;
    double closest_dist = 10000;
    for(IGMesh::face_iterator itf = m_mesh->faces_begin();
        !itf.isDone(); itf.next()){
        Face f = itf.value();
        if((*color)[f.getID()]==ASurfID){
            std::vector<Node> f_nodes = f.get<Node>();
            math::Triangle t(f_nodes[0].getPoint(),
                             f_nodes[1].getPoint(),
                             f_nodes[2].getPoint());
            math::Point p = t.project(AP);
            if(p.distance2(AP)<closest_dist){
                closest_dist=p.distance2(AP);
                closest_face=f;
            }
        }
    }
    // The point is moved too!!
    
    std::vector<Node> closest_nodes = closest_face.get<Node>();
    math::Triangle t(closest_nodes[0].getPoint(),
                     closest_nodes[1].getPoint(),
                     closest_nodes[2].getPoint());
    
    AP = t.project(AP);
    return closest_face;
    
}

/*---------------------------------------------------------------------------*/
Edge HDMeshGenerator::closestEdge(math::Point& AP, const int ACurvID)
{
    Variable<int>* color = m_mesh->getVariable<int>(GMDS_EDGE,"BND_CURVE_COLOR");
    Edge closest_edge;
    double closest_dist = 10000;
    for(IGMesh::edge_iterator ite = m_mesh->edges_begin();
        !ite.isDone(); ite.next()){
        Edge e = ite.value();
        if((*color)[e.getID()]==ACurvID){
            std::vector<Node> e_nodes = e.get<Node>();
            math::Segment s(e_nodes[0].getPoint(),
                            e_nodes[1].getPoint());
            math::Point p = s.project(AP);
            if(p.distance2(AP)<closest_dist){
                closest_dist=p.distance2(AP);
                closest_edge=e;
            }
        }
    }
    // The point is moved too!!
    std::vector<Node> closest_nodes = closest_edge.get<Node>();
    math::Segment s(closest_nodes[0].getPoint(),
                    closest_nodes[1].getPoint());
    
    AP = s.project(AP);
    return closest_edge;
    
}
/*---------------------------------------------------------------------------*/
bool HDMeshGenerator::getOppositeFace(const TCellID ANodeID,
                                      const Region&  AR,
                                      Face & AOut)
{
    std::vector<Face> fs = AR.get<Face>();
    for(auto f:fs){
        std::vector<TCellID> f_nids = f.getIDs<Node>();
        
        if(f_nids[0]!=ANodeID && f_nids[1]!=ANodeID && f_nids[2]!=ANodeID){
            AOut = f;
            return true;
        }
        
    }
    return false;
}
/*---------------------------------------------------------------------------*/
bool HDMeshGenerator::getOppositeRegion(const TCellID ANodeID,
                                        const Region&  AR,
                                        Region & AOut)
{
    Face opp_face;
    if(getOppositeFace(ANodeID, AR, opp_face)){
        std::vector<Region> f_r = opp_face.get<Region>();
        if(f_r.size()==1){
            return false;
        }
        if(f_r[0].getID()==AR.getID()){
            AOut =f_r[1];
            return true;
        }
        else{
            AOut = f_r[0];
            return true;
        }
    }
    
    return false;
}
/*---------------------------------------------------------------------------*/
bool HDMeshGenerator::getOppositeBndFace(const gmds::TCellID ANodeID,
                                         const gmds::Face&   AF,
                                         gmds::Face&         AOut)
{
    std::vector<Edge> es = AF.get<Edge>();
    for(auto e:es){
        std::vector<TCellID> e_nids = e.getIDs<Node>();
        
        if(e_nids[0]!=ANodeID && e_nids[1]!=ANodeID ){
            
            std::vector<Face> e_f = e.get<Face>();
            std::vector<Face> bnd_faces;
            for(auto f:e_f){
                //we do not keep the face we come from
                if(f.getID()==AF.getID())
                    continue;
                
                if(m_mesh->isMarked(f,m_bm.mark_face_on_surf))
                    bnd_faces.push_back(f);
            }
            
            if(bnd_faces.empty()){
                return false;
            }
            
            AOut =bnd_faces[0];
            return true;
        }
        
    }
    return false;
}
/*---------------------------------------------------------------------------*/
Cell::Data HDMeshGenerator::getRegionContaining(const math::Point& APnt,
                                                const TCellID      ARegionID)
{
    gmds::Region current_r = m_mesh->get<Region>(ARegionID);
    
    while(true){
        std::vector<Node> n = current_r.get<Node>();
        math::Point p[4] ={
            n[0].getPoint(), n[1].getPoint(),
            n[2].getPoint(), n[3].getPoint()
        };
        
        double coeff[4]={0, 0, 0, 0};
        
        math::Point::computeBarycentric(p[0], p[1], p[2], p[3], APnt,
                                        coeff[0],coeff[1],
                                        coeff[2],coeff[3]);
        
        if(coeff[0]>=0 && coeff[1]>=0 && coeff[2]>=0 &&  coeff[3]>=0)
            return Cell::Data(3,current_r.getID());
        
        Region next_r;
        bool on_bnd;
        if(coeff[0]<0)
            on_bnd=getOppositeRegion(n[0].getID(), current_r, next_r);
        else if(coeff[1]<0)
            on_bnd=getOppositeRegion(n[1].getID(), current_r, next_r);
        else if(coeff[2]<0)
            on_bnd=getOppositeRegion(n[2].getID(), current_r, next_r);
        else //(coeff[3]<0)
            on_bnd=getOppositeRegion(n[3].getID(), current_r, next_r);
        
        if(!on_bnd){
            current_r=next_r;
        }
        else{
            return Cell::Data(3,current_r.getID());
        }
    }
}
/*---------------------------------------------------------------------------*/
Cell::Data HDMeshGenerator::
getBoundaryFaceContaining(math::Point&  AP,
                          const TCellID ATetID,
                          const double  ADistance,
                          const int     ASurfID)
{
    gmds::Region current_r = m_mesh->get<Region>(ATetID);
    
    std::set<TCellID> close_reg = getCloseRegionsFrom(AP,current_r,ADistance);
    std::set<TCellID> close_bnd_faces;
    
    Variable<int>* color = m_mesh->getVariable<int>(GMDS_FACE,"BND_SURFACE_COLOR");
    for(auto i:close_reg){
        Region r = m_mesh->get<Region>(i);
        std::vector<Face> r_faces = r.get<Face>();
        for(auto f:r_faces){
            if((*color)[f.getID()]==ASurfID){
                close_bnd_faces.insert(f.getID());
            }
        }
    }
    
    Face closest_face;
    double closest_dist = 10000;
    for(auto i:close_bnd_faces){
        Face f = m_mesh->get<Face>(i);
        std::vector<Node> f_nodes = f.get<Node>();
        math::Triangle t(f_nodes[0].getPoint(),
                         f_nodes[1].getPoint(),
                         f_nodes[2].getPoint());
        math::Point p = t.project(AP);
        if(p.distance2(AP)<closest_dist){
            closest_dist=p.distance2(AP);
            closest_face=f;
        }
    }
    // The point is moved too!!
    
    std::vector<Node> closest_nodes = closest_face.get<Node>();
    math::Triangle t(closest_nodes[0].getPoint(),
                     closest_nodes[1].getPoint(),
                     closest_nodes[2].getPoint());
    
    AP = t.project(AP);
    return Cell::Data(2,closest_face.getID());
}

/*---------------------------------------------------------------------------*/
Cell::Data HDMeshGenerator::
getBoundaryEdgeContaining(math::Point&  AP,
                          const TCellID ATetID,
                          const double  ADistance,
                          const int     ACurvID)
{
    gmds::Region current_r = m_mesh->get<Region>(ATetID);
    
    std::set<TCellID> close_reg = getCloseRegionsFrom(AP,current_r,ADistance);
    std::set<TCellID> close_bnd_edges;
    
    Variable<int>* color = m_mesh->getVariable<int>(GMDS_EDGE,"BND_CURVE_COLOR");
    for(auto i:close_reg){
        Region r = m_mesh->get<Region>(i);
        std::vector<Edge> r_edges = r.get<Edge>();
        for(auto e:r_edges){
            if((*color)[e.getID()]==ACurvID){
                close_bnd_edges.insert(e.getID());
            }
        }
    }
    
    Edge closest_edge;
    double closest_dist = 10000;
    for(auto i:close_bnd_edges){
        Edge e = m_mesh->get<Edge>(i);
        std::vector<Node> e_nodes = e.get<Node>();
        math::Segment s(e_nodes[0].getPoint(),
                        e_nodes[1].getPoint());
        math::Point p = s.project(AP);
        if(p.distance2(AP)<closest_dist){
            closest_dist=p.distance2(AP);
            closest_edge=e;
        }
    }
    
    std::vector<Node> closest_nodes = closest_edge.get<Node>();
    math::Segment s(closest_nodes[0].getPoint(),
                    closest_nodes[1].getPoint());
    
    // The point is moved too!!
    AP = s.project(AP);

    return Cell::Data(1,closest_edge.getID());
}
/*---------------------------------------------------------------------------*/
std::set<TCellID>
HDMeshGenerator::getCloseRegionsFrom(const math::Point& AFromPnt,
                                     const Region& AFromTet,
                                     const double AEpsilon)
{
    std::set<TCellID> region_ids;
    
    //starting from AFromTet, we look for all tet t such that the distance
    //to one face of t is < to AEpsilon
    std::vector<Region> candidates;
    candidates.push_back(AFromTet);
    std::set<TCellID> done;
    
    
    while(!candidates.empty()){
        Region c = candidates.back();
        candidates.pop_back();
        
        done.insert(c.getID());
        
        std::vector<Face> fs = c.get<Face>();
        for(auto f:fs){
            std::vector<Node> n = f.get<Node>();
            math::Plane pl(n[0].getPoint(),
                           n[1].getPoint(),
                           n[2].getPoint());
            if(pl.project(AFromPnt).distance(AFromPnt)<AEpsilon){
                //means close to this face
                
                //Add the region
                region_ids.insert(c.getID());
                
                //test to add the opposite region
                std::vector<TCellID> f_regs = f.getIDs<Region>();
                if(f_regs.size()>1){
                    //otherwise nothing to do we are on the bnd
                    TCellID opp_reg_id = (c.getID()==f_regs[0])?f_regs[1]:f_regs[0];
                    if(region_ids.find(opp_reg_id)==region_ids.end() &&
                       done.find(opp_reg_id)==done.end()){
                        //we add it as a candidate since it is not
                        candidates.push_back(m_mesh->get<Region>(opp_reg_id));
                    }
                }
            }
        }
    }
    return region_ids;
}

/*---------------------------------------------------------------------------*/
void HDMeshGenerator::flagCavityFaces(std::vector<Face>& AInFaces,
                                      std::vector<std::vector<Face> >& ACavFaces)
{
    // We build a N2F connectivity between cavity faces and their
    // nodes only
    std::map<TCellID,std::vector<Face> > local_N2F;
    
    for (auto f:AInFaces){
        std::vector<TCellID> f_node_ids = f.getIDs<Node>();
        
        for(auto n_id:f_node_ids){
            local_N2F[n_id].push_back(f);
        }
    }
    
    //We store a cavity number for each cavity face, with
    // -1 meaning no cavity number assigned
    std::map<TCellID,int> cavity_id;
    for (auto f:AInFaces){
        cavity_id[f.getID()]=-1;
    }
    
    //Now we go through all the cavity faces in order to build
    // individual cavities
    int current_cavity_id=1;
    
    for (auto f:AInFaces){
        
        if (cavity_id[f.getID()]!=-1){
            //already assigned
            continue;
        }
        
        std::vector<Face> heap;
        heap.push_back(f);
        while(!heap.empty()){
            //We pick the last element
            Face current = heap.back();
            heap.pop_back();
            //we put it in the current cavity
            cavity_id[current.getID()]=current_cavity_id;
            
            //We put in the heap all the faces sharing an edge with
            //current and not still assigned to a cavity
            
            
            std::vector<TCellID> current_nodes = current.getIDs<Node>();
            auto nb_nodes = current_nodes.size();
            for(auto i_n=0; i_n<current_nodes.size(); i_n++){
                TCellID ni = current_nodes[i_n];
                TCellID nj = current_nodes[(i_n+1)%nb_nodes];
                
                std::vector<Face> ni_fs = local_N2F[ni];
                std::vector<Face> nj_fs = local_N2F[nj];
                std::vector<Face> nij_fs;
                
                for(auto fi:ni_fs){
                    for(auto fj:nj_fs){
                        if(fi.getID()==fj.getID() && //adj to both points
                           cavity_id[fi.getID()]==-1 && //not already assigned
                           fi.getID()!=current.getID()){ //not the current face
                            nij_fs.push_back(fi);
                        }
                    }
                }
                if(nij_fs.size()==1){
                    heap.push_back(nij_fs[0]);
                }
                else if(nij_fs.size()==2){
                    //We look for the face which is not reachable from current
                    //by a series of regions
                    
                    Face from = current;
                    //just one adj region by definition a cavity face
                    std::vector<Region> from_regs = from.get<Region>();
                    if(from_regs.empty()){
                        std::cout<<"From   : "<<from.getID()<<std::endl;
                        std::cout<<"Other 0: "<<nij_fs[0].getID()<<std::endl;
                        std::cout<<"Other 1: "<<nij_fs[1].getID()<<std::endl;
                        std::cout<<"Node ni: "<<ni<<" -> ";
                        for(auto fi:ni_fs){
                            std::cout<<fi.getID()<<" ";
                        }
                        std::cout<<endl;
                        std::cout<<"Node nj: "<<nj<<" -> ";
                        for(auto fi:nj_fs){
                            std::cout<<fi.getID()<<" ";
                        }
                        std::cout<<endl;
                    }
                    Region from_reg = from.get<Region>()[0];
                    auto in_the_volume = true;
                    Face next;
                    while(in_the_volume){
                        next = getFace(from, from_reg, ni, nj);
                        std::vector<Region> next_r = next.get<Region>();
                        if(next_r.size()==1){
                            in_the_volume=false;
                        }
                        else{
                            if(next_r[0].getID()==from_reg.getID()){
                                from_reg=next_r[1];
                            }
                            else{
                                from_reg=next_r[0];
                            }
                            from = next;
                        }
                    }
                    //next is the face we looked for
                    if(next.getID()==nij_fs[0].getID()){
                        heap.push_back(nij_fs[1]);
                    }
                    else{
                        heap.push_back(nij_fs[0]);
                    }
                    
                }
                else if(nij_fs.size()>2){
                    //WE LOOK FOR THE TWO CLOSEST FACES OF F AROUND [ni,nj]
                    Face candidates[2];
                    std::vector<double> dot;
                    std::vector<math::Vector3d> cross;
                    dot.resize(nij_fs.size());
                    cross.resize(nij_fs.size());
                    
                    math::Point c = 0.5* (m_hexdom.get<Node>(ni).getPoint() +
                                          m_hexdom.get<Node>(nj).getPoint() );
                    
                    math::Vector3d v_edge(c,m_hexdom.get<Node>(ni).getPoint());
                    v_edge.normalize();
                    //                    std::cout<<"FACE: "<<current.getID()<<std::endl;
                    //                    std::cout<<"Nodes: "<<ni<<" "<<nj<<std::endl;
                    //                    std::cout<<"c= "<<c<<std::endl;
                    math::Vector3d normal_ref = current.normal();
                    math::Vector3d v_ref = v_edge.cross(normal_ref);
                    math::Vector3d w_ref(c,current.center());
                    if(v_ref.dot(w_ref)<0)
                        v_ref = -v_ref;
                    v_ref.normalize();
                    //                    std::cout<<"\t vref "<<v_ref<<std::endl;
                    //                    math::Point ci (c.X()+v_ref.X(),
                    //                                    c.Y()+v_ref.Y(),
                    //                                    c.Z()+v_ref.Z());
                    //                  std::cout<<"point ref:"<<ci<<std::endl;
                    for(auto i=0; i<nij_fs.size(); i++){
                        
                        math::Vector3d normal_i = nij_fs[i].normal();
                        math::Vector3d vi = v_edge.cross(normal_i);
                        math::Vector3d wi(c,nij_fs[i].center());
                        if(vi.dot(wi)<0)
                            vi = -vi;
                        
                        vi.normalize();
                        dot[i]= v_ref.dot(vi);
                        cross[i]=v_ref.cross(vi);
                        //                        std::cout<<" - Face "<<nij_fs[i].getID()<<std::endl;
                        //                        std::cout<<"\t v"<<i<<": "<<vi<<std::endl;
                        //                        std::cout<<"\t d"<<i<<": "<<dot[i]<<std::endl;
                        //                        std::cout<<"\t c"<<i<<": "<<cross[i]<<std::endl;
                        //                        math::Point ci (c.X()+vi.X(),
                        //                                        c.Y()+vi.Y(),
                        //                                        c.Z()+vi.Z());
                        //                        std::cout<<"\point "<<i<<": "<<ci<<std::endl;
                        
                    }
                    
                    //We keep the ones with the 2 highest dot products
                    std::vector<int> orient[2];
                    orient[0].push_back(0);
                    for(auto i=1; i<nij_fs.size(); i++){
                        if(cross[0].dot(cross[i])>0){
                            //                            std::cout<<"In 0: "<<i<<std::endl;
                            orient[0].push_back(i);
                        }
                        else{
                            //                          std::cout<<"In 1: "<<i<<std::endl;
                            orient[1].push_back(i);
                        }
                    }
                    
                    //                    for(auto fi:ni_fs){
                    //                        std::cout<<" "<<fi.getID();
                    //                    }
                    //                    std::cout<<std::endl;
                    //                    for(auto fi:nj_fs){
                    //                        std::cout<<" "<<fi.getID();
                    //                    }
                    //                    std::cout<<std::endl;
                    //                    for(auto fi:nij_fs){
                    //                        std::cout<<" "<<fi.getID();
                    //                    }
                    //                    std::cout<<std::endl;
                    if(orient[0].empty() || orient[1].empty())
                        throw GMDSException("Orientation issue for cavity non manifold");
                    
                    for(auto i=0;i<2;i++){
                        auto dot_max = dot[orient[i][0]];
                        auto id = 0;
                        for(auto j=1;j<orient[i].size();j++){
                            if(dot[orient[i][j]]>dot_max){
                                dot_max=dot[orient[i][j]];
                                id = j;
                            }
                        }
                        candidates[i]=nij_fs[orient[i][id]];
                    }
                    
                    //First candidate is the true-oriented with the highest dot product
                    
                    //AND NOW WE HAVE ONLY TWO CANDIDATES
                    Face from = current;
                    //just one adj region by definition a cavity face
                    Region from_reg = from.get<Region>()[0];
                    auto in_the_volume = true;
                    Face next;
                    while(in_the_volume){
                        next = getFace(from, from_reg, ni, nj);
                        std::vector<Region> next_r = next.get<Region>();
                        if(next_r.size()==1){
                            in_the_volume=false;
                        }
                        else{
                            if(next_r[0].getID()==from_reg.getID()){
                                from_reg=next_r[1];
                            }
                            else{
                                from_reg=next_r[0];
                            }
                            from = next;
                        }
                    }
                    //next is the face we looked for
                    if(next.getID()==candidates[0].getID()){
                        heap.push_back(candidates[1]);
                    }
                    else{
                        heap.push_back(candidates[0]);
                    }
                    
                }
                
            }
        }
        current_cavity_id++;
    }
    
    Variable<int>* v_cavity = m_hexdom.getVariable<int>(GMDS_FACE, "Cavity");

    // cavities are now created as a vector of faces
    ACavFaces.resize(current_cavity_id);
    for (auto f:AInFaces){
        (*v_cavity)[f.getID()]=cavity_id[f.getID()];
        ACavFaces[cavity_id[f.getID()]].push_back(f);
    }

}

/*---------------------------------------------------------------------------*/
void HDMeshGenerator::buildCavities()
{
    //==============================================================
    // (1) Add faces in the model for computing boundary data
    //==============================================================
    IGMeshDoctor doctor(&m_hexdom);
    doctor.buildFacesAndR2F();
    // F2R and N2F are now built
    doctor.updateUpwardConnectivity();
    
    //==============================================================
    // (2) Extract cavity faces
    //==============================================================
    Variable<int>* var_cl   = m_hexdom.getVariable<int>(GMDS_NODE,"classification");
    Variable<int>* var_surf = m_hexdom.getVariable<int>(GMDS_NODE,"surface_id");

    // v_bnd[i] = 0 means in-between hex face
    // v_bnd[i] = 1 means surface face
    // v_bnd[i] = 2 means cavity face
    Variable<int>* v_bnd    = m_hexdom.newVariable<int>(GMDS_FACE, "Bnd");

    std::vector<Face> cavity_faces;
    
    int in=0, on=0, between=0;
    for(IGMesh::face_iterator itf = m_hexdom.faces_begin();
        !itf.isDone(); itf.next()){
        Face f = itf.value();
        if(f.get<Region>().size()==1){
            //On the boundary of the hex mesh
            std::vector<TCellID> nodes = f.getIDs<Node>();
            int nb_nodes_in_vol=0;
            std::set<int> surf_ids;
            for(auto n:nodes){
                if((*var_cl)[n]==IN_VOLUME)
                    nb_nodes_in_vol++;
                else if((*var_cl)[n]==ON_SURFACE){
                    surf_ids.insert((*var_surf)[n]);
                }
            }
            if(nb_nodes_in_vol!=0 || surf_ids.size()>1){
                (*v_bnd)[f.getID()]=2;
                between++;
                cavity_faces.push_back(f);
            }
            else {
                //It doesn't mean that the node is not part
                //of a cavity
                (*v_bnd)[f.getID()]=1;
                on++;
            }
        }
        else {
            (*v_bnd)[f.getID()]=0;
            in++;
        }
    }
    
    //===============================================================
    // (3) cavity must be now created. A cavity is a set of faces
    //     sharing edges and forming a manifold surface.
    //===============================================================
    std::vector<std::vector<Face> > cavities;

    flagCavityFaces(cavity_faces, cavities);
    
    if(m_param_gl.with_debug_files)
        writeOutput();

    //===============================================================
    //FIRST WE extract for each cavity some bnd loops
    //===============================================================

    std::vector<std::vector<Node> > all_loops;
    double point_tolerance = 1e-2;
    for(auto cavity:cavities){
        std::vector<std::vector<Node> > local_loops;
        bool succeed = extractBndLoops(cavity, local_loops);
        if(succeed){
            std::vector<std::vector<Node> > clean_loops;
            clean_loops.resize(local_loops.size());
            for(auto i=0; i<local_loops.size(); i++){
                clean_loops[i].clear();
                std::vector<Node> li = local_loops[i];
                for(auto j=0;j<li.size();j++){
                    Node nj = li[j];
                    Node nk = li[(j+1)%li.size()];
                    if(nj.getPoint().distance2(nk.getPoint())>point_tolerance){
                        clean_loops[i].push_back(nj);
                    }
                }
            }

            all_loops.insert(all_loops.end(),
                             clean_loops.begin(),
                             clean_loops.end());
        }
        
    }
    
    //===============================================================
    //All loops are imprinted together on the surface mesh
    //===============================================================
    imprintLoop(all_loops);
    
    if(m_param_gl.with_debug_files){
        writeOutput();
    }
    //===============================================================
    //All the new faces have been now added in the hex dom mesh,
    //we can so extract the final boundaries of cavities which
    // totally enclose cavities so.
    std::vector<std::vector<Face> > final_cavities;

    extractFinalCavityBoundaries(final_cavities);

    //===============================================================
    //Some cavities can have remaining holes, we fill them now
    for(auto i_c=0; i_c<final_cavities.size(); i_c++){
        if(final_cavities[i_c].empty())
            continue;
        closeCavity(final_cavities[i_c]);
    }
    
    if(m_param_gl.with_debug_files){
        writeOutput();
    }
    //===============================================================
    // Finally we fill in each cavity considering algorithm options
    if(m_param_hexdom.with_whisker_weaving){
        for(auto i_c=0; i_c<final_cavities.size(); i_c++){
            if(final_cavities[i_c].empty())
                continue;
            applyWhiskerWeaving(final_cavities[i_c]);
        }
    }
    else if(m_param_hexdom.with_pyramids){
        throw GMDSException("Pyramid transition not yet implemented");
    }
    else{
        for(auto i_c=0; i_c<final_cavities.size(); i_c++){
            if(final_cavities[i_c].empty())
                continue;
            buildTetMesh(final_cavities[i_c]);
        }
    }
    
    if(m_param_gl.with_debug_files){
        writeOutput();
    }

}

/*---------------------------------------------------------------------------*/
void HDMeshGenerator::applyWhiskerWeaving(std::vector<Face>& AC)
{
    //the boundary mesh we provide to the WW algorithm
    IGMesh bnd_mesh(MeshModel(DIM3|F|N|F2N));
    //==================================================================
    // We keep in mind original nodes to project back the generated mesh
    //==================================================================
    std::set<TCellID> from_nodes;
   
    for(auto f:AC){
        std::vector<TCellID> ns = f.getIDs<Node>();
        from_nodes.insert(ns.begin(),ns.end());
    }
    std::vector<TCellID> to_nodes;
    to_nodes.reserve(from_nodes.size());
    std::map<TCellID, TCellID> node_map_1, node_map_2;
    //==================================================================
    // We copy AC nodes into the boundary mesh
    //==================================================================
    for(auto from_id:from_nodes){
        math::Point p= m_hexdom.get<Node>(from_id).getPoint();
        Node to = bnd_mesh.newNode(p);
        node_map_1[to.getID()]=from_id;
        node_map_2[from_id]=to.getID();
        to_nodes.push_back(to.getID());
    }
    //==================================================================
    // We copy AC faces into the boundary mesh and we keek an inward
    // reference for having the direction of propagation
    //==================================================================
    math::Vector3d ref_dir;
    TCellID        ref_face;
    bool           has_ref =false;
    for(auto f:AC){
        std::vector<TCellID> ns = f.getIDs<Node>();
        if(ns.size()==3){
            bnd_mesh.newTriangle(node_map_2[ns[0]],
                                 node_map_2[ns[1]],
                                 node_map_2[ns[2]]);
        }
        else if (ns.size()==4){
            Face q = bnd_mesh.newQuad(node_map_2[ns[0]],
                                     node_map_2[ns[1]],
                                     node_map_2[ns[2]],
                                     node_map_2[ns[3]]);
            
            if(!has_ref){
                std::vector<Region> fr = f.get<Region>();
                if(!fr.empty()){
                    //means we are along the boundary of already
                    //created hex. elements
                    ref_dir  = getOutputNormal(f, fr[0]);
                    ref_face = q.getID();
                    has_ref  = true;
                }
            }
        }
    }
    static int cpt=0;
    WhiskerWeaving ww;
    ww.setDebugFile(m_param_gl.output_dir+"/ww_"+to_string(cpt++)+"_");
    ww.execute(&bnd_mesh, ref_face, ref_dir);
    //We build faces of the mesh given to the WW algorithm
    // A inward normal is computed to as a reference
     math::Vector3d ref_normal;
}
/*---------------------------------------------------------------------------*/
void HDMeshGenerator::buildTetMesh(std::vector<Face>& AC)
{
    IGMesh bnd_mesh(MeshModel(DIM3|F|N|F2N));
    //We get the nodes of the boundary
    std::set<TCellID> from_nodes;
    for(auto f:AC){
        std::vector<TCellID> ns = f.getIDs<Node>();
        from_nodes.insert(ns.begin(),ns.end());
    }
    std::vector<TCellID> to_nodes;
    to_nodes.reserve(from_nodes.size());
    std::map<TCellID, TCellID> node_map_1, node_map_2;
    for(auto from_id:from_nodes){
        math::Point p= m_hexdom.get<Node>(from_id).getPoint();
        Node to = bnd_mesh.newNode(p);
        node_map_1[to.getID()]=from_id;
        node_map_2[from_id]=to.getID();
        to_nodes.push_back(to.getID());
        
    }
    
    for(auto f:AC){
        std::vector<TCellID> ns = f.getIDs<Node>();
        if(ns.size()==3){
            bnd_mesh.newTriangle(node_map_2[ns[0]],
                                 node_map_2[ns[1]],
                                 node_map_2[ns[2]]);
        }
        else if (ns.size()==4){
            bnd_mesh.newTriangle(node_map_2[ns[0]],
                                 node_map_2[ns[1]],
                                 node_map_2[ns[2]]);
            bnd_mesh.newTriangle(node_map_2[ns[0]],
                                 node_map_2[ns[2]],
                                 node_map_2[ns[3]]);
        }
    }
    IGMesh tet_mesh(MeshModel(DIM3|R|N|R2N));

    TetMeshManipulator tm(0);
    tm.tetrahedralizeTetgen(&bnd_mesh, &tet_mesh);
    
    //Now the tetrahedra stored in tet_mesh must be added in m_hexdom
    //To be better afterward
    for(IGMesh::region_iterator itr =tet_mesh.regions_begin();
        !itr.isDone(); itr.next()){
        Region r = itr.value();
        std::vector<Node> r_nodes = r.get<Node>();
        Node n0 = m_hexdom.newNode(r_nodes[0].getPoint());
        Node n1 = m_hexdom.newNode(r_nodes[1].getPoint());
        Node n2 = m_hexdom.newNode(r_nodes[2].getPoint());
        Node n3 = m_hexdom.newNode(r_nodes[3].getPoint());
        m_hexdom.newTet(n0,n1,n2,n3);
    }
    
}
/*---------------------------------------------------------------------------*/
void HDMeshGenerator::closeCavity(std::vector<Face>& AC)
{
    
    // We build a N2F connectivity between cavity faces and their
    // nodes only
    std::map<TCellID,std::vector<Face> > local_N2F;
    
    for (auto f:AC){
        std::vector<TCellID> f_node_ids = f.getIDs<Node>();
        
        for(auto n_id:f_node_ids){
            local_N2F[n_id].push_back(f);
        }
    }
    
    //a cavity must be a closed surface, if not we close it.
    std::vector<FakeEdge> bnd_edges;
    for(auto f:AC){
        std::vector<TCellID> n_ids = f.getIDs<Node>();
        for(auto i=0; i<n_ids.size();i++){
            TCellID id_i = n_ids[i];
            TCellID id_j = n_ids[(i+1)%n_ids.size()];
            std::vector<Face> ni_fs = local_N2F[id_i];
            std::vector<Face> nj_fs = local_N2F[id_j];
            std::vector<Face> nij_fs;
            
            for(auto fi:ni_fs){
                for(auto fj:nj_fs){
                    if(fi.getID()==fj.getID() && //adj to both points
                       fi.getID()!=f.getID()){ //not the current face
                        nij_fs.push_back(fi);
                    }
                }
            }//for(auto fi:ni_fs)
            
            if(nij_fs.size()<1){
                //we have a bnd_edge!!!
                bnd_edges.push_back(FakeEdge(id_i,id_j));
            }
        }//for(auto i=0; i<n_ids.size();i++)
    }
    if(bnd_edges.empty()){
        //closed cavity, nothing to do
        return;
    }
    
    Variable<int>* v_cavity = m_hexdom.getVariable<int>(GMDS_FACE, "Cavity");
    int cavity_color = (*v_cavity)[AC[0].getID()];

    //look for loops of bnd edges (only 3-sided loop right now)
    std::vector<bool> edge_done;
    edge_done.resize(bnd_edges.size());
    for(auto i=0; i<edge_done.size(); i++){
        edge_done[i]=false;
    }

    
    FakeEdge cur_edge = bnd_edges[0];
    int cur_index = 0;
    bool all_done = false;
    while(!all_done){
        TCellID e1 = cur_edge.getID().getID1();
        TCellID e2 = cur_edge.getID().getID2();
        FakeEdge edge_1, edge_2;
        int index_1, index_2;
        TCellID other_e1, other_e2;
        for(auto i=0; i<bnd_edges.size();i++){
            FakeEdge ei=bnd_edges[i];
            if(ei==cur_edge)
                continue;
            if(edge_done[i])
                continue;
            if(ei.getID().getID1()==e1){
                edge_1=ei;
                index_1=i;
                other_e1=ei.getID().getID2();
            }
            else if(ei.getID().getID2()==e1){
                edge_1=ei;
                index_1=i;
                other_e1=ei.getID().getID1();
            }
            else if(ei.getID().getID1()==e2){
                edge_2=ei;
                index_2=i;
                other_e2=ei.getID().getID2();
            }
            else if(ei.getID().getID2()==e2){
                edge_2=ei;
                index_2=i;
                other_e2=ei.getID().getID1();
            }
        }
        
        if(other_e1==other_e2){
//            std::cout<<"CLOSE WITH "<<e1<<", "<<e2<<", "<<other_e1<<std::endl;
            Face t = m_hexdom.newTriangle(e1,e2,other_e1);
            (*v_cavity)[t.getID()]=cavity_color;
            AC.push_back(t);
        }
        else{
//            std::cout<<"CLOSE ONLY TRIANGULAR FACES"<<std::endl;
        }
        //look for the next edge to work on!
        edge_done[cur_index]=true;
        edge_done[index_1]  =true;
        edge_done[index_2]  =true;
        bool find_next = false;
        for(auto idx=0; idx<edge_done.size() && !find_next; idx++){
            if(edge_done[idx]==false){
                cur_index=idx;
                cur_edge = bnd_edges[idx];
                find_next=true;
            }
        }
        if(!find_next){
            all_done=true;
        }
    }
    

}
/*---------------------------------------------------------------------------*/
void HDMeshGenerator::
extractFinalCavityBoundaries(std::vector<std::vector<Face> >& AFinalCavities)
{
    Variable<int>* v_cavity = m_hexdom.getVariable<int>(GMDS_FACE, "Cavity");
    
    std::vector<Face> cavity_faces;

    for(IGMesh::face_iterator itf = m_hexdom.faces_begin();
        !itf.isDone(); itf.next()){
        Face f = itf.value();
        if(f.get<Region>().size()==0 || (*v_cavity)[f.getID()]!=0){
            //means we have a cavity face
            cavity_faces.push_back(f);
        }
         
    }
    
    //Now we reinitialize the  cavity variable before extracting cavities again
    for(IGMesh::face_iterator itf = m_hexdom.faces_begin();
        !itf.isDone(); itf.next()){
        Face f = itf.value();
        (*v_cavity)[f.getID()]=0;
    }
    
    flagCavityFaces(cavity_faces,AFinalCavities);
    
    
    if(m_param_gl.with_debug_files){
        IGMesh cavity_mesh(MeshModel(DIM3|F|F2N));
        for(auto f:cavity_faces){
            std::vector<Node> f_nodes = f.get<Node>();
            std::vector<Node> cav_nodes;
            cav_nodes.reserve(f_nodes.size());
            for(auto n:f_nodes){
                cav_nodes.push_back(cavity_mesh.newNode(n.getPoint()));
            }
            cavity_mesh.newFace(cav_nodes);
        }
        
        VTKWriter<IGMesh> nws(cavity_mesh);
        nws.write(m_param_gl.output_dir+"/CAVITY_BOUNDARY",F|N);

    }
    
}
/*---------------------------------------------------------------------------*/
bool HDMeshGenerator::extractBndLoops(std::vector<Face>& AFaces,
                                      std::vector<std::vector<Node> >& ALoops)
{
    
    if(AFaces.empty())
        return true;
    //===============================================================
    // (1) First, we orient faces locally to the cavity. Local
    //     numbering will be such that it is going counter-clockwise
    //     the outward normal
    //===============================================================
    for(auto f:AFaces){
        //This is a boundary face, so it has only one adjacent region
        math::Vector3d n = getOutputNormalOfABoundaryFace(f);
        std::vector<Node> fn = f.get<Node>();
        math::Vector3d v01(fn[0].getPoint(),fn[1].getPoint());
        math::Vector3d v12(fn[1].getPoint(),fn[2].getPoint());
        if(n.dot(v01.cross(v12))<0){
            //we reverse the numerotation
            std::reverse(std::begin(fn), std::end(fn));
            f.set(fn);
        }
    }
    //===============================================================
    // (2) We build boundary loops
    //===============================================================
    //First we build temporary edges
    std::map<FakeEdge, int> edge_count;
    //A fake edge e consits in 2 node ids with e.first < e.second
    // orientation is true if it is the real edge goes in the same
    // direction, false otherwise
    std::map<FakeEdge, bool> edge_orientation;
    //For each boundary edge, we will keep the out normal of the face it
    // comes from
    std::map<FakeEdge, math::Vector3d> edge_vec;

    for(auto f:AFaces) {
        std::vector<Node> fn = f.get<Node>();
        //we goes through f nodes following the true face orientation
        // thanks to step (1)
        for(auto i=0; i<fn.size(); i++){
            Node ni = fn[i];
            Node nj = fn[(i+1)%fn.size()];
            FakeEdge eij(ni.getID(),nj.getID());
            
            if(edge_count.find(eij)==edge_count.end()){
                //New edge, including bnd_edge
                edge_count[eij]=1;
                if(eij.getID().getID1()==ni.getID())
                    edge_orientation[eij]=true;//same orientation
                else
                    edge_orientation[eij]=false;//opp. orientation
                
                edge_vec[eij]=getOutputNormalOfABoundaryFace(f);

            }
            else{
                edge_count[eij]++;
            }
        }
    }
    // Second, we only keep boundary edges and we build a graph
    // on it
    std::vector<FakeEdge> bnd_E;
    std::set<TCellID> bnd_N;
    //We get the boundary edges, and keep the edges
    for(auto eij : edge_count){
        if(eij.second==1){
            bnd_E.push_back(eij.first);
            bnd_N.insert(eij.first.first());
            bnd_N.insert(eij.first.second());
        }
    }
    //the graph is initialized from nodes
    OrientedGraph g(bnd_N);
    
    //then oriented edges are added.
    std::vector<bool> done;
    done.resize(bnd_E.size());

    //we transfer inward cavity vectors from fake edges
    //to  oriented graph edges
    std::map<int,math::Vector3d> in_cavity_vec;
    for(auto i_e=0; i_e<bnd_E.size(); i_e++){
        FakeEdge ei = bnd_E[i_e];
        done[i_e]=false;
        if(edge_orientation[ei]==true){
            g.addEdge(i_e, ei.first(),ei.second());
            //std::cout<<"edge "<<ei.first()<<" -> "<<ei.second()<<std::endl;
        }
        else{
            g.addEdge(i_e, ei.second(),ei.first());
            //std::cout<<"edge "<<ei.second()<<" -> "<<ei.first()<<std::endl;
        }
        in_cavity_vec[i_e]= edge_vec[ei];
    }
    g.updateNodes();

    std::vector<std::vector<GraphEdge*> > loops;
    for(int i=0; i<g.nbEdges();i++){
        
        if(done[i]==true)
            continue;
        
        std::vector<GraphEdge*> loop;
        GraphEdge* first = g.edge(i);
//        std::cout<<"*****************************"<<std::endl;
//        std::cout<<"first= "
//        <<first->tail()->id()<<" -> "
//        <<first->head()->id()<<std::endl;
        
        //we count the number of times we go through bnd nodes
        // for creating this loop
        std::map<TCellID, int> nb_traversed;
        for(auto n:bnd_N){
            nb_traversed[n]=0;
        }
        done[i]=true;
        
        loop.push_back(first);
        nb_traversed[first->head()->id()]=1;
        nb_traversed[first->tail()->id()]=1;
        GraphEdge* prev =first;
        
        GraphEdge* current    = NULL;
        GraphNode* prev_node  = NULL;
        GraphNode* first_node = NULL;
        
        GraphNode* head = prev->head();
        GraphNode* tail = prev->tail();
        //We get all the edge going out of head
        std::vector<GraphEdge*> nexts = prev->getEdgesStartingFrom(head);
        bool find_next=false;
        for(auto j=0; j<nexts.size() && !find_next;j++){
            if(done[nexts[j]->id()]==false){
                find_next=true;
                current =nexts[j];
                prev_node = head;
                first_node= tail;
            }
        }
        if(!find_next){
            throw GMDSException("Error in the edge boundary cavity building process");
        }
        GraphNode* next_node  = current->head();
        done[current->id()]=true;
        loop.push_back(current);

        while(next_node!=first_node){
//            std::cout<<"current edge "<<"["
//            <<current->head()->id()<<", "
//            <<current->tail()->id()<<"]\n ";
//            std::cout<<"prev node= "<<prev_node->id()<<std::endl;
//            std::cout<<"next node= "<<next_node->id()<<std::endl;
            //We take the next potential edges
            nexts = current->getEdgesStartingFrom(next_node);
            GraphEdge* next = NULL;
            
            //Among all of them, we pick the first available one
            bool find_next = false;
            for(auto j=0; j<nexts.size() && !find_next;j++){
                if(done[nexts[j]->id()]==false){
                    find_next=true;
                    next = nexts[j];
                }
            }
            if(!find_next)
                throw GMDSException("Error in the edge boundary cavity building process(2)");
            
            prev_node=next_node;
            next_node = next->head();
            
            nb_traversed[prev_node->id()]++;
            
            if(nb_traversed[prev_node->id()]>1){
                //we just close an inner loop
                //going backward in the list, we look for the last
                //edge adjacent to this node
                std::vector<GraphEdge*> inner_loop;
                bool found_first_edge = false;
                GraphEdge* last_edge = loop.back();
                loop.pop_back();
                inner_loop.push_back(last_edge);
                while(!found_first_edge){
                    GraphEdge* last_edge = loop.back();
                    loop.pop_back();
                    inner_loop.push_back(last_edge);
                    if(last_edge->tail()==prev_node ||
                       last_edge->head()==prev_node){
                        found_first_edge=true;
                    }
                };
                //necessary to keep the right orientation
                std::reverse(inner_loop.begin(),inner_loop.end());
                loops.push_back(inner_loop);
            }

            //we traverse the node
            //next edge becomes current edge
            current=next;
            prev_node = next_node;
            done[current->id()]=true;
            loop.push_back(current);

        }//while(next_node!=first_node)
        
        loops.push_back(loop);
        
    }
    
    //===============================================================
    // Loops do not self-intersect now
    // So just work with list of nodes
    //===============================================================
    std::vector<std::vector<Node> > loop_nodes;
    //keep the classification of each point
    std::vector<std::vector<int> > loop_classification;
    std::vector<std::vector<int> > loop_geom_id;
    loop_nodes.reserve(loops.size());
    loop_classification.reserve(loops.size());
    loop_geom_id.reserve(loops.size());
    
    for(auto l:loops){
        std::vector<Node> l_nodes;
        l_nodes.reserve(l.size());
        for(auto e:l){
            l_nodes.push_back(m_hexdom.get<Node>(e->head()->id()));
        }
        loop_nodes.push_back(l_nodes);
    }
    
    if(m_param_gl.with_debug_files){
        writeLoops(loop_nodes);
    }

    ALoops = loop_nodes;
    return true;
}
/*---------------------------------------------------------------------------*/
bool HDMeshGenerator::imprintLoop(std::vector<std::vector<Node> >& ALoops)
{
    Variable<int>* v_cl   = m_hexdom.getVariable<int>(GMDS_NODE,"classification");
    Variable<int>* v_surf = m_hexdom.getVariable<int>(GMDS_NODE,"surface_id");
    Variable<int>* v_curv = m_hexdom.getVariable<int>(GMDS_NODE,"curve_id");
    //===============================================================
    // We keep in mind the points that are on curve or surface and
    // which are unused
    std::vector<math::Point> free_point_on_curve;
    std::vector<math::Point> free_point_on_vertex;
    
    
    for(auto i=0; i<m_pnt.size();i++){
        if(m_used[i])
            continue;
        
        //So we have an unused point
        if(m_classification[i]==1){
            free_point_on_curve.push_back(m_pnt[i]);
        }
        else if(m_classification[i]==0){
            free_point_on_vertex.push_back(m_pnt[i]);
        }
    }

    
    //===============================================================
    // STEP 1 - WE EXTRACT THE TETRA MESH BOUNDARY WITH THE RIGHT
    // ORIENTATION
    //===============================================================
    BoundaryOperator boundaryOp(m_mesh);
    IGMesh bnd_mesh(MeshModel(DIM3|F|F2N|N2F));
    Variable<int>* bnd_curve_color = bnd_mesh.newVariable<int>(GMDS_NODE, "BND_CURVE_COLOR");
    Variable<int>* bnd_surf_color  = bnd_mesh.newVariable<int>(GMDS_FACE, "BND_SURFACE_COLOR");
    Variable<int>* bnd_vertex_color= bnd_mesh.newVariable<int>(GMDS_NODE, "BND_VERTEX_COLOR");
    Variable<int>* surf_color      = m_mesh->getVariable<int>(GMDS_FACE, "BND_SURFACE_COLOR");
    Variable<int>* curve_color     = m_mesh->getVariable<int>(GMDS_NODE, "BND_CURVE_COLOR");
    Variable<int>* vertex_color    = m_mesh->getVariable<int>(GMDS_NODE,"BND_VERTEX_COLOR");
    std::vector<TCellID> bnd_faces;
    std::set<TCellID> bnd_nodes;
    
    Variable<math::Vector>* var_normal =
    bnd_mesh.newVariable<math::Vector>(GMDS_FACE,"NORMAL");
    
    for(IGMesh::face_iterator itf = m_mesh->faces_begin();
        !itf.isDone(); itf.next()){
        if(itf.value().getIDs<Region>().size()==1){
            bnd_faces.push_back(itf.value().getID());
            std::vector<TCellID> bnd_n = itf.value().getIDs<Node>();
            bnd_nodes.insert(bnd_n.begin(),bnd_n.end());
        }
    }
    std::map<TCellID,TCellID> vol_2_surf;
    for(auto n:bnd_nodes){
        Node n_surf = bnd_mesh.newNode(m_mesh->get<Node>(n).getPoint());
        vol_2_surf[n]=n_surf.getID();
        (*bnd_curve_color)[n_surf.getID()]= (*curve_color)[n];
        (*bnd_vertex_color)[n_surf.getID()]= (*vertex_color)[n];
    }
    for(auto surf:bnd_faces){
        Face from = m_mesh->get<Face>(surf);
        std::vector<TCellID> surf_nodes = from.getIDs<Node>();
        
        // We create the triangle in the right direction
        math::Vector3d n = getOutputNormalOfABoundaryFace(from);
        std::vector<Node> fn = from.get<Node>();
        math::Vector3d v01(fn[0].getPoint(),fn[1].getPoint());
        math::Vector3d v02(fn[0].getPoint(),fn[2].getPoint());
        Face to;
        if(n.dot(v01.cross(v02))<0){
            //we reverse the numerotation
            to = bnd_mesh.newTriangle(vol_2_surf[surf_nodes[0]],
                                      vol_2_surf[surf_nodes[2]],
                                      vol_2_surf[surf_nodes[1]]);
            
        }
        else{
            to = bnd_mesh.newTriangle(vol_2_surf[surf_nodes[0]],
                                      vol_2_surf[surf_nodes[1]],
                                      vol_2_surf[surf_nodes[2]]);
            
        }
        math::Vector v1(to.get<Node>()[0].getPoint(),
                        to.get<Node>()[1].getPoint());
        math::Vector v2(to.get<Node>()[0].getPoint(),
                        to.get<Node>()[2].getPoint());
        math::Vector vn=v1.cross(v2);
        vn.normalize();
        (*var_normal)[to.getID()]=vn;
        (*bnd_surf_color)[to.getID()]= (*surf_color)[from.getID()];
    }
    
    
    
    IGMeshDoctor doc(&bnd_mesh);
    doc.updateUpwardConnectivity();

    Variable<int>* var_loop_faces = bnd_mesh.newVariable<int>(GMDS_FACE,"loop_faces");
    static int file_nb=0;
    if(m_param_gl.with_debug_files){

        VTKWriter<IGMesh> bnd_writer(bnd_mesh);
        std::string file_name = m_param_gl.output_dir+"/SURFACE_MESH_"+std::to_string(file_nb);
        bnd_writer.write(file_name,F|N);
        //===============================================================
        {
            IGMesh normal_mesh(MeshModel(DIM3|F|F2N));
            for(IGMesh::face_iterator itf = bnd_mesh.faces_begin();
                !itf.isDone(); itf.next()){
                Face f = itf.value();
                math::Point c = f.center();
                math::Vector v = (*var_normal)[f.getID()];
                math::Point c2 = c-v;
                Node n1 = normal_mesh.newNode(c);
                Node n2 = normal_mesh.newNode(c2);
                normal_mesh.newTriangle(n1,n1,n2);
            }
            
            VTKWriter<IGMesh> nw(normal_mesh);
            file_name = m_param_gl.output_dir+"/NORMAL_SURF_MESH_"+std::to_string(file_nb);
            nw.write(file_name,F|N);
        }
    }
    //===============================================================
    // STEP 2 - We insert each loop into the surface boundary
    //===============================================================
    //When we create a patch some boundary nodes are those coming
    //from the hex mesh, and some other are on surface boundary
    //we put last ones into this vector for comparison
    
    std::vector<Node> new_bnd_patch_nodes;
    
    for(auto l:ALoops){
        //=============================================================
        //we prepare the data to provide to the insertion algorithm
        std::vector<math::Point> p;
        std::vector<int> classification;
        std::vector<int> geom_id; //meaningful for curve and surface only
        p.resize(l.size());
        classification.resize(l.size());
        geom_id.resize(l.size());
        for(auto j=0; j<l.size();j++){
            p[j]             = l[j].getPoint();
            classification[j]= (*v_cl)[l[j].getID()];
            
            if(classification[j]==2){
                geom_id[j]= (*v_surf)[l[j].getID()];
            }
            else if(classification[j]==1){
                geom_id[j]= (*v_curv)[l[j].getID()];
            }
        }
        
        //=============================================================
        //We perform the insertion
        TriangularSurfaceManipulator tm(&bnd_mesh);
        std::vector<TCellID> loop_nodes;
        std::vector<TCellID> loop_faces;
        
        std::vector<bool> loop_nodes_from_pnt;
        
        tm.insertLoop(p,classification,geom_id,
                      loop_faces, loop_nodes,
                      loop_nodes_from_pnt);
        
        //=============================================================
        // We keep in mind which point of the inserted loop is coming
        // from the loop point given as an input (indeed new points
        // resulting from intersection with triangular edges were
        // added)
        loop_nodes_from_pnt.clear();
        loop_nodes_from_pnt.resize(false,loop_nodes.size());
        for(auto in=0; in<loop_nodes.size();in++){
            Node ni =bnd_mesh.get<Node>(loop_nodes[in]);
            math::Point pi = ni.getPoint();
            double tol =1e-4;
            bool found_from=false;
            for(auto ip=0; ip<p.size()&& !found_from;ip++){
                math::Point pp = p[ip];
                if(pp.distance2(pi)<tol)
                    found_from=true;
            }
            loop_nodes_from_pnt[in]=found_from;
        }
        
        //=============================================================
        //Now we add points coming from unused points if they are in the
        // vicinity of loop_faces
        double min_xyz[3]={DBL_MAX, DBL_MAX, DBL_MAX};
        double max_xyz[3]={DBL_MIN, DBL_MIN, DBL_MIN};
        for(auto cur_id:loop_nodes){
            math::Point cur_pnt = bnd_mesh.get<Node>(cur_id).getPoint();
            if(cur_pnt.X()<min_xyz[0]){
                min_xyz[0]= cur_pnt.X();
            }
            else if(cur_pnt.X()>max_xyz[0]){
                max_xyz[0]= cur_pnt.X();
            }
            
            if(cur_pnt.Y()<min_xyz[1]){
                min_xyz[1]= cur_pnt.Y();
            }
            else if(cur_pnt.Y()>max_xyz[1]){
                max_xyz[1]= cur_pnt.Y();
            }
            
            if(cur_pnt.Z()<min_xyz[2]){
                min_xyz[2]= cur_pnt.Z();
            }
            else if(cur_pnt.Z()>max_xyz[2]){
                max_xyz[2]= cur_pnt.Z();
            }
        }
        min_xyz[0]=min_xyz[0]-0.1;
        min_xyz[1]=min_xyz[1]-0.1;
        min_xyz[2]=min_xyz[2]-0.1;
        
        max_xyz[0]=max_xyz[0]+0.1;
        max_xyz[1]=max_xyz[1]+0.1;
        max_xyz[2]=max_xyz[2]+0.1;
        
        std::vector<math::Point> to_insert;
        for(auto pnt:free_point_on_curve){
            if(pnt.X()>min_xyz[0] && pnt.X()<max_xyz[0] &&
               pnt.Y()>min_xyz[1] && pnt.Y()<max_xyz[1] &&
               pnt.Z()>min_xyz[2] && pnt.Z()<max_xyz[2]){
                to_insert.push_back(pnt);
            }
        }
        for(auto pnt:free_point_on_vertex){
            if(pnt.X()>min_xyz[0] && pnt.X()<max_xyz[0] &&
               pnt.Y()>min_xyz[1] && pnt.Y()<max_xyz[1] &&
               pnt.Z()>min_xyz[2] && pnt.Z()<max_xyz[2]){
                to_insert.push_back(pnt);
            }
        }
        
        if(m_param_gl.with_debug_files){
            VTKWriter<IGMesh> bnd_writer(bnd_mesh);
            bnd_writer.write(m_param_gl.output_dir+
                             "/SURFACE_MESH_INSERTION_BEF_FREE_PNTS"+
                             std::to_string(file_nb),F|N);
        }
        //insertion of the free points into the mesh
        std::vector<Node> free_nodes = tm.insertFreePoints(to_insert,
                                                           loop_faces);
        loop_faces = tm.extractEnclosedFaces(loop_nodes);

        for(auto id:loop_faces){
            (*var_loop_faces)[id]=1;
        }
        
        if(m_param_gl.with_debug_files){
            VTKWriter<IGMesh> bnd_writer(bnd_mesh);
            bnd_writer.write(m_param_gl.output_dir+
                             "/SURFACE_MESH_INSERTION_"+
                             std::to_string(file_nb),F|N);
            
            IGMesh normal_mesh(MeshModel(DIM3|F|F2N));
            for(IGMesh::face_iterator itf = bnd_mesh.faces_begin();
                !itf.isDone(); itf.next()){
                Face f = itf.value();
                math::Point c = f.center();
                math::Vector v = (*var_normal)[f.getID()];
                math::Point c2 = c-v;
                Node n1 = normal_mesh.newNode(c);
                Node n2 = normal_mesh.newNode(c2);
                normal_mesh.newTriangle(n1,n1,n2);
            }
            VTKWriter<IGMesh> nw(normal_mesh);
            nw.write(m_param_gl.output_dir+
                     "/NORMAL_SURF_MESH_INSERTION_"+std::to_string(file_nb),F|N);
        }
        //===============================================================
        // We split the loop faces depending on the surface they belong
        // to
        std::map<int,std::vector<TCellID> > surf_to_faces_full;
        for(auto id:loop_faces){
            surf_to_faces_full[(*bnd_surf_color)[id]].push_back(id);
        }
        
        //===============================================================
        //now for faces associated to a surface, we split them into
        //connex components
        std::map<int,std::vector<TCellID> > surf_to_faces;

        int surf_index=1;
        for(auto surf:surf_to_faces_full){
            std::vector<TCellID> full_faces = surf.second;
            //starting from the full set of faces, we generate patchs
            //via a growing process
            
            // Faces are marked as soon as added in a patch
            std::map<TCellID,bool> already_done;
            for(auto id:full_faces){
                already_done[id]=false;
            }
            TCellID free_id = foundFreeIndex(already_done);
            std::set<TCellID> cur_faces;
            while(free_id!=NullID){
                //we have our starting face;
                std::vector<TCellID> heap;
                heap.push_back(free_id);
                while(!heap.empty()){
                    TCellID cur_id = heap.back();
                    heap.pop_back();
                    already_done[cur_id]=true;
                
                    Face current = bnd_mesh.get<Face>(cur_id);
                    cur_faces.insert(cur_id);
                    
                    //Get adjacent faces which are in the same patch too
                    std::vector<Node> cur_nodes = current.get<Node>();
                    for(auto i=0; i<3; i++){
                        Node ni = cur_nodes[i];
                        Node nj = cur_nodes[(i+1)%3];
                        std::vector<TCellID> fij =getFaces(ni,nj);
                        for(auto f:fij){
                            if(f==cur_id)
                                continue;
                            
                            if(!isIn(f,full_faces))
                                continue;
                            
                            if(already_done[f]==true)
                                continue;
                            
                            //f is the id of a face of this patch and
                            //is not the current face, and not yet
                            //treated
                            heap.push_back(f);
                        }
                    }
                }
                std::vector<TCellID> vec_faces;
                vec_faces.insert(vec_faces.end(),cur_faces.begin(),cur_faces.end());
                
                surf_to_faces[surf_index]=vec_faces;
                cur_faces.clear();
                surf_index++;
                free_id=foundFreeIndex(already_done);
            }

        }
        
        //===============================================================
        // Each node of loop_nodes with a flag "true" means it correspond
        // to a node of the hex. mesh
        std::map<TCellID,TCellID> hex2tri_nodes;
        std::map<TCellID,TCellID> tri2hex_nodes;
        std::set<TCellID> bnd_nodes_for_paving;
        std::set<TCellID> bnd_nodes_for_paving_on_loop;
        //Nodes of l belongs to the hex_dom mesh
        int i_hex=0;

        for(auto i_tri=0; i_tri<loop_nodes.size();i_tri++){
            if(loop_nodes_from_pnt[i_tri]==false)
                continue;
            
            //the current node comes from a hex point so
            Node hex_node = l[i_hex];
            hex2tri_nodes[hex_node.getID() ] = loop_nodes[i_tri];
            tri2hex_nodes[loop_nodes[i_tri]] = hex_node.getID();
            i_hex++;
            bnd_nodes_for_paving.insert(loop_nodes[i_tri]);
            bnd_nodes_for_paving_on_loop.insert(loop_nodes[i_tri]);

        }
        //and free nodes too
        for(auto n:free_nodes){
            bnd_nodes_for_paving.insert(n.getID());
        }
        
        //===============================================================
        // We split the loop of nodes according to the change of surfaces
        std::map<int, std::vector<TCellID> > bnd_patch_per_surf;
        for(auto surf:surf_to_faces){
            std::vector<TCellID> faces = surf.second;
//            std::cout<<"("<<surf.first<<"): ";
//            for(auto fi:faces){
//                std::cout<<fi<<" ";
//            }
//            std::cout<<std::endl;
            //We get the boundary nodes of faces. To do it we get boundary
            // fake edges
            std::map<FakeEdge::EdgeID, int> edge_count;
            for(auto f_id:faces){
                Face f = bnd_mesh.get<Face>(f_id);
                std::vector<TCellID> fn = f.getIDs<Node>();
                FakeEdge e01(fn[0],fn[1]);
                FakeEdge e12(fn[1],fn[2]);
                FakeEdge e20(fn[2],fn[0]);
                edge_count[e01.getID()]++;
                edge_count[e12.getID()]++;
                edge_count[e20.getID()]++;
            }
            
            //We keep only the edges with a count of 1 (boundary so)
            std::vector<FakeEdge> bnd_edges;
//            std::cout<<"======= BND EDGES ======="<<std::endl;
            for(auto ec:edge_count){
                if(ec.second==1){
                    bnd_edges.push_back(FakeEdge(ec.first.getID1(),
                                                 ec.first.getID2()));
//                    std::cout<<"("<<ec.first.getID1()<<", "<<
//                    ec.first.getID2()<<")\n";
                }
            }
            
            //Now we order the edges as a loop starting from the first one
            FakeEdge first_edge = bnd_edges[0];
            FakeEdge cur_edge = first_edge;
            std::vector<TCellID> ordered_nodes;
            ordered_nodes.reserve(bnd_edges.size());
            TCellID first_id = first_edge.getID().getID1();
            TCellID cur_id   = first_edge.getID().getID2();
            ordered_nodes.push_back(first_edge.getID().getID1());
            while(cur_id!=first_id){
                ordered_nodes.push_back(cur_id);
                bool find_next=false;
                for(auto i=0;i<bnd_edges.size() && !find_next; i++){
                    FakeEdge ei = bnd_edges[i];
                    if(ei==cur_edge)
                        continue;
                    
                    if(ei.getID().getID1()==cur_id){
                        cur_id=ei.getID().getID2();
                        cur_edge= ei;
                        find_next=true;
                    }
                    else if(ei.getID().getID2()==cur_id){
                        cur_id=ei.getID().getID1();
                        cur_edge= ei;
                        find_next=true;
                    }
                }
                
            }
            //the orientation can be the wrong one. We use the normal
            //of adjacent enclosed faces to check that
            Node n0  = bnd_mesh.get<Node>(ordered_nodes[0]);
            Node n1 = bnd_mesh.get<Node>(ordered_nodes[1]);
            std::vector<TCellID> f01, f0, f1;
            f0 = n0.getIDs<Face>();
            f1 = n1.getIDs<Face>();
            for(auto id0:f0){
                for(auto id1:f1){
                    if(id0==id1){
                        f01.push_back(id0);
                    }
                }
            }
            TCellID adj_enclosed_face = NullID;
            for(auto id01:f01){
                if(std::find(faces.begin(),faces.end(),id01)!=faces.end()){
                    adj_enclosed_face=id01;
                }
            }
            if(adj_enclosed_face==NullID){
                throw GMDSException("Orientation issue");
            }
            
            std::vector<TCellID> ref = bnd_mesh.get<Face>(adj_enclosed_face).getIDs<Node>();
            if((ref[0]==n1.getID() && ref[1]==n0.getID()) ||
               (ref[1]==n1.getID() && ref[2]==n0.getID()) ||
               (ref[2]==n1.getID() && ref[0]==n0.getID()) ){
                //opposite direction
                std::reverse(ordered_nodes.begin(),ordered_nodes.end());
//                std::cout<<"\t REVERSE BND"<<std::endl;
               }
            // Now we have an ordered list of nodes enclosing the surface
            // patch to be meshed.
            // We keep now the nodes corresponding to hex nodes
            std::vector<TCellID> bnd_paving;
            for(auto id:ordered_nodes){
                if(bnd_nodes_for_paving.find(id)!=bnd_nodes_for_paving.end()){
                    bnd_paving.push_back(id);
//                    std::cout<<"\t in "<<id<<std::endl;
                }
            }
            bnd_patch_per_surf[surf.first]=bnd_paving;
        }

        //===============================================================
        // EACH PATCH/SURFACE IS NOW MESHED VIA PAVING
        //===============================================================
        for(auto bnd_nodes:bnd_patch_per_surf){
            //loop of nodes
            std::vector<TCellID> loop_bnd_node_ids = bnd_nodes.second;
           
            //surface rep
            std::vector<TCellID> faces = surf_to_faces[bnd_nodes.first];

            std::vector<Face> surf_faces;
            std::vector<math::Vector> surf_normals;

            surf_faces.reserve(faces.size());
            surf_normals.reserve(faces.size());
            for(auto idx:faces){
                surf_faces.push_back(bnd_mesh.get<Face>(idx));
                surf_normals.push_back((*var_normal)[idx]);
            }
            std::vector<math::Point> bnd_points;
            bnd_points.resize(loop_bnd_node_ids.size());
            for(auto i=0; i<loop_bnd_node_ids.size();i++){
                bnd_points[i] = bnd_mesh.get<Node>(loop_bnd_node_ids[i]).getPoint();
            }
            
            
            IGMesh paving_mesh(MeshModel(DIM3|F|F2N));

            //Call to the paving algorithm
            meshSurfacePatch(surf_faces,surf_normals, bnd_points,
                             paving_mesh);
            
            //paving_mesh contains al local view of the patch mesh.
            //Now we look for nodes of paving_mesh that corresponds to nodes of the initial mesh
            
            //nodes are new or correspond to nodes in l or in new_bnd_patch_nodes
            std::map<TCellID, TCellID> paving_to_hex_mesh_nodes;
            double tolerance = 1e-4;
            for(IGMesh::node_iterator itn = paving_mesh.nodes_begin();
                !itn.isDone(); itn.next()){
                Node ni = itn.value();
                math::Point pi = ni.getPoint();
                
                //First we compare to the location of hex nodes
                bool come_from_hex=false;
                for(auto i_l = 0; i_l<l.size() && !come_from_hex; i_l++){
                    math::Point pl = l[i_l].getPoint();
                    if(pi.distance2(pl)<tolerance){
                        come_from_hex=true;
                        paving_to_hex_mesh_nodes[ni.getID()]=l[i_l].getID();
                    }
                }
                
                //Second we compare to the location of previous added nodes in the
                //hex dom mesh
                bool come_from_prev_add=false;
                for(auto i_l = 0;
                    i_l<new_bnd_patch_nodes.size() && !come_from_prev_add; i_l++){
                    math::Point pl = new_bnd_patch_nodes[i_l].getPoint();
                    if(pi.distance2(pl)<tolerance){
                        come_from_prev_add=true;
                        paving_to_hex_mesh_nodes[ni.getID()]=new_bnd_patch_nodes[i_l].getID();
                    }
                }
                
                //Finally we add the point if we didnt encounter it
                if(!come_from_hex && ! come_from_prev_add){
                    Node new_ni = m_hexdom.newNode(pi);
                    paving_to_hex_mesh_nodes[ni.getID()]= new_ni.getID();
                    new_bnd_patch_nodes.push_back(new_ni);
                }
            }
            //Now we add faces!!
            for(IGMesh::face_iterator itf = paving_mesh.faces_begin();
                !itf.isDone(); itf.next()){
                Face fi = itf.value();
                std::vector<TCellID> node_ids = fi.getIDs<Node>();
                
                if(node_ids.size()==3){
                    //new triangle
                    Face new_fi = m_hexdom.newTriangle(paving_to_hex_mesh_nodes[node_ids[0]],
                                                       paving_to_hex_mesh_nodes[node_ids[1]],
                                                       paving_to_hex_mesh_nodes[node_ids[2]]);
                }
                else if(node_ids.size()==4){
                    //new quad
                    Face new_fi = m_hexdom.newQuad(paving_to_hex_mesh_nodes[node_ids[0]],
                                                   paving_to_hex_mesh_nodes[node_ids[1]],
                                                   paving_to_hex_mesh_nodes[node_ids[2]],
                                                   paving_to_hex_mesh_nodes[node_ids[3]]);
                }
            }
            
        }
        file_nb++;
        
    }

    return false;
}
/*---------------------------------------------------------------------------*/
void HDMeshGenerator::
meshSurfacePatch(std::vector<Face>& ASurfFaceRep,
                 std::vector<math::Vector>& ASurfNormalRep,
                 std::vector<math::Point>& ABndPoints,
                 IGMesh& APatch)
{
    static int file_nb=0;

    CMLSurfEvalImpl surf_rep(ASurfFaceRep, ASurfNormalRep);
    CMLPaver paver(&surf_rep);
    paver.set_mixed_elements(true);
    paver.set_sizing_function(CML::LINEAR_SIZING);
    //===============================================================
    // Then, we build a quad mesh using only the boundary nodes coming
    // from the hexahadral elements
    
    int nb_points = ABndPoints.size();
    double* points = new double[3*nb_points];
    
    int* loop = new int[nb_points];
    for(auto i_p=0; i_p<nb_points; i_p++){

        math::Point pi =ABndPoints[i_p];
        points[3*i_p]=pi.X();
        points[3*i_p+1]=pi.Y();
        points[3*i_p+2]=pi.Z();
        loop[i_p]=i_p;
    }
    int nb_loops=1;
    int loop_size[1] ={nb_points};
  //  std::cout<<"Paver set boundary mesh"<<std::endl;
    
    paver.set_boundary_mesh(nb_points, points,
                            nb_loops,loop_size,loop);
    int nb_gen_pnts;
    int nb_gen_quads;
    paver.generate_mesh(nb_gen_pnts,nb_gen_quads);
   // std::cout<<"PAVER (P,Q)= ("<<nb_gen_pnts<<", "<<nb_gen_quads<<")"<<std::endl;
    
    double* gen_pnts = new double[3*nb_gen_pnts];
    int* gen_quads = new int[4*nb_gen_quads];
    paver.get_mesh(nb_gen_pnts, gen_pnts, nb_gen_quads, gen_quads);
    
    for(auto i=0;i<nb_gen_pnts;i++){
        APatch.newNode(gen_pnts[3*i  ],
                       gen_pnts[3*i+1],
                       gen_pnts[3*i+2]);
    }
    for(auto i=0; i<nb_gen_quads;i++){
        // The algorith can create degenerate quads and not triangles
        //We check that now
        std::set<int> points_index;
        points_index.insert(gen_quads[4*i]);
        points_index.insert(gen_quads[4*i+1]);
        points_index.insert(gen_quads[4*i+2]);
        points_index.insert(gen_quads[4*i+3]);
        if(points_index.size()==4){
        APatch.newQuad(gen_quads[4*i],
                       gen_quads[4*i+1],
                       gen_quads[4*i+2],
                       gen_quads[4*i+3]);
        }
        else  if(points_index.size()==3){
            int index[3];
            int pi=0;
            for(auto p:points_index){
                index[pi++]=p;
            }
            APatch.newTriangle(index[0], index[1], index[2]);
        }
        else{
            throw GMDSException("Wrong number of nodes for a face generated by paving");
        }
            
    }
    
    if(m_param_gl.with_debug_files){
        VTKWriter<IGMesh> nw(APatch);
        nw.write(m_param_gl.output_dir+"/PAVER"+std::to_string(file_nb),F|N);
        //===============================================================
        
        IGMesh paving_surf_mesh(MeshModel(DIM3|F|F2N));
        for(auto origin:ASurfFaceRep){
            std::vector<Node> origin_nodes = origin.get<Node>();
            Node n0 = paving_surf_mesh.newNode(origin_nodes[0].getPoint());
            Node n1 = paving_surf_mesh.newNode(origin_nodes[1].getPoint());
            Node n2 = paving_surf_mesh.newNode(origin_nodes[2].getPoint());
            paving_surf_mesh.newTriangle(n0,n1,n2);
        }
        VTKWriter<IGMesh> nws(paving_surf_mesh);
        nws.write(m_param_gl.output_dir+
                  "/PAVER_SURF"+std::to_string(file_nb),F|N);
        
        //===============================================================
        
        IGMesh normal_mesh(MeshModel(DIM3|F|F2N));
        for(auto i=0;i<ASurfFaceRep.size();i++){
            Face f = ASurfFaceRep[i];
            math::Point c = f.center();
            math::Vector v = ASurfNormalRep[i];
            math::Point c2 = c+v;
            Node n1 = normal_mesh.newNode(c);
            Node n2 = normal_mesh.newNode(c2);
            normal_mesh.newTriangle(n1,n1,n2);
        }
        VTKWriter<IGMesh> nwn(normal_mesh);
        nwn.write(m_param_gl.output_dir+
                  "/PAVER_NORMAL_"+std::to_string(file_nb),F|N);
        
        file_nb++;
    }

    
    delete[] gen_pnts;
    delete[] gen_quads;
    delete[] points;
    delete[] loop;
}
/*---------------------------------------------------------------------------*/
std::vector<TCellID> HDMeshGenerator::
getFaces(const gmds::Node& ANI, const gmds::Node& ANJ) const
{
    std::vector<TCellID> fij, fi, fj;
    fi = ANI.getIDs<Face>();
    fj = ANJ.getIDs<Face>();
    for(auto i:fi){
        for(auto j:fj){
            if(i==j){
                fij.push_back(i);
            }
        }
    }
    return fij;
}
/*---------------------------------------------------------------------------*/
bool HDMeshGenerator::isIn(const TCellID AI,
                           const std::vector<TCellID>& AV) const
{
    return std::find(AV.begin(),AV.end(),AI)!=AV.end();
}

/*---------------------------------------------------------------------------*/
TCellID HDMeshGenerator::
foundFreeIndex(const std::map<TCellID, bool>& AVec) const
{
    for(auto elt:AVec){
        if(elt.second==false)
            return elt.first;
    }
    return NullID;
}
/*---------------------------------------------------------------------------*/
Face HDMeshGenerator::getFace(const Face& AFrom,
                              const Region& AR,
                              const TCellID ANI,
                              const TCellID ANJ)
{
    std::vector<Face> faces = AR.get<Face>();
    for(auto f:faces){
        std::vector<TCellID> fn = f.getIDs<Node>();
        auto found_ni = false, found_nj=false;
        for(auto n_id:fn){
            if(n_id==ANI)
                found_ni=true;
            else if(n_id==ANJ)
                found_nj=true;
        }
        if(found_ni && found_nj && f.getID()!=AFrom.getID())
            return f;
    }
    throw GMDSException("HDMeshGenerator::getFace(..) - No next face found");
}
/*---------------------------------------------------------------------------*/
bool HDMeshGenerator::isIn(const int AFrom,
                           const int ATo,
                           const std::vector<OrientedEdge>& AEdgeSet,
                           OrientedEdge& AOutEdge) const
{
    OrientedEdge witness(AFrom,ATo);
    for(unsigned int i=0; i<AEdgeSet.size(); i++){
        if(witness  == AEdgeSet[i]){
            AOutEdge = AEdgeSet[i];
            return true;
        }
    }
    return false;
}
/*---------------------------------------------------------------------------*/
void HDMeshGenerator::
buildEdges(std::vector<std::vector<OrientedEdge> >& AInEdges,
           std::vector<std::vector<OrientedEdge> >& AOutEdges)
{
    int nb_points = m_pnt.size();
    
    AOutEdges.clear();
    AOutEdges.resize(nb_points);
    std::vector<int> on_curves;
    
    //======================================================================
    // STEP 1 - BUILD GLOBAL EDGES
    //======================================================================
    for(int i=0; i<nb_points; i++){
        
        if(m_classification[i]==ON_CURVE)
            on_curves.push_back(i);
        
        std::vector<OrientedEdge> orient_edges = AInEdges[i];
        
        for(auto e:orient_edges){
            int from = e.first;
            int to   = e.second;
            OrientedEdge e_out;
            //===============================================================
            //We check if ei has already been put in the out set of edges
            if(isIn(from,to,AOutEdges[from], e_out)){
                //YES, nothing to do so
                continue;
            }
            if(m_classification[from]>m_classification[to]){
                // Wse build the edge systematically
                // And the edge can be kept as the opposite
                AOutEdges[from].push_back(e);
                
                OrientedEdge e_inv(to,from);
                AOutEdges[to].push_back(e_inv);
            }
            else if(m_classification[from]==m_classification[to]){
                //We connect only if the connection exists in both sides
                OrientedEdge e_inv;
                if(isIn(to, from,
                        AInEdges[to], e_inv)){
                    //YES, We connect
                    AOutEdges[from].push_back(e);
                    AOutEdges[to  ].push_back(e_inv);
                }
            }
        }//for(unsigned int j=0; j<oriented_edges_i.size(); j++)
        
    }//for(int i=0; i<nb_points; i++)
    //======================================================================
    // STEP 2 - CLEAN SOME EDGES
    //======================================================================
    for(auto i:on_curves){
        //i= index of a point
        std::vector<OrientedEdge> edges = AOutEdges[i];
        std::vector<OrientedEdge> final_edges;
        for(auto j=0; j<edges.size(); j++){
            OrientedEdge ej = edges[j];
            
            if(m_classification[ej.second]==ON_CURVE ||
               m_classification[ej.second]==ON_VERTEX ){
                //Another edge connected to a surface can be a better choice
                bool found_best=false;
                for(auto k=0; k<edges.size(); k++){
                    OrientedEdge ek = edges[k];
                    if(k==j ||
                       m_classification[ek.second]==ON_CURVE ||
                       m_classification[ek.second]==ON_VERTEX){
                        continue;
                    }
                    //We compare length and direction of ej and ek
                    math::Vector3d vj(m_pnt[ej.first], m_pnt[ej.second]);
                    math::Vector3d vk(m_pnt[ek.first], m_pnt[ek.second]);
                    if(vk.norm2()<vj.norm2()){
                        vj.normalize();
                        vk.normalize();
                        if(vj.dot(vk)>0.8)
                            found_best=true;
                    }
                }//for(auto k=0; k<edges.size(); k++)
                if(!found_best){
                    final_edges.push_back(ej);
                }
            }
            else{
                final_edges.push_back(ej);
            }
        }//for(auto j=0; j<edges.size(); j++)
        
        //The cleaned set of edges is assigned to i
        AOutEdges[i]=final_edges;
        
    }
}
/*---------------------------------------------------------------------------*/
int HDMeshGenerator::
filterPointsForBuildingOrientedEdge(const std::vector<int>& AID,
                                    const std::vector<double>& ADot,
                                    const std::vector<double>& ADist){
    //======================================================================
    //We get the best aligned candidate
    //======================================================================

    int    best_candidate_id =-1;
    double best_dot          = 0;
    for(auto i_c=0; i_c<AID.size();i_c++){
        if(ADot[i_c]>best_dot){
            best_dot=ADot[i_c];
            best_candidate_id=i_c;
        }
    }
    //Then we check if a closer one is in a reasonable angle tolerance
    int            ref_id     = best_candidate_id;
    double         ref_dist   = ADist[ref_id];
    
   // std::cout<<"\t best aligned: "<<best_candidate_id<<std::endl;
    
    int    new_best_id   = ref_id;
    double new_best_dist = ref_dist;
    
    for(unsigned int i_c=0; i_c<AID.size();i_c++){
        //Only closer point deserve to be checked
        if(i_c==ref_id || ADist[i_c]>=new_best_dist)
            continue;
        
        double dot_prod = abs(best_dot-ADot[i_c]);
        if((dot_prod<0.045) ||//5 degree of diff, we switch to this one,
           //we are quite close in angle
           (dot_prod<0.33 && ADist[i_c]<0.7*ref_dist) )// only switch if distance is really dimished
        {
            new_best_dist=ADist[i_c];
            new_best_id = i_c;
        }
    }//for(unsigned int i_c=0; i_c<final_candidates.size();i_c++)
    
    return new_best_id;
}
/*---------------------------------------------------------------------------*/
void HDMeshGenerator::
buildOrientedEdgesInVolume(const int            APntID,
                           const math::Vector3d AVec[][2],
                           math::Point          APnt[][2],
                           int                  AIndex[][2],
                           bool                 AFound[][2])
{
    int i = APntID;
    math::Point pi = m_pnt[i];
    
    // Means the point is in a stable area. So it is connected to 6 other
    // points if it is far way from a singularity area.
    double tol = 0.8;
    
    //======================================================================
    // STEP 1 - WE DETECT THE POSSIBLE POINTS
    //======================================================================
    for(int axis=0; axis<3; axis++){
        
        for(int dir=0; dir<2; dir++){
            
            math::Vector3d v_axis=AVec[axis][dir];
            v_axis.normalize(); //likely useless
            
            math::Point next_pi;
            if(!computeVolumePointFrom(i, v_axis,next_pi)){
                continue; //means we reached a FF-singular area
            }
            math::Vector3d v(pi,next_pi);
            std::vector<int> candidates = m_filter[i];
            math::Point pw(pi.X()+v.X(),
                           pi.Y()+v.Y(),
                           pi.Z()+v.Z());
            
            v.normalize();
            
            std::vector<int>            final_candidates;
            std::vector<double>         final_candidates_dot;
            std::vector<double>         final_candidates_dist;
            
            for(unsigned int j=0; j<candidates.size(); j++){
                
                math::Point pj = m_pnt[candidates[j]];
                math::Vector3d vij(pi,pj);
                vij.normalize();
                if(v.dot(vij)>tol) {
                    //pj is candidate so
                    final_candidates.push_back(candidates[j]);
                    final_candidates_dot.push_back(v.dot(vij));
                    final_candidates_dist.push_back(pi.distance(pj));
                    
                }//if(v.dot(vij)>m_dot_tolerance)
                
            }//for(unsigned int j=0; j<candidates.size(); j++)
            
            //=========================================================
            // Now, we have to parse candidates to find the best one.
            // We start from the best aligned one, and we decrease while satisfying
            // the m_spacing criterion.
            // But if two candidates are only 10° different, we select via
            // distance.
            if(final_candidates.empty()){
                //We don't have found a candidate point
                continue;
            }
            
            int j = filterPointsForBuildingOrientedEdge(final_candidates,
                                                        final_candidates_dot,
                                                        final_candidates_dist);
            
            
            APnt  [axis][dir] = m_pnt[final_candidates[j]];
            AFound[axis][dir] = true;
            AIndex[axis][dir] = final_candidates[j];
            
        }//for(int dir=0;dir<2;dir++)
        
    }//for(int axis=0; axis<3; axis++)
}
/*---------------------------------------------------------------------------*/
void HDMeshGenerator::
buildOrientedEdgesOnSurface(const int            APntID,
                            const math::Vector3d AVec[][2],
                            math::Point          APnt[][2],
                            int                  AIndex[][2],
                            bool                 AFound[][2])
{
    int i = APntID;
    math::Point pi = m_pnt[i];

    // Means the point is in a stable area. So it is connected to 6 other
    // points if it is far way from a singularity area.
    double tol = 0.8;
    
    //We are on a surface, we must take care of points outside of the domain
    math::Vector3d n=m_normal[i];
    n.normalize();
    int surf_id = m_surface[i];
    
    // as the point is on the boundary, the classification process assigned
    // it to a starting triangle
    if(m_mesh_data[i].dim!=2)
        throw GMDSException("No face in buildOrientedEdgesOnSurface");
    gmds::Face t = m_mesh->get<Face>(m_mesh_data[i].id);
    
    //======================================================================
    // STEP 1 - WE DETECT THE POSSIBLE POINTS, WHICH ARE RESTRICTED TO BE ON
    //          THE BOUNDARY
    //======================================================================
    for(auto axis=0; axis<3; axis++){
        
        for(auto dir=0; dir<2; dir++){
            
            math::Vector3d v_axis=AVec[axis][dir];
            v_axis.normalize(); //likely useless
            math::Vector3d v;
            if((v_axis.dot(n))>tol){
                //We go out of the domain
                continue;
            }
            else if((v_axis.dot(n))<-tol){
                //we go inside the domain
                math::Point next_pi;
                if(computeVolumePointFrom(i, v_axis,next_pi))
                    v=math::Vector3d(pi,next_pi);
                else
                    continue; //it means we reached a FF-singular area
            }
            else{
                //we are on the surface
                math::Point next_pi;
                if(computeSurfacePointFrom(i, v_axis,next_pi))
                    v=math::Vector3d(pi,next_pi);
                else
                    continue; //it means we reached a FF-singular area

            }
            math::Point pw(pi.X()+v.X(),
                           pi.Y()+v.Y(),
                           pi.Z()+v.Z());
            std::vector<int> candidates = m_filter[i];
            
            v.normalize();
            
            std::vector<int>            final_candidates;
            std::vector<double>         final_candidates_dot;
            std::vector<double>         final_candidates_dist;
            
            
            for(unsigned int j=0; j<candidates.size(); j++){
                math::Point pj = m_pnt[candidates[j]];
                math::Vector3d vij(pi,pj);
                vij.normalize();
                
                //only volume and boundary well-aligned points are added
                if(v.dot(vij)>tol){
                    
                    if(m_classification[candidates[j]]==ON_SURFACE){
                        //if the point is on a surface it must be
                        // on the same surface! No it can be a thin layer
                        
                        if ( m_surface[candidates[j]]==surf_id){
                            //If they are on the same surface, we must discard
                            //configuration with cylinder shapes.
                            math::Vector3d nj=m_normal[candidates[j]];
                            if(nj.dot(n)<=0.0)
                                continue;
                        }
                    }
                    
                    //pj is candidate so
                    final_candidates.push_back(candidates[j]);
                    final_candidates_dot.push_back(v.dot(vij));
                    final_candidates_dist.push_back(pi.distance(pj));
                }//if(v.dot(vij)>m_dot_tolerance)
                
            }//for(unsigned int j=0; j<candidates.size(); j++)
            //=========================================================
            // Now, we have to parse candidates to find the best one.
            // We start from the best aligned one, and we decrease while satisfying
            // the m_spacing criterion.
            // But if two candidates are only 10° different, we select via
            // distance.
            if(final_candidates.empty()){
                //We don't have found a candidate point
                continue;
            }
            
            int j = filterPointsForBuildingOrientedEdge(final_candidates,
                                                        final_candidates_dot,
                                                        final_candidates_dist);
            
            APnt  [axis][dir] = m_pnt[final_candidates[j]];
            AFound[axis][dir] = true;
            AIndex[axis][dir] = final_candidates[j];
        }//for(int dir=0;dir<2;dir++)
        
    }//for(int axis=0; axis<3; axis++)
    
}
/*---------------------------------------------------------------------------*/
void HDMeshGenerator::
buildOrientedEdgesOnCurve(const int            APntID,
                          const math::Vector3d AVec[][2],
                          math::Point          APnt[][2],
                          int                  AIndex[][2],
                          bool                 AFound[][2])
{
    int i = APntID;
    math::Point pi = m_pnt[i];

//    std::cout<<"ORIENTED EDGES FROM "<<APntID<<" ON CURVE"<<std::endl;
    // Means the point is in a stable area. So it is connected to 6 other
    // points if it is far way from a singularity area.
    double tol = 0.9;
    
    //======================================================================
    // STEP 1 - WE DETECT THE POSSIBLE POINTS, WHICH ARE RESTRICTED TO BE ON
    //          THE BOUNDARY
    //======================================================================
    for(int axis=0; axis<3; axis++){
        
        for(int dir=0; dir<2; dir++){
            
            math::Vector3d v=AVec[axis][dir];
            v.normalize(); //likely useless
            
            std::vector<int>            candidates = m_filter[i];
            std::vector<int>            final_candidates;
            std::vector<double>         final_candidates_dot;
            std::vector<double>         final_candidates_dist;

            for(unsigned int j=0; j<candidates.size(); j++){
                math::Point pj = m_pnt[candidates[j]];
                math::Vector3d vij(pi,pj);
                vij.normalize();
                
                //only point well-aligned and on the boundary are taken
                //into account
                if(v.dot(vij)>tol &&
                   (m_classification[candidates[j]]==ON_CURVE  ||
                    m_classification[candidates[j]]==ON_VERTEX)){
                       
                       
                       //pj is candidate so
                       final_candidates.push_back(candidates[j]);
                       final_candidates_dot.push_back(v.dot(vij));
                       final_candidates_dist.push_back(pi.distance(pj));
                   }//if(v.dot(vij)>m_dot_tolerance)
                
            }//for(unsigned int j=0; j<candidates.size(); j++)
            //=========================================================
            // Now, we have to parse candidates to find the best one.
            // We start from the best aligned one, and we decrease while satisfying
            // the m_spacing criterion.
            // But if two candidates are only 10° different, we select via
            // distance.
            if(final_candidates.empty()){
                //We don't have found a candidate point
                continue;
            }

            int j = filterPointsForBuildingOrientedEdge(final_candidates,
                                                        final_candidates_dot,
                                                        final_candidates_dist);
            APnt  [axis][dir] = m_pnt[final_candidates[j]];
            AFound[axis][dir] = true;
            AIndex[axis][dir] = final_candidates[j];

        }//for(int dir=0;dir<2;dir++)
        
    }//for(int axis=0; axis<3; axis++)
}
/*---------------------------------------------------------------------------*/
void HDMeshGenerator::
buildOrientedEdges(std::vector<std::vector<OrientedEdge> >& AEdges)
{
    //One vector of edges per point. This vector can be empty so.
    AEdges.clear();
    AEdges.resize(m_pnt.size());
    
    //We go through all the points and we store for each of them the list of
    // oriented edges
    for(auto i=0; i<m_pnt.size(); i++){

        math::Point pi = m_pnt[i];
        
        if(m_type[i]==FRAME_SING){
            //only regular and param sing points are considered for
            //hex formation
            continue;
        }
        //==================================================================
        //We get the 6 vectors starting from pi following its chart
        math::Chart ci = m_chart[i];
        math::Vector3d ci_vectors[3][2] = {
            {ci.VX(), -ci.VX()},
            {ci.VY(), -ci.VY()},
            {ci.VZ(), -ci.VZ()}
        };
        
        //==================================================================
        // 6 corresponding points in directions x [0], y[1] and z[2]
        // the point in [0] is in positive direction
        // the point in [1] is in negative direction
        // For instance p[1][1] is the pnt starting from pi following -ci[1]
        math::Point p[3][2];
        int         p_index[3][2];
        bool        p_found[3][2]={{false,false},{false,false},{false,false}};

        
        if(m_classification[i]==IN_VOLUME){
            buildOrientedEdgesInVolume(i, ci_vectors, p, p_index, p_found);
        }
        else if(m_classification[i]==ON_SURFACE){
            buildOrientedEdgesOnSurface(i, ci_vectors, p, p_index, p_found);;
        }
        else if(m_classification[i]==ON_CURVE){
            buildOrientedEdgesOnCurve(i, ci_vectors, p, p_index, p_found);
        }
        
        //======================================================================
        // STEP 2 - ORIENTED EDGES ARE NOW CREATED
        //======================================================================
        std::vector<OrientedEdge> edges_i;
        for(auto axis=0; axis<3; axis++){
            for(auto dir=0; dir<2; dir++){
                if(p_found[axis][dir]==false)
                    continue;
                
                //We have an edge to create
                OrientedEdge e(i,                 // first point index
                               p_index[axis][dir],// second point index
                               axis,              // used chart axis in i
                               dir                // used chart direction in i
                               );
                
                edges_i.push_back(e);
            
            } //for(int dir=0;dir<2;dir++)
            
        } //for(auto axis=0; axis<3; axis++)
        AEdges[i]=edges_i;
        
    }
}
/*---------------------------------------------------------------------------*/
void HDMeshGenerator::
buildHexCorners(std::vector<std::vector<OrientedEdge> >& AEdges)
{
    for(unsigned int i=0; i<m_pnt.size(); i++){
        math::Point pi = m_pnt[i];
        std::vector<OrientedEdge> edges = AEdges[i];
        //We get the 6 vectors starting from pi following its chart
        math::Chart ci = m_chart[i];
        math::Vector3d ci_vectors[3][2] = {
            {ci.VX(), -ci.VX()},
            {ci.VY(), -ci.VY()},
            {ci.VZ(), -ci.VZ()}
        };
        
        int  p_index[3][2];
        bool p_found[3][2]= {
            {false,false},
            {false,false},
            {false,false}
        };

        bool found_undef_axis =false;
        int nb_undef_axi=0;
        for(unsigned int j=0; j<edges.size();j++){
            OrientedEdge ej = edges[j];
            if(ej.axis==-1 || ej.dir==-1){
                found_undef_axis=true;
                nb_undef_axi++;
            }
            else{
                p_found[ej.axis][ej.dir]=true;
                p_index[ej.axis][ej.dir]=ej.second;
            }
        }
        if(found_undef_axis){
            buildCornersAsSolidAngles(i, edges);
        }
        else{
            //now hex corner structure are created from 1 to 8 for pi
            //CORNER 1 - (X,Y,Z)
            if(p_found[0][0] && p_found[1][0] && p_found[2][0]){
                addCorner(i,
                          p_index[0][0], ci_vectors[0][0],
                          p_index[1][0], ci_vectors[1][0],
                          p_index[2][0], ci_vectors[2][0]);
            }
            //CORNER 2 - (X,Y,-Z)
            if(p_found[0][0] && p_found[1][0] && p_found[2][1]){
                addCorner(i,
                          p_index[0][0], ci_vectors[0][0],
                          p_index[1][0], ci_vectors[1][0],
                          p_index[2][1], ci_vectors[2][1]);
            }
            //CORNER 3 - (X,Y,-Z)
            if(p_found[0][0] && p_found[1][1] && p_found[2][0]){
                addCorner(i,
                          p_index[0][0], ci_vectors[0][0],
                          p_index[1][1], ci_vectors[1][1],
                          p_index[2][0], ci_vectors[2][0]);
            }
            //CORNER 4 - (X,-Y,-Z)
            if(p_found[0][0] && p_found[1][1] && p_found[2][1]){
                addCorner(i,
                          p_index[0][0], ci_vectors[0][0],
                          p_index[1][1], ci_vectors[1][1],
                          p_index[2][1], ci_vectors[2][1]);
            }
            //CORNER 5 - (-X,Y,Z)
            if(p_found[0][1] && p_found[1][0] && p_found[2][0]){
                addCorner(i,
                          p_index[0][1], ci_vectors[0][1],
                          p_index[1][0], ci_vectors[1][0],
                          p_index[2][0], ci_vectors[2][0]);
            }
            //CORNER 6 - (-X,Y,-Z)
            if(p_found[0][1] && p_found[1][0] && p_found[2][1]){
                addCorner(i,
                          p_index[0][1], ci_vectors[0][1],
                          p_index[1][0], ci_vectors[1][0],
                          p_index[2][1], ci_vectors[2][1]);
            }
            //CORNER 7 - (-X,Y,-Z)
            if(p_found[0][1] && p_found[1][1] && p_found[2][0]){
                addCorner(i,
                          p_index[0][1], ci_vectors[0][1],
                          p_index[1][1], ci_vectors[1][1],
                          p_index[2][0], ci_vectors[2][0]);
            }
            //CORNER 8 - (-X,-Y,-Z)
            if(p_found[0][1] && p_found[1][1] && p_found[2][1]){
                addCorner(i,
                          p_index[0][1], ci_vectors[0][1],
                          p_index[1][1], ci_vectors[1][1],
                          p_index[2][1], ci_vectors[2][1]);
            }
        }
    }//for(unsigned int i=0; i<m_pnt.size(); i++)...
}

/*---------------------------------------------------------------------------*/
void HDMeshGenerator::
buildCornersAsSolidAngles(const int AOrigin,
                          std::vector<OrientedEdge> & AEdges)
{
    
    int origin_id = AOrigin;

    //=====================================================================
    // STEP 1 - Compute combination of 3 integers among [0;AEdges.size()]
    //=====================================================================
    std::vector<std::vector<int> > triplets;
    computeCombinations(AEdges.size(),3, triplets);

    //All edges come from the same points
    math::Point origin_pnt(m_pnt[origin_id]);
    //=====================================================================
    // STEP 2 - Filter valid combination, that is vector triplets such that
    // no other vector is inside their positive space quadrant
    //=====================================================================
    for(unsigned int i=0; i<triplets.size(); i++){
        std::vector<int> ti = triplets[i];
        math::Point p[3] = {
            m_pnt[AEdges[ti[0]].second],
            m_pnt[AEdges[ti[1]].second],
            m_pnt[AEdges[ti[2]].second]
        };
        
        // Points in p were only selected on a combinatorial way. If they
        // are geometrically lying on the same plan, they do not define a
        // valid corner
        
        math::Vector3d v12(origin_pnt,p[0]);
        math::Vector3d v13(origin_pnt,p[1]);
        math::Vector3d v14(origin_pnt,p[2]);

        
//      if (origin_pnt.areCoplanar(p[0], p[1], p[2])){
      if (fabs(v12.dot(v13.cross(v14)))<1e-11){
            continue;
        }
        
        bool valid = true;
        
        for(unsigned int j=0; j<AEdges.size() && valid; j++){
            
            int j_index = AEdges[j].second;
            
            //We don't check the vectors of the current basis
            if(j_index==AEdges[ti[0]].second ||
               j_index==AEdges[ti[1]].second ||
               j_index==AEdges[ti[2]].second)
                continue;
            
            math::Point pj =m_pnt[j_index];
            
            double coord[4];
            try {
                math::Point::computeBarycentric(origin_pnt,
                                                p[0], p[1], p[2], pj,
                                                coord[0], coord[1], coord[2],
                                                coord[3]);
            } catch (GMDSException& e) {
                valid=false;
            }
            if(coord[1]>0 && coord[2]>0 && coord[3]>0)
                valid=false;
            
            //Find where is vj comparing to v[0], v[1] and v[2]
        }
        if(valid){
            
            //So we can build a hex corner from this set of edges
            math::Vector3d v[3] = {
                math::Vector3d(origin_pnt,p[0]),
                math::Vector3d(origin_pnt,p[1]),
                math::Vector3d(origin_pnt,p[2])
            };
            addCorner(origin_id,
                      AEdges[ti[0]].second, v[0],
                      AEdges[ti[1]].second, v[1],
                      AEdges[ti[2]].second, v[2]);
        }
    }
    
    
    
}
/*---------------------------------------------------------------------------*/
void HDMeshGenerator::
addCorner(const int AIndex,
          const int AIndex1, const math::Vector3d& AV1,
          const int AIndex2, const math::Vector3d& AV2,
          const int AIndex3, const math::Vector3d& AV3)
{
    HexCorner corner;
    corner.p = AIndex;
    corner.index =m_hc_mapping[AIndex].size();
    corner.free=true;
    
    corner.adj[0] = AIndex1;
    corner.vec[0] = AV1;
    
    corner.adj[1] = AIndex2;
    corner.vec[1] = AV2;
    
    corner.adj[2] = AIndex3;
    corner.vec[2] = AV3;
    
    m_hc_mapping[AIndex].push_back(corner);
}
/*---------------------------------------------------------------------------*/
int HDMeshGenerator::findFreeCorners(const std::vector<HexCorner>& AIn,
                                     const int AFrom,
                                     std::vector<HexCorner>& AOut)
{
    AOut.clear();
    for(auto const& c:AIn){
        if(!c.free) {
            //this corner is already used for a building another hex
            continue;
        }
        //c must point to the point AFrom
        bool found = false;
        for(int i=0; i<3; i++){
            if(c.adj[i]==AFrom){
                found=true;
            }
        }
        if(found)
            AOut.push_back(c);
    }
    return AOut.size();
}
/*---------------------------------------------------------------------------*/
std::vector<int> HDMeshGenerator::findCommonPoints(const int AFrom,
                                                   const int AI,
                                                   const int AJ)
{
    std::vector<HexCorner> corners_i, corners_j;
    findFreeCorners(m_hc_mapping[AI], AFrom, corners_i);
    findFreeCorners(m_hc_mapping[AJ], AFrom, corners_j);
    
    std::vector<int> common;
    
    for(auto ci:corners_i){
        for(auto cj:corners_j){
            
            for(int i=0; i<3; i++) {
                if(ci.adj[i]==AFrom)
                    continue;
                for(int j=0; j<3; j++) {
                    if(cj.adj[j]==AFrom)
                        continue;
                    
                    if(ci.adj[i]==cj.adj[j])
                        common.push_back(ci.adj[i]);
                }
            }
        }
    }
    return common;
}

/*---------------------------------------------------------------------------*/
std::vector<HDMeshGenerator::HexCorner>
HDMeshGenerator::findFreeCorners(const std::vector<int>& AOrigin,
                                 const int AI,
                                 const int AJ)
{
    std::vector<HexCorner>   out;
    
    for(auto& org: AOrigin){
        std::vector<HexCorner> corners_i;
        findFreeCorners(m_hc_mapping[org], AI, corners_i);
        for(auto ci:corners_i) {
            
            bool found = false;
            
            for(int i=0; i<3; i++) {
                if(ci.adj[i]==AJ)
                    found =true;
            }
            if(found)
                out.push_back(ci);
        }
    }
    return out;
}

/*---------------------------------------------------------------------------*/
bool HDMeshGenerator::getCorner(const int AOrigin,
                                const int AI,
                                const int AJ,
                                const int AK,
                                HexCorner& AOut)
{
    for(auto& c: m_hc_mapping[AOrigin]){
        if((c.adj[0]==AI && c.adj[1]==AJ && c.adj[2]==AK) ||
           (c.adj[0]==AI && c.adj[2]==AJ && c.adj[1]==AK) ||
           (c.adj[1]==AI && c.adj[0]==AJ && c.adj[2]==AK) ||
           (c.adj[1]==AI && c.adj[2]==AJ && c.adj[0]==AK) ||
           (c.adj[2]==AI && c.adj[0]==AJ && c.adj[1]==AK) ||
           (c.adj[2]==AI && c.adj[1]==AJ && c.adj[0]==AK) ){
            AOut = c;
            return true;
        }
    }
    return false;
}

/*---------------------------------------------------------------------------*/
bool HDMeshGenerator::isCorner(const HexCorner& AC,
                               const int AOrigin,
                               const int AI,
                               const int AJ,
                               const int AK )
{
    if(AC.p!=AOrigin)
        return false;
    
    if((AC.adj[0]==AI && AC.adj[1]==AJ && AC.adj[2]==AK) ||
       (AC.adj[0]==AI && AC.adj[2]==AJ && AC.adj[1]==AK) ||
       (AC.adj[1]==AI && AC.adj[0]==AJ && AC.adj[2]==AK) ||
       (AC.adj[1]==AI && AC.adj[2]==AJ && AC.adj[0]==AK) ||
       (AC.adj[2]==AI && AC.adj[0]==AJ && AC.adj[1]==AK) ||
       (AC.adj[2]==AI && AC.adj[1]==AJ && AC.adj[0]==AK) ){
        return true;
    }
    return false;
}
/*---------------------------------------------------------------------------*/
bool HDMeshGenerator::findCommmonLastCorner(const HexCorner& ACorner1,
                                            const HexCorner& ACorner2,
                                            const HexCorner& ACorner3,
                                            HexCorner&       ACornerOut)
{
    //We look for the point shared bu ACorner1, ACorner2 and ACorner3
    int common_pnt_12[2]={-1,-1};
    int i_12=0;
    for(int i1=0; i1<3; i1++){
        int adj1= ACorner1.adj[i1];
        for(int i2=0; i2<3; i2++){
            if(adj1 == ACorner2.adj[i2]){
                common_pnt_12[i_12++]=adj1;
            }
        }
    }//for(int i1=0; i1<3; i1++){
    int common_pnt=-1;
    for(int i3=0; i3<3; i3++){
        int adj3= ACorner3.adj[i3];
        for(int i12=0; i12<3; i12++){
            if(adj3 == common_pnt_12[i12]){
                common_pnt =adj3;
            }
        }
    }//for(int i1=0; i1<3; i1++){
    
    if(common_pnt==-1){
        return false;
    }
    
    // Among all corner, we need the one pointing to ACorner1, ACorner2 and
    // ACorner3
    std::vector<HexCorner> candidates = m_hc_mapping[common_pnt];
    for(unsigned int i = 0; i<candidates.size(); i++){
        HexCorner ci= candidates[i];
        if(!ci.free)
            continue;
        int nb_same=0;
        for(int j=0; j<3; j++){
            if(ci.adj[j]==ACorner1.p ||
               ci.adj[j]==ACorner2.p ||
               ci.adj[j]==ACorner3.p )
                nb_same++;
            
        }
        if(nb_same==3){
            ACornerOut = ci;
            return true;
        }
    }
    return false;
}
/*---------------------------------------------------------------------------*/
void HDMeshGenerator::buildHexahedral()
{
    Variable<int>* var_id = m_hexdom.newVariable<int>(GMDS_NODE,"HD_ID");
    Variable<int>* var_type = m_hexdom.newVariable<int>(GMDS_REGION,"TYPE");
    Variable<int>* var_cl = m_hexdom.newVariable<int>(GMDS_NODE,"classification");
    Variable<int>* var_surf = m_hexdom.newVariable<int>(GMDS_NODE,"surface_id");
    Variable<int>* var_curv = m_hexdom.newVariable<int>(GMDS_NODE,"curve_id");
    
    /* We go trought all the points and we build hexahedral elements
     * for each associated free corner */
    for(unsigned int i=0; i<m_pnt.size(); i++){
        
        math::Point pi = m_pnt[i];
        std::vector<HexCorner> corners_i = m_hc_mapping[i];
        //Now we go throught the associated corners
        
        for(unsigned int j=0; j<corners_i.size();j++){
            HexCorner cj = corners_i[j];
            
            if(!cj.free) //this corner is already used for a hex
                continue;
            //=====================================================
            // STEP 1 - We look for 3 adjacent compatible corners
            // to form a hexahedral element
            //=====================================================
            std::vector<int> common_01 =findCommonPoints(cj.p,cj.adj[0],cj.adj[1]);
            std::vector<int> common_02 =findCommonPoints(cj.p,cj.adj[0],cj.adj[2]);
            std::vector<int> common_12 =findCommonPoints(cj.p,cj.adj[1],cj.adj[2]);
            
            if(common_01.empty() || common_02.empty() || common_12.empty())
                continue;
            
            std::vector<HexCorner> corners_01 = findFreeCorners(common_01, cj.adj[0],cj.adj[1]);
            std::vector<HexCorner> corners_02 = findFreeCorners(common_02, cj.adj[0],cj.adj[2]);
            std::vector<HexCorner> corners_12 = findFreeCorners(common_12, cj.adj[1],cj.adj[2]);
            
            if(corners_01.empty()|| corners_02.empty()|| corners_12.empty())
                continue;
            
            bool found_last_corner = false;
            
            HexCorner final_01, final_02, final_12, final_012;
            for(auto c_01:corners_01){
                for(auto c_02:corners_02){
                    for(auto c_12:corners_12){
                        if(c_01.p==c_02.p || c_01.p==c_12.p ||c_02.p==c_12.p)
                            continue;
                        
                        HexCorner c_012;
                        if(findCommmonLastCorner(c_01,c_02,c_12,c_012)){
                            found_last_corner=true;
                            final_01=c_01;
                            final_02=c_02;
                            final_12=c_12;
                            final_012=c_012;
                        }
                    }
                }
            }
            if(!found_last_corner)
                continue;

            HexCorner final_0, final_1, final_2;

            if(!getCorner(cj.adj[0], cj.p,final_01.p,final_02.p, final_0) ||
               !getCorner(cj.adj[1], cj.p,final_01.p,final_12.p, final_1) ||
               !getCorner(cj.adj[2], cj.p,final_02.p,final_12.p, final_2))
                continue;
          
            //Now we have the 8 corners, we can build the hexahedron
            int idx[8] = {i,
                final_0.p,
                final_01.p,
                final_1.p,
                final_2.p,
                final_02.p,
                final_012.p,
                final_12.p
            };
            
            Node n[8];
            for(int i_n=0; i_n<8; i_n++){
                int index=idx[i_n];
                std::map<int, Node>::iterator it=m_node_mapping.find(index);
                if(it==m_node_mapping.end()){
                    //New node to add
                    n[i_n]=m_hexdom.newNode(m_pnt[index]);
                    (*var_id)[n[i_n].getID()] = index;
                    (*var_cl)[n[i_n].getID()] = m_classification[index];
                    (*var_surf)[n[i_n].getID()] = m_surface[index];
                    (*var_curv)[n[i_n].getID()] = m_curve[index];
                    m_node_mapping[index]=n[i_n];
                    m_used[index] = true;
                }
                else {
                    //existing node
                    n[i_n] = it->second;
                }
            }
            
            Region r = m_hexdom.newHex(n[0], n[1], n[2], n[3],
                                       n[4], n[5], n[6], n[7]);
            (*var_type)[r.getID()]=GMDS_HEX;
            //Eventually, used corners are no more free!
            m_hc_mapping[cj.p][cj.index].free=false;
            m_hc_mapping[final_0.p][final_0.index].free=false;
            m_hc_mapping[final_1.p][final_1.index].free=false;
            m_hc_mapping[final_2.p][final_2.index].free=false;
            m_hc_mapping[final_12.p][final_12.index].free=false;
            m_hc_mapping[final_01.p][final_01.index].free=false;
            m_hc_mapping[final_02.p][final_02.index].free=false;
            m_hc_mapping[final_012.p][final_012.index].free=false;
            
        }//for(unsigned j=0; j<corners_i.size();j++)
        
    }//for(unsigned int i=0; i<m_pnt.size(); i++)
}

/*---------------------------------------------------------------------------*/
void HDMeshGenerator::
buildPrism3(const std::vector<std::vector<OrientedEdge> >& AEdges)
{
    Variable<int>* var_id = m_hexdom.getVariable<int>(GMDS_NODE,"HD_ID");
    Variable<int>* var_type = m_hexdom.getVariable<int>(GMDS_REGION,"TYPE");

    /* We go trought all the points and we build prism elements
     * from free corners */
    for(unsigned int i=0; i<m_pnt.size(); i++){
        
        math::Point pi = m_pnt[i];
        std::vector<HexCorner> corners_i = m_hc_mapping[i];
        
        //Now we go throught the associated corners
        
        for(unsigned int j=0; j<corners_i.size();j++){
            HexCorner cj = corners_i[j];
            
            if(!cj.free) //this corner is already used
                continue;
            //=====================================================
            // STEP 1 - We look for an edge between two adj pnt
            //=====================================================
            std::vector<OrientedEdge> e[3] = {
                AEdges[cj.adj[0]],
                AEdges[cj.adj[1]],
                AEdges[cj.adj[2]]
            };
            
            bool connect[3][3] = {
                {false,false,false},
                {false,false,false},
                {false,false,false}
            };
            for(int i_e=0; i_e<3; i_e++){
                for(auto current_edge:e[i_e]){
                    for(int i_pnt=0; i_pnt<3; i_pnt++){
                        if(i_e!=i_pnt && current_edge.second==cj.adj[i_pnt]){
                            connect[i_e][i_pnt]=1;
                        }
                    }
                }
            }
            // some redundancies in the previous loop since edges
            // are going in both direction
            
            //if we have one edge, is ok, if more, impossible to
            // build a 3-side prism
            int cpt=0;
            int index[3] = {-1,-1,-1};
            for(int cpt_i=0; cpt_i<3; cpt_i++){
                for(int cpt_j=0; cpt_j<cpt_i; cpt_j++){
                    if(connect[cpt_i][cpt_j]){
                        cpt++;
                        index[0]=cpt_i;
                        index[1]=cpt_j;
                    }
                }
            }
            
            if(cpt!=1)
                continue;
            
            //So we have a corner with exactly two end points connected
            //by an edge
            if((index[0]==0 && index[1]==1) ||
               (index[1]==0 && index[0]==1))
                index[2]=2;
            else if((index[0]==0 && index[1]==2) ||
                    (index[1]==0 && index[0]==2))
                index[2]=1;
            else if((index[0]==2 && index[1]==1) ||
                    (index[1]==2 && index[0]==1))
                index[2]=0;
            
            
            //=====================================================
            // STEP 1 - We look for 3 adjacent compatible corners
            // to form a hexahedral element
            //=====================================================
            std::vector<int> common_02 =findCommonPoints(cj.p,
                                                         cj.adj[index[0]],
                                                         cj.adj[index[2]]);
            std::vector<int> common_12 =findCommonPoints(cj.p,
                                                         cj.adj[index[1]],
                                                         cj.adj[index[2]]);
            
            if(common_02.empty() || common_12.empty())
                continue;
            
            std::vector<HexCorner> corners_02 = findFreeCorners(common_02,
                                                                cj.adj[index[0]],
                                                                cj.adj[index[2]]);
            
            std::vector<HexCorner> corners_12 = findFreeCorners(common_12,
                                                                cj.adj[index[1]],
                                                                cj.adj[index[2]]);

            
            if(corners_02.empty()|| corners_12.empty())
                continue;
            
            
            HexCorner final_02, final_12;
            bool found_final_pnts=false;
            
            for(auto c_02:corners_02){
            
                for(auto c_12:corners_12){
                
                    if(isCorner(c_02,c_02.p, c_12.p,
                                cj.adj[index[0]],
                                cj.adj[index[2]]) &&
                       isCorner(c_12,c_12.p, c_02.p,
                                cj.adj[index[1]],
                                cj.adj[index[2]]) ){
                           found_final_pnts=true;
                           final_02=c_02;
                           final_12=c_12;
                    }
                }
            }
            if(!found_final_pnts)
                continue;
            
            HexCorner final_0, final_1, final_2;
            if(!getCorner(cj.adj[index[0]], cj.p, cj.adj[index[1]], final_02.p, final_0) ||
               !getCorner(cj.adj[index[1]], cj.p, cj.adj[index[0]], final_12.p, final_1) ||
               !getCorner(cj.adj[index[2]], cj.p, final_02.p,final_12.p, final_2))
                continue;
            
            //Now we have the 6 corners, we can build the prism
            int idx[6] = {i,
                final_0.p,
                final_1.p,
                final_2.p,
                final_02.p,
                final_12.p
            };
            
            Node n[6];
            for(int i_n=0; i_n<6; i_n++){
                int index=idx[i_n];
                std::map<int, Node>::iterator it=m_node_mapping.find(index);
                if(it==m_node_mapping.end()){
                    //New node to add
                    n[i_n]=m_hexdom.newNode(m_pnt[index]);
                    (*var_id)[n[i_n].getID()] = index;
                    m_node_mapping[index]=n[i_n];

                }
                else {
                    //existing node
                    n[i_n] = it->second;
                }
            }
            
            Region r = m_hexdom.newPrism3(n[0], n[1], n[2], n[3],
                                          n[4], n[5]);
            (*var_type)[r.getID()]=GMDS_PRISM3;

            //Eventually, used corners are no more free!
            m_hc_mapping[cj.p][cj.index].free=false;
            m_hc_mapping[final_0.p][final_0.index].free=false;
            m_hc_mapping[final_1.p][final_1.index].free=false;
            m_hc_mapping[final_2.p][final_2.index].free=false;
            m_hc_mapping[final_12.p][final_12.index].free=false;
            m_hc_mapping[final_02.p][final_02.index].free=false;
            
        }//for(unsigned j=0; j<corners_i.size();j++)
        
    }//for(unsigned int i=0; i<m_pnt.size(); i++)
}
/*---------------------------------------------------------------------------*/
bool HDMeshGenerator::computeVolumePointFrom(const int APntIndex,
                                             const math::Vector3d& AV,
                                             math::Point& APnt)
{
//    std::cout<<"Volume point: "<<APntIndex<<" ("<<m_mesh_data[APntIndex].dim
//    <<", "<<m_classification[APntIndex]<<")"<<std::endl;
    math::Point p = m_pnt[APntIndex];
    double max_dist = m_spacing;
    Region t;
    if(m_param_hexdom.with_edge_interpolation){
        
        if(m_mesh_data[APntIndex].dim==3){
            t = m_mesh->get<Region>(m_mesh_data[APntIndex].id);
        }
        else if(m_mesh_data[APntIndex].dim==2){
            //inner face, we start from the adjacent tet in which
            // AV is going through
            Face f = m_mesh->get<Face>(m_mesh_data[APntIndex].id);
            std::vector<TCellID> fn_ids = f.getIDs<Node>();
            std::vector<Region> fr = f.get<Region>();
            if(fr.size()==1){
                t = fr[0];
            }
            else{
                gmds::Node n[2];
                for(int i=0; i<2; i++){
                    std::vector<Node> reg_nodes =fr[i].get<Node>();
                    for(auto ni:reg_nodes){
                        if(ni.getID()!=fn_ids[0] &&
                           ni.getID()!=fn_ids[1] &&
                           ni.getID()!=fn_ids[2] ){
                            n[i]=ni;
                        }
                    }
                }
                if(AV.dot(math::Vector3d(p,n[0].getPoint()))>0){
                    t=fr[0];
                }
                else{
                    t=fr[1];
                }
            }
        }
        else if(m_mesh_data[APntIndex].dim==1){
            //inner edge, we start from the first adjacent tet
            // by slightly moving p into it
            Edge e = m_mesh->get<Edge>(m_mesh_data[APntIndex].id);
            std::vector<Region> fe = e.get<Region>();
            t = fe[0];
        }
        else{
            //inner point, we start from the first adjacent tet
            // by slightly moving p into it
            Node n = m_mesh->get<Node>(m_mesh_data[APntIndex].id);
            t = n.get<Region>()[0];
        }
        
        // Now starting from p in t, we let p be transported by the
        // frame field in a maximum distance of max_dist
        
        Tools toolkit(m_mesh,0, m_rot_field);
        Tools::PointVolumetricData pvd(p,AV,t);
        
        return (toolkit.followFlow(pvd,max_dist,APnt));
    }
    else{
        APnt.X() = p.X()+AV.X();
        APnt.Y() = p.Y()+AV.Y();
        APnt.Z() = p.Z()+AV.Z();
        return true;
    }
    
}
/*---------------------------------------------------------------------------*/
bool HDMeshGenerator::computeSurfacePointFrom(const int APntIndex,
                                              const math::Vector3d& AV,
                                              math::Point& APnt)
{
    math::Point p = m_pnt[APntIndex];
    
    if(m_param_hexdom.with_edge_interpolation){
        
        double max_dist = m_spacing;
        Face t;
        if(m_mesh_data[APntIndex].dim==2){
            t = m_mesh->get<Face>(m_mesh_data[APntIndex].id);
        }
        else if(m_mesh_data[APntIndex].dim==1){
            //inner edge, we start from the first adjacent boundary
            //triangle
            Edge e = m_mesh->get<Edge>(m_mesh_data[APntIndex].id);
            std::vector<Face> fe = e.get<Face>();
            bool found_bnd_face=false;
            for(auto i=0; i<fe.size(); i++){
                Face fi = fe[i];
                if(fi.getIDs<Region>().size()==1){
                    found_bnd_face=true;
                    t=fi;
                }
            }
            if(!found_bnd_face)
                throw GMDSException("Error no bnd face in HDMeshGenerator::computeSurfacePointFrom for an edge");
        }
        else{//inner point, we start from the first adjacent tet
            // by slightly moving p into it
            Node n = m_mesh->get<Node>(m_mesh_data[APntIndex].id);
            std::vector<Face> fn = n.get<Face>();
            bool found_bnd_face=false;
            for(auto i=0; i<fn.size(); i++){
                Face fi = fn[i];
                if(fi.getIDs<Region>().size()==1){
                    found_bnd_face=true;
                    t=fi;
                }
            }
            if(!found_bnd_face)
                throw GMDSException("Error no bnd face in HDMeshGenerator::computeSurfacePointFrom for an edge");
        }
        
        
        // Now starting from p in t, we let p be transported by the
        // frame field in a maximum distance of max_dist on the surface
        Tools toolkit(m_mesh,0, m_rot_field);
        Tools::PointSurfacicData pvd(p,AV,t);
        
        return (toolkit.followFlow(pvd,max_dist,m_bm.mark_edge_on_curv,APnt));
    }
    else{
        APnt.X() = p.X()+AV.X();
        APnt.Y() = p.Y()+AV.Y();
        APnt.Z() = p.Z()+AV.Z();
        return true;
    }
}

/*---------------------------------------------------------------------------*/
void HDMeshGenerator::updateTetMesh()
{
    
    Variable<int>* v_cavity = m_hexdom.getVariable<int>(GMDS_FACE, "Cavity");
    Variable<int>* v_class  = m_hexdom.getVariable<int>(GMDS_NODE,"classification");

    std::map<int, std::set<Node> > cavity_n;
    for(auto itf = m_hexdom.faces_begin(); !itf.isDone(); itf.next()){
        Face f = itf.value();
        if((*v_cavity)[f.getID()]!=0){
            int cav_id = (*v_cavity)[f.getID()];
            std::vector<Node> fn=f.get<Node>();
            for(auto ni:fn)
                cavity_n[cav_id].insert(ni);
        }
    }
    
    for(auto cav_nodes:cavity_n){
        std::vector<math::Point> cav_pnts;
        std::vector<int> cav_classification;
        cav_pnts.resize(cav_nodes.second.size());
        cav_classification.resize(cav_nodes.second.size());
        int k=0;
        for(auto ni:cav_nodes.second){
            cav_pnts[k]=ni.getPoint();
            cav_classification[k] = (*v_class)[ni.getID()];
            k++;
        }
        TetMeshManipulator tmm(m_mesh);
//        IGMesh temp(MeshModel(DIM3|R|N|R2N));
        tmm.insertPoints(m_mesh,cav_pnts,cav_classification);
    }
   
}
/*---------------------------------------------------------------------------*/
void HDMeshGenerator::writeInput()
{
    Log::mng() << "Write point and chart ...";
    MeshModel model(DIM3 | R  | N | R2N );
    IGMesh mesh_pnts  (model);
    IGMesh mesh_charts(model);
    Variable<int>* var_pnts   = mesh_pnts.newVariable<int>(GMDS_REGION, "type");
    Variable<int>* var_charts_type = mesh_charts.newVariable<int>(GMDS_REGION, "type");
    Variable<int>* var_charts_idx = mesh_charts.newVariable<int>(GMDS_REGION, "index");
    Variable<int>* var_charts_cl = mesh_charts.newVariable<int>(GMDS_REGION, "classification");
    Variable<int>* var_curv = mesh_charts.newVariable<int>(GMDS_REGION, "curv_number");
    Variable<int>* var_surf = mesh_charts.newVariable<int>(GMDS_REGION, "surf_number");
    Variable<int>* var_used = mesh_charts.newVariable<int>(GMDS_REGION, "used");

    
    double cube_size = m_spacing/5.0;
    
    for(unsigned int i=0; i<m_pnt.size();i++){
        math::Point pi = m_pnt[i];
        math::Chart ci = m_chart[i];
        int ti = m_type[i];
        
        //POINT OUTPUT
        Node   ni  = mesh_pnts.newNode(pi);
        Region ri = mesh_pnts.newTet(ni,ni,ni,ni);
        (*var_pnts)[ri.getID()]= ti;
        
        //CHART OUTPUT
        math::Vector vx = ci.X();
        math::Vector vy = ci.Y();
        math::Vector vz = ci.Z();
        math::Point center = pi;
        math::Point p1 = center + (vx + vy - vz)*cube_size;
        Node n1 = mesh_charts.newNode(p1);
        math::Point p2 = center + (vx - vy - vz)*cube_size;
        Node n2 = mesh_charts.newNode(p2);
        math::Point p3 = center + (vx + vy + vz).opp()*cube_size;
        Node n3 = mesh_charts.newNode(p3);
        math::Point p4 = center + (vy - vx - vz)*cube_size;
        Node n4 = mesh_charts.newNode(p4);
        
        math::Point p5 = center + (vx + vy + vz)*cube_size;
        Node n5 = mesh_charts.newNode(p5);
        math::Point p6 = center + (vx - vy + vz)*cube_size;
        Node n6 = mesh_charts.newNode(p6);
        math::Point p7 = center + (vx + vy - vz).opp()*cube_size;
        Node n7 = mesh_charts.newNode(p7);
        math::Point p8 = center + (vy - vx + vz)*cube_size;
        Node n8 = mesh_charts.newNode(p8);
        Region r = mesh_charts.newHex(n1, n2, n3, n4, n5, n6, n7, n8);
        (*var_charts_type)[r.getID()]= ti;
        (*var_charts_idx)[r.getID()] = i;
        (*var_charts_cl)[r.getID()] = m_classification[i];
        (*var_curv)[r.getID()] = m_curve[i];
        (*var_surf)[r.getID()] = m_surface[i];
        (*var_used)[r.getID()] = m_used[i];
    }
    static int nb_file=0;
    std::string file_pnts, file_charts;
    file_pnts=m_param_gl.output_dir+"/HDM_INPUT_PNTS_"+to_string(nb_file);
    file_charts=m_param_gl.output_dir+"/HDM_INPUT_CHARTS_"+to_string(nb_file);
    nb_file++;
    VTKWriter<IGMesh> writer_pnts(mesh_pnts);
    writer_pnts.write(file_pnts, DIM3 | R | N);
    VTKWriter<IGMesh> writer_charts(mesh_charts);
    writer_charts.write(file_charts, DIM3 | R | N);
    
    Log::mng()<<" done!\n";
    Log::mng().flush();
}

/*---------------------------------------------------------------------------*/
void HDMeshGenerator::
writeEdges(std::vector<std::vector<OrientedEdge> >& AEdges,
           const std::string& AFileName)
{
    Log::mng()<< "Write  edges ...\n";
    MeshModel model(DIM3 | F  | N | F2N );
    IGMesh mesh_edges  (model);
    
    for(unsigned int i=0; i<m_pnt.size();i++){

        std::vector<OrientedEdge> edges_i = AEdges[i];
        for(unsigned int j=0; j<edges_i.size(); j++){
            math::Point p0 = m_pnt[edges_i[j].first];
            math::Point p1 = m_pnt[edges_i[j].second];
            Node   n0 = mesh_edges.newNode(p0);
            Node   n1 = mesh_edges.newNode(p1);
            mesh_edges.newTriangle(n0,n1,n1);
            
        }
    }
    VTKWriter<IGMesh> writer(mesh_edges);
    writer.write(AFileName, DIM3 | F | N);
    Log::mng()<<" done!\n";
    Log::mng().flush();
}


/*---------------------------------------------------------------------------*/
void HDMeshGenerator::writeOutput()
{
    static int nb_file=0;
    Log::mng()<< "Write hex-dom mesh (";
    Log::mng()<< m_hexdom.getNbRegions()<<" regions) ...";
    VTKWriter<IGMesh> writer(m_hexdom);
    std::string file_name=m_param_gl.output_dir+"/HDM_MESH_"+to_string(nb_file++);
    writer.write(file_name, DIM3 | R| F| N);
    Log::mng()<<" done!\n";
    
    Log::mng()<<"Initial spacing: "<<m_spacing<<"\n";
    Log::mng().flush();
    
    double l=0;
    for(IGMesh::region_iterator itr = m_hexdom.regions_begin();
        !itr.isDone(); itr.next()){
        Region r = itr.value();
        if(r.getType()!=GMDS_HEX)
            continue;
        std::vector<Node> n = r.get<Node>();
        l+= n[0].getPoint().distance(n[1].getPoint());
        l+= n[1].getPoint().distance(n[2].getPoint());
        l+= n[2].getPoint().distance(n[3].getPoint());
        l+= n[3].getPoint().distance(n[0].getPoint());
        
        l+= n[4].getPoint().distance(n[5].getPoint());
        l+= n[5].getPoint().distance(n[6].getPoint());
        l+= n[6].getPoint().distance(n[7].getPoint());
        l+= n[7].getPoint().distance(n[4].getPoint());
        
        l+= n[0].getPoint().distance(n[4].getPoint());
        l+= n[1].getPoint().distance(n[5].getPoint());
        l+= n[2].getPoint().distance(n[6].getPoint());
        l+= n[3].getPoint().distance(n[7].getPoint());

    }
    Log::mng()<<"Final average length: "<<l/(8*m_hexdom.getNbRegions())<<"\n";
    Log::mng().flush();
}

/*---------------------------------------------------------------------------*/
void HDMeshGenerator::writeLoops(std::vector<std::vector<Node> >& ALoops)
{
    Log::mng()<< "Write loops (";
    Log::mng()<< ALoops.size()<<") ...";
    
    IGMesh mesh_loops(MeshModel(DIM3 | F | N | F2N));
    Variable<int>* v = mesh_loops.newVariable<int>(GMDS_FACE, "LOOP");

    for(auto i=0; i<ALoops.size(); i++){
        std::vector<Node> loop= ALoops[i];
        for(auto j=0; j<loop.size();j++){
            Node n1 = mesh_loops.newNode(loop[j].getPoint());
            Node n2 = mesh_loops.newNode(loop[(j+1)%loop.size()].getPoint());
            Face f = mesh_loops.newTriangle(n1,n1,n2);
            (*v)[f.getID()]=i;
        }
    }
    
    static int nb_file=0;
    std::string file_name=m_param_gl.output_dir+"/LOOPS_CAVITY_"+to_string(nb_file++);
    VTKWriter<IGMesh> writer(mesh_loops);
    writer.write(file_name, DIM3 | F| N);
    Log::mng()<<" done!\n";
    Log::mng().flush();
}
/*----------------------------------------------------------------------------*/
math::Vector3d HDMeshGenerator::
getOutputNormal(Face& AFace, Region& ARegion)
{
    std::vector<Node> region_nodes = ARegion.get<Node>();
    std::vector<Node> face_nodes = AFace.get<Node>();
    
    math::Point  reg_center  = ARegion.center();
    math::Point  face_center = AFace.center();
    math::Vector face_normal = AFace.normal();
    math::Vector v_ref(face_center, reg_center);
    
    if(face_normal.dot(v_ref)>0)
        return math::Vector3d(-face_normal.X(), -face_normal.Y(), -face_normal.Z());
    
    return math::Vector3d(face_normal.X(), face_normal.Y(), face_normal.Z());

}
/*----------------------------------------------------------------------------*/
math::Vector3d HDMeshGenerator::
getOutputNormalOfABoundaryFace(Face& AFace)
{
    std::vector<Region> adj_regions = AFace.get<Region>();
    if (adj_regions.size() != 1)
        throw GMDSException("A boundary face must be adjacent to only 1 region!!!");
    
    return getOutputNormal(AFace, adj_regions[0]);
}
/*---------------------------------------------------------------------------*/


