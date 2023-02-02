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
#include <GMDS/IO/VTKWriter.h>
#include <GMDS/Utils/Log.h>
#include "GMDS/Math/Constants.h"
#include "GMDS/Math/Numerics.h"
/*---------------------------------------------------------------------------*/
//GMDS File Header
#include <GMDS/Math/Numerics.h>
#include <GMDS/Math/Plane.h>
#include <GMDS/Utils/Log.h>
/*---------------------------------------------------------------------------*/
// FHeDo File Headers
#include <Tools.h>
#include "SingularityRegionMesher.h"
#include <algorithm>
/*---------------------------------------------------------------------------*/
using namespace gmds;
using namespace fhedo;
/*---------------------------------------------------------------------------*/
SingularityRegionMesher::
SingularityRegionMesher(IGMesh* AMesh,
                        ParamsGlobal& AGlobal,
                        double AEdgeSize,
                        ParamsMark& AMarks):
m_mesh(AMesh),
m_params_gl(AGlobal),
m_edge_size(AEdgeSize),
m_marks(AMarks),
m_graph(AMesh)
{
    m_rot = m_mesh->getVariable<math::AxisAngleRotation>(GMDS_NODE,
                                                         "rotation_field");
    
    MeshModel model(DIM3|F|F2N);
    m_line_mesh = new IGMesh(model);
}
/*---------------------------------------------------------------------------*/
SingularityRegionMesher::~SingularityRegionMesher()
{
    if(m_line_mesh)
        delete m_line_mesh;
}
/*---------------------------------------------------------------------------*/
void SingularityRegionMesher::execute()
{
    initializeBooleanMarks();
    
    //=========================================================================
    // DETECTION OF SINGULAR ELEMENTS
    //=========================================================================
    detectSingularCells();
    
    //========================================================================
    // STEP 1 - Singularity points inside the volume
    //========================================================================
    Log::mng()<<"=========================================================\n";
    Log::mng()<<"Creation of volume singularity points\n";
    Log::mng().flush();
    createVolumeSingularityPoints();
    Log::mng()<< "\t DONE\n";
    Log::mng().flush();
    
    //========================================================================
    // STEP 2 - Singularity points on the surface
    //========================================================================
    Log::mng()<< "========================================================\n";
    Log::mng()<< "Creation of boundary singularity points\n";
    createBoundarySingularityPoints();
    Log::mng()<< "\t DONE\n";
    Log::mng().flush();
    
    //========================================================================
    // STEP 3 - Singularity lines inside the volume
    //========================================================================
    //at this point we have made all cluster of sing. points
    //now we need to compute the lines of sing into separatrices
    Log::mng()<< "========================================================\n";
    Log::mng()<< "Creation of volume singularity lines\n";
    createVolumeSingularityLines();
    Log::mng()<< "\t DONE\n";
    Log::mng().flush();
    
    //========================================================================
    // STEP 4 - Build surface quad patch around surface sing point
    //========================================================================
    Log::mng()<< "========================================================\n";
    Log::mng()<< "Creation of surface patch\n";
    createSurfaceQuadPatches();
    createLineTubes();
    Log::mng()<< "\t DONE\n";
    Log::mng().flush();
    writeOutputSingle(m_params_gl.output_dir+"/final_lines");
    finalizeBooleanMarks();
}/*---------------------------------------------------------------------------*/
void SingularityRegionMesher::initializeBooleanMarks()
{
    m_mark_line_sing = m_mesh->getNewMark<Region>();
    m_mark_volume_pnt_sing = m_mesh->getNewMark<Region>();
    m_mark_face_sing = m_mesh->getNewMark<Face>();
    m_mark_edge_sing = m_mesh->getNewMark<Edge>();
    m_mark_node_sing = m_mesh->getNewMark<Node>();
}
/*---------------------------------------------------------------------------*/
void SingularityRegionMesher::finalizeBooleanMarks()
{
    m_mesh->unmarkAll<Region>(m_mark_volume_pnt_sing);
    m_mesh->freeMark<Region>(m_mark_volume_pnt_sing);
    
    m_mesh->unmarkAll<Region>(m_mark_line_sing);
    m_mesh->freeMark<Region>(m_mark_line_sing);
    
    m_mesh->unmarkAll<Face>(m_mark_face_sing);
    m_mesh->freeMark<Face>(m_mark_face_sing);
    
    m_mesh->unmarkAll<Edge>(m_mark_edge_sing);
    m_mesh->freeMark<Edge>(m_mark_edge_sing);
    
    m_mesh->unmarkAll<Node>(m_mark_node_sing);
    m_mesh->freeMark<Node>(m_mark_node_sing);
}
/*----------------------------------------------------------------------------*/
void SingularityRegionMesher::detectSingularCells()
{
    Tools tool(m_mesh,0,m_rot);
    IGMesh::region_iterator itr = m_mesh->regions_begin();
    for (; !itr.isDone(); itr.next()){
        Region r = itr.value();
        if(tool.isFFSingular(r)){
            m_sing_tet.push_back(r.getID());
            m_mesh->mark(r, m_mark_line_sing);
            int type = 0;
            std::vector<Face> fs = r.get<Face>();
            for(auto f:fs){
                if (tool.isFFSingular(f)){
                    type++;
                    if(m_mesh->isMarked(f, m_marks.mark_face_on_surf))
                        m_sing_bnd_tri.push_back(f.getID());
                }
            }
            m_sing_tet_type[r.getID()]=type;
        }
        
    }
    Log::mng()<<"Singular tets    : "<<m_sing_tet.size()<<")\n";
    Log::mng()<<"Singular bnd tris: "<<m_sing_bnd_tri.size()<<")\n";
    Log::mng().flush();
}
/*----------------------------------------------------------------------------*/
void SingularityRegionMesher::createBoundarySingularityPoints()
{
    //=========================================================================
    // For each boudary region, we create a singularity point with
    // slots oriented inside the volume and along the surface.
    //=========================================================================
    Tools tool(m_mesh,0,0);
    
    //angles approach gives worst results
    bool with_angles=false;
    std::vector<bool> already_done(m_sing_bnd_tri.size(),false);
    for (auto i=0; i<m_sing_bnd_tri.size();i++) {
        
        //already done as adjacent to another bnd singular triangle
        if(already_done[i])
            continue;
        
        TCellID id = m_sing_bnd_tri[i];
        //f is a singular boundary face
        Face f = m_mesh->get<Face>(id);
        std::vector<TCellID> other_f;
        std::vector<std::pair<TCellID, TCellID> > other_e;
        std::vector<Node> nodes = f.get<Node>();
        //We look for a singular boundary faces sharing an edge with f;
        for(auto i_n=0; i_n<nodes.size(); i_n++){
            Node ni = nodes[i_n];
            Node nj = nodes[(i_n+1)%nodes.size()];
            std::vector<TCellID> faces_i = ni.getIDs<Face>();
            std::vector<TCellID> faces_j = nj.getIDs<Face>();
            std::vector<TCellID> common_ij;
            for(auto fi:faces_i){
                for(auto fj:faces_j){
                    if(fi==fj){
                        common_ij.push_back(fi);
                    }
                }
            }
            for(auto c_ij:common_ij){
                if(m_mesh->isMarked<Face>(c_ij,m_mark_face_sing) &&
                   m_mesh->isMarked<Face>(c_ij,m_marks.mark_face_on_surf)){
                    other_f.push_back(c_ij);
                    other_e.push_back(std::pair<TCellID,TCellID>(ni.getID(),nj.getID()));
                }
            }
        }
        std::vector<math::Vector3d> sepa;

        math::Point sing_pnt(0,0,0);
        math::Vector n(0,0,0);
        if(other_f.empty()){
            //only f has a triangular face
            sing_pnt = f.center();
            math::Vector3d n3d = getOutputNormal(f,f.get<Region>()[0]);
            n =math::Vector(n3d.X(),n3d.Y(),n3d.Z());

        }
        else if(other_f.size()==1){
            Face f2 = m_mesh->get<Face>(other_f[0]);
            TCellID i1 = other_e[0].first;
            TCellID i2 = NullID;
            TCellID i3 = other_e[0].second;
            TCellID i4 = NullID;
            std::vector<TCellID> n1_ids = f.getIDs<Node>();
            std::vector<TCellID> n2_ids = f2.getIDs<Node>();
            for(auto id_1:n1_ids){
                if(id_1!=i1 && id_1!=i3){
                    i2=id_1;
                }
            }
            for(auto id_2:n2_ids){
                if(id_2!=i1 && id_2!=i3){
                    i4=id_2;
                }
            }
            nodes.resize(3);
            nodes[0]= m_mesh->get<Node>(i1);
            nodes[1]= m_mesh->get<Node>(i2);
            nodes[2]= m_mesh->get<Node>(i3);
            //nodes[3]= m_mesh->get<Node>(i4);
            
            sing_pnt = f.center();
            math::Vector3d n3d = getOutputNormal(f,f.get<Region>()[0]);
            n =math::Vector(n3d.X(),n3d.Y(),n3d.Z());
        }
        else{
            throw GMDSException("Impossible to get more than 2 adjacent boundary singular faces");
        }
        
        if(with_angles){
            sepa = defineBoundarySlotsViaAngles(nodes, sing_pnt, n);
        }
        else{
            sepa = defineBoundarySlotsViaVectors(nodes, sing_pnt, n);
        }
        std::cout<<"NB POINTS: "<<sepa.size()<<std::endl;
        SurfaceSingularityPoint* new_sing_pnt = m_graph.newSurfacePoint();
        new_sing_pnt->setLocation(sing_pnt);
        new_sing_pnt->addMeshFace(f);
        for(auto s :sepa){
            math::Vector vs(s.X(),s.Y(),s.Z());
            vs.normalize();
            math::Point p = sing_pnt+m_edge_size*vs;
            new_sing_pnt->newSlot(p, vs,
                                  f.getID(), 2,
                                  true);

        }
        //slot inward
        new_sing_pnt->newSlot(sing_pnt, n.opp(),
                              f.getID(), 2, false);
    }
    
    
}
/*----------------------------------------------------------------------------*/
std::vector<math::Vector3d> SingularityRegionMesher::
defineBoundarySlotsViaAngles(const std::vector<Node>& ANodes,
                             const math::Point&       ASingLoc,
                             const math::Vector&      ANormal)
{
    std::vector<math::Vector3d> sep;
    
    //=====================================================================
    //we got the rotation assigned to each vertex of f
    math::AxisAngleRotation rot[3];
    for(auto i=0; i<ANodes.size(); i++){
        rot[i] = (*m_rot)[ANodes[i].getID()];
    }
    
    //=====================================================================
    //We pick on chart vector living in the plane of f
    math::Chart c[3] = {rot[0].toChart(), rot[1].toChart(), rot[2].toChart()};
    math::Vector3d v[3];
    math::Vector3d n3d(ANormal);
    for(auto i=0; i<3; i++){
        for(auto j=0; j<3; j++){
            if(abs(c[i][j].dot(n3d))<0.1){
                v[i]=c[i][j];
            }
        }
    }

    //=====================================================================
    // We compute the singularity index (+1 or -1)
    //=====================================================================
    // We compute reference angle in the face plane using Palacio technique
    math::Vector ref(ANodes[0].getPoint(),ANodes[1].getPoint());
    
    ref.normalize();
    
    math::Vector vc[3] = {
        math::Vector(v[0].X(), v[0].Y(), v[0].Z()),
        math::Vector(v[1].X(), v[1].Y(), v[1].Z()),
        math::Vector(v[2].X(), v[2].Y(), v[2].Z()),
    };
    double angle[3]={0,0,0};
    for(auto i=0;i<3;i++){
        TCoord a = vc[i].orientedAngle(ref,ANormal);
        if(a<0){
            a = math::Constants::PI2+a;
        }
        angle[i] = math::modulo2PI(4*a);
    }
    
    //=====================================================================
    // Then the reference vector
    math::Vector ref0(cos(angle[0]),sin(angle[0]),0.0);
    math::Vector ref1(cos(angle[1]),sin(angle[1]),0.0);
    math::Vector ref2(cos(angle[2]),sin(angle[2]),0.0);
    
    //=====================================================================
    // And so the face index k = 1 or -1 depending on the sing type
    double w01 =ref0.orientedAngle(ref1);
    double w12 =ref1.orientedAngle(ref2);
    double w20 =ref2.orientedAngle(ref0);
    double index_d = (w01+w12+w20)/math::Constants::PI2;
    int index = round(index_d);
    

    
    
    
    //=====================================================================
    // As we have the index, we can deduce singularity direction as done
    // in H. Fogg PHD mansucrit (see Bunin's papers too)
    
    //we need to find a first direction
    bool found_first_separatrix = false;
    math::Vector3d first_sepa;
    for(auto i=0; i<3 && !found_first_separatrix; i++){
        const auto j=(i+1)%3;
        math::Point pi = ANodes[i].getPoint();
        math::Point pj = ANodes[j].getPoint();
//        std::cout<<"From "<<nodes[i].getID()<<" to "<<nodes[j].getID()<<std::endl;
        //we work along the edge [i,j]
        math::Vector vij(pi,pj);
        
        //all ref angles are recomputed for this edge
        for(auto k=0;k<3;k++){
            TCoord a = vij.orientedAngle(vc[k],ANormal);
            if(a<0){
                a = math::Constants::PI2+a;
            }
            angle[k] = math::modulo2PI(4*a);
        }
        
        //We get the four angles of the cross in point i and j
        //cross in the plane of f
        
        double cross_angle_i[4];
        TCoord a = angle[i]/4.0;
        for(auto k=0; k<4; k++){
            cross_angle_i[k] = math::modulo2PI(a+k*math::Constants::PIDIV2);
        }
        
        double cross_angle_j[4];
        a = angle[j]/4.0;
        for(auto k=0; k<4; k++){
            cross_angle_j[k] =math::modulo2PI(a+k*math::Constants::PIDIV2);
        }
        //For each vector, we check if we find a point where to
        //go out along edge [i,j]
        math::Vector vi(pi,ASingLoc);
        math::Vector vj(pj,ASingLoc);
        double alpha_i = vij.orientedAngle(vi,ANormal);
        double alpha_j = vij.orientedAngle(vj,ANormal);
        
        for(auto k=0; k<4 && !found_first_separatrix;k++){
            double ai = cross_angle_i[k];
            double aj = cross_angle_j[0];
            
            double dij = abs(aj-ai);
            if(dij>math::Constants::PI)
                dij=math::Constants::PI2-dij;
            int match_j=0;
            for(auto l=1;l<4;l++){
                double dij_l = abs(ai-cross_angle_j[l]);
                if(dij_l>math::Constants::PI)
                    dij_l=math::Constants::PI2-dij_l;
                
                
                if(dij_l<dij){
                    dij=dij_l;
                    aj = cross_angle_j[l];
                    match_j=l;
                }
            }
            double delta = aj-ai;
            double t= (alpha_i-ai)/(delta-(alpha_j-alpha_i));
            if(t>=0 && t<=1){
                found_first_separatrix=true;
                math::Point pnt_sepa = (1-t)*pi + t*pj;
                first_sepa = math::Vector3d(ASingLoc,pnt_sepa);
                first_sepa.normalize();
            }
        }
        
    }//for(auto i=0; i<3 && !found_first_separatrix; i++)
    
    if(!found_first_separatrix){
        throw GMDSException("Impossible to compute a separatrix");
    }
    
    double separatrix_angle = math::Constants::PI2/(4+index);
    
    int nb_sing = 4+index;
    sep.resize(nb_sing);
    sep[0]=first_sepa;
    math::AxisAngleRotation rot_sepa(ANormal, separatrix_angle);
    for(auto s= 1; s<nb_sing;s++){
        sep[s] = rot_sepa*sep[s-1];
        
    }

    return sep;
}
/*----------------------------------------------------------------------------*/
std::vector<math::Vector3d> SingularityRegionMesher::
defineBoundarySlotsViaVectors(const std::vector<Node>& ANodes,
                              const math::Point&       ASingLoc,
                              const math::Vector&      ANormal)
{
    std::vector<math::Vector3d> sep;
    
    double nx = ANormal.X();
    double ny = ANormal.Y();
    double nz = ANormal.Z();
    double cx = ASingLoc.X();
    double cy = ASingLoc.Y();
    double cz = ASingLoc.Z();

    
    //=====================================================================
    //we got the rotation assigned to each vertex of f
    std::vector<math::AxisAngleRotation> rot;
    rot.resize(ANodes.size());
    for(auto i=0; i<ANodes.size(); i++){
        rot[i] = (*m_rot)[ANodes[i].getID()];
    }
    
    //=====================================================================
    //We pick on chart vector living closed to the plane of f
    math::Vector3d v[3];
    math::Vector3d n3d(ANormal);
    for(auto i=0; i<ANodes.size(); i++){
        math::Chart c = rot[i].toChart();
        for(auto j=0; j<3; j++){
            if(abs(c[j].dot(n3d))<0.1){
                v[i]=c[j];
            }
        }
    }
    //=====================================================================
    //We project charts into the plane to avoid numerical issues
    math::Plane pl(ASingLoc,ANormal);
    for(auto i=0; i<ANodes.size(); i++){
        math::Vector vi(v[i].X(),v[i].Y(),v[i].Z());
        math::Point p = ASingLoc+vi;
        p = pl.project(p);
        v[i]= math::Vector3d(ASingLoc,p);
        v[i].normalize();
    }
    //=====================================================================
    //We get the four vector defining each cross in each node
    std::vector<std::vector<math::Vector3d> > cross_vec;
    cross_vec.resize(ANodes.size());
    for(auto i=0; i<ANodes.size(); i++){
        cross_vec[i].resize(4);
        cross_vec[i][0] = v[i];
        math::AxisAngleRotation rot_pi2(ANormal, math::Constants::PIDIV2);
        for(auto j= 1; j<4; j++){
            cross_vec[i][j] = rot_pi2*cross_vec[i][j-1];
            
        }

    }
   
    
    for(auto i=0; i<ANodes.size();i++){
        std::vector<double> param_ij;
        int j = (i+1)%ANodes.size();
        Node ni = ANodes[i];
        Node nj = ANodes[j];
        math::Point pi = ni.getPoint();
        math::Point pj = nj.getPoint();
        double pix = pi.X();
        double piy = pi.Y();
        double piz = pi.Z();
        double pjx = pj.X();
        double pjy = pj.Y();
        double pjz = pj.Z();

        std::vector<math::Vector3d> vis = cross_vec[i];
        std::vector<math::Vector3d> vjs = cross_vec[j];
     
        bool found=false;
        //We check with the 4 couple of vectors between i and j
        for(auto k=0; k<4 && !found; k++){
            math::Vector3d vi = vis[k];
            math::Vector3d vj;
            double d=-2;
            for(auto l=0; l<4; l++){
                double dot_il = vi.dot(vjs[l]);
                if(dot_il>d){
                    vj = vjs[l];
                    d=dot_il;
                }
            }
            
            double vix = vi.X();
            double viy = vi.Y();
            double viz = vi.Z();
            double vjx = vj.X();
            double vjy = vj.Y();
            double vjz = vj.Z();
            //2nd order equation in alpha the paramrter
            double a =
            -( (nz*piy - ny*piz - nz*pjy + ny*pjz)*vix
              -(nz*pix - nx*piz - nz*pjx + nx*pjz)*viy
              +(ny*pix - nx*piy - ny*pjx + nx*pjy)*viz
              -(nz*piy - ny*piz - nz*pjy + ny*pjz)*vjx
              +(nz*pix - nx*piz - nz*pjx + nx*pjz)*vjy
              -(ny*pix - nx*piy - ny*pjx + nx*pjy)*vjz);
            double c =
            -(cz*ny - cy*nz + nz*piy - ny*piz)*vix
            +(cz*nx - cx*nz + nz*pix - nx*piz)*viy
            -(cy*nx - cx*ny + ny*pix - nx*piy)*viz;
            
            double b =
            (cz*ny - cy*nz + 2*nz*piy - 2*ny*piz - nz*pjy + ny*pjz)*vix
            -(cz*nx - cx*nz + 2*nz*pix - 2*nx*piz - nz*pjx + nx*pjz)*viy
            +(cy*nx - cx*ny + 2*ny*pix - 2*nx*piy - ny*pjx + nx*pjy)*viz
            -(cz*ny - cy*nz + nz*piy - ny*piz)*vjx
            +(cz*nx - cx*nz + nz*pix - nx*piz)*vjy
            -(cy*nx - cx*ny + ny*pix - nx*piy)*vjz;
            std::vector<double> x;
            int nb_x = math::solve2ndDegreePolynomial(a,b,c,x);
            //we can have 0,1 or 2 solutions
            for(auto sol:x){
                
                if(sol>=0 && sol <1){
                    param_ij.push_back(sol);
                }
            }//for(auto sol:x)
            
        }
        
        for(auto i=0;i<param_ij.size();i++){
            bool close = false;
            for(auto j=i+1;j<param_ij.size()&& !close;j++){
                if(abs(param_ij[i]-param_ij[j])<1e-2)
                    close=true;
            }
            if(!close){
                double s =param_ij[i];
                math::Point p_sol = (1-s)*pi + s*pj;
                math::Vector3d v_sol(ASingLoc, p_sol);
                v_sol.normalize();
                sep.push_back(v_sol);
            }
        }
    }
    
    //vectors of sep are not necessarily well sorted. We do it now
    std::vector<math::Vector3d> final_sep;
    final_sep.resize(sep.size());

    final_sep[0]=sep[0];
    math::Vector ref(sep[0].X(),sep[0].Y(),sep[0].Z());
    std::vector<double> angle(sep.size(),0);
    angle[0]=0;
    for(auto i=1;i<sep.size();i++){
        math::Vector vi(sep[i].X(),sep[i].Y(),sep[i].Z());
        angle[i] = ref.angleIn02PI(vi,ANormal);
        std::cout<<"angle "<<i<<": "<<angle[i]<<std::endl;
    }
    for(auto i=1;i<sep.size();i++){
        int index =1;
        for(auto j=1;j<sep.size();j++){
            if(i==j)
                continue;
            if(angle[i]>angle[j])
                index++;
        }
        final_sep[index]=sep[i];
        std::cout<<index<<": "<<final_sep[index]<<std::endl;
    }
    
    
    return final_sep;
    

}
/*----------------------------------------------------------------------------*/
void SingularityRegionMesher::createVolumeSingularityPoints()
{
    
    auto nb_clusters = 0;
    
    //=========================================================================
    // FIRST LOOP ON REGIONS TO GET ALL THE 3-SING. TETS
    //=========================================================================
    
    for (auto i:m_sing_tet){
        
        Region r = m_mesh->get<Region>(i);
        int s_type = m_sing_tet_type[i];
        
        if (s_type >= 3 && (!m_mesh->isMarked(r, m_mark_volume_pnt_sing)))
        {
            nb_clusters++;
            //we are currently in a region that is part of a cluster and has
            //not been used yet, we are going to use it now
            
            m_mesh->mark(r, m_mark_volume_pnt_sing);
            std::vector<Region> cluster_region;
            
            cluster_region.push_back(r);
            bool stop = false;
            while (!stop) {
                stop=true;
                for (auto ri:cluster_region){
                    std::vector<Face> ri_faces = ri.get<Face>();
                    for (auto f:ri_faces) {
                        
                        std::vector<Region> f_regs = f.get<Region>();
                        //we get the opposite regions or ri
                        Region r_opp =(f_regs[0].getID()==ri.getID())?f_regs[1]:f_regs[0];
                        
                        
                        if (  m_mesh->isMarked(r_opp, m_mark_line_sing) &&
                            (!m_mesh->isMarked(r_opp, m_mark_volume_pnt_sing))) {
                            cluster_region.push_back(r_opp);
                            m_mesh->mark(r_opp, m_mark_volume_pnt_sing);
                            stop=false;
                        }
                        
                        
                    }//for (auto f:ri_faces)
                    
                }//for (auto ri:cluster_region)
                
            }//while (!stop)
            //==================================================================
            // We create the volume singularity point
            //==================================================================
            createOneVolumeSingularityPoint(cluster_region);
        }//if (singTypeTmp == 3 && (!m_mesh->isMarked(Rtmp, m_markClusterSingDone)))
    
    }//for (; !itr.isDone(); itr.next())
    
    //=========================================================================
    Log::mng()<<"We have " << nb_clusters << " clusters of 3-sing. tetrahedra, ";
    Log::mng()<< "and so " << m_graph.getNbVolumePoints() << " sing. graph points\n";
    Log::mng().flush();
}
/*----------------------------------------------------------------------------*/
void SingularityRegionMesher::
createOneVolumeSingularityPoint(std::vector<gmds::Region>& ACluster)
{
    VolumeSingularityPoint* singularity = m_graph.newVolumePoint();
    //==================================================================
    // First, we compute the center of mass of the cluster and we
    // associate the region to the singularity point
    //==================================================================
    auto nb_tets = ACluster.size();
    auto x = 0.0, y = 0.0, z = 0.0;
    for(auto r: ACluster) {
        //association
        singularity->addMeshRegion(r);
        
        //contribution to the point location
        math::Point c = r.center();
        x+=c.X();
        y+=c.Y();
        z+=c.Z();
    }
    //==================================================================
    //location
    math::Point sing_pnt (x/nb_tets, y/nb_tets, z/nb_tets);
    singularity->setLocation(sing_pnt);
    
    //==================================================================
    //slots definition
    
    //We get the outer faces of the cluster
    std::vector<Face> cluster_bnd_faces;
    for(auto r: ACluster) {
        std::vector<Face> local_faces = r.get<Face>();
        for (auto f:local_faces) {
            std::vector<Region> local_regions = f.get<Region>();
            if (local_regions.size() == 1)
                cluster_bnd_faces.push_back(f);
            else {
                Region r0 = local_regions[0];
                Region r1 = local_regions[1];
                if (( m_mesh->isMarked(r0, m_mark_volume_pnt_sing) &&
                     !m_mesh->isMarked(r1, m_mark_volume_pnt_sing)) ||
                    (!m_mesh->isMarked(r0, m_mark_volume_pnt_sing) &&
                     m_mesh->isMarked(r1, m_mark_volume_pnt_sing)))
                    cluster_bnd_faces.push_back(f);
            }
        }
    }
    Tools tool(m_mesh,0,m_rot);
    //and now we look for slots
    for (auto f:cluster_bnd_faces){
        
        m_mesh->mark(f, m_mark_face_sing);
        for (auto e:f.get<Edge>())
            m_mesh->mark(e, m_mark_edge_sing);
        
        
        for (auto n:f.get<Node>())
            m_mesh->mark(n, m_mark_node_sing);
        
        if (tool.isFFSingular(f)){
            //A Singularty line goes through this face
            math::Point face_center = f.center();
            math::Vector line_vec(sing_pnt, face_center);
            singularity->newSlot(face_center, line_vec, f.getID(), 2, false);
        }
    }
    
    Log::mng()<< "Volume sing. point created at "
    << singularity->getLocation()
    << " in a cluster of " << singularity->getNbMeshCells()
    << " with " << singularity->getSlots().size() << " slots \n";
    Log::mng().flush();
    
}
/*---------------------------------------------------------------------------*/
void SingularityRegionMesher::createSurfaceQuadPatches()
{
    std::vector<SingularityPoint*> singularity_points = m_graph.getPoints();
    
    for (auto i=0; i<singularity_points.size(); i++){
        SingularityPoint* p =singularity_points[i];
        if(p->getGeomType()!=SingularityPoint::SURFACE)
            continue;
        
        
        math::Point loc = p->getLocation();
        
        std::vector<Node> q;
        for (auto s:p->getSlots()){
            if (!s->isOnSurface)
                continue;
            
            math::Vector v = s->direction;
            v.normalize();
            q.push_back(m_line_mesh->newNode(loc+m_edge_size*v));
            
        }// for (auto s:p->getSlots())
        Node origin = m_line_mesh->newNode(loc);
        
        m_sing_pnt_to_patch_center[p] = origin.getID();

        std::vector<TCellID> faces;
        for(auto i=0; i<q.size(); i++){
            const auto j = (i+1)%q.size();
            math::Point qi = q[i].getPoint();
            math::Point qj = q[j].getPoint();
            math::Vector vi(loc,qi);
            math::Vector vj(loc,qj);
            
            Node r = m_line_mesh->newNode(loc+vi+vj);
            Face f= m_line_mesh->newQuad(origin, q[i], r, q[j]);
            
            faces.push_back(f.getID());
            
            m_sing_pnt_to_patch_circle[p].push_back(q[i].getID());
            m_sing_pnt_to_patch_circle[p].push_back(r.getID());
        }
        
        
    }// for (auto p:singularity_points)
    
    
}

/*---------------------------------------------------------------------------*/
void SingularityRegionMesher::createLineTubes()
{
    
    for (auto l:m_graph.getVolumeLines()){
        
        double len = l->length();
        
        int nb_layers = round(len/m_edge_size);
        //============================================================
        // STEP 1 - Initialization of origin and target date
        //============================================================
        //Origin (0) and target (1) patch data
        SingularityPoint*    pnt[2]    = {
            l->getEndPoints()[0],
            l->getEndPoints()[1]
        };
        math::Vector         dir[2]    = {
            l->getTangent(0),
            l->getTangent(1)}
        ;
        TCellID              center[2] = {
            m_sing_pnt_to_patch_center[pnt[0]],
            m_sing_pnt_to_patch_center[pnt[1]]
        };
        std::vector<TCellID> circle[2] = {
            m_sing_pnt_to_patch_circle[pnt[0]],
            m_sing_pnt_to_patch_circle[pnt[1]]
        };
        
        std::vector<math::Point> circle_pnt[2];
        math::Chart chart[2];
        std::vector<math::Vector3d> coords[2];
        for(auto i=0; i<2; i++){
            circle_pnt[i].reserve(circle[i].size());
            for(auto n:circle[i]){
                circle_pnt[i].push_back(m_line_mesh->get<Node>(n).getPoint());
            }
            
            //We build a chart of the origin patch
            math::Point oc = m_line_mesh->get<Node>(center[i]).getPoint();
            math::Point o0 = m_line_mesh->get<Node>(circle[i][0]).getPoint();
            math::Vector axis(oc,o0);
            dir[i].normalize();
            axis.normalize();
            chart[i] = math::Chart(axis,dir[i]);
            math::Vector3d v_x = chart[i].VX();
            math::Vector3d v_y = chart[i].VY();
            math::Vector3d v_z = chart[i].VZ();
            
            //We compute the cooordinates of circle points
            // into the chart representation in the plane.
            coords[i].resize(circle[i].size());
            for(auto ic=0; ic<circle[i].size(); ic++){
                
                Node ni = m_line_mesh->get<Node>(circle[i][ic]);
                math::Vector3d vi(oc, ni.getPoint());
                coords[i][ic]=math::Vector3d(vi.dot(v_x),
                                             vi.dot(v_y),
                                             vi.dot(v_z));
            }
            
        }
        
        //============================================================
        // STEP 2 - We create inner layer from origin to target
        //============================================================
        math::Chart prev_chart(chart[0]);
        
        //to store charts from layer 1 to N-1 viewed from origin
        std::vector<math::Chart> layer_chart;

        std::vector<std::vector<math::Point> > layer_circle_pnt;
        std::vector<std::vector<TCellID> >     layer_circle;
        layer_circle_pnt.resize(nb_layers+1); //+1 for the copy of the target surf
        layer_circle.resize(nb_layers+1); //+1 for the copy of the target surf
        for(auto& lc:layer_circle_pnt){
            lc.resize(circle[0].size());
        }
        for(auto& lc:layer_circle){
            lc.resize(circle[0].size());
        }
        layer_circle_pnt[0]= circle_pnt[0];
        for(auto i=1;i<nb_layers;i++){
            
            double param = ((double)i)/((double)nb_layers);
            math::Vector loc_dir = l->getTangent(param);
            math::Point  cur_loc = l->getPoint  (param);
            
            //Direction is ponderate by Bi2 Bernstein polynomial in order to
            //avoid boundary misalignement at the end points.
            //B02(t) = (1-t)^2
            //B12(t) = 2(1-t)t
            //B22(t) = t^2
            double b02 = (1-param)*(1-param);
            double b12 = 2*(1-param)*param;
            double b22 = param*param;
            math::Vector cur_dir = b02*dir[0] + b12*loc_dir + b22*dir[1];

            // We rotate the previous chart to be aligned as most as
            // possible with the current direction
            math::Chart cur_chart(prev_chart);
            cur_chart.align(cur_dir);
            math::Vector3d v_x = cur_chart.VX();
            math::Vector3d v_y = cur_chart.VY();
            math::Vector3d v_z = cur_chart.VZ();
            //We create the points
            for(auto ic=0; ic<circle[0].size(); ic++){
                math::Vector3d c = coords[0][ic];
                math::Vector3d vi =c.X()*v_x + c.Y()*v_y + c.Z()*v_z;
                math::Point pi = cur_loc+math::Vector(vi.X(),vi.Y(),vi.Z());
                layer_circle_pnt[i][ic]=pi;
            }
            //go to the next layer
            prev_chart=cur_chart;
        }
        
        //We deal with the last circle, the target one
        math::Vector cur_dir = dir[1];
        math::Point  cur_loc = l->getPoint(1);
        // We rotate the previous chart to be aligned as most as
        // possible with the current direction. For that, we look
        // for the correspondance between first point only
        math::Chart cur_chart(prev_chart);
        cur_chart.align(dir[1]);
        math::Vector3d v_x = cur_chart.VX();
        math::Vector3d v_y = cur_chart.VY();
        math::Vector3d v_z = cur_chart.VZ();
        math::Vector3d c = coords[0][0];
        math::Vector3d v0 =c.X()*v_x + c.Y()*v_y + c.Z()*v_z;
        math::Point p0 = cur_loc+math::Vector(v0.X(),v0.Y(),v0.Z());
        //Among all the target points, we look the one which is the
        //closest to p0
        int n_idx=0;
        double min_dist = p0.distance2(circle_pnt[1][0]);
        for(auto in = 1; in<circle_pnt[1].size(); in++){
            double dist_i = p0.distance2(circle_pnt[1][in]);
            if(dist_i<min_dist){
                min_dist = dist_i;
                n_idx=in;
            }
        }
        //We look circle orientation now
        c = coords[0][1];
        math::Vector3d v1 =c.X()*v_x + c.Y()*v_y + c.Z()*v_z;
        math::Point    p1 = cur_loc+math::Vector(v1.X(),v1.Y(),v1.Z());
        math::Point    t0 = circle_pnt[1][n_idx                        ];
        math::Point    t1 = circle_pnt[1][(n_idx+1)%circle_pnt[1].size()];
        math::Vector p01(p0,p1);
        math::Vector t01(t0,t1);
        bool same_orient = (p01.dot(t01)>0);
        //now we match with the nodes of the target patch
        //We create the points
        if(!same_orient){
            std::reverse(circle[1].begin(),circle[1].end());
            n_idx=circle[1].size()-1-n_idx;
        }
        std::vector<TCellID> tmp(circle[1]);
        for(auto ic=0; ic<circle[1].size(); ic++){
            circle[1][ic]=tmp[(n_idx+ic)%circle[1].size()];
        }
        //Now circle 1 and circle 0 numerotation are the same
        for(auto ic=0; ic<circle[0].size(); ic++){
            //we have a - since origin and target patchs are supposed to be
            // oriented in an opposite fashion
            Node t_ni =m_line_mesh->get<Node>(circle[1][ic]);
           
            layer_circle_pnt[nb_layers][ic]=t_ni.getPoint();
            //target circle is also modified
            layer_circle[nb_layers][ic]=t_ni.getID();
        }
        
        //============================================================
        // STEP 3 - We receed the traversal of layers and create pnts
        // using interpolation from data coming from 0->1 and 1->0
        // traversal
        //============================================================
        //Now we create tubular side nodes really as linear interpolation from
        // projected origin nodes and target ones on each layer
        
        //As the target circle has been reoriented according to the
        //original one, we update some quantities.
        math::Point oc = m_line_mesh->get<Node>(center[1]).getPoint();
        math::Point o0 = m_line_mesh->get<Node>(circle[1][0]).getPoint();
        math::Vector axis(oc,o0);
        dir[1].normalize();
        axis.normalize();
        chart[1] = math::Chart(axis,dir[1]);
        v_x = chart[1].VX();
        v_y = chart[1].VY();
        v_z = chart[1].VZ();
        
        //We compute the cooordinates of circle points
        // into the chart representation in the plane.
        coords[1].resize(circle[1].size());
        for(auto ic=0; ic<circle[1].size(); ic++){
            
            Node ni = m_line_mesh->get<Node>(circle[1][ic]);
            math::Vector3d vi(oc, ni.getPoint());
            coords[1][ic]=math::Vector3d(vi.dot(v_x),
                                         vi.dot(v_y),
                                         vi.dot(v_z));
        }
        
         prev_chart =chart[1];
        
//        std::vector<std::vector<math::Point> > layer_circle_pnt;
//        std::vector<std::vector<TCellID> >     layer_circle;

        for(auto i=nb_layers-1;i>0;i--){
            
            double param = ((double)i)/((double)nb_layers);
            math::Vector loc_dir = l->getTangent(param);
            math::Point  cur_loc = l->getPoint  (param);
            
            //Direction is ponderate by Bi2 Bernstein polynomial in order to
            //avoid boundary misalignement at the end points.
            //B02(t) = (1-t)^2
            //B12(t) = 2(1-t)t
            //B22(t) = t^2
            double b02 = (1-param)*(1-param);
            double b12 = 2*(1-param)*param;
            double b22 = param*param;
            math::Vector cur_dir = b02*dir[0] + b12*loc_dir + b22*dir[1];
            
            // We rotate the previous chart to be aligned as most as
            // possible with the current direction
            math::Chart cur_chart(prev_chart);
            cur_chart.align(cur_dir);
            math::Vector3d v_x = cur_chart.VX();
            math::Vector3d v_y = cur_chart.VY();
            math::Vector3d v_z = cur_chart.VZ();
            //We create the points and we interpolate the layer nodes
            //using value coming from the 0->1 traversal and this one
            for(auto ic=0; ic<circle[1].size(); ic++){
                math::Vector3d c = coords[1][ic];
                math::Vector3d vi =c.X()*v_x + c.Y()*v_y + c.Z()*v_z;
                math::Point pi = cur_loc+math::Vector(vi.X(),vi.Y(),vi.Z());
                math::Point oi = layer_circle_pnt[i][ic];
                
                layer_circle_pnt[i][ic]=(1-param)*oi+param*pi;
            }
            //Now
            
            //go to the next layer
            prev_chart=cur_chart;
        }
        
        //============================================================
        // STEP 4 - Nodes are created
        //============================================================
        for(auto i=1; i<nb_layers;i++){

            //current points where created from the original surface
            std::vector<math::Point> pnts = layer_circle_pnt[i];
            layer_circle[i].resize(pnts.size());
            for(auto j=0; j<pnts.size(); j++){
//                math::Point pnt_target =layer_circle_pnt[nb_layers][j];
//                math::Point p =(1-t)*pnts[j]+t*pnt_target;
                Node nj = m_line_mesh->newNode(pnts[j]);
                layer_circle[i][j]=nj.getID();

            }
        }
        //first layer is original nodes
        layer_circle[0]=circle[0];
        //last layer was done previously
        //============================================================
        // STEP 5 - And eventually, we create tubular side faces
        //============================================================
        for(auto i=1; i<nb_layers+1;i++){
            std::vector<TCellID> l1 = layer_circle[i-1];
            std::vector<TCellID> l2 = layer_circle[i  ];
            for(auto j=0; j<l1.size(); j++){
                Node n0 = m_line_mesh->get<Node>(l1[j]);
                Node n1 = m_line_mesh->get<Node>(l1[(j+1)%l1.size()]);
                Node n3 = m_line_mesh->get<Node>(l2[j]);
                Node n2 = m_line_mesh->get<Node>(l2[(j+1)%l1.size()]);
                m_line_mesh->newTriangle(n0,n1,n2);
                m_line_mesh->newTriangle(n0,n2,n3);
            }
        }
    }// for (auto p:singularity_points)
    
    
}
/*----------------------------------------------------------------------------*/
void SingularityRegionMesher::createVolumeSingularityLines()
{
    //at this point we have made all cluster of sing into Singularity3D
    //now we need to compute the lines of sing into separatrices
    int m_face_sepa = m_mesh->getNewMark<Face>();
    
    
    std::vector<SingularityPoint*> singularity_points = m_graph.getPoints();
    
    for (unsigned int i = 0; i < singularity_points.size(); i++){
        SingularityPoint* pnt_i = singularity_points[i];
        
        //We have a  singularity point
        std::vector<SingularityPoint::Slot*> slots_i = pnt_i->getSlots();
        std::cout << "==> We create sing line from sing. pnt " << pnt_i->getLocation()
        << ", having " << slots_i.size() << std::endl;
        for (unsigned int j = 0; j < slots_i.size(); j++){
            
            SingularityPoint::Slot* s_j = slots_i[j];
            if (!s_j->isOnSurface){
                if (!s_j->isLaunched){
                    createOneSingularityLineFrom(pnt_i, s_j, m_face_sepa);
                }
            }
        }
    }
    
    //=========================================================================
    //SMOOTHING OF SEPARATRICES
    //=========================================================================
    std::vector<SingularityLine*> sing_lines = m_graph.getLines();
    for (unsigned int i = 0; i < sing_lines.size(); i++){
        SingularityLine* line_i = sing_lines[i];
        line_i->smooth(500);
    }
    writeOutput("lines");
    m_mesh->unmarkAll<Face>(m_face_sepa);
    m_mesh->freeMark<Face>(m_face_sepa);
    
}

/*----------------------------------------------------------------------------*/
void SingularityRegionMesher::
createOneSingularityLineFrom(SingularityPoint* ASingPoint,
                             SingularityPoint::Slot* ASlot,
                             const int AMark)
{
    //we have a new sep to create here starting with Ftmp
    //first we look if we have to create a new sing, on the boundary,
    //or to retrieve the sing inside the mesh volume
    
    SingularityPoint* starting_singularity_pnt = ASingPoint;
    
    Face starting_face = m_mesh->get<Face>(ASlot->starting_cell_id);
    
    // This way to get the adjacent region is okay both for a volume or a
    // surface singularity point
    std::vector<Region> adj_regions = starting_face.get<Region>();
    Region firstReg;
    for (unsigned int i_reg = 0; i_reg < adj_regions.size(); i_reg++){
        if (!m_mesh->isMarked(adj_regions[i_reg], m_mark_volume_pnt_sing))
            firstReg = adj_regions[i_reg];
    }
    
    
    SingularityLine* new_line = m_graph.newVolumeLine();
    new_line->setNumber(m_graph.getNbLines());
    new_line->addSingularityPoint(starting_singularity_pnt);
    new_line->addDiscretizationPoint(starting_singularity_pnt->getLocation());
    ASlot->isLaunched = true;
    ASlot->line = new_line;
    
    
    Face current_face = starting_face;
    //we get the face center
    math::Point current_pnt = ASlot->location;
    
    // tetrahedral element in which we currently propagate/build the sing. line
    Region current_tet = firstReg;
    
    
    //=====================================================================
    //at this point we have a sep, a point to insert, and a region to
    // continue the sep with we will now iteratively continue until reaching
    // either a boundary or a singularity point
    //=====================================================================
    bool line_done = false;

    Tools tool(m_mesh,0,m_rot);
    
    while (!line_done) {
        new_line->addDiscretizationPoint(current_pnt);
        m_mesh->mark(current_face, AMark);
        //Now look for the opposite face
        int nb_sing_faces = 0;
        
        std::vector<Face> current_tet_faces = current_tet.get<Face>();
        int stop_loop = 0;
        
        for (unsigned int i = 0; i < 4; i++){
            
            Face current_tet_face = current_tet.get<Face>()[i];
            if (!stop_loop && current_tet_face.getID() != current_face.getID()){
                
                if (tool.isFFSingular(current_tet_face)){
                    stop_loop = 1;
                    //We have the new face
                    nb_sing_faces++;
                    current_face = current_tet_face;
                    std::vector<Node> current_face_nodes = current_face.get<Node>();
                    
                    //we get the center of the face as next discretization point
                    current_pnt = current_face.center();
                    
                    //Now, we distinguish the case of a boundary face and an inner face
                    if (current_face.get<Region>().size() == 1){
                        //BOUNDARY
                        //We reached the domain boundary, so we will have to stop after this point creation
                        line_done = true;
                        std::vector<SurfaceSingularityPoint*> surf_sing = m_graph.getSurfacePoints();
                        
                        for (unsigned int i_sing = 0; i_sing < surf_sing.size(); i_sing++){
                            
                            SurfaceSingularityPoint* surf_sing_i = surf_sing[i_sing];
                            Face sing_face = surf_sing_i->getMeshFace();
                            if (sing_face.getID() == current_face.getID()) {
                                
                                //IMPORTANT Fill the singularity slot !!!!
                                bool can_be_added = surf_sing_i->addLineFromVolume(new_line);
                                if (can_be_added){
                                    //to finish the new_line on this sing. point, we must have a compatible free slot
                                    new_line->addSingularityPoint(surf_sing_i);
                                    new_line->addDiscretizationPoint(surf_sing_i->getLocation());
                                    //the line is finished
                                    line_done = true;
                                }
                            }
                        }
                        
                        
                        m_mesh->mark(current_face, AMark);
                    } //if (current_face.get<Region>().size() == 1)
                    else if (m_mesh->isMarked(current_tet_face, m_mark_face_sing)){
                        
                        //The current_tet_face is incident to one region in a cluster and one out
                        //We reached a sing. cluster, so we will have to stop after this point creation
                        
                        new_line->addDiscretizationPoint(current_pnt);
                        m_mesh->mark(current_face, AMark);

                        //need to find the corresponding cluster
                        TCellID id_tet=NullID;
                        if (current_face.get<Region>()[0].getID() == current_tet.getID())
                            current_tet = current_face.get<Region>()[1];
                        else
                            current_tet = current_face.get<Region>()[0];
                        id_tet = current_tet.getID();
                        int point_index = 0;
                        std::vector<VolumeSingularityPoint*> singularities = m_graph.getVolumePoints();
                        for (unsigned int i = 0; i < singularities.size(); i++) {
                            //for each singularity point, we check its mesh cells
                            std::vector<Region> current_cells = singularities[i]->getMeshRegions();
                            for (unsigned int j = 0; j < current_cells.size(); j++){
                                if (current_cells[j].getType() == GMDS_TETRA &&
                                    current_cells[j].getID() == id_tet) {
                                    //we got the right cluster
                                    point_index = i;
                                }
                            }
                        }

                        //new_line->addDiscretizationPoint(singularities[point_index]->getLocation());
                        
                        //IMPORTANT Fill the singularity slot !!!!
                        bool can_be_added = singularities[point_index]->addLine(new_line, current_face);
                        if (can_be_added){
                            //to finish the new_line on this sing. point, we must have a compatible free slot
                            new_line->addSingularityPoint(singularities[point_index]);
                            new_line->addDiscretizationPoint(singularities[point_index]->getLocation());
                            //the line is finished
                            line_done = true;
                        }
                        
                    }
                    else{
                        //continue business as usual
                        if (current_face.get<Region>()[0] == current_tet){
                            current_tet = current_face.get<Region>()[1];
                        }
                        else {
                            current_tet = current_face.get<Region>()[0];
                        }
                    }
                }
            }//if (!stop_loop && current_tet_face.getID() != current_face.getID())
        }
        
    } //while (!line_done)
    Log::mng() << "A sing. line has been built\n";
    Log::mng().flush();
    writeOutput("singularity_line");
}


/*----------------------------------------------------------------------------*/
math::Vector SingularityRegionMesher::
getOutputNormal(const Face& AF, const Region& AR)
{
    std::vector<Node> r_nodes = AR.get<Node>();
    std::vector<Node> f_nodes = AF.get<Node>();
    
    if (r_nodes.size() != 4)
        throw GMDSException("SingularityGraphBuilder::getOutputNormal can only be used on tetrahedral regions");
    if (f_nodes.size() != 3)
        throw GMDSException("SingularityGraphBuilder::getOutputNormal can only be used on triangular faces");
    
    //we go through all the nodes of ARegion to find the one that do not belong
    //to AF
    for (auto n:r_nodes){
        if (n != f_nodes[0] && n != f_nodes[1] && n != f_nodes[2]){
            //n is the node opposite to the face AF
            math::Vector normal = AF.normal();
            math::Vector in_vector(f_nodes[0].getPoint(), n.getPoint());
            if (normal.dot(in_vector)>0.0) {
                return math::Vector(-normal.get(0),
                                    -normal.get(1),
                                    -normal.get(2));
            }
            else {
                return normal;
            }
        }
    }
    throw GMDSException("SingularityGraphBuilder::getOutputNormal unexpected behaviour");
}/*----------------------------------------------------------------------------*/
math::Vector SingularityRegionMesher::
getInputNormal(const Face& AF, const Region& AR)
{
    math::Vector outVec = getOutputNormal(AF, AR);
    return outVec.opp();
}

/*----------------------------------------------------------------------------*/
void SingularityRegionMesher::writeOutput(const std::string& AFileName)
{
    static int out = 0;
    std::string file_name = m_params_gl.output_dir+"/"+AFileName+"_"+to_string(out);
    writeOutputSingle(file_name);
    out++;
}
/*----------------------------------------------------------------------------*/
void SingularityRegionMesher::writeOutputSingle(const std::string& AFileName)
{
    m_graph.createVTKOutputFile(AFileName);
    VTKWriter<IGMesh> writer(*m_line_mesh);
    std::cout<<"NB QUAD: "<<m_line_mesh->getNbFaces()<<std::endl;
    writer.write(m_params_gl.output_dir+"/line_quad", DIM3 | F | N);
}
/*----------------------------------------------------------------------------*/

