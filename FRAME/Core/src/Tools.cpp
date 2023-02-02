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
/*----------------------------------------------------------------------------*/
/*
 * Tools.cpp
 *
 *  Created on: Sept. 18, 2015
 *      Author: F. Ledoux
 */
/*----------------------------------------------------------------------------*/
#include "Tools.h"
#include <GMDS/Math/Line.h>
#include <GMDS/Math/Ray.h>
#include <GMDS/Math/Numerics.h>
#include <GMDS/Math/Plane.h>
#include <GMDS/Math/Segment.h>
#include <GMDS/Math/Quaternion.h>
/*----------------------------------------------------------------------------*/
#include <sstream>
#include <set>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
Tools::Tools(IGMesh* AMesh,
             Variable<math::Cross2D>* AField,
             Variable<math::AxisAngleRotation>* ARotField)
: m_mesh(AMesh), m_field(AField), m_rot_field(ARotField)
{}
/*----------------------------------------------------------------------------*/
math::Chart Tools::computeChartIn(const math::Point& APnt,
                                  const Face&        AFace)
{
    std::vector<Node> n =AFace.get<Node>();
    //===================================================================
    //STEP 1 - We compute the location of APnt into AFace
    //===================================================================
    double coeff[3]={0, 0, 0};

    math::Point::computeBarycentric(n[0].getPoint(),n[1].getPoint(),
                                    n[2].getPoint(), APnt,
                                    coeff[0],coeff[1], coeff[2]);

    //===================================================================
    //STEP 2 - We extract the quaternion representation of each frame  at
    //         the face corners
    //===================================================================
    std::vector<math::AxisAngleRotation> r;
    r.resize(3);
    for(int i=0;i<3;i++){
        r[i]=(*m_rot_field)[n[i].getID()];
    }

    std::vector<math::Quaternion> qs;
    qs.resize(3);
    for(int i=0;i<3;i++){
        qs[i]=math::Quaternion(r[i].toChart());
    }

    std::vector<TCoord> ws;
    ws.resize(3);
    ws[0]=coeff[0];
    ws[1]=coeff[1];
    ws[2]=coeff[2];
    
    //===================================================================
    //STEP 3 - We compute the mean quaternion and return the corresponding
    //         chart
    //===================================================================
    math::Quaternion q = math::Quaternion::mean(qs,
                                                ws);
    
    return math::Chart(q);
    
}
/*---------------------------------------------------------------------------*/
char Tools::orient3d(const gmds::math::Point& AP0,
                     const gmds::math::Point& AP1,
                     const gmds::math::Point& AP2,
                     const gmds::math::Point& AP3)
{
    double p0[3]={AP0.X(),AP0.Y(),AP0.Z()};
    double p1[3]={AP1.X(),AP1.Y(),AP1.Z()};
    double p2[3]={AP2.X(),AP2.Y(),AP2.Z()};
    double p3[3]={AP3.X(),AP3.Y(),AP3.Z()};
    return GEO::PCK::orient_3d(p0, p1, p2, p3);
}
/*---------------------------------------------------------------------------*/
bool Tools::isIn(const math::Point& AP,
                 const Face& ATri,
                 bool& AOnEdge0,
                 bool& AOnEdge1,
                 bool& AOnEdge2)
{
    std::vector<Node> n = ATri.get<Node>();
    
    math::Point pnt[3] = {
        n[0].getPoint(),
        n[1].getPoint(),
        n[2].getPoint()
    };
    math::Plane plane(pnt[0],pnt[1],pnt[2]);
    math::Point p = plane.project(AP);
    
    // We look if p is insie or outside of the triangle defined by ATri
    math::Vector normal = plane.getNormal();
    normal.normalize();
    char ori[3] ={
        orient3d(p, pnt[1] , pnt[2] , p+normal),
        orient3d(p, pnt[2] , pnt[0] , p+normal),
        orient3d(p, pnt[0] , pnt[1] , p+normal)
    };
    if ((ori[0] >= 0 && ori[1] >= 0 && ori[2] >= 0 ) ||
        (ori[0] <= 0 && ori[1] <= 0 && ori[2] <= 0 ) ) {
        if(ori[0]==GEO::ZERO){
            AOnEdge0=true;
        }
        if(ori[1]==GEO::ZERO){
            AOnEdge1=true;
        }
        if(ori[2]==GEO::ZERO){
            AOnEdge2=true;
        }
        return true;
    }
    return false;
}

/*----------------------------------------------------------------------------*/
void Tools::computeFuzzyHeuns(const math::Point&                AFromPnt,
                              const math::Vector3d&             AFromDir,
                              const std::vector<Face>&          AFaces,
                              const std::vector<math::Triangle>&ATri,
                              math::Point&                      AToPnt,
                              math::Vector3d&                   AToDir,
                              int&                              AToFaceId)
{
    
    // check whether the line intersects the triangle
    math::Line ray(AFromPnt,math::Vector(AFromDir.X(),
                                         AFromDir.Y(),
                                         AFromDir.Z()));
    
//    std::cout<<"Line "<<ray.getFirstPoint()<<" - "
//    <<ray.getSecondPoint()<<std::endl;
    //===================================================================
    double param[4]     = {-1, -1,-1,-1};
    math::Point p[4];
    for(auto i=0; i<4; i++){
        math::Plane pli = ATri[i].getPlaneIncluding();
        if(ray.intersect3D(pli, p[i], param[i])){
            //intersect the plane, but the triangle??
            //before computing the bar coordinate, we eliminate intersection
            //at the infinity almost due to almost parallel ray and plane
            std::vector<Node> n = AFaces[i].get<Node>();
            math::Point pn[3] ={
                n[0].getPoint(),
                n[1].getPoint(),
                n[2].getPoint()};
            //We compute the barycentric coords.
            bool on_edge[3]={false,false,false};
            if(!isIn(p[i],AFaces[i], on_edge[0],on_edge[1],on_edge[2])){
                param[i]=-1;

            };
//            double coeff[3]={0, 0, 0};
//            std::cout<<"HERE IN FACE"<<AFaces[i]<<" with point "
//            <<p[i]<<" and param "<<param[i]<<std::endl;
//            p[i]=pli.project(p[i]);
//            math::Point::computeBarycentric(pn[0], pn[1], pn[2], p[i],
//                                            coeff[0],coeff[1],
//                                            coeff[2]);
//            
//            if(coeff[0]<0 || coeff[1]<0 || coeff[2]<0 ){
//                param[i]=-1;
//            }
            
        }
        else{
            param[i]=-1;
        }
//          std::cout<<"\t intersection with face "<<AFaces[i].getID()<<" -> "<<param[i]<<std::endl;
    }
    double best_param =param[0];
    auto out_index = 0;

    for(auto i=1; i<4; i++){
        if(param[i]>best_param){
            out_index=i;
            best_param=param[i];
        }
    }
    if(best_param<=1e-8)
        throw GMDSException("Tools::computeFuzzyHeuns: No out face (1)");
    AToPnt=p[out_index];
//    std::cout<<"\t OUT PNT: "<<AToPnt<<std::endl;
    //===================================================================
    //compute the frame in out_pnt
    math::Chart ci = computeChartIn(AToPnt, AFaces[out_index]);
    
    //===================================================================
    //among the 6 vectors of ci, we take the one which is the
    // best aligned with dir and we start the process a second
    // time
    math::Vector3d ci_vectors[6] = {
        ci.VX(), -ci.VX(),  ci.VY(),
        -ci.VY(),  ci.VZ(), -ci.VZ()
    };
    math::Vector3d heuns_corr = ci_vectors[0];
    double best_align_dot = AFromDir.dot(ci_vectors[0]);
    
    for(int i=0; i<6;i++){
        if(best_align_dot<AFromDir.dot(ci_vectors[i])){
            heuns_corr = ci_vectors[i];
            best_align_dot = AFromDir.dot(ci_vectors[i]);
        }
    }
    
    AToDir= heuns_corr;
    AToFaceId = out_index;
}


/*----------------------------------------------------------------------------*/
void Tools::computeFuzzyHeuns(const math::Point&                AFromPnt,
                              const math::Vector3d&             AFromDir,
                              const std::vector<Edge>&          AEdges,
                              const std::vector<math::Segment>& ASeg,
                              math::Point&                      AToPnt,
                              math::Vector3d&                   AToDir,
                              int&                              AToFaceId)
{
    //===================================================================
    math::Point ray_from = AFromPnt;
    math::Point ray_to = AFromPnt+math::Vector(AFromDir.X(),
                                               AFromDir.Y(),
                                               AFromDir.Z());
    double param[3] = {-1, -1, -1};
    double param_seg[3] = {-1, -1, -1};
    math::Point p[3];
    for(auto i=0; i<3; i++){
//        std::cout<<"Edge "<<AEdges[i].getIDs<Node>()[0]<<" - "
//        <<AEdges[i].getIDs<Node>()[1]<<std::endl;
        //bool ok =
        
        math::Point s0 =ASeg[i].getPoint(0);
        math::Point s1 =ASeg[i].getPoint(1);
        math::Point ray_to_i = math::Plane(ray_from,s0,s1).project(ray_to);
        math::Ray ray(ray_from,ray_to_i);

        ray.intersect3D(ASeg[i], p[i], param_seg[i], param[i]);
//        if(ok){
//            std::cout<<"INTERSECT Edge "<<AEdges[i].getIDs<Node>()[0]<<" - "
//            <<AEdges[i].getIDs<Node>()[1]<<": "<<param_seg[i]<<"(s), "
//            <<param[i]<<"(r)"<<std::endl;
//        }
//        else{
//            std::cout<<"NOT INTERSECT Edge "<<AEdges[i].getIDs<Node>()[0]<<" - "
//            <<AEdges[i].getIDs<Node>()[1]<<": "<<param_seg[i]<<"(s), "
//            <<param[i]<<"(r)"<<std::endl;
//
//        }
    }
    double best_param =param[0];
    auto out_index = 0;
    
    for(auto i=1; i<3; i++){
        if(param[i]>best_param){
            out_index=i;
            best_param=param[i];
        }
    }
    if(best_param<=1e-8){

        for(int i=0; i<3;i++){
            std::cout<<"Edge "<<AEdges[i].getIDs<Node>()[0]<<" - "
            <<AEdges[i].getIDs<Node>()[1]<<": "<<param_seg[i]<<"(s), "
            <<param[i]<<"(r)"<<std::endl;
        }
        throw GMDSException("Tools::computeFuzzyHeuns: No out face (2)");
    }
    double out_param = param_seg[out_index];
    AToPnt = p[out_index];
    //===================================================================
    //compute the frame in out_pnt
    Edge e = AEdges[out_index];
    std::vector<TCellID> ne = e.getIDs<Node>();
    math::AxisAngleRotation r[2] = {
        (*m_rot_field)[ne[0]],
        (*m_rot_field)[ne[1]]
    };
    
    std::vector<math::Quaternion> qs;
    qs.resize(2);
    for(int i=0;i<2;i++){
        qs[i]=math::Quaternion(r[i].toChart());
    }
    
    std::vector<TCoord> ws;
    ws.resize(2);
    ws[0]=out_param;
    ws[1]=1-out_param;

    math::Quaternion q = math::Quaternion::mean(qs, ws);
    
    math::Chart ci = math::Chart(q);
    //===================================================================
    //among the 6 vectors of ci, we take the one which is the
    // best aligned with dir and we start the process a second
    // time
    math::Vector3d ci_vectors[6] = {
        ci.VX(), -ci.VX(),  ci.VY(),
        -ci.VY(),  ci.VZ(), -ci.VZ()
    };
    math::Vector3d heuns_corr = ci_vectors[0];
    double best_align_dot = AFromDir.dot(ci_vectors[0]);
    
    for(int i=0; i<6;i++){
        if(best_align_dot<AFromDir.dot(ci_vectors[i])){
            heuns_corr = ci_vectors[i];
            best_align_dot = AFromDir.dot(ci_vectors[i]);
        }
    }
    
    AToDir= heuns_corr;
    AToFaceId = out_index;
}
/*----------------------------------------------------------------------------*/
bool Tools::followFlow(const PointVolumetricData& AData,
                       const double               AMaxDist,
                       math::Point&               APnt)
{
//   std::cout<<"================FLOW VOLUME ================"<<std::endl;
    double current_dist=0;
    math::Vector3d dir(AData.dir);
    math::Point from = AData.pnt;
    gmds::Region current_cell = AData.tet;
    
    while (current_dist<AMaxDist) {
//        if(isFFSingular(current_cell))
//            return false;

        //to avoid starting in the wrong being on an edge
        from  = 0.95*from+ 0.05*current_cell.center();

        //================================================
        // We build the triangular shapes corresponding to each
        // tet face
        std::vector<Face> cf = current_cell.get<Face>();
//        std::cout<<"------------------"<<std::endl;
//        std::cout<<"From point "<<from<<std::endl;
//        math::Vector dir_v(dir.X(),dir.Y(), dir.Z());
//        std::cout<<"        to "<<from+dir_v<<std::endl;
//        std::cout<<"Region "<<current_cell.getID()<<std::endl;
        
        std::vector<math::Triangle> t;
        t.resize(4);
        for(unsigned int i=0; i<4; i++){
            std::vector<Node> cn = cf[i].get<Node>();
            t[i]= math::Triangle(cn[0].getPoint(),
                                 cn[1].getPoint(),
                                 cn[2].getPoint());
            
        }
        //================================================
        // PHASE 1
        //================================================
        math::Point out_pnt;
        math::Vector3d out_dir;
        int out_index=-1;
        computeFuzzyHeuns(from, dir, cf, t,
                          out_pnt, out_dir, out_index);
        //================================================
        // PHASE 2
        //================================================
        dir = 0.5*dir + 0.5*out_dir;
        computeFuzzyHeuns(from, dir, cf, t,
                          out_pnt, out_dir, out_index);
        
        //================================================
        // FINALIZATION
        //================================================
        // Now we have our point to get out of current_cell
        math::Vector v(from,out_pnt);
        if(current_dist+v.norm() >=AMaxDist){
            v.normalize();
            double remaining_dist = AMaxDist-current_dist;
            
            //we stop in this tet
            APnt = from+remaining_dist*v;
            return true;
        }
        else{
            Face out_face = cf[out_index];
     //       std::cout<<"Out face: "<<out_face.getID()<<std::endl;
            //Fuzzy approach to avoid topological cases
            math::Point new_from =0.95*out_pnt+ 0.05*out_face.center();
            dir = out_dir;
            current_dist +=math::Vector(from,new_from).norm();
            from  = new_from;
            //We compute the next current cell using the face we go
            //from
            std::vector<Region> adj_out_face = out_face.get<Region>();
            if(adj_out_face.size()==1){
                //we are on the boudary even if we have not reach the
                //expected distance
                APnt = from;
                return true;
            }
            else{
                if(adj_out_face[0].getID()==current_cell.getID())
                    current_cell=adj_out_face[1];
                else
                    current_cell=adj_out_face[0];
            }
   //         std::cout<<"\t -->next region: "<<current_cell.getID()<<std::endl;
            //not the same current dist due to the fuzzy displacement
            if(current_dist >=AMaxDist){
                //we stop here so
                APnt = from;
                return true;
            }
        }//else
    };
    
    throw GMDSException("Error in Tools::followFlow");
    
}
/*---------------------------------------------------------------------------*/
math::Chart::Mapping Tools::getRij(const TCellID AFrom,
                                   const TCellID ATo) const
{
    math::AxisAngleRotation rotation_from = (*m_rot_field)[AFrom];
    math::AxisAngleRotation rotation_to   = (*m_rot_field)[ATo];
    
    math::Chart chart_from  = rotation_from.toChart();
    math::Chart chart_to    = rotation_to.toChart();
    
    return math::Chart::Mapping(chart_from,chart_to);
}
/*---------------------------------------------------------------------------*/
bool Tools::isFFSingular(const Face& AF)
{
    std::vector<TCellID> n = AF.getIDs<Node>();
    math::Chart::Mapping m01 = getRij(n[0], n[1]);
    math::Chart::Mapping m12 = getRij(n[1], n[2]);
    math::Chart::Mapping m20 = getRij(n[2], n[0]);
    
    return !(m20*m12*m01).isIdentity();
}
/*---------------------------------------------------------------------------*/
bool Tools::isFFSingular(const Region& AR)
{
    std::vector<Face> f = AR.get<Face>();
    return (isFFSingular(f[0]) ||
            isFFSingular(f[1]) ||
            isFFSingular(f[2]) ||
            isFFSingular(f[3]));
}
/*----------------------------------------------------------------------------*/
bool Tools::followFlow(const PointSurfacicData& AData,
                       const double             AMaxDist,
                       const int                AMarkEdgeOnCurve,
                       math::Point&             APnt)
{
//    std::cout<<"==================== FLOW SURF ===================="<<std::endl;

    double current_dist=0;
    math::Vector3d dir(AData.dir);
    math::Point from = AData.pnt;
    gmds::Face current_cell = AData.tri;

    
    while (current_dist<AMaxDist) {
//        std::cout<<"************************"<<std::endl;
//        if(isFFSingular(current_cell))
//            return false;
        //to avoid starting in the wrong being on an edge
        from  = 0.95*from+ 0.05*current_cell.center();
        
        // out_dir lives in the plan of the face we come from, while
        // dir must live in the plan of the next face, now current_cell
        math::Point p (from.X()+dir.X(),
                       from.Y()+dir.Y(),
                       from.Z()+dir.Z());
        std::vector<Node> cur_nodes = current_cell.get<Node>();
        math::Plane pl(cur_nodes[0].getPoint(),
                       cur_nodes[1].getPoint(),
                       cur_nodes[2].getPoint());
        math::Point proj =pl.project(p);
        dir = math::Vector3d(from,proj);
        

        
//        std::cout<<"From point "<<from<<std::endl;
//        math::Vector dir_v(dir.X(),dir.Y(), dir.Z());
//        std::cout<<"        to "<<from+dir_v<<std::endl;
//        std::cout<<"Face "<<current_cell.getID()<<std::endl;
        //================================================
        // We build the triangular shapes corresponding to each
        // tet face
        std::vector<Edge> cf = current_cell.get<Edge>();
        std::vector<math::Segment> s;
        s.resize(3);
        for(unsigned int i=0; i<3; i++){
            std::vector<Node> cn = cf[i].get<Node>();
            s[i]= math::Segment(cn[0].getPoint(), cn[1].getPoint());
        }
        //================================================
        // PHASE 1
        //================================================
        math::Point out_pnt;
        math::Vector3d out_dir;
        int out_index=-1;
        computeFuzzyHeuns(from, dir, cf, s,
                          out_pnt, out_dir, out_index);
        //================================================
        // PHASE 2
        //================================================
        dir = 0.5*dir + 0.5*out_dir;
        math::Point p2 (from.X()+dir.X(),
                        from.Y()+dir.Y(),
                        from.Z()+dir.Z());
        dir = math::Vector3d(from,pl.project(p2));

        computeFuzzyHeuns(from, dir, cf, s,
                          out_pnt, out_dir, out_index);
        
        //================================================
        // FINALIZATION
        //================================================
        // Now we have our point to get out of current_cell
        math::Vector v(from,out_pnt);
//        std::cout<<"OUT PNT: "<<out_pnt<<std::endl;
        if(current_dist+v.norm() >=AMaxDist){
            v.normalize();
            double remaining_dist = AMaxDist-current_dist;
            
            //we stop in this tet
            APnt = from+remaining_dist*v;
            return true;
        }
        else{
            Edge out_edge = cf[out_index];
            //Fuzzy approach to avoid topological cases
            from  = 0.95*out_pnt+ 0.05*out_edge.center();
            current_dist +=v.norm();
            //We compute the next current cell using the face we go
            //from
            std::vector<Face> adj_out_edge = out_edge.get<Face>();
            std::vector<Face> adj_bnd;
            for(auto f:adj_out_edge){
                std::vector<Region> adj_f = f.get<Region>();
                if(adj_f.size()==1)
                    adj_bnd.push_back(f);
            }
            if(adj_bnd.size()==1 ||
               m_mesh->isMarked(out_edge, AMarkEdgeOnCurve)){
                //we are on the boudary even if we have not reach the
                //expected distance
                APnt = from;
                return true;
            }
            else{
                if(adj_bnd[0].getID()==current_cell.getID())
                    current_cell=adj_bnd[1];
                else
                    current_cell=adj_bnd[0];
            }
            dir=out_dir;

        }
    };
    
    throw GMDSException("Error in Tools::followFlow");
    
}
/*----------------------------------------------------------------------------*/
void Tools::traverseTriangle(const Face&         AFace,
                             const math::Point&  AInPnt,
                             const math::Vector& AInVec,
                             const int                 AInCellDim,
                             const TCellID       AInCellID,
                             math::Point&        AOutPnt,
                             math::Vector&       AOutVec,
                             int&                      AOutCellDim,
                             TCellID&            AOutCellID)
{
    if(AInCellDim==0){
        Node from_node = m_mesh->get<Node>(AInCellID);
        traverseTriangle(AFace, from_node, AInPnt, AInVec,
                         AOutPnt, AOutVec, AOutCellDim, AOutCellID);
    }
    else if(AInCellDim==1){
        Edge from_edge = m_mesh->get<Edge>(AInCellID);
        traverseTriangle(AFace, from_edge, AInPnt, AInVec,
                         AOutPnt, AOutVec, AOutCellDim, AOutCellID);
    }
    else
        throw GMDSException("findNextCell: we can only comes from a node or an edge");
}

/*----------------------------------------------------------------------------*/
void Tools::
traverseTriangle(const Face&         AFace,
                 const Node&         ANode,
                 const math::Point&  AInPnt,
                 const math::Vector& AInVec,
                 math::Point&        AOutPnt,
                 math::Vector&       AOutVec,
                 int&                      AOutCellDim,
                 TCellID&            AOutCellID)
{
    Node node_from = ANode;
    //=====================================================================
    // We look for the opposite edge
    //=====================================================================
    Edge opposite_edge;
    std::vector<Edge> current_edges = AFace.get<Edge>();
    for (unsigned int i = 0; i < current_edges.size(); i++) {
        Edge ei = current_edges[i];
        std::vector<Node> ei_nodes = ei.get<Node>();
        if (node_from != ei_nodes[0] && node_from != ei_nodes[1])
            opposite_edge = ei;
    }
    
    //=====================================================================
    // Other nodes are the nodes of the opposite edge
    //=====================================================================
    std::vector<Node> other_nodes = opposite_edge.get<Node>();
    
    //=====================================================================
    // STEP 2: We compute a first out pnt and vector
    //=====================================================================
    //  std::cout << "========================= Heun's 1 ===================="
    //	    << std::endl;
    math::Point  out_pnt;
    math::Vector out_vec;
    math::Vector in_vec = AInVec;
    in_vec.normalize();
    int out_id = heunsComputation(node_from,
                                  AInPnt, in_vec,
                                  other_nodes[0],
                                  other_nodes[1],
                                  opposite_edge,
                                  out_pnt, out_vec);
    
    //  std::cout << "First computed out point: "<< out_pnt << std::endl;
    //  std::cout << "First computed out vect: " << out_vec << std::endl;
    //=====================================================================
    // STEP 3: We recompute out pnt and vector using Heun's method
    //=====================================================================
    //  std::cout << "========================= Heun's 2 ====================" << std::endl;
    //  std::cout << "Recomputation via Heun's method" << std::endl;
    //math::Vector out_vec = cross(out_q, ACurrentFace.normal()).closestVector(AINVec);
    
    //We compute the approximation only if we are always in the face
    if (out_id== 3){
        math::Vector v_Heun = 0.5*AInVec + 0.5*out_vec;
        v_Heun.normalize();
        
        out_id = heunsComputation(node_from,
                                  AInPnt, v_Heun,
                                  other_nodes[0],
                                  other_nodes[1],
                                  opposite_edge,
                                  out_pnt, out_vec);
    }
    
    //  std::cout << "Second computed out point: " << out_pnt << std::endl;
    
    AOutPnt = out_pnt;
    AOutVec = out_vec;
    
    if (out_id == 1 || out_id == 2)
        AOutCellDim = 0;
    else if (out_id == 3)
        AOutCellDim = 1;
    else {
        std::cout<<"OUT ID: "<<out_id<<std::endl;
        throw GMDSException("Wrong output cell dimension in Heun's method");
    }
    
    if (out_id == 1) {
        AOutCellID = other_nodes[0].getID();
    }
    else if (out_id == 2) {
        AOutCellID = other_nodes[1].getID();
    }
    else if (out_id == 3) {
        AOutCellID = opposite_edge.getID();
    }
    
}
/*----------------------------------------------------------------------------*/
void Tools::
traverseTriangle(const Face&         AFace,
                 const Edge&         AEdge,
                 const math::Point&  AInPnt,
                 const math::Vector& AInVec,
                 math::Point&        AOutPnt,
                 math::Vector&       AOutVec,
                 int&                      AOutCellDim,
                 TCellID&            AOutCellID)
{
    
    Edge edge_in = AEdge;
    Edge other_edges[2];
    int nb_edges = 0;
    std::vector<Edge> current_edges = AFace.get<Edge>();
    for (unsigned int i = 0; i < current_edges.size(); i++)
    {
        if (current_edges[i] != edge_in)
            other_edges[nb_edges++] = current_edges[i];
    }
    
    
    Node other_node;
    std::vector<Node> current_nodes = AFace.get<Node>();
    std::vector<Node> nodes_in = edge_in.get<Node>();
    for (unsigned int i = 0; i < current_nodes.size(); i++){
        Node current_node = current_nodes[i];
        if (current_node != nodes_in[0] && current_node != nodes_in[1])
            other_node = current_node;
    }
    
    //  std::cout << "EDGE IN  : " << nodes_in[0].getID()
    //	    << ", " << nodes_in[1].getID() << std::endl;
    
    std::vector<Node> nodes_out0 = other_edges[0].get<Node>();
    std::vector<Node> nodes_out1 = other_edges[1].get<Node>();
    //  std::cout << "OTHER TRIANGLE NODE: " << other_node.getID() << std::endl;
    
    //=====================================================================
    // STEP 2: We compute a first out pnt and vector
    //=====================================================================
    //  std::cout << "========================= Heun's 1 ===================="
    //	    << std::endl;
    math::Point  out_pnt;
    math::Vector out_vec;
    math::Vector in_vec = AInVec;
    in_vec.normalize();
    int out_id = heunsComputation(edge_in,
                                  AInPnt, in_vec,
                                  other_node,
                                  edge_in.get<Node>()[0],
                                  edge_in.get<Node>()[1],
                                  other_edges[0],
                                  other_edges[1],
                                  out_pnt, out_vec);
    
    //  std::cout << "First computed out point: " << out_pnt << std::endl;
    //  std::cout << "First computed out vect: " << out_vec << std::endl;
    //=====================================================================
    // STEP 3: We recompute out pnt and vector using Heun's method
    //=====================================================================
    //  std::cout << "========================= Heun's 2 ====================" << std::endl;
    //  std::cout << "Recomputation via Heun's method" << std::endl;
    
    if (out_id == 1 || out_id == 4 || out_id == 5){
        
        math::Vector v_Heun = 0.5*AInVec + 0.5*out_vec;
        v_Heun.normalize();
        
        out_id = heunsComputation(edge_in,
                                  AInPnt, v_Heun,
                                  other_node,
                                  edge_in.get<Node>()[0],
                                  edge_in.get<Node>()[1],
                                  other_edges[0],
                                  other_edges[1],
                                  out_pnt, out_vec);
    }
    //std::cout << "Second computed out point: " << out_pnt << std::endl;
    
    
    AOutPnt = out_pnt;
    AOutVec = out_vec;
    
    if (out_id == 1 || out_id == 2 || out_id == 3)
        AOutCellDim = 0;
    else if (out_id == 4 || out_id == 5)
        AOutCellDim = 1;
    else
        throw GMDSException("Wrong output cell dimension in Heun's method");
    
    if (out_id == 1) {
        AOutCellID = other_node.getID();
    }
    else if (out_id == 2) {
        AOutCellID = edge_in.get<Node>()[0].getID();
    }
    else if (out_id == 3) {
        AOutCellID = edge_in.get<Node>()[1].getID();
    }
    else if (out_id == 4) {
        AOutCellID = other_edges[0].getID();
    }
    else {
        AOutCellID = other_edges[1].getID();
    }
    
}
/*----------------------------------------------------------------------------*/
int Tools::
heunsComputation(const Edge&         AINEdge,
                 const math::Point&  AINPnt,
                 const math::Vector& AINVec,
                 const Node&         AOPPNode,
                 const Node&         AINNode1,
                 const Node&         AINNode2,
                 const Edge&         AOUTEdge1,
                 const Edge&         AOUTEdge2,
                 math::Point&        AOUTPnt,
                 math::Vector&       AOUTVec)

{
    math::Point  pout = AINPnt + AINVec;
    
    // std::cout << "===================================" << std::endl;
    // std::cout << "compute OUT cross from edge "<<AINEdge.getID() << std::endl;
    // std::cout << "\t AIN Pnt: " << AINPnt << std::endl;
    // std::cout << "\t AIN Vec: " << AINVec << std::endl;
    // std::cout << "\t \t out pnt: " << pout << std::endl;
    // std::cout << "\t Node candidate : " << AOPPNode.getID() << std::endl;
    // std::cout << "\t Edge candidate1: "
    // 	    << AOUTEdge1.get<Node>()[0].getID() << ",  "
    // 	    << AOUTEdge1.get<Node>()[1].getID() << std::endl;
    // std::cout << "\t Edge candidate2: "
    // 	    << AOUTEdge2.get<Node>()[0].getID() << ",  "
    // 	    << AOUTEdge2.get<Node>()[1].getID() << std::endl;
    // std::cout << "Node " << AOPPNode.getID() << " - " << AOPPNode.getPoint() << std::endl;
    // std::cout << "Node " << AINNode1.getID() << " - " << AINNode1.getPoint() << std::endl;
    // std::cout << "Node " << AINNode2.getID() << " - " << AINNode2.getPoint() << std::endl;
    
    math::Vector v_in = AINVec;
    v_in.normalize();
    
    //================================================
    // Go through the opposite node AOPPNode???
    //================================================
    math::Point opp_node_loc = AOPPNode.getPoint();
    math::Vector v_opp(AINPnt, opp_node_loc);
    v_opp.normalize();
    
    if (math::near(v_opp.dot(v_in) - 1,0)) {
        //   std::cout << "INTERSECT NODE --> out in a node" << std::endl;
        // intersected point = node
        AOUTPnt = opp_node_loc;
        computeOutVectorAtPoint(AOPPNode,AINVec, AOUTVec);
        return 1;
    }
    //================================================
    // Go through AINNode1???
    //================================================
    math::Point node1_loc = AINNode1.getPoint();
    math::Vector v1(AINPnt, node1_loc);
    v1.normalize();
    
    if (math::near(v1.dot(v_in) - 1,0)) {
        //    std::cout << "INTERSECT NODE --> out in a node" << std::endl;
        // intersected point = node
        AOUTPnt = node1_loc;
        computeOutVectorAtPoint(AINNode1,AINVec, AOUTVec);
        return 2;
    }
    //================================================
    // Go through AINNode2???
    //================================================
    math::Point node2_loc = AINNode2.getPoint();
    math::Vector v2(AINPnt, node2_loc);
    v2.normalize();
    
    if (math::near(v2.dot(v_in) - 1,0)) {
        //    std::cout << "INTERSECT NODE --> out in a node" << std::endl;
        // intersected point = node
        AOUTPnt = node2_loc;
        computeOutVectorAtPoint(AINNode2,AINVec, AOUTVec);
        return 3;
    }
    
    
    //================================================
    // Go through the first edge ???
    // And not through one of its end points due to
    // previous tests.
    //================================================
    bool intersectEdge1 = computeOutVectorFromRayAndEdge(AOUTEdge1,
                                                         AINPnt,
                                                         AINVec,
                                                         AOUTPnt,
                                                         AOUTVec);
    
    // std::cout<<"First intersection: "<<intersectEdge1<<std::endl;
    if (intersectEdge1)
        return 4;
    //================================================
    // Go through the second edge ???
    // And not through one of its end points due to
    // previous tests.
    //================================================
    bool intersectEdge2 = computeOutVectorFromRayAndEdge(AOUTEdge2,
                                                         AINPnt,
                                                         AINVec,
                                                         AOUTPnt,
                                                         AOUTVec);
    //  std::cout<<"Second intersection: "<<intersectEdge2<<std::endl;
    if (intersectEdge2)
        return 5;
    
    
    //================================================
    // If we arrive here, it means that the intersection
    // of our ray is out of the triangle
    //
    // Our choice is to take the closest end point
    //================================================
    double value_opp = v_opp.dot(v_in);
    double value_1 = v1.dot(v_in);
    double value_2 = v2.dot(v_in);
    
    if(value_opp>= value_1 && value_opp>=value_2 ) {
        //the outpoint is the first one
        AOUTPnt = opp_node_loc;
        computeOutVectorAtPoint(AOPPNode,AINVec, AOUTVec);
        return 1;
    }
    else if(value_1>= value_opp && value_1 >= value_2 ) {
        //the outpoint is the second one
        AOUTPnt = node1_loc;
        computeOutVectorAtPoint(AINNode1,AINVec, AOUTVec);
        return 2;
    }
    else {
        AOUTPnt = node2_loc;
        computeOutVectorAtPoint(AINNode2,AINVec, AOUTVec);
        return 3;
    }
    throw GMDSException("ERROR: No out point in Tools::heunsComputation");
}
/*----------------------------------------------------------------------------*/
int Tools::
heunsComputation(const Node&         AFROMNode,
                 const math::Point&  AINPnt,
                 const math::Vector& AINVec,
                 const Node&         AOPPNode1,
                 const Node&         AOPPNode2,
                 const Edge&         AOPPEdge,
                 math::Point&        AOUTPnt,
                 math::Vector&       AOUTVec)

{
    math::Point  pout = AINPnt + AINVec;
    
    // std::cout << "===================================" << std::endl;
    // std::cout << "compute OUT cross from node"<<AFROMNode.getID() << std::endl;
    // std::cout << "\t AIN Pnt: " << AINPnt << std::endl;
    // std::cout << "\t AIN Vec: " << AINVec << std::endl;
    // std::cout << "\t \t out pnt: " << pout << std::endl;
    // std::cout << "\t Node candidate 1: " << AOPPNode1.getID() << std::endl;
    // std::cout << "\t Node candidate 2: " << AOPPNode2.getID() << std::endl;
    // std::cout << "\t Edge candidate: "
    // 	    << AOPPEdge.get<Node>()[0].getID() << ",  "
    // 	    << AOPPEdge.get<Node>()[1].getID() << std::endl;
    
    math::Vector v_in = AINVec;
    v_in.normalize();
    
    //================================================
    // Go through the opposite node AOPPNode1???
    //================================================
    math::Point opp_node_loc1 = AOPPNode1.getPoint();
    math::Vector v_opp1(AINPnt, opp_node_loc1);
    v_opp1.normalize();
    
    if (math::near(v_opp1.dot(v_in) - 1,0)) {
        //   std::cout << "INTERSECT NODE --> out in a node" << std::endl;
        // intersected point = node
        AOUTPnt = opp_node_loc1;
        computeOutVectorAtPoint(AOPPNode1,AINVec, AOUTVec);
        return 1;
    }
    //================================================
    // Go through the opposite node AOPPNode1???
    //================================================
    math::Point opp_node_loc2 = AOPPNode2.getPoint();
    math::Vector v_opp2(AINPnt, opp_node_loc2);
    v_opp1.normalize();
    
    if (math::near(v_opp2.dot(v_in) - 1,0)) {
        //  std::cout << "INTERSECT NODE --> out in a node" << std::endl;
        // intersected point = node
        AOUTPnt = opp_node_loc2;    
        computeOutVectorAtPoint(AOPPNode2,AINVec, AOUTVec);
        return 2;
    }
    
    //================================================
    // Go through the opposite edge
    // And not through one of its end points due to 
    // previous tests.
    //================================================
    bool intersectEdge = computeOutVectorFromRayAndEdge(AOPPEdge,
                                                        AINPnt,
                                                        AINVec,
                                                        AOUTPnt,		
                                                        AOUTVec);
    
    //  std::cout<<"Edge intersection: "<<intersectEdge<<std::endl;
    if (intersectEdge)
        return 3;
    
    
    //================================================
    // If we arrive here, it means that the intersection 
    // of our ray with the opposite edge is out of the 
    // segment defined by AOppEdge
    // Our choice is to take the closest end point
    //================================================
    double value1 = v_opp1.dot(v_in);
    double value2 = v_opp2.dot(v_in);
    if(value1>value2) {
        //the outpoint is the first one  
        AOUTPnt = opp_node_loc1;    
        computeOutVectorAtPoint(AOPPNode1,AINVec, AOUTVec);
        return 1;
    }
    else {
        //the outpoint is the second one  
        AOUTPnt = opp_node_loc2;    
        computeOutVectorAtPoint(AOPPNode2,AINVec, AOUTVec);
        return 2;
    }
    
    return 0;
}
/*----------------------------------------------------------------------------*/
void Tools::
computeOutVectorAtPoint(const Node&         ANode,
                        const math::Vector& AInVec,
                        math::Vector&       AOutVec)
{
    math::Cross2D c =(*m_field)[ANode.getID()]; 
    std::vector<math::Vector> c_vectors = c.componentVectors();
    
    int ref_index = 0;
    double ref_dot_product = AInVec.dot(c_vectors[0]);
    for(int i=1; i<4;i++){
        double dot_product_i =  AInVec.dot(c_vectors[i]);
        if(dot_product_i>ref_dot_product){
            ref_dot_product = dot_product_i;
            ref_index = i;
        }
        
    }
    AOutVec = c_vectors[ref_index];
}

/*----------------------------------------------------------------------------*/
bool Tools::
computeOutVectorFromRayAndEdge(const Edge&         AEdge,
                               const math::Point&  AINPnt,
                               const math::Vector& AINVec,
                               math::Point&        AOUTPnt,
                               math::Vector&       AOUTVec)
{
    math::Ray ray(AINPnt, AINVec);
    
    math::Segment seg(AEdge.get<Node>()[0].getPoint(),
                      AEdge.get<Node>()[1].getPoint());
    
    math::Point pnt_intersection;
    double param_intersection = 0.0;
    bool intersectEdge = ray.intersect2D(seg,
                                         pnt_intersection,
                                         param_intersection);
    
    
    if(!intersectEdge)
        return false;
    
    AOUTPnt = pnt_intersection;
    
    Node n0 = AEdge.get<Node>()[0];
    Node n1 = AEdge.get<Node>()[1];
    
    math::Cross2D c0 =(*m_field)[n0.getID()];
    math::Cross2D c1 =(*m_field)[n1.getID()];
    
    math::Cross2D c = math::Cross2D::mean(c0, 1-param_intersection,
                                          c1, param_intersection  );
    
    std::vector<math::Vector> c_vectors = c.componentVectors();
    
    int ref_index = 0;
    double ref_dot_product = AINVec.dot(c_vectors[0]);
    for(int i=1; i<4;i++){
        double dot_product_i =  AINVec.dot(c_vectors[i]);
        if(dot_product_i>ref_dot_product){
            ref_dot_product = dot_product_i;
            ref_index = i;
        }
        
    }
    AOUTVec = c_vectors[ref_index];
    
    return true;
}

/*----------------------------------------------------------------------------*/
void Tools::findNextCell(const math::Point&  AFromPnt,
                         const math::Vector& AFromVec,
                         const int AFromCellDim,
                         const TCellID AFromCellID,
                         int& AToCellDim,
                         TCellID& AToCellID)
{
    if(AFromCellDim==0){
        Node from_node = m_mesh->get<Node>(AFromCellID);
        findNextCell(AFromPnt,AFromVec,from_node, AToCellDim,AToCellID);
    }
    else if(AFromCellDim==1){
        Edge from_edge = m_mesh->get<Edge>(AFromCellID);
        findNextCell(AFromPnt,AFromVec,from_edge, AToCellDim,AToCellID);
    }
    else
        throw GMDSException("findNextCell: we can only comes from a node or an edge");
}
/*----------------------------------------------------------------------------*/
void Tools::findNextCell(const math::Point&  AFromPnt,
                         const math::Vector& AFromVec,
                         const Node& AFromNode,
                         int& AToCellDim,
                         TCellID& AToCellID)
{
    std::vector<Face> adj_faces = AFromNode.get<Face>();
    
    //==============================================================
    // LOOK FOR A FACE
    //==============================================================
    bool find_next_cell = false;
    Face next_face;
    for(unsigned int i=0; i<adj_faces.size() && !find_next_cell; i++){
        Face current_face = adj_faces[i];
        
        if(isGoingInto(AFromPnt, AFromVec, AFromNode, current_face)) {
            next_face = current_face;
            find_next_cell = true;
        }
    } //for(unsigned int i=0; i<adj_faces.size(); i++){
    
    if(!find_next_cell)
        throw GMDSException("findNextCell: no right face for a node");
    
    //==============================================================
    // LOOK FOR AN EDGE
    //==============================================================
    // We have found a candidate face, now we check if edges are not
    // a better answer
    bool find_better_edge= false;
    Edge next_edge;
    std::vector<Edge> next_edges = next_face.get<Edge>();
    for(unsigned int i=0; i<next_edges.size() && !find_better_edge; i++){
        Edge current_edge = next_edges[i];
        if(isAlong(AFromVec, AFromNode, current_edge)) {
            next_edge = current_edge;
            find_better_edge = true;
        }
    }
    if(find_better_edge){
        AToCellDim = 1;
        AToCellID = next_edge.getID();
    }
    else if(find_next_cell){
        AToCellDim = 2;
        AToCellID = next_face.getID();
    }
}
/*----------------------------------------------------------------------------*/
void Tools::findNextCell(const math::Point&  AFromPnt,
                         const math::Vector& AFromVec,
                         const Edge& AFromEdge,
                         int& AToCellDim,
                         TCellID& AToCellID)
{
    std::vector<Face> adj_faces = AFromEdge.get<Face>();
    bool find_next_cell = false;
    
    for(unsigned int i=0; i<adj_faces.size() && !find_next_cell; i++){
        Face current_face = adj_faces[i];
        
        if(isGoingInto(AFromPnt, AFromVec, AFromEdge, current_face)) {
            find_next_cell = true;
            AToCellDim = 2;
            AToCellID = current_face.getID();
        }
    } //for(unsigned int i=0; i<adj_faces.size(); i++){
}
/*----------------------------------------------------------------------------*/
bool Tools::isGoingInto(const math::Point& APnt,
                        const math::Vector& AVec,
                        const Node& AFromNode,
                        const Face& AFace)
{
    // APnt is located in AFromNode
    
    // We take the nodes of AFace opposite to AFromNode
    std::vector<Node> current_nodes = AFace.get<Node>();
    std::vector<Node> other_nodes;
    for (unsigned int j = 0; j < current_nodes.size(); j++) {
        if (current_nodes[j].getID() != AFromNode.getID())
            other_nodes.push_back(current_nodes[j]);
    }
    if (other_nodes.size()!=2)
        throw GMDSException("isGoingInto(): Error, I don't find 2 opposite nodes");
    
    math::Ray from_ray(APnt, AVec);
    math::Segment opp_seg(other_nodes[0].getPoint(), other_nodes[1].getPoint());
    math::Point intersection_pnt;
    double intersection_param;
    return from_ray.intersect2D(opp_seg, intersection_pnt, intersection_param);
    
}
/*----------------------------------------------------------------------------*/
bool Tools::isGoingInto(const math::Point& APnt,
                        const math::Vector& AVec,
                        const Edge& AFromEdge,
                        const Face& AFace)
{
    // APnt is located on AFromEdge
    
    // We take the node of AFace that is not incident to AFromEdge
    std::vector<Node> current_nodes = AFace.get<Node>();
    std::vector<Node> from_nodes    = AFromEdge.get<Node>();
    Node opp_node;
    for (unsigned int j = 0; j < current_nodes.size(); j++)
    {
        if (current_nodes[j] != from_nodes[0] &&
            current_nodes[j] != from_nodes[1])
            opp_node = current_nodes[j];
    }
    if (opp_node.getID() == NullID)
        throw GMDSException("isGoingInto(): Error, no opposite node");
    
    math::Point  opp_pnt = opp_node.getPoint();
    math::Vector edge_vector (from_nodes[0].getPoint(), from_nodes[1].getPoint());
    math::Vector edge_witness(from_nodes[0].getPoint(), opp_pnt);
    math::Vector edge_ortho = edge_vector.cross(math::Vector(0,0,1));
    
    if(edge_ortho.dot(edge_witness)<0)
        edge_ortho=edge_ortho.opp();
    
    return (AVec.dot(edge_ortho) >= 0.0);
}
/*----------------------------------------------------------------------------*/
bool Tools::isAlong(const math::Vector& AVec,
                                        const Node& AFromNode,
                                        Edge& AEdge)
{
    std::vector<Node> n = AEdge.get<Node>();
    
    Node other_node;
    bool found_origin = false;
    for(unsigned int i=0;i<n.size();i++){
        Node current_node = n[i];
        if(current_node.getID()==AFromNode.getID())
            found_origin = true;
        else
            other_node = current_node;
    }
    if(!found_origin)
        return false;
    
    
    math::Vector v_edge(AFromNode.getPoint(),other_node.getPoint());
    v_edge.normalize();
    
    math::Vector v_dir = AVec;
    v_dir.normalize();
    
    return math::near(v_dir.dot(v_edge)-1, 0);
}
/*----------------------------------------------------------------------------*/


