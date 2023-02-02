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
#include <GMDS/Math/Triangle.h>
#include <GMDS/Math/Segment.h>
#include <GMDS/Math/Plane.h>
#include <GMDS/Math/Tetrahedron.h>
#include <GMDS/Math/Numerics.h>
#include <GMDS/Math/Ray.h>
#include <GMDS/Utils/Log.h>
/*---------------------------------------------------------------------------*/
// STL File Headers
#include <random>
#include <set>
#include <algorithm>
/*---------------------------------------------------------------------------*/
// FRAME File Headers
#include "TriangularSurfaceManipulator.h"
/*---------------------------------------------------------------------------*/
using namespace gmds;
#define DEBUG_GMDS
/*---------------------------------------------------------------------------*/
/** Local node numbering so that the triangle formed with:
*   - local node  i
*   - m_local_tri_node2edge[i][0]
*   - m_local_tri_node2edge[i][1]
*   has the same orientation as the original triangle for any vertex i.
*/
const int TriangularSurfaceManipulator::m_local_tri_node2edge[3][2] = {
    {1, 2},
    {2, 0},
    {0, 1}
};
/*---------------------------------------------------------------------------*/
TriangularSurfaceManipulator::TriangularSurfaceManipulator(IGMesh* AMesh)
:m_mesh(AMesh)
{
    GEO::PCK::initialize();
}
/*---------------------------------------------------------------------------*/
TriangularSurfaceManipulator::~TriangularSurfaceManipulator()
{
    GEO::PCK::terminate();
}
/*---------------------------------------------------------------------------*/
std::vector<Node> TriangularSurfaceManipulator::
insertFreePoints(std::vector<math::Point>& APntToInsert,
                 std::vector<TCellID>& AInFaces)
{
    
    // for each face, we keep in mind the index of the points inserted
    // inside it
    std::map<TCellID,std::vector<int> >  faces_to_new_nodes;

    std::vector<Face> faces;
    faces.reserve(AInFaces.size());
    for(auto id:AInFaces){
        faces.push_back(m_mesh->get<Face>(id));
    }
    std::vector<Node> new_nodes;
    bool insert_one=false;
    for(auto i=0; i<APntToInsert.size(); i++){
        try{
            math::Point p = APntToInsert[i];
            Node new_node;
            std::vector<Face> fs = findFaces(p,faces, new_node);
            if(fs.size()<2)
                throw GMDSException("A free point must be on an edge or a vertex");
             insert_one = true;
            std::cout<<"Free point:"<<p<<std::endl;
            new_nodes.push_back(new_node);
            for(auto f:fs){
                std::cout<<"\t put in "<<f.getID()<<std::endl;
                faces_to_new_nodes[f.getID()].push_back(new_nodes.size()-1);
            }

        }
        catch(GMDSException& e){
            //We jump this point so
            continue;
        }
        
    }
    if(!insert_one){
        std::cout<<"No point to insert"<<std::endl;
        return new_nodes;
    }
    
    
    std::cout<<"------- FREE POINT INSERTION --------"<<std::endl;
    std::vector<TCellID> faces_to_delete;
    faces_to_delete.reserve(faces_to_new_nodes.size());
    std::cout<<"here"<<std::endl;
    for(auto d:faces_to_new_nodes){
        
        Face f = m_mesh->get<Face>(d.first);
        std::vector<int> to_insert =d.second;
        std::vector<Node> ns;
        ns.resize(to_insert.size());
        std::cout<<"face "<<f.getID()<<std::endl;
        for(auto i=0;i<to_insert.size();i++){
            std::cout<<" insert "<<to_insert[i]<<std::endl;
            ns[i]=new_nodes[to_insert[i]];
        }
        
        if(split(f,ns))
            faces_to_delete.push_back(f.getID());
    }
    std::cout<<"here 2"<<std::endl;
    for(auto f_id:faces_to_delete){
        Face f = m_mesh->get<Face>(f_id);
        m_mesh->deleteFace(f);
    }
     std::cout<<"here 3"<<std::endl;
    return new_nodes;
}

/*---------------------------------------------------------------------------*/
void TriangularSurfaceManipulator::
insertLoop(std::vector<math::Point>& ALoop,
           std::vector<int>& AClassification,
           std::vector<int>& AGeomID,
           std::vector<gmds::TCellID>& AEnclosedFaces,
           std::vector<gmds::TCellID>& ALoopNodes,
           std::vector<bool>& ALoopNodesFromPnt)
{
    
    //we get the points forming the loop. Those points will be enriched by
    // others. Nothe that contrary to ALoop, points will not defined a cycle
    std::vector<Node> new_nodes;
    std::vector<bool> new_from_pnt;
    // for each face, we keep in mind the index of the points inserted
    // inside it
    std::map<TCellID,std::vector<int> >  datas;

    //=======================================================================
    // (1) We detect the triangle containing the first point
    //=======================================================================
    //The detection of the first face is not optimized
    std::cout<<"=========FIND FACE==========="<<std::endl;
    Node first_node;
    std::vector<Face> first_faces = findFaces(ALoop[0],first_node);
    
    // we keep in mind the point and face association
    new_nodes.push_back(first_node);
    new_from_pnt.push_back(true); //it is an initial node
   
    
    for(auto f:first_faces){
        datas[f.getID()].push_back(0);
    }
    //=======================================================================
    // (2) We goes throug the complete loop to find all the intersected faces
    //     and computing intersection points that will be stored in points
    //=======================================================================
    std::vector<Face>  current_faces=first_faces, next_faces;
    for(auto i=0; i<ALoop.size(); i++){
        int current = i;
        int next = (i+1)%ALoop.size();
        int cur_cl  = AClassification[current];
        int next_cl = AClassification[next];
        bool is_last= (i==ALoop.size()-1);
        
        bool on_same_curve = false;
        int curve_id=-1;
        if(cur_cl<2 && next_cl<2){
            //on curve or point
            if(cur_cl==1 && next_cl==1){
                //both on curve
                if(AGeomID[current]==AGeomID[next]){
                    //same curve
                    on_same_curve=true;
                    curve_id=AGeomID[current];
                }
            }
            //We cannot conclude with current data
            else if(cur_cl==1 && next_cl==0){
                on_same_curve=true;
            }
            else if(cur_cl==0 && next_cl==1){
                on_same_curve=true;
            }
        }
        if(on_same_curve){
            walkOnCurve(ALoop[current],
                        ALoop[next],
                        cur_cl,
                        next_cl,
                        AGeomID[current],
                        AGeomID[next],
                        current_faces,
                        new_nodes,
                        new_from_pnt,
                        datas,
                        next_faces,
                        is_last);
            
        }
        else{
            walk(ALoop[current],
                 ALoop[next],
                 current_faces,
                 new_nodes,
                 new_from_pnt,
                 datas,
                 next_faces,
                 is_last);
        }
        //the last inserted point correspond to next. So we get one the
        // relevant faces taking the first face of the last entry
        current_faces = next_faces;
    }
    std::vector<Node> clean_nodes;
    std::vector<bool> clean_from_pnt;
    clean_nodes.reserve(new_nodes.size());
    clean_from_pnt.reserve(new_from_pnt.size());

    int same_id=0, same_geo=0;
    for(auto i=0;i<new_nodes.size();i++){
        Node ni = new_nodes[i];
        Node nj = new_nodes[(i+1)%new_nodes.size()];
        if(ni.getID()==nj.getID()){
            same_id++;
            continue;
        }
        if(ni.getPoint().distance2(nj.getPoint())<1e-5){
            std::cout<<"TOO CLOSE: "<<ni.getID()<<", "<<nj.getID()<<std::endl;
            same_geo++;
            continue;
        
        }
        clean_nodes.push_back(new_nodes[i]);
        clean_from_pnt.push_back(new_from_pnt[i]);
    }
    std::vector<math::Point> new_points;
    new_points.resize(clean_nodes.size());
    for(auto i=0;i<clean_nodes.size();i++){
        new_points[i] = clean_nodes[i].getPoint();
    }
    writePoints(new_points);
    VTKWriter<IGMesh> bnd_writer(*m_mesh);
    bnd_writer.write("BEFORE_SPLIT",F|N);

    //=======================================================================
    // (3) We imprint the loop into the mesh
    //=======================================================================
    // We split each encountred triangle using the list of nodes we created
    // for each face, we keep in mind the index of the points inserted
    // inside it
    std::vector<TCellID> faces_to_delete;
    faces_to_delete.reserve(datas.size());
    for(auto d:datas){
        Face f = m_mesh->get<Face>(d.first);
        std::vector<int> to_insert =d.second;
        std::vector<Node> ns;
        ns.resize(to_insert.size());
        for(auto i=0;i<to_insert.size();i++){
            ns[i]=new_nodes[to_insert[i]];
        }
        
        if(split(f,ns))
            faces_to_delete.push_back(f.getID());
    }
    for(auto f_id:faces_to_delete){
        Face f = m_mesh->get<Face>(f_id);
        m_mesh->deleteFace(f);
    }
    bnd_writer.write("BEFORE_BALANCE",F|N);

    //=======================================================================
    // (4) Edge reconstruction. Missing edges are rebuilt using edge balance
    //=======================================================================
    for(auto i=0; i<clean_nodes.size();i++){
        Node ni = clean_nodes[i];
        Node nj = clean_nodes[(i+1)%clean_nodes.size()];
        if(ni.getID()==nj.getID())
            continue;
        //build the edge if it does not exist
      buildEdge(ni,nj);
    }
    
    //=======================================================================
    // (5) We store the nodes of the loop
    //=======================================================================
    ALoopNodes.reserve(clean_nodes.size());
    for(auto n:clean_nodes){
        ALoopNodes.push_back(n.getID());
    }
    ALoopNodesFromPnt.reserve(clean_from_pnt.size());
    for(auto b:clean_from_pnt){
        ALoopNodesFromPnt.push_back(b);
    }
    //=======================================================================
    // (5) And extract the inner faces
    //=======================================================================
    AEnclosedFaces = extractEnclosedFaces(clean_nodes);
}
/*---------------------------------------------------------------------------*/
std::vector<TCellID> TriangularSurfaceManipulator::
extractEnclosedFaces(const std::vector<TCellID>& ALoopNodes)
{
    std::vector<TCellID> enclosed_faces;
    
    //=======================================================================
    // (1) We build a vector of oriented edges that are used to advance into
    // the area to extract
    //=======================================================================
    std::vector<std::pair<TCellID,TCellID> > edge_loop, init_loop;
    edge_loop.reserve(ALoopNodes.size());
    auto nb_nodes = ALoopNodes.size();
    for(auto i=0; i<nb_nodes; i++){
        TCellID id1 = ALoopNodes[i];
        TCellID id2 = ALoopNodes[(i+1)%nb_nodes];
        if(id1==id2)
            continue;
        edge_loop.push_back(std::pair<TCellID,TCellID>(id1,id2));
    }
    init_loop = edge_loop;
    while(!edge_loop.empty()){
        //We pick an edge of the loop
        std::pair<TCellID,TCellID> e = edge_loop.back();
        edge_loop.pop_back();
        Node from = m_mesh->get<Node>(e.first);
        Node to   = m_mesh->get<Node>(e.second);
        
        std::vector<Face> fs = getFaces(from,to);
        //among the two faces sharing e, one is oriented in the same direction
        // as e. We pick it
        Face in_face;
        for(auto f:fs){
            std::vector<TCellID> fn_ids = f.getIDs<Node>();
            if((fn_ids[0]==e.first && fn_ids[1]==e.second) ||
               (fn_ids[1]==e.first && fn_ids[2]==e.second) ||
               (fn_ids[2]==e.first && fn_ids[0]==e.second) ){
                in_face=f;
            }
        }
        
        if(std::find(enclosed_faces.begin(),enclosed_faces.end(),
                     in_face.getID())== enclosed_faces.end()){
            //means this face is not yet added
//            std::cout<<" ADD FACE= "<<in_face.getID()<<" from "<<
//            e.first <<" "<< e.second<<std::endl;
            enclosed_faces.push_back(in_face.getID());
            std::vector<TCellID> n_ids = in_face.getIDs<Node>();
            for(auto i=0; i<3; i++){
                TCellID id1 = n_ids[i];
                TCellID id2 = n_ids[(i+1)%3];
                if(id1==e.first && id2==e.second)
                    continue;
                std::pair<TCellID,TCellID> e1(id2,id1);
                if(std::find(init_loop.begin(),init_loop.end(),
                             e1)!= init_loop.end())
                    continue;
                std::pair<TCellID,TCellID> e2(id1,id2);
                if(std::find(init_loop.begin(),init_loop.end(),
                             e2)!= init_loop.end())
                    continue;
                
                //other edge, we add it
                edge_loop.push_back(e1);
            }
            
        }
    }
    return enclosed_faces;
}

/*---------------------------------------------------------------------------*/
std::vector<TCellID> TriangularSurfaceManipulator::
extractEnclosedFaces(const std::vector<Node>& ALoopNodes)
{
    std::vector<TCellID> loop_ids;
    loop_ids.reserve(ALoopNodes.size());
    for(auto n:ALoopNodes){
        loop_ids.push_back(n.getID());
    }
    return extractEnclosedFaces(loop_ids);
    
}
/*---------------------------------------------------------------------------*/
void TriangularSurfaceManipulator::buildEdge(Node& ANI, Node& ANJ)
{
    //this function is called for triangles obtained after splitting a single
    //triangle. As a consequence, we assume they are almost coplanar
    Variable<math::Vector>* var_normal =
    m_mesh->getVariable<math::Vector>(GMDS_FACE,"NORMAL");
    std::vector<Face> fij = getFaces(ANI,ANJ);
    
    //the edge is already here
    if(fij.size()==2)
        return;
    
    std::cout<<"Build edge "<<ANI.getID()<<" - "<<ANI.getPoint()<<std::endl
             <<"           "<<ANJ.getID()<<" - "<<ANJ.getPoint()<<std::endl;
    //Among all the faces of ANI, we look for the one sij is going through
    Face current_face = getOutFace(ANI, ANJ.getPoint());
    
    std::vector<Node> current_nodes = current_face.get<Node>();
    //the way opp_nodes is filled in is important to keep the right orientation
    Node opp_nodes[2];
    if(current_nodes[0]==ANI){
        opp_nodes[0] = current_nodes[1];
        opp_nodes[1] = current_nodes[2];
    }
    else if(current_nodes[1]==ANI){
        opp_nodes[0] = current_nodes[2];
        opp_nodes[1] = current_nodes[0];
    }
    else{//[2]
        opp_nodes[0] = current_nodes[0];
        opp_nodes[1] = current_nodes[1];
    }
    
    bool terminate=false;
    Node current_node = ANI;
    while(!terminate){
        
        std::vector<Face> f_opp = getFaces(opp_nodes[0], opp_nodes[1]);
        Face next_face = (f_opp[0]==current_face)?f_opp[1]:f_opp[0];
        std::vector<Node> next_nodes = next_face.get<Node>();
        Node next_node;
        for(auto n:next_nodes){
            if(n!=opp_nodes[0] && n!=opp_nodes[1]){
                next_node = n;
            }
        }
        std::cout<<"-------------------"<<std::endl;
        std::cout<<"Current node: "<<current_node.getID()<<std::endl;
        std::cout<<"Edge node 0 : "<<opp_nodes[0].getID()<<std::endl;
        std::cout<<"Edge node 1 : "<<opp_nodes[1].getID()<<std::endl;
        std::cout<<"Next    node: "<<next_node.getID()<<std::endl;
        std::cout<<"Current face: "<<current_face.getID()<<std::endl;
        std::cout<<"Next    face: "<<next_face.getID()<<std::endl;
        std::cout<<"Edge of 0,1: ";
        for(auto k:f_opp){
            std::cout<<k.getID()<<" ";
        }
        std::cout<<std::endl;
        //we balance the edge
        //NEW FACE 1
        std::vector<Node> n;
        n.resize(3);
        n[0]=current_node;
        n[1]=opp_nodes[0];
        n[2]=next_node;
        current_face.set(n);
        //NEW FACE 2
        n[1]=next_node;
        n[2]=opp_nodes[1];
        next_face.set(n);
        
        
        {
            math::Vector v1(current_face.get<Node>()[0].getPoint(),
                            current_face.get<Node>()[1].getPoint());
            math::Vector v2(current_face.get<Node>()[0].getPoint(),
                            current_face.get<Node>()[2].getPoint());
            math::Vector vn=v1.cross(v2);
            vn.normalize();
            (*var_normal)[current_face.getID()]=vn;
        }

        
        {
            math::Vector v1(next_face.get<Node>()[0].getPoint(),
                            next_face.get<Node>()[1].getPoint());
            math::Vector v2(next_face.get<Node>()[0].getPoint(),
                            next_face.get<Node>()[2].getPoint());
            math::Vector vn=v1.cross(v2);
            vn.normalize();
            (*var_normal)[next_face.getID()]=vn;
        }

        //Update N->F
        current_node.add(next_face);
        next_node.add(current_face);
        opp_nodes[0].remove(next_face);
        opp_nodes[1].remove(current_face);

        if(next_node==ANJ){
            terminate=true;
        }
        else{
            //We need to find the next triangle to go through
            math::Plane cur_plane(current_node.getPoint(),
                                  opp_nodes[0].getPoint(),
                                  opp_nodes[1].getPoint());
            math::Point next_pnt = cur_plane.project(next_node.getPoint());
            math::Segment s[2] = {
                math::Segment(opp_nodes[0].getPoint(),next_pnt),
                math::Segment(opp_nodes[1].getPoint(),next_pnt)
            };
            math::Ray r(ANI.getPoint(),ANJ.getPoint());
            math::Point intersection_pnt;
            double param_seg=0, param_ray=0;
            int index=-1;
            for(auto i=0; i<2; i++){
                if(r.intersect3D(s[i], intersection_pnt,
                                 param_seg, param_ray)){
                    if(param_ray>0){
                        index=i;
                    }
                }
            }//for(auto i=0;...
            if(index==0){
                //current node remains ANI
                //opp_nodes[0] remains the same too
                opp_nodes[1] = next_node;
                //current_face remains the same
            }
            else{
                //current node remains ANI
                opp_nodes[0] = next_node;
                //opp_nodes[1] remains the same too
                current_face=next_face;
                
            }
        }
    }
}
/*---------------------------------------------------------------------------*/
bool TriangularSurfaceManipulator::split(Face& AF,
                                         std::vector<Node>& AN)
{

    Variable<int>* surf_color = m_mesh->getVariable<int>(GMDS_FACE, "BND_SURFACE_COLOR");

    Variable<math::Vector>* var_normal =
    m_mesh->getVariable<math::Vector>(GMDS_FACE,"NORMAL");

    std::vector<Face> faces(1,AF);
    
    std::vector<Node> f_nodes = AF.get<Node>();
    std::cout<<"=================================="<<std::endl;
    std::cout<<"SPLIT FACE "<<AF.getID()<<": ";
    std::cout<<f_nodes[0].getID()<<" ";
    std::cout<<f_nodes[1].getID()<<" ";
    std::cout<<f_nodes[2].getID()<<std::endl;
    std::cout<<" with nodes= "<<std::endl;;
    for(auto n:AN){
        std::cout<<"\t"<<n.getID()<<" "<<n.getPoint()<<std::endl;;
    }
    
    std::vector<TCellID> nodes_on_edge[3];
    std::vector<TCellID> nodes_inside;

    math::Point p[3] ={
        f_nodes[0].getPoint(),
        f_nodes[1].getPoint(),
        f_nodes[2].getPoint()
    };
    
    math::Segment s01(p[0],p[1]);
    math::Segment s12(p[1],p[2]);
    math::Segment s20(p[2],p[0]);
    
    //=======================================================================
    // (1) Nodes are dispatched along if they belon to an edge of AF or not
    //=======================================================================
    for(auto n:AN){
        std::cout<<"Node "<<n.getID();
        if(n==f_nodes[0] || n==f_nodes[1] || n==f_nodes[2])
            continue;
        
        math::Point pnt = n.getPoint();
        bool on_an_edge =false;
        std::vector<TCoord> coords(3,0.0);
        bool on_edge[3] = {false,false,false};

    //    if(isIn(pnt,AF, on_edge[0], on_edge[1], on_edge[2])){
            if(isIn(pnt,AF, coords)){

            
            for(auto i=0; i<3;i++){
                if(coords[i]<0.0001)
                    on_edge[i]=true;
            }
            if((on_edge[0] && on_edge[1]) ||
               (on_edge[1] && on_edge[2]) ||
               (on_edge[2] && on_edge[0])){
                //means that the point is very close to one triangle
                //vertex, but not under the tolerance
                //we put the point on the edge it is the closest
                int keep_on = 0;
                double c = coords[0];
                for(auto i=1;i<3;i++){
                    if(coords[i]<c){
                        c=coords[i];
                        keep_on=i;
                    }
                }
                on_an_edge=true;
                on_edge[keep_on]=true;
                on_edge[(keep_on+1)%3]=false;
                on_edge[(keep_on+2)%3]=false;
            }
            
            if(on_edge[0]){
                //should be on edge 12
                std::cout<<" on edge 1"<<std::endl;
                nodes_on_edge[1].push_back(n.getID());
                on_an_edge=true;
            }
            else if(on_edge[1]){
                //should be on edge 02
                std::cout<<" on edge 2"<<std::endl;
                nodes_on_edge[2].push_back(n.getID());
                on_an_edge=true;
            }
            else if(on_edge[2]){
                //should be on edge 01
                std::cout<<" on edge 0"<<std::endl;
                nodes_on_edge[0].push_back(n.getID());
                on_an_edge=true;

            }
            if(!on_an_edge){
                std::cout<<" inside"<<std::endl;
                nodes_inside.push_back(n.getID());
            }
            
            
        }
    }
    std::cout<<"INSIDE: ";
    for(auto x:nodes_inside){
        std::cout<<x<<" ";
    }
    std::cout<<"\n Edge 0: ";
    for(auto x:nodes_on_edge[0]){
        std::cout<<x<<" ";
    }
    
    std::cout<<"\n Edge 1: ";
    for(auto x:nodes_on_edge[1]){
        std::cout<<x<<" ";
    }
    
    std::cout<<"\n Edge 2: ";
    for(auto x:nodes_on_edge[2]){
        std::cout<<x<<" ";
    }
    std::cout<<std::endl;
    
    //=======================================================================
    // (2) Nodes on edges are inserted first
    //=======================================================================
    bool split_done =false;
    for(auto i=0; i<3; i++){
        if(nodes_on_edge[i].empty())
            continue;
        
        split_done = true;
        //extremities nodes
        Node ni = f_nodes[i];
        Node nj = f_nodes[(i+1)%3];
        math::Point pi = ni.getPoint();
        // We order the points to insert from ni to nj but points on ni or
        // nj are discarded
        std::map<double,TCellID> distance_to_ni;
        for(auto n_id:nodes_on_edge[i]){
            if(n_id==ni.getID() || n_id==nj.getID())
                continue;
            
            math::Point p = m_mesh->get<Node>(n_id).getPoint();
            distance_to_ni[pi.distance2(p)]=n_id;
        }
        //We sort nodes now
  //      std::sort(distance_to_ni.begin(),distance_to_ni.end());
        //and build the final series of nodes
        std::vector<Node> ordered_nodes;
        ordered_nodes.push_back(ni);
        for(auto in_node: distance_to_ni){
            ordered_nodes.push_back(m_mesh->get<Node>(in_node.second));
        }
        ordered_nodes.push_back(nj);
        
        // And now we split the face sharing edge ni,nj along this series
        // of new nodes.
        Face f_ij;
        int f_ij_index;
        for(auto i_f=0; i_f<faces.size(); i_f++){
            std::vector<TCellID> f_ids= faces[i_f].getIDs<Node>();
            if((f_ids[0]==ni.getID() && f_ids[1]==nj.getID()) ||
               (f_ids[1]==ni.getID() && f_ids[0]==nj.getID()) ||
               (f_ids[1]==ni.getID() && f_ids[2]==nj.getID()) ||
               (f_ids[2]==ni.getID() && f_ids[1]==nj.getID()) ||
               (f_ids[0]==ni.getID() && f_ids[2]==nj.getID()) ||
               (f_ids[2]==ni.getID() && f_ids[0]==nj.getID()) ){
                f_ij = faces[i_f];
                f_ij_index=i_f;
            }
        }
        
        std::vector<Node> f_ij_nodes = f_ij.get<Node>();
        Node nk;
        for(auto n:f_ij_nodes){
            if(n!=ni && n!=nj){
                nk=n;
            }
        }
        //We split the triangle
        for(auto i_n=0; i_n<ordered_nodes.size()-1;i_n++){
            Node n1 = ordered_nodes[i_n];
            Node n2 = ordered_nodes[i_n+1];
            Face new_f = m_mesh->newTriangle(n1,n2,nk);
          
            math::Vector v1(new_f.get<Node>()[0].getPoint(),
                            new_f.get<Node>()[1].getPoint());
            math::Vector v2(new_f.get<Node>()[0].getPoint(),
                            new_f.get<Node>()[2].getPoint());
            math::Vector vn=v1.cross(v2);
            vn.normalize();
            (*var_normal)[new_f.getID()]=vn;

            (*surf_color)[new_f.getID()] = (*surf_color)[AF.getID()];
            std::cout<<"\t new face "<<new_f.getID()<<": "<<n1.getID()<<" "<<n2.getID()
            <<" "<<nk.getID()<<std::endl;
            n1.add(new_f);
            n2.add(new_f);
            nk.add(new_f);
            faces.push_back(new_f);
        }
        ni.remove(f_ij);
        nj.remove(f_ij);
        nk.remove(f_ij);
        
        //the splitted face is removed now from the set of candidates
        //but only if different from AF because of the hole-filling policy
        //in GMDS containers
        if(f_ij!=AF){
            std::cout<<"Remove face "<<f_ij.getID()<<std::endl;
            faces[f_ij_index]= faces.back();
            faces.pop_back();
            m_mesh->deleteFace(f_ij);
        }
    }
    //=======================================================================
    // (3) Inner nodes
    //=======================================================================

    
    for(auto n_id:nodes_inside){
        Node n = m_mesh->get<Node>(n_id);
        math::Point p = n.getPoint();
        std::cout<<"-------------------------"<<std::endl;
        std::cout<<" INNER NODE "<<n.getID()<<": "<<p<<std::endl;
        bool found_f = false;
        Face f;
        int f_index =-1;
        bool on_edge[3] = {false,false,false};
        if(faces.size()==1){
            f= faces[0];
            f_index=0;
            found_f=true;
            std::cout<<"\t one face only"<<std::endl;
        }
        else{
            for(int i_f=0;i_f<faces.size() && !found_f; i_f++){
                Face cur = faces[i_f];
                if(cur==AF)
                    continue;
                std::vector<TCoord> coords(3,0.0);
                if(isIn(p,cur, coords)){
                    found_f=true;
//                    for(auto i=0; i<3;i++){
//                        if(coords[i]<0)
//                            found_f=false;
//                    }
                    if(found_f){
                        
                        f_index = i_f;
                        f=cur;
                        for(auto i=0; i<3;i++){
                            if(coords[i]<0.0001)
                                on_edge[i]=true;
                        }
                    }
                }
            }
        }
        if(!found_f){
            throw GMDSException("Error in splitting");
        }
        std::cout<<"\t in face "<<f.getID()<<std::endl;
        std::cout<<"\t with on_edge = "<<on_edge[0]<<" "<<on_edge[1]
        <<" "<<on_edge[2]<<std::endl;
        std::vector<Node> f_nodes = f.get<Node>();
        Node common0, common1, other;
        bool on_bnd=false;
        //we have the face to split
        if((on_edge[0] && on_edge[1]) ||
           (on_edge[1] && on_edge[2]) ||
           (on_edge[2] && on_edge[0])){
            //nothing to do, the point is already in the mesh
            ;
        }
        else if(on_edge[0]){
            //means the point must be inserted on edge [1,2]
            common0 = f_nodes[1];
            common1 = f_nodes[2];
            other   = f_nodes[0];
            on_bnd=true;
        }
        else if(on_edge[1]){
            //means the point must be inserted on edge [0,2]
            common0 = f_nodes[2];
            common1 = f_nodes[0];
            other   = f_nodes[1];
            on_bnd=true;
        }
        else if(on_edge[2]){
            //means the point must be inserted on edge [0,1]
            common0 = f_nodes[0];
            common1 = f_nodes[1];
            other   = f_nodes[2];
            on_bnd=true;
        }
        if(on_bnd){
            std::vector<Face> fs = getFaces(common0,common1);
            std::cout<<"On edge and nb faces = "<<fs.size()<<std::endl;
            if(fs.size()==1){
                throw GMDSException("IMPOSSIBLE HERE");
//                //split one triangle
//                Face cur = fs[0];
//                std::vector<Node> cur_nodes = cur.get<Node>();
//                Face t0 = m_mesh->newTriangle(cur_nodes[0], cur_nodes[1],n);
//                Face t1 = m_mesh->newTriangle(cur_nodes[1], cur_nodes[2],n);
//                Face t2 = m_mesh->newTriangle(cur_nodes[2], cur_nodes[0],n);
//                (*surf_color)[t0.getID()] = (*surf_color)[AF.getID()];
//                (*surf_color)[t1.getID()] = (*surf_color)[AF.getID()];
//                (*surf_color)[t2.getID()] = (*surf_color)[AF.getID()];
//
//                cur_nodes[0].add(t0);
//                cur_nodes[1].add(t0);
//                n.add(t0);
//                cur_nodes[1].add(t1);
//                cur_nodes[2].add(t1);
//                n.add(t1);
//                cur_nodes[2].add(t2);
//                cur_nodes[0].add(t2);
//                n.add(t2);
//                
//                cur_nodes[0].remove(f);
//                cur_nodes[1].remove(f);
//                cur_nodes[2].remove(f);
//                
//                {
//                    math::Vector v1(t0.get<Node>()[0].getPoint(),
//                                    t0.get<Node>()[1].getPoint());
//                    math::Vector v2(t0.get<Node>()[0].getPoint(),
//                                    t0.get<Node>()[2].getPoint());
//                    math::Vector vn=v1.cross(v2);
//                    vn.normalize();
//                    (*var_normal)[t0.getID()]=vn;
//                }
//                
//                {
//                    math::Vector v1(t1.get<Node>()[0].getPoint(),
//                                    t1.get<Node>()[1].getPoint());
//                    math::Vector v2(t1.get<Node>()[0].getPoint(),
//                                    t1.get<Node>()[2].getPoint());
//                    math::Vector vn=v1.cross(v2);
//                    vn.normalize();
//                    (*var_normal)[t1.getID()]=vn;
//                }
//                
//                {
//                    math::Vector v1(t2.get<Node>()[0].getPoint(),
//                                    t2.get<Node>()[1].getPoint());
//                    math::Vector v2(t2.get<Node>()[0].getPoint(),
//                                    t2.get<Node>()[2].getPoint());
//                    math::Vector vn=v1.cross(v2);
//                    vn.normalize();
//                    (*var_normal)[t2.getID()]=vn;
//                }
            }
            else {
                //We have the two common nodes, but we need to preserve
                //the orientation. fs[0] is our reference here
                Node n1 = other;
                Node n0 = common0;
                Node n2 = common1;
                //n0, n1 and n2 are ordered in the right way so
                
                Face cur = (f==fs[0])?fs[1]:fs[0];
                std::vector<Node> cur_n = cur.get<Node>();
                Node n3;
                for(auto n_f:cur_n){
                    if(n_f.getID()!=n0.getID() &&
                       n_f.getID()!=n2.getID() )
                        n3=n_f;
                }
                
                std::cout<<"\n so: "<<n0.getID()<<" "<<
                n1.getID()<<" "<<n2.getID()<<" "<<n3.getID()<<std::endl;
                Face t0 = m_mesh->newTriangle(n1,n0,n);
                Face t1 = m_mesh->newTriangle(n2,n1,n);
                Face t2 = m_mesh->newTriangle(n3,n2,n);
                Face t3 = m_mesh->newTriangle(n0,n3,n);
                (*surf_color)[t0.getID()] = (*surf_color)[AF.getID()];
                (*surf_color)[t1.getID()] = (*surf_color)[AF.getID()];
                (*surf_color)[t2.getID()] = (*surf_color)[AF.getID()];
                (*surf_color)[t3.getID()] = (*surf_color)[AF.getID()];

                n0.add(t0);
                n1.add(t0);
                n.add(t0);
                n1.add(t1);
                n2.add(t1);
                n.add(t1);
                n2.add(t2);
                n3.add(t2);
                n.add(t2);
                n3.add(t3);
                n0.add(t3);
                n.add(t3);
                n0.remove(fs[0]);
                n0.remove(fs[1]);
                n1.remove(fs[0]);

                n2.remove(fs[0]);
                n2.remove(fs[1]);
                n3.remove(fs[1]);

                faces.push_back(t0);
                faces.push_back(t1);
                faces.push_back(t2);
                faces.push_back(t3);
                
                {
                    math::Vector v1(t0.get<Node>()[0].getPoint(),
                                    t0.get<Node>()[1].getPoint());
                    math::Vector v2(t0.get<Node>()[0].getPoint(),
                                    t0.get<Node>()[2].getPoint());
                    math::Vector vn=v1.cross(v2);
                    vn.normalize();
                    (*var_normal)[t0.getID()]=vn;
                }
                
                {
                    math::Vector v1(t1.get<Node>()[0].getPoint(),
                                    t1.get<Node>()[1].getPoint());
                    math::Vector v2(t1.get<Node>()[0].getPoint(),
                                    t1.get<Node>()[2].getPoint());
                    math::Vector vn=v1.cross(v2);
                    vn.normalize();
                    (*var_normal)[t1.getID()]=vn;
                }
                
                {
                    math::Vector v1(t2.get<Node>()[0].getPoint(),
                                    t2.get<Node>()[1].getPoint());
                    math::Vector v2(t2.get<Node>()[0].getPoint(),
                                    t2.get<Node>()[2].getPoint());
                    math::Vector vn=v1.cross(v2);
                    vn.normalize();
                    (*var_normal)[t2.getID()]=vn;
                }
                
                {
                    math::Vector v1(t3.get<Node>()[0].getPoint(),
                                    t3.get<Node>()[1].getPoint());
                    math::Vector v2(t3.get<Node>()[0].getPoint(),
                                    t3.get<Node>()[2].getPoint());
                    math::Vector vn=v1.cross(v2);
                    vn.normalize();
                    (*var_normal)[t3.getID()]=vn;
                }

            }
        }
        else{
            //inside f
            Face t0 = m_mesh->newTriangle(f_nodes[0], f_nodes[1],n);
            Face t1 = m_mesh->newTriangle(f_nodes[1], f_nodes[2],n);
            Face t2 = m_mesh->newTriangle(f_nodes[2], f_nodes[0],n);
            (*surf_color)[t0.getID()] = (*surf_color)[AF.getID()];
            (*surf_color)[t1.getID()] = (*surf_color)[AF.getID()];
            (*surf_color)[t2.getID()] = (*surf_color)[AF.getID()];

            f_nodes[0].add(t0);
            f_nodes[1].add(t0);
            n.add(t0);
            f_nodes[1].add(t1);
            f_nodes[2].add(t1);
            n.add(t1);
            f_nodes[2].add(t2);
            f_nodes[0].add(t2);
            n.add(t2);
            
            f_nodes[0].remove(f);
            f_nodes[1].remove(f);
            f_nodes[2].remove(f);
            faces.push_back(t0);
            faces.push_back(t1);
            faces.push_back(t2);
            
            {
                math::Vector v1(t0.get<Node>()[0].getPoint(),
                                t0.get<Node>()[1].getPoint());
                math::Vector v2(t0.get<Node>()[0].getPoint(),
                                t0.get<Node>()[2].getPoint());
                math::Vector vn=v1.cross(v2);
                vn.normalize();
                (*var_normal)[t0.getID()]=vn;
            }
            
            {
                math::Vector v1(t1.get<Node>()[0].getPoint(),
                                t1.get<Node>()[1].getPoint());
                math::Vector v2(t1.get<Node>()[0].getPoint(),
                                t1.get<Node>()[2].getPoint());
                math::Vector vn=v1.cross(v2);
                vn.normalize();
                (*var_normal)[t1.getID()]=vn;
            }
            
            {
                math::Vector v1(t2.get<Node>()[0].getPoint(),
                                t2.get<Node>()[1].getPoint());
                math::Vector v2(t2.get<Node>()[0].getPoint(),
                                t2.get<Node>()[2].getPoint());
                math::Vector vn=v1.cross(v2);
                vn.normalize();
                (*var_normal)[t2.getID()]=vn;
            }

        }
        
        //the splitted face is removed now from the set of candidates
        //but only if different from AF because of the hole-filling policy
        //in GMDS containers
        if(f!=AF){
            faces[f_index]= faces.back();
            faces.pop_back();
            m_mesh->deleteFace(f);
        }

    }
    return split_done;
}
/*---------------------------------------------------------------------------*/
gmds::Face TriangularSurfaceManipulator::
getOutFace(const gmds::Node& ANI, const gmds::math::Point ATo)
{
    std::vector<Face> faces = ANI.get<Face>();

    math::Point from_pnt = ANI.getPoint();
    math::Point to_pnt = ATo;
    math::Vector direction(from_pnt, to_pnt);
    
    math::Point intersection_pnt;
    double param_seg=0;
    double param_ray=0;
    for (auto i=0; i<faces.size(); i++){
        Face fi = faces[i];
        
        
        math::Vector normal = fi.normal();
        //means we are on a wrong face
        if(std::abs(normal.dot(direction))>0.5)
            continue;
        
        std::vector<Node> nodes = fi.get<Node>();
        math::Plane pl(nodes[0].getPoint(),
                       nodes[1].getPoint(),
                       nodes[2].getPoint());
        
        math::Vector local_direction(from_pnt, pl.project(ATo));
        math::Ray r(from_pnt,local_direction);
        
        math::Segment opposite_seg;
        if(nodes[0]==ANI){
            opposite_seg=math::Segment(nodes[1].getPoint(),
                                       nodes[2].getPoint());
        }
        else if(nodes[1]==ANI){
            opposite_seg=math::Segment(nodes[0].getPoint(),
                                       nodes[2].getPoint());
        }
        else {// nodes[2]==from_node
            opposite_seg=math::Segment(nodes[0].getPoint(),
                                       nodes[1].getPoint());
        }
        
        if(r.intersect3D(opposite_seg, intersection_pnt,
                         param_seg, param_ray)){
            if(param_ray>0){
                return fi;
            }
        }
        
    }//for (auto i=0; i<AFromFaces.size() && !found_first; i++)
    
    throw GMDSException("Found no face to go out");
}
/*---------------------------------------------------------------------------*/
std::vector<Face> TriangularSurfaceManipulator::
findFaces(const math::Point& APnt, Node &AOutNode)
{
    std::cout<<"Look for face containing : "<<APnt<<std::endl;
    std::vector<Face> faces;
    double tolerance2 = 0.001;
    //TO IMPROVE
    bool not_found = true;
    for(IGMesh::face_iterator itf = m_mesh->faces_begin();
        !itf.isDone() && not_found ; itf.next()){
        //we eliminate the triangles wich are too far (it is important since
        //after we work with the projected point
        Face current = itf.value();
        std::vector<Node> current_nodes = current.get<Node>();
        math::Plane pl(current_nodes[0].getPoint(),
                       current_nodes[1].getPoint(),
                       current_nodes[2].getPoint());
        math::Point proj = pl.project(APnt);
        if(proj.distance2(APnt)>tolerance2){
            //this point is too far from the plane of the triangle we test
            // and so the next computation could be non significative
            continue;
        }
        std::vector<TCoord> coords(3,0.0);
        bool on_edge[3] = {false,false,false};
 //       if(isIn(proj,current, on_edge[0], on_edge[1], on_edge[2])){

        if(isIn(proj,current, coords)){
            for(auto i=0; i<3;i++){
                if(coords[i]<0.001)
                    on_edge[i]=true;
            }
            not_found=false;

            std::cout<<"Found face "<<itf.value().getID()<<" with edge values "
            <<on_edge[0]<<" "
            <<on_edge[1]<<" "
            <<on_edge[2]<<std::endl;

            faces = getFaces(itf.value(),on_edge[0],on_edge[1], on_edge[2]);
            if(on_edge[0] && on_edge[1]){
                //means on node 2
                AOutNode = current_nodes[2];
                
            }
            else if(on_edge[1] && on_edge[2]){
                //means on node 0
                AOutNode = current_nodes[0];
            }
            else if(on_edge[2] && on_edge[0]){
                //means on node 1
                AOutNode = current_nodes[1];
            }
            else{
                Node new_node =m_mesh->newNode(proj);
                //node on an edge or inside the face
                std::cout<<"FIND FACE CREATE NODE: "<<new_node.getID()<<std::endl;;
                AOutNode = new_node;
            }
        }
    }
    
    if(not_found){
        throw GMDSException("Impossible to find a triangle containing a pnt");
    }
    

    return faces;
}
/*---------------------------------------------------------------------------*/
std::vector<Face> TriangularSurfaceManipulator::
findFaces(const math::Point& APnt,
          const std::vector<Face>& AAmongFaces,
          Node &AOutNode)
{
    std::cout<<"Look for face containing : "<<APnt<<std::endl;
    std::vector<Face> faces;
    double tolerance2 = 0.001;
    //TO IMPROVE
    bool not_found = true;
    for(auto current:AAmongFaces){
        //we eliminate the triangles wich are too far (it is important since
        //after we work with the projected point
        std::vector<Node> current_nodes = current.get<Node>();
        math::Plane pl(current_nodes[0].getPoint(),
                       current_nodes[1].getPoint(),
                       current_nodes[2].getPoint());
        math::Point proj = pl.project(APnt);
        if(proj.distance2(APnt)>tolerance2){
            //this point is too far from the plane of the triangle we test
            // and so the next computation could be non significative
            continue;
        }
        std::vector<TCoord> coords(3,0.0);
        bool on_edge[3] = {false,false,false};
        //       if(isIn(proj,current, on_edge[0], on_edge[1], on_edge[2])){
        
        if(isIn(proj,current, coords)){
            for(auto i=0; i<3;i++){
                if(coords[i]<0.001)
                    on_edge[i]=true;
            }
            not_found=false;
            
            std::cout<<"Found face "<<current.getID()<<" with edge values "
            <<on_edge[0]<<" "
            <<on_edge[1]<<" "
            <<on_edge[2]<<std::endl;
            
            faces = getFaces(current,on_edge[0],on_edge[1], on_edge[2]);
            if(on_edge[0] && on_edge[1]){
                //means on node 2
                AOutNode = current_nodes[2];
                
            }
            else if(on_edge[1] && on_edge[2]){
                //means on node 0
                AOutNode = current_nodes[0];
            }
            else if(on_edge[2] && on_edge[0]){
                //means on node 1
                AOutNode = current_nodes[1];
            }
            else{
                Node new_node =m_mesh->newNode(proj);
                //node on an edge or inside the face
                std::cout<<"FIND FACE CREATE NODE: "<<new_node.getID()<<std::endl;;
                AOutNode = new_node;
            }
        }
    }
    
    if(not_found){
        throw GMDSException("Impossible to find a triangle containing a pnt");
    }
    
    
    return faces;
}

/*---------------------------------------------------------------------------*/
void TriangularSurfaceManipulator::
walk(const math::Point& AFrom,
     math::Point& ATo,
     const std::vector<Face>& AFromFaces,
     std::vector<Node>& ANewNodes,
     std::vector<bool>& ANewFromPnt,
     std::map<TCellID, std::vector<int> >&AData,
     std::vector<Face>& AOutFaces,
     const bool AIsLast)
{
    double point_tolerance=1e-4;
    
#ifdef DEBUG_GMDS
    Log::mng()<<"------------------- NEW WALK --------------------\n";
    Log::mng()<<"From "<<AFrom<<" to "<<ATo<<"\n";
    Log::mng()<<"From faces: ";
    for(auto f:AFromFaces){
        Log::mng()<<" "<<f.getID();
    }
    Log::mng().flush();
#endif //DEBUG_GMDS
    
    math::Point cur_pnt = AFrom;
    Face        cur_tri;
    bool come_from_node = false;
    Node incoming_node;
    //=========================================================================
    // INITIALIZATION - FIND THE FIRST FACE
    //=========================================================================
    // A starting point can be on a curve edge, and start from a face on the
    // wrong side of the curve, it means so that the projection is totally
    // false
    // We correct that during an initialization phase
   
    std::cout<<std::endl;
    if(AFromFaces.size()==1){
        cur_tri=AFromFaces[0];
        
        std::vector<Node> current_nodes = cur_tri.get<Node>();
        math::Plane pl(current_nodes[0].getPoint(),
                       current_nodes[1].getPoint(),
                       current_nodes[2].getPoint());
        cur_pnt = pl.project(AFrom);
    }
    else  if(AFromFaces.size()==2){
        math::Vector init_dir(AFrom,ATo);
        init_dir.normalize();
        //We start from an edge shared by two faces
        Face f0 = AFromFaces[0];
        Face f1 = AFromFaces[1];
        
        std::vector<Node> f0_nodes = f0.get<Node>();
        std::vector<Node> f1_nodes = f1.get<Node>();
        std::vector<Node> f01_nodes;
        for(auto n0: f0_nodes){
            for(auto n1: f1_nodes){
                if(n0==n1)
                    f01_nodes.push_back(n0);
            }
        }
        math::Vector edge_vec(f01_nodes[0].getPoint(),
                              f01_nodes[1].getPoint());
        edge_vec.normalize();
        
        math::Vector normal0 = f0.normal();
        math::Vector normal1 = f1.normal();
        normal0.normalize();
        normal1.normalize();
        
        math::Point c0 = f0.center();
        math::Point c1 = f1.center();
        math::Vector d0   = normal0.cross(edge_vec);
        math::Vector d1   = normal1.cross(edge_vec);

        math::Vector w0(f01_nodes[0].getPoint(),c0);
        math::Vector w1(f01_nodes[0].getPoint(),c1);
        
        if(d0.dot(w0)<0)
            d0 = d0.opp();
        if(d1.dot(w1)<0)
            d1 = d1.opp();
        
        //the normal which is the most orthogonal with
        //the direction is the one to be kept
        if(d0.dot(init_dir) > d1.dot(init_dir)){
            cur_tri=f0;
        }
        else{
            cur_tri=f1;
        }
    }
    else{
        //We  are on a point. We look for the first triangle containing this point
        //First of all, we get the common node. For that, we go through the 3 first
        // faces
        std::vector<TCellID> n0 = AFromFaces[0].getIDs<Node>();
        std::vector<TCellID> n1 = AFromFaces[1].getIDs<Node>();
        std::vector<TCellID> n2 = AFromFaces[2].getIDs<Node>();
        TCellID common_node_id =-1;
        for(auto i0:n0){
            for(auto i1:n1){
                for(auto i2:n2){
                    if(i0==i1 && i0==i2){
                        common_node_id=i0;
                    }
                }
            }
        }
        if(common_node_id==-1){
            throw GMDSException("No common node between 3 triangles");
        }
        

        Node from_node = m_mesh->get<Node>(common_node_id);
        cur_tri = getOutFace(from_node, ATo);
    }
    
#ifdef DEBUG_GMDS
    std::cout<<"INIT FACE CORRECTION. New Face "<<cur_tri.getID()<<std::endl;
#endif //DEBUG_GMDS
    
    bool finished = false;
    Face prev_tri = cur_tri;
    while(!finished){
        
#ifdef DEBUG_GMDS
        std::cout<<"Walk in "<<cur_tri.getID()<<std::endl;
        std::cout<<" \t from "<<cur_pnt<<std::endl;
        std::cout<<" \t to   "<<ATo<<std::endl;
#endif //DEBUG_GMDS
        
        math::Vector cur_normal = cur_tri.normal();
        std::vector<Node> cur_nodes = cur_tri.get<Node>();
        math::Plane cur_plane(cur_nodes[0].getPoint(),
                              cur_nodes[1].getPoint(),
                              cur_nodes[2].getPoint());
        math::Point dir_pnt= cur_plane.project(ATo);
        math::Vector cur_dir(cur_pnt,dir_pnt);
        
        if(cur_dir.norm2()<point_tolerance){
#ifdef DEBUG_GMDS
            std::cout<<"termination case on a point or curve"<<std::endl;
#endif //DEBUG_GMDS
            //===========================================================
            // TERMINATION CASE ON A POINT OR A CURVE
            //===========================================================
            //means we have already reach the point ATo which was so on a
            //node or an edge.
            
            if(come_from_node){
                //We return all the faces adj to the incoming node
                AOutFaces = incoming_node.get<Face>();
            }
            else{
                //On the edge between the current and the previous
                // face
                AOutFaces.clear();
                AOutFaces.resize(2);
                AOutFaces[0]=cur_tri;
                AOutFaces[1]=prev_tri;
#ifdef DEBUG_GMDS
                std::cout<<"Faces: "<<AOutFaces[0].getID()<<
                ", "<<AOutFaces[0].getID()<<std::endl;
#endif //DEBUG_GMDS
            }
            
            return;
        }
        
        math::Ray ray(cur_pnt,dir_pnt);
        
        math::Segment side[3] = {
            math::Segment(cur_nodes[0].getPoint(), cur_nodes[1].getPoint()),
            math::Segment(cur_nodes[1].getPoint(), cur_nodes[2].getPoint()),
            math::Segment(cur_nodes[2].getPoint(), cur_nodes[0].getPoint())
        };
        std::pair<TCellID,TCellID> side_ids[3]={
            std::pair<TCellID,TCellID>(cur_nodes[0].getID(),cur_nodes[1].getID()),
            std::pair<TCellID,TCellID>(cur_nodes[1].getID(),cur_nodes[2].getID()),
            std::pair<TCellID,TCellID>(cur_nodes[2].getID(),cur_nodes[0].getID())
        };
#ifdef DEBUG_GMDS
        std::cout<<"Edge 0 "<<cur_nodes[0].getID()<<", "<<cur_nodes[1].getID()<<std::endl;
        std::cout<<"Edge 1 "<<cur_nodes[1].getID()<<", "<<cur_nodes[2].getID()<<std::endl;
        std::cout<<"Edge 2 "<<cur_nodes[2].getID()<<", "<<cur_nodes[0].getID()<<std::endl;
        if(come_from_node)
            std::cout<<"Incoming node: "<<incoming_node.getID()<<std::endl;

#endif// DEBUG_GMDS

        //==================================================================
        // ALIGN WITH AN EDGE OF THE CURRENT FACE
        //==================================================================
        //we check first colinearities

        int align_with = -1;
        math::Vector normalized_dir = cur_dir;
        
        normalized_dir.normalize();
        for(auto i=0; i<3; i++){
            math::Point p0 = side[i].getPoint(0);
            math::Point p1 = side[i].getPoint(1);
            if(cur_pnt.distance2(p0)<cur_pnt.distance2(p1)){
                math::Point p_tmp = p0;
                p0 = p1;
                p1 = p_tmp;
            }
            math::Vector vi(p0,p1);
            vi.normalize();
            math::Vector vip(p0, cur_pnt);
            vip.normalize();
            std::cout<<"Side "<<side_ids[i].first<<", "<<side_ids[i].second<<std::endl;
            std::cout<<vi<<", "<<vip<<std::endl;
            std::cout<<"\t "<<(1-abs(normalized_dir.dot(vi)))<<std::endl;
            std::cout<<"\t "<<1-abs(vip.dot(vi))<<std::endl;
            bool is_parallel = (1-abs(normalized_dir.dot(vi))<point_tolerance);
            bool is_in = (1-abs(vip.dot(vi))<point_tolerance);
            std::cout<<"\t -> "<<is_parallel<<", "<<is_in<<std::endl;
            if(is_in && is_parallel){
                align_with=i;
            }
        }
        if(align_with!=-1){
#ifdef DEBUG_GMDS
            std::cout<<"Aligned with an edge"<<std::endl;
#endif// DEBUG_GMDS

            //we are along an edge
            math::Vector to0(cur_pnt,side[align_with].getPoint(0));
            Node next_node;
            Node prev_node;
            if(cur_dir.dot(to0)>0){
                //we go toward point 0
                next_node = m_mesh->get<Node>(side_ids[align_with].first);
                prev_node = m_mesh->get<Node>(side_ids[align_with].second);
            }
            else{
                //we go toward point 1
                next_node = m_mesh->get<Node>(side_ids[align_with].second);
                prev_node = m_mesh->get<Node>(side_ids[align_with].first);

            }
            
#ifdef DEBUG_GMDS
            std::cout<<"\t go towards node "<<next_node.getID()<<std::endl;
#endif// DEBUG_GMDS

            //We have the next node and the current point.
            //Question is ATO betweeen those? If yes we are over
            math::Vector to_next_node(cur_pnt,next_node.getPoint());
            math::Vector to_destination(cur_pnt,ATo);
            std::cout<<"Edge length: "<<to_next_node.norm2()<<std::endl;
            std::cout<<"Remainding dist: "<<to_destination.norm2()<<std::endl;
            if(to_next_node.norm2()>to_destination.norm2()){
                //stop on this edge!!!
                math::Segment seg(prev_node.getPoint(),next_node.getPoint());
                ATo = seg.project(ATo);
                Node new_node = m_mesh->newNode(ATo);
                ANewNodes.push_back(new_node);
                ANewFromPnt.push_back(true);
#ifdef DEBUG_GMDS
                std::cout<<"(11) FINISH ON A CURVE EDGE "<<new_node.getID()<<std::endl;
#endif //DEBUG_GMDS
                std::vector<Face> adj_f = getFaces(prev_node, next_node);
                for(auto f:adj_f){
                    std::cout<<f.getID()<<" ";
                    AData[f.getID()].push_back(ANewNodes.size()-1);
                }
                AOutFaces = adj_f;
                return;
            }
            else if(ATo.distance2(next_node.getPoint())<point_tolerance){
                //We stop  on the next point
                ATo = next_node.getPoint();
                ANewNodes.push_back(next_node);
                ANewFromPnt.push_back(true);
#ifdef DEBUG_GMDS
                std::cout<<"(12) FINISH ON A CURVE NODE "<<next_node.getID()<<std::endl;
#endif// DEBUG_GMDS

                //All the faces containing next_node will be modified
                for(auto f:next_node.get<Face>()){
                    AData[f.getID()].push_back(ANewNodes.size()-1);
                }
                AOutFaces = next_node.get<Face>();
                return;
            }
            else{
                Face next_face;
                
                findOutFace(next_node, ATo, cur_tri, next_face);
                
                prev_tri=cur_tri;
                cur_tri = next_face;
                cur_pnt = next_node.getPoint();
                
                come_from_node=true;
                incoming_node=next_node;
                
                ANewNodes.push_back(next_node);
                ANewFromPnt.push_back(false);

                //All the faces containing next_node will be modified
                for(auto f:next_node.get<Face>()){
                    AData[f.getID()].push_back(ANewNodes.size()-1);
                }
            }
        }
        //==================================================================
        // NOT ALIGNED, SO INTERSECT IN AN EDGE OR A VERTEX
        //==================================================================
        else{
            std::cout<<"Not aligned with an edge"<<std::endl;

            if(come_from_node){
                std::cout<<"come from a node"<<std::endl;
                Node n1, n2;
                if(cur_nodes[0]==incoming_node){
                    n1 = cur_nodes[1];
                    n2 = cur_nodes[2];
                }
                else if(cur_nodes[1]==incoming_node){
                    n1 = cur_nodes[0];
                    n2 = cur_nodes[2];
                }
                else {// nodes[2]==from_node
                    n1 = cur_nodes[0];
                    n2 = cur_nodes[1];
                }
                std::cout<<"Incoming node: "<<incoming_node.getID()<<std::endl;
                std::cout<<"Other nodes: "<<n1.getID()<<", "<<n2.getID()<<std::endl;
                math::Segment opposite_seg=math::Segment(n1.getPoint(),
                                                         n2.getPoint());
                math::Point intersection;
                double param_seg=0, param_ray=0;
                ray.intersect3D(opposite_seg, intersection,param_seg,param_ray);
                std::cout<<"Intersection: "<<intersection<<std::endl;
                std::cout<<"Param Seg: "<<param_seg<<std::endl;
                std::cout<<"Param Ray: "<<param_ray<<std::endl;
                std::cout<<"1.0-param_seg"<<abs(1.0-param_seg)<<std::endl;
                if(abs(param_seg)<0.05){
                    std::cout<<"param_seg = "<<param_seg<<std::endl;
                    std::vector<Face> fs = getFaces(n1,incoming_node);
                    prev_tri=cur_tri;

                    math::Vector edge_vec(incoming_node.getPoint(), n1.getPoint());
                    math::Vector to_vec  (incoming_node.getPoint(), ATo);
                    
                    if(ATo.distance2(n1.getPoint())<point_tolerance){
                        //We stop in n1
                        finished = true;
                        
                        if(!AIsLast){
                            ATo=n1.getPoint();
                            ANewNodes.push_back(n1);
                            ANewFromPnt.push_back(true);

                            //All the faces containing the current node will be modified
                            for(auto f:n1.get<Face>()){
                                AData[f.getID()].push_back(ANewNodes.size()-1);
                            }
                        }
                        AOutFaces = n1.get<Face>();

                    }
                    if(edge_vec.norm2()>to_vec.norm()){
                        //means we stop on this edge
                        finished = true;
                       
                        //we add the final point
                        //before creating a node, we look if it was not created
                        // by another touching loop
                        std::vector<Face> faces =getFaces(incoming_node,n1);
                        bool found = false;
                        Node found_node;
                        for(auto i_f=0; i_f<faces.size()&&!found; i_f++){
                            Face fi = faces[i_f];
                            std::vector<int> vicinity_index= AData[fi.getID()];
                            for(auto i_p=0;i_p<vicinity_index.size()&&!found; i_p++){
                                math::Point pv = ANewNodes[vicinity_index[i_p]].getPoint();
                                if(ATo.distance2(pv)<point_tolerance){
                                    found=true;
                                    found_node = ANewNodes[vicinity_index[i_p]];
                                }
                                
                            }
                            
                        }
                        if(found){
                            if(!AIsLast){
                                ANewNodes.push_back(found_node);
                                ANewFromPnt.push_back(true);
                            }
                        }
                        else if(!AIsLast){
                            math::Segment edge_seg(incoming_node.getPoint(),n1.getPoint());
                            ATo = edge_seg.project(ATo);

                            Node new_node = m_mesh->newNode(ATo);
                            std::cout<<"(1) New Node: "<<new_node.getID()<<std::endl;
                            ANewNodes.push_back(new_node);
                            ANewFromPnt.push_back(true);

                            //All the faces containing the current node will be modified
                            for(auto f:faces){
                                AData[f.getID()].push_back(ANewNodes.size()-1);
                            }
                        }
                        AOutFaces = faces;
                        
                    }
                    else{
                        come_from_node=true;
                        incoming_node=n1;
                        
                        Face next_face;
                        findOutFace(incoming_node, ATo, cur_tri, next_face);
                        prev_tri= cur_tri;
                        cur_tri = next_face;
                        
                        ANewNodes.push_back(n1);
                        ANewFromPnt.push_back(false);

                        //All the faces containing n1 will be modified
                        for(auto f:incoming_node.get<Face>()){
                            AData[f.getID()].push_back(ANewNodes.size()-1);
                        }
                    }
                }
                else if(abs(1.0-param_seg)<0.05){
                    std::cout<<"param_seg = "<<param_seg<<std::endl;
                    std::vector<Face> fs = getFaces(n2,incoming_node);
                    prev_tri=cur_tri;
                    
                    math::Vector edge_vec(incoming_node.getPoint(), n2.getPoint());
                    math::Vector to_vec  (incoming_node.getPoint(), ATo);
                    std::cout<<"IN 2 with dist= "<<ATo.distance2(n2.getPoint())<<std::endl;
                    if(ATo.distance2(n2.getPoint())<point_tolerance){
                        //We stop in n2
                        finished = true;
                        
                        ATo=n2.getPoint();
                        if(!AIsLast){
                            ANewNodes.push_back(n2);
                            ANewFromPnt.push_back(true);

                            //All the faces containing the current node will be modified
                            for(auto f:n2.get<Face>()){
                                AData[f.getID()].push_back(ANewNodes.size()-1);
                            }
                        }
                        AOutFaces = n2.get<Face>();
                        
                    }
                    if(edge_vec.norm2()>to_vec.norm()){
                        //means we stop on this edge
                        finished = true;
                       
                        //we add the final point
                        //before creating a node, we look if it was not created
                        // by another touching loop
                        std::vector<Face> faces =getFaces(incoming_node,n2);
                        std::cout<<"Use faces between "<<incoming_node.getID()<<" and "
                        <<n2.getID()<<std::endl;
                        
                        bool found = false;
                        Node found_node;
                        for(auto i_f=0; i_f<faces.size()&&!found; i_f++){
                            Face fi = faces[i_f];
                            std::vector<int> vicinity_index= AData[fi.getID()];
                            for(auto i_p=0;i_p<vicinity_index.size()&&!found; i_p++){
                                math::Point pv = ANewNodes[vicinity_index[i_p]].getPoint();
                                if(ATo.distance2(pv)<point_tolerance){
                                    found=true;
                                    found_node = ANewNodes[vicinity_index[i_p]];
                                }
                                
                            }
                            
                        }
                        if(found){
                            if(!AIsLast){
                                ANewNodes.push_back(found_node);
                                ANewFromPnt.push_back(true);
                            }
                            
                        }
                        else  if(!AIsLast){
                            math::Segment edge_seg(incoming_node.getPoint(),n2.getPoint());
                            ATo = edge_seg.project(ATo);
                            Node new_node = m_mesh->newNode(ATo);
                            std::cout<<"(2) New Node: "<<new_node.getID()<<std::endl;

                            ANewNodes.push_back(new_node);
                            ANewFromPnt.push_back(true);

                            //All the faces containing the current node will be modified
                            for(auto f:faces){
                                std::cout<<"link from face "<<f.getID()<<" to "
                                <<new_node.getID()<<std::endl;
                                AData[f.getID()].push_back(ANewNodes.size()-1);
                            }
                        }
                        AOutFaces = faces;

                    }
                    else{
                        std::cout<<" So we are here"<<std::endl;
                        come_from_node=true;
                        incoming_node=n2;
                        
                        Face next_face;
                        findOutFace(incoming_node, ATo, cur_tri, next_face);
                        prev_tri=cur_tri;
                        cur_tri = next_face;
                        
                        ANewNodes.push_back(n2);
                        ANewFromPnt.push_back(false);

                        //All the faces containing n2 will be modified
                        for(auto f:incoming_node.get<Face>()){
                            AData[f.getID()].push_back(ANewNodes.size()-1);
                        }
                    }
                    
                }
                else{
                    std::cout<<" in edge"<<std::endl;

                    //on the edge
                    come_from_node=false;
                    std::vector<Face> fs = getFaces(n1,n2);
                    prev_tri=cur_tri;
                    cur_tri = (fs[0]==cur_tri)?fs[1]:fs[0];
                    math::Point next_pnt =( (1-param_seg)*n1.getPoint()
                                           +   param_seg *n2.getPoint());
                    math::Vector v_next(cur_pnt,next_pnt);
                    
                    if(ATo.distance2(next_pnt)<point_tolerance){
                        std::cout<<"Stop 1"<<std::endl;
                        //We stop we reached the point on the edge made of
                        //point n1 and n2
                        finished = true;
                        ATo=next_pnt;
                        
                        //before creating a node, we look if it was not created
                        // by another touching loop
                        bool found = false;
                        Node found_node;
                        for(auto i_f=0; i_f<fs.size()&&!found; i_f++){
                            Face fi = fs[i_f];
                            std::vector<int> vicinity_index= AData[fi.getID()];
                            for(auto i_p=0;i_p<vicinity_index.size()&&!found; i_p++){
                                math::Point pv = ANewNodes[vicinity_index[i_p]].getPoint();
                                if(ATo.distance2(pv)<point_tolerance){
                                    found=true;
                                    found_node = ANewNodes[vicinity_index[i_p]];
                                }
                                
                            }
                            
                        }
                        if(found){
                            if(!AIsLast){
                                ANewNodes.push_back(found_node);
                                ANewFromPnt.push_back(true);

                            }
                            
                        }
                        else  if(!AIsLast){
                            Node new_node = m_mesh->newNode(ATo);
                            std::cout<<"(3) New Node: "<<new_node.getID()<<std::endl;

                            ANewNodes.push_back(new_node);
                            ANewFromPnt.push_back(true);

                            //All the faces containing the current node will be modified
                            for(auto f:fs){
                                AData[f.getID()].push_back(ANewNodes.size()-1);
                            }
                        }
                        AOutFaces = fs;

                    }
                    if(v_next.norm()>cur_dir.norm()){
                        std::cout<<"Stop 2"<<std::endl;
                        //means we stop in this face;
                        finished = true;
                        
                        std::vector<Face> faces;
                        faces.resize(1);
                        faces[0]=prev_tri;
                        
                        //before creating a node, we look if it was not created
                        // by another touching loop
                        bool found = false;
                        Node found_node;
                        Face fi = prev_tri;
                        std::vector<int> vicinity_index= AData[fi.getID()];
                        for(auto i_p=0;i_p<vicinity_index.size()&&!found; i_p++){
                            math::Point pv = ANewNodes[vicinity_index[i_p]].getPoint();
                            if(ATo.distance2(pv)<point_tolerance){
                                found=true;
                                found_node = ANewNodes[vicinity_index[i_p]];
                            }
                            
                            
                        }
                        if(found){
                            if(!AIsLast){
                                 ANewNodes.push_back(found_node);
                                ANewFromPnt.push_back(true);

                            }
                            
                        }
                        else  if(!AIsLast) {
                            Node new_node = m_mesh->newNode(ATo);

                            std::cout<<"(4) New Node: "<<new_node.getID()<<std::endl;

                            ANewNodes.push_back(new_node);
                            ANewFromPnt.push_back(true);

                            AData[prev_tri.getID()].push_back(ANewNodes.size()-1);
                        }
                        AOutFaces = faces;
                    }
                    else{
                        cur_pnt=next_pnt;
                        //We keep going
                        
                        //before creating a node, we look if it was not created
                        // by another touching loop
                        bool found = false;
                        Node found_node;
                        for(auto i_f=0; i_f<fs.size()&&!found; i_f++){
                            Face fi = fs[i_f];
                            std::vector<int> vicinity_index= AData[fi.getID()];
                            for(auto i_p=0;i_p<vicinity_index.size()&&!found; i_p++){
                                math::Point pv = ANewNodes[vicinity_index[i_p]].getPoint();
                                if(cur_pnt.distance2(pv)<point_tolerance){
                                    found=true;
                                    found_node = ANewNodes[vicinity_index[i_p]];
                                }
                                
                            }
                            
                        }
                        if(found){
                            ANewNodes.push_back(found_node);
                            ANewFromPnt.push_back(false);

                            
                        }
                        else {
                            Node new_node = m_mesh->newNode(cur_pnt);
                            std::cout<<"(5) New Node: "<<new_node.getID()<<std::endl;

                            ANewNodes.push_back(new_node);
                            ANewFromPnt.push_back(false);

                            //All the faces containing the current node will be modified
                            for(auto f:fs){
                                AData[f.getID()].push_back(ANewNodes.size()-1);
                            }
                        }
                        std::cout<<"keep going"<<std::endl;
                    }
                    
                }
                
            }
            else{
                std::cout<<"general case"<<std::endl;

                //We arrive in the face from an edge
                math::Point p[3];
                double param_seg[3]={-1,-1,-1};
                double param_ray[3]={-1,-1,-1};
                bool intersect[3]={false,false,false};
                for(auto i=0; i<3; i++){
                    intersect[i]=ray.intersect3D(side[i], p[i], param_seg[i],
                                                 param_ray[i]);
                    if(param_ray[i]<0.001){
                        intersect[i]=false;
                        //it means that we remain on the edge, otherwise we
                        // would go and back from this face and the previous
                        // one
                    }
                }
                std::cout<<"Seg: "<<param_seg[0]<<", "<<param_seg[1]<<", "
                <<param_seg[2]<<std::endl;
                std::cout<<"Ray: "<<param_ray[0]<<", "<<param_ray[1]<<", "
                <<param_ray[2]<<std::endl;
                std::cout<<"Int: "<<intersect[0]<<", "<<intersect[1]<<", "
                <<intersect[2]<<std::endl;
                if(intersect[0] && intersect[1] && intersect[2]){
                    throw GMDSException("Cannot intersect the 3 side of a triangle");
                }
                int on_edge=-1;
                int on_node=-1;
                math::Point out_pnt;
                if(!intersect[0] && !intersect[1] && !intersect[2]){
                    //We encounter a numerical issue where a the direction we go
                    //to must be aligned with the edge we arrive from
                    //We rectify it so
                    throw GMDSException("Numerical issue in the ray/segment intersection");

                }
                if(intersect[0] && intersect[1]){
                    std::cout<<"H 1"<<std::endl;
                    //intersect the edge [0,1] and [1,2]
                    // so either we get out in node 1,
                    //or it means that one of the edge is intersected
                    //out of the triangle.
                    if(p[0].distance2(cur_nodes[1].getPoint())<point_tolerance){
                        on_node = 1;
                        out_pnt=cur_nodes[1].getPoint();
                    }
                    else if(cur_pnt.distance2(p[0])<cur_pnt.distance2(p[1])){
                        if(cur_pnt.distance2(p[0])>point_tolerance){
                            //means we intersect edge 0 first
                            on_edge = 0;
                            out_pnt = p[0];
                        }
                        else{
                            //means we intersect edge 1 first
                            on_edge = 1;
                            out_pnt = p[1];

                        }
                    }
                    else{
                        if(cur_pnt.distance2(p[1])>point_tolerance){
                            //means we intersect edge 1 first
                            on_edge = 1;
                            out_pnt = p[1];
                        }
                        else{
                            //means we intersect edge 0 first
                            on_edge = 0;
                            out_pnt = p[0];
                            
                        }
                    }
                    
                }
                else if(intersect[0] && intersect[2]){
                    std::cout<<"H 2"<<std::endl;
                    //intersect the edge [0,1] and [0,2] so get out in node 0
                    // so either we get out in node 0,
                    //or it means that one of the edge is intersected
                    //out of the triangle.
                    if(p[0].distance2(cur_nodes[0].getPoint())<point_tolerance){
                        on_node = 0;
                        out_pnt=cur_nodes[0].getPoint();
                    }
                    else if(cur_pnt.distance2(p[0])<cur_pnt.distance2(p[2])){
                        if(cur_pnt.distance2(p[0])>point_tolerance){
                            //means we intersect edge 0 first
                            on_edge = 0;
                            out_pnt = p[0];
                        }
                        else{
                            //means we intersect edge 2 first
                            on_edge = 2;
                            out_pnt = p[2];
                            
                        }
                    }
                    else{
                        if(cur_pnt.distance2(p[2])>point_tolerance){
                            //means we intersect edge 2 first
                            on_edge = 2;
                            out_pnt = p[2];
                        }
                        else{
                            //means we intersect edge 0 first
                            on_edge = 0;
                            out_pnt = p[0];
                            
                        }
                    }
                }
                else if(intersect[1] && intersect[2]){
                    std::cout<<"H 3"<<std::endl;
                    //intersect the edge [0,2] and [1,2] so get out in node 2
                    // so either we get out in node 1,
                    //or it means that one of the edge is intersected
                    //out of the triangle.
                    std::cout<<p[1]<<std::endl;
                    std::cout<<cur_nodes[2].getPoint()<<std::endl;
                    std::cout<<p[1].distance2(cur_nodes[2].getPoint())<<std::endl;
                    std::cout<<p[1].distance2(cur_nodes[2].getPoint())<<std::endl;
                    if(p[1].distance2(cur_nodes[2].getPoint())<point_tolerance){
                        on_node = 2;
                        out_pnt=cur_nodes[2].getPoint();
                    }
                    else if(cur_pnt.distance2(p[1])<cur_pnt.distance2(p[2])){
                        if(cur_pnt.distance2(p[1])>point_tolerance){
                            //means we intersect edge 1 first
                            on_edge = 1;
                            out_pnt = p[1];
                        }
                        else{
                            //means we intersect edge 2 first
                            on_edge = 2;
                            out_pnt = p[2];
                            
                        }
                    }
                    else{
                        if(cur_pnt.distance2(p[2])>point_tolerance){
                            //means we intersect edge 2 first
                            on_edge = 2;
                            out_pnt = p[2];
                        }
                        else{
                            //means we intersect edge 1 first
                            on_edge = 1;
                            out_pnt = p[1];
                            
                        }
                    }
                }
                else if(intersect[0]){
                    std::cout<<"H 4"<<std::endl;
                    //intersect the edge [0,1]
                    on_edge = 0;
                    out_pnt = p[0];
                    //but very close to one extremity??
                    if(out_pnt.distance2(cur_nodes[1].getPoint())<point_tolerance){
                        on_node = 1;
                        out_pnt=cur_nodes[1].getPoint();
                    }
                    if(out_pnt.distance2(cur_nodes[0].getPoint())<point_tolerance){
                        on_node = 0;
                        out_pnt=cur_nodes[0].getPoint();
                    }
                
                }
                else if(intersect[1]){
                    std::cout<<"H 5"<<std::endl;
                    //intersect the edge [1,2]
                    on_edge = 1;
                    out_pnt = p[1];
                    //but very close to one extremity??
                    if(out_pnt.distance2(cur_nodes[1].getPoint())<point_tolerance){
                        on_node = 1;
                        out_pnt=cur_nodes[1].getPoint();
                    }
                    if(out_pnt.distance2(cur_nodes[2].getPoint())<point_tolerance){
                        on_node = 2;
                        out_pnt=cur_nodes[2].getPoint();
                    }
                }
                else if(intersect[2]){
                    std::cout<<"H 6"<<std::endl;
                    //intersect the edge [0,2]
                    on_edge = 2;
                    out_pnt = p[2];
                    //but very close to one extremity??
                    if(out_pnt.distance2(cur_nodes[2].getPoint())<point_tolerance){
                        std::cout<<"\t and out in node  "<<cur_nodes[2]<<std::endl;
                        on_node = 2;
                        out_pnt=cur_nodes[2].getPoint();
                    }
                    if(out_pnt.distance2(cur_nodes[0].getPoint())<point_tolerance){
                        std::cout<<"\t and out in node  "<<cur_nodes[0]<<std::endl;
                        on_node = 0;
                        out_pnt=cur_nodes[0].getPoint();
                    }
                }
                std::cout<<"DISTANCE: "<<out_pnt.distance2(ATo)<<std::endl;
                std::cout<<"ON EDGE:  "<<on_edge<<std::endl;
                std::cout<<"Distance to arrrival: "<<out_pnt.distance2(ATo)<<std::endl;
                if(out_pnt.distance2(ATo)<point_tolerance){
                    std::cout<<"Arrived 1"<<std::endl;
                    //it means we finish the algorithm on this point
                    finished = true;
                    
                    //we add the final point
                    ATo=out_pnt;
                    if(on_node!=-1){
                        //we finish on a point
                        Node out_node = cur_nodes[on_node];
                        std::vector<Face> fs = out_node.get<Face>();
                        if(!AIsLast){
                            ANewNodes.push_back(out_node);
                            ANewFromPnt.push_back(true);

                            //All the faces containing the current node will
                            //be modified
                            for(auto f:fs){
                                AData[f.getID()].push_back(ANewNodes.size()-1);
                            }
                        }
                        AOutFaces = fs;

                    }
                    else if(on_edge!=-1){
                        //we finish on a curve
                        int idx[2]={-1,-1};
                        if(on_edge==0){
                            idx[0]=0;
                            idx[1]=1;
                        }
                        else if(on_edge==1){
                            idx[0]=1;
                            idx[1]=2;
                        }
                        else{//2
                            idx[0]=0;
                            idx[1]=2;
                        }

                        //intersect the edge [0,1]
                        std::vector<Face> fs = getFaces(cur_nodes[idx[0]],
                                                        cur_nodes[idx[1]]);
                        
                        
                        //before creating a node, we look if it was not created
                        // by another touching loop
                        bool found = false;
                        Node found_node;
                        for(auto i_f=0; i_f<fs.size()&&!found; i_f++){
                            Face fi = fs[i_f];
                            std::vector<int> vicinity_index= AData[fi.getID()];
                            for(auto i_p=0;i_p<vicinity_index.size()&&!found; i_p++){
                                math::Point pv = ANewNodes[vicinity_index[i_p]].getPoint();
                                if(out_pnt.distance2(pv)<point_tolerance){
                                    found=true;
                                    found_node = ANewNodes[vicinity_index[i_p]];
                                }
                                
                            }
                            
                        }
                        if(found){
                            if(!AIsLast){
                                ANewNodes.push_back(found_node);
                                ANewFromPnt.push_back(true);

                            }
                            
                        }
                        else  if(!AIsLast){

                            Node new_node = m_mesh->newNode(out_pnt);
                            std::cout<<"(6) New Node: "<<new_node.getID()<<std::endl;

                            ANewNodes.push_back(new_node);
                            ANewFromPnt.push_back(true);
                            //All the faces containing the current node will be modified
                            for(auto f:fs){
                                AData[f.getID()].push_back(ANewNodes.size()-1);
                            }
                        }

                        AOutFaces = fs;
                    }
                    else{
                        //we are in the current face
                        std::cout<<"*********************"<<std::endl;
                        std::cout<<"For point "<<out_pnt<<" -> stop in face "
                        <<cur_tri.getID()<<std::endl;
                        
                        std::vector<Face> faces(1,cur_tri);
                        
                        
                        
                        //before creating a node, we look if it was not created
                        // by another touching loop
                        bool found = false;
                        Node found_node;
                        Face fi = cur_tri;
                        std::vector<int> vicinity_index= AData[fi.getID()];
                        for(auto i_p=0;i_p<vicinity_index.size()&&!found; i_p++){
                            math::Point pv = ANewNodes[vicinity_index[i_p]].getPoint();
                            if(out_pnt.distance2(pv)<point_tolerance){
                                found=true;
                                found_node = ANewNodes[vicinity_index[i_p]];
                            }
                            
                        }
                        if(found){
                            if(!AIsLast){
                                
                                ANewNodes.push_back(found_node);
                                ANewFromPnt.push_back(true);
                                
                            }
                        }
                        else  {
                            Node new_node = m_mesh->newNode(out_pnt);
                            std::cout<<"(7) New Node: "<<new_node.getID()<<std::endl;

                            ANewNodes.push_back(new_node);
                            ANewFromPnt.push_back(true);
                            //All the faces containing the current node will be modified
                            for(auto f:faces){
                                AData[f.getID()].push_back(ANewNodes.size()-1);
                            }
                        }
                        
                        AOutFaces = faces;
                    }
                }
                else if(cur_pnt.distance2(out_pnt)>cur_pnt.distance2(ATo)){
                    std::cout<<"Arrived in a face"<<std::endl;
                    //We are in the final face
                    
                    finished = true;
                    
                    //we add the final point
                    //All the faces containing the current node will be modified
                    std::vector<Face> faces(1,cur_tri);
                    
                    
                    //before creating a node, we look if it was not created
                    // by another touching loop
                    bool found = false;
                    Node found_node;
                    Face fi = cur_tri;
                    std::vector<int> vicinity_index= AData[fi.getID()];
                    for(auto i_p=0;i_p<vicinity_index.size()&&!found; i_p++){
                        math::Point pv = ANewNodes[vicinity_index[i_p]].getPoint();
                        if(ATo.distance2(pv)<point_tolerance){
                            found=true;
                            found_node = ANewNodes[vicinity_index[i_p]];
                        }
                        
                    }
                    if(found){
                        if(!AIsLast){
                            ANewNodes.push_back(found_node);
                            ANewFromPnt.push_back(true);
                        }
                        
                    }
                    else if(!AIsLast) {
                        Node new_node = m_mesh->newNode(ATo);
                        std::cout<<"(8) New Node: "<<new_node.getID()<<std::endl;

                        ANewNodes.push_back(new_node);
                        ANewFromPnt.push_back(true);
                        //All the faces containing the current node will be modified
                        for(auto f:faces){
                            AData[f.getID()].push_back(ANewNodes.size()-1);
                        }
                    }
                    
                    AOutFaces = faces;
                    
                }
                else if(on_node!=-1){
                    std::cout<<"Get through node"<<std::endl;
                    Node out_node = cur_nodes[on_node];
                    //We go out by a node
                    Face next_face;
                    math::Vector next_dir;
                    findOutFace(out_node, ATo, cur_tri, next_face);
                    prev_tri=cur_tri;
                    cur_tri = next_face;
                    std::cout<<"cur  face "<<prev_tri.getID()<<std::endl;
                    std::cout<<"next face "<<cur_tri.getID()<<std::endl;
                    cur_pnt = out_node.getPoint();
                    come_from_node =true;
                    incoming_node=out_node;
                    ANewNodes.push_back(out_node);
                    ANewFromPnt.push_back(false);
                    //All the faces containing the current node will be modified
                    std::cout<<"Rink of faces= ";
                    for(auto f:out_node.get<Face>()){
                        AData[f.getID()].push_back(ANewNodes.size()-1);
                        std::cout<<f.getID()<<" ";
                    }
                    std::cout<<std::endl;
                }
                else if (on_edge!=-1){
                    std::cout<<"Get through edge"<<std::endl;
                    int idx[2]={-1,-1};
                    if(on_edge==0){
                        idx[0]=0;
                        idx[1]=1;
                    }
                    else if(on_edge==1){
                        idx[0]=1;
                        idx[1]=2;
                    }
                    else{//2
                        idx[0]=0;
                        idx[1]=2;
                    }
                    //intersect the edge
                    std::vector<Face> fs = getFaces(cur_nodes[idx[0]],
                                                    cur_nodes[idx[1]]);
                    prev_tri=cur_tri;
                    cur_tri = (fs[0]==cur_tri)?fs[1]:fs[0];
                    cur_pnt = out_pnt;
                    
                    come_from_node=false;
                    
                    
                    //before creating a node, we look if it was not created
                    // by another touching loop
                    bool found = false;
                    Node found_node;
                    for(auto i_f=0; i_f<fs.size()&&!found; i_f++){
                        Face fi = fs[i_f];
                        std::vector<int> vicinity_index= AData[fi.getID()];
                        for(auto i_p=0;i_p<vicinity_index.size()&&!found; i_p++){
                            math::Point pv = ANewNodes[vicinity_index[i_p]].getPoint();
                            if(cur_pnt.distance2(pv)<point_tolerance){
                                found=true;
                                found_node = ANewNodes[vicinity_index[i_p]];
                            }
                            
                        }
                        
                    }
                    if(found){
                        ANewNodes.push_back(found_node);
                        ANewFromPnt.push_back(false);
                        
                    }
                    else  {
                        Node new_node = m_mesh->newNode(cur_pnt);
                        std::cout<<"(9) New Node: "<<new_node.getID()<<std::endl;

                        ANewNodes.push_back(new_node);
                        ANewFromPnt.push_back(false);
                        
                        //All the faces containing the current node will be modified
                        for(auto f:fs){
                            AData[f.getID()].push_back(ANewNodes.size()-1);
                        }
                    }
                    
                    AOutFaces = fs;
                }
                
            }
        }
        
        
    }
}

/*---------------------------------------------------------------------------*/
void TriangularSurfaceManipulator::
walkOnCurve(const math::Point& AFrom,
            math::Point& ATo,
            const int AFromClassif,
            const int AToClassif,
            const int AFromGeomID,
            const int AToGeomID,
            const std::vector<Face>& AFromFaces,
            std::vector<Node>& ANewNodes,
            std::vector<bool>& ANewFromPnt,
            std::map<TCellID, std::vector<int> >&AData,
            std::vector<Face>& AOutFaces,
            const bool AIsLast)
{
    std::cout<<"---------------- NEW WALK ON CURVE -----------------"<<std::endl;
    //We only walk through nodes and edges. To do that, we use the common curve
    //id to walk through edges having this curve classification only.
 
    Variable<int>* var_color = 0;
    int surf_colors[2]={-1,-1};
    //found adjacent surfaces
    var_color = m_mesh->getVariable<int>(GMDS_FACE, "BND_SURFACE_COLOR");
    std::cout<<"From faces: ";
    for(auto f:AFromFaces)
        std::cout<<f.getID()<<" ";
    std::cout<<std::endl;
    std::cout<<"From ("<<AFromClassif<<", "<<AFromGeomID<<"): "<<AFrom<<std::endl;
    std::cout<<"To   ("<<AToClassif<<", "<<AToGeomID<<"): "<<ATo<<std::endl;;
    math::Vector init_dir(AFrom,ATo);

    math::Point cur_pnt = AFrom;
    Node cur_node;

    bool finished = false;
    bool come_from_vertex = false;
    Node next_node_from_vertex;
    
    std::vector<Face> from_faces = AFromFaces;

    //=========================================================================
    // INITIALIZATION - FIND THE FIRST NODE
    //=========================================================================
    // the starting point can be on a curve edge, or on a vertex
    // We correct that during an initialization phase
    if(from_faces.size()==1){
        throw GMDSException(" a curve walking can not start inside a face");
    }
    else  if(from_faces.size()==2){
        //We start from an edge shared by two faces
        Face f0 = from_faces[0];
        Face f1 = from_faces[1];
        
        surf_colors[0]= (*var_color)[f0.getID()];
        surf_colors[1]= (*var_color)[f1.getID()];
        
        std::vector<Node> f0_nodes = f0.get<Node>();
        std::vector<Node> f1_nodes = f1.get<Node>();
        std::vector<Node> f01_nodes;
        for(auto n0: f0_nodes){
            for(auto n1: f1_nodes){
                if(n0==n1){
                    f01_nodes.push_back(n0);
                    std::cout<<"Common node "<<n0<<std::endl;
                }
            }
        }

        if(surf_colors[0]==surf_colors[1]){
            //means we had an issue with the geometric tolerance.
            //the current point is likely very close to one point of
            //the edge, we must so get the faces of this node as new reference
            Node n0=f01_nodes[0];
            Node n1=f01_nodes[1];
    
            if(AFrom.distance2(n0.getPoint())<AFrom.distance2(n1.getPoint())){
                from_faces = n0.get<Face>();
            }
            else{
                from_faces = n1.get<Face>();

            }
        }
        
    }
    if(from_faces.size()==2){
        //We start from an edge shared by two faces
        Face f0 = from_faces[0];
        Face f1 = from_faces[1];
        
        surf_colors[0]= (*var_color)[f0.getID()];
        surf_colors[1]= (*var_color)[f1.getID()];
        
        std::vector<Node> f0_nodes = f0.get<Node>();
        std::vector<Node> f1_nodes = f1.get<Node>();
        std::vector<Node> f01_nodes;
        for(auto n0: f0_nodes){
            for(auto n1: f1_nodes){
                if(n0==n1)
                    f01_nodes.push_back(n0);
            }
        }
        
        if(surf_colors[0]==surf_colors[1]){
            //means we had an issue with the geometric tolerance.
              throw GMDSException(" a curve walking with color issue");
        }

        //We get this edge and look towards which point to walk
        math::Vector v0(cur_pnt,f01_nodes[0].getPoint());
        math::Vector v1(cur_pnt,f01_nodes[1].getPoint());
        
        if(v0.dot(init_dir)>v1.dot(init_dir)){
            //walk toward the node 0
            if(v0.norm()>init_dir.norm()){
                //the end is already reached
                finished = true;
            }
            cur_node=f01_nodes[0];
            
        }
        else{
            //walk toward the node 1
            if(v1.norm()>init_dir.norm()){
                //the end is already reached
                finished = true;
            }
            cur_node=f01_nodes[1];
        }
        
        if(finished){
            if(!AIsLast){
                //We ensure to stay on the edge
                math::Segment seg(f01_nodes[0].getPoint(),f01_nodes[1].getPoint());
                ATo = seg.project(ATo);
                Node to_insert = m_mesh->newNode(ATo);
                ANewNodes.push_back(to_insert);
                ANewFromPnt.push_back(true);
                AData[f0.getID()].push_back(ANewNodes.size()-1);
                AData[f1.getID()].push_back(ANewNodes.size()-1);
                AOutFaces.clear();
                AOutFaces.push_back(f0);
                AOutFaces.push_back(f1);
            }
        }
        else{

            ANewNodes.push_back(cur_node);
            ANewFromPnt.push_back(false);
            //All the faces containing next_node will be modified
            for(auto f:cur_node.get<Face>()){
                AData[f.getID()].push_back(ANewNodes.size()-1);
            }
        }
    }
    else{
        //We  are on a vertex, we look for  at most two edges classified
        //on the same curve and we keep the one getting in the right direction
        
        std::vector<TCellID> n0 = from_faces[0].getIDs<Node>();
        std::vector<TCellID> n1 = from_faces[1].getIDs<Node>();
        std::vector<TCellID> n2 = from_faces[2].getIDs<Node>();
        TCellID common_node_id =-1;
        for(auto i0:n0){
            for(auto i1:n1){
                for(auto i2:n2){
                    if(i0==i1 && i0==i2){
                        common_node_id=i0;
                    }
                }
            }
        }
        if(common_node_id==-1){
            throw GMDSException("No common node between 3 triangles");
        }
        Node from_node = m_mesh->get<Node>(common_node_id);
        std::cout<<"Start from node: "<<from_node.getID()<<std::endl;
        
        if(AFromClassif==0){
            come_from_vertex=true;
            //Get surface colors
            std::set<int> colors;
            std::vector<TCellID> from_faces = from_node.getIDs<Face>();
            for(auto f:from_faces) {
                colors.insert((*var_color)[f]);
            }
            std::cout<<"Colors= ";
            for(auto c:colors){
                std::cout<<c<<" ";
            }
            std::cout<<std::endl;
            if(colors.size()<2){
                throw GMDSException("Color error");
            }
            else{
                //Look for the node separing two colors and the most
                //alignd with the direction to go to
                std::vector<Node> candidates;
                std::vector<int> color1;
                std::vector<int> color2;
                std::vector<Face> adj_faces=from_node.get<Face>();
                std::set<TCellID> adj_nodes;
                for(auto f:adj_faces){
                    std::vector<TCellID> nf = f.getIDs<Node>();
                    adj_nodes.insert(nf.begin(),nf.end());
                }
                for(auto adj_n:adj_nodes){
                    if(adj_n==from_node.getID())
                        continue;
                    
                    Node adj_node = m_mesh->get<Node>(adj_n);
                    std::vector<Face> neigh_faces = getFaces(from_node,adj_node);
                    if((*var_color)[neigh_faces[0].getID()]!=
                       (*var_color)[neigh_faces[1].getID()]){
                        candidates.push_back(adj_node);
                    }
                }
                
                math::Vector ref_dir = init_dir;
                ref_dir.normalize();
                int ref=0;
                math::Vector cur_dir(from_node.getPoint(),candidates[0].getPoint());
                cur_dir.normalize();
                double ref_dot=ref_dir.dot(cur_dir);
                for(auto i=1; i<candidates.size();i++){
                    
                    math::Vector cur_dir(from_node.getPoint(),candidates[i].getPoint());
                    cur_dir.normalize();
                    if(ref_dir.dot(cur_dir)>ref_dot){
                        ref=i;
                        ref_dot=ref_dir.dot(cur_dir);
                    }
                }
                
                next_node_from_vertex=candidates[ref];
                //we have the most aligned point, we get the adj colors
                std::vector<Face> neigh_faces = getFaces(from_node,next_node_from_vertex);
                surf_colors[0]= (*var_color)[neigh_faces[0].getID()];
                surf_colors[1]= (*var_color)[neigh_faces[1].getID()];
                

            }
        }
        else{
            //Get surface colors
            std::set<int> colors;
            std::vector<TCellID> from_faces = from_node.getIDs<Face>();
            for(auto f:from_faces) {
                colors.insert((*var_color)[f]);
            }
            if(colors.size()!=2){
                throw GMDSException("Color error");
            }
            
            surf_colors[0]=*(colors.begin());
            surf_colors[1]=*(++colors.begin());
        }
        cur_node = from_node;
    }
    
    while(!finished){
        std::cout<<"Walk from node "<<cur_node.getID()<<std::endl;
        std::cout<<" \t to   "<<ATo<<std::endl;

        cur_pnt=cur_node.getPoint();
        math::Vector cur_dir(cur_pnt,ATo);
        
        
        if(cur_dir.norm()<1e-5){
            //arrived on the point
            ATo = cur_node.getPoint();
            ANewNodes.push_back(cur_node);
            ANewFromPnt.push_back(true);
            std::cout<<"(5) PUSH NODE ON CURVE "<<cur_node.getID()<<std::endl;

            //All the faces containing cur_node will be modified
            for(auto f:cur_node.get<Face>()){
                AData[f.getID()].push_back(ANewNodes.size()-1);
            }
            AOutFaces=cur_node.get<Face>();
            return;
        }
        Node next_node;
        if(come_from_vertex){
            //only possible on the first step from a geom vertex
            come_from_vertex=false;
            next_node = next_node_from_vertex;
        }
        else{
            // We look for the edge incident to cur_node and on the same
            // curve
            std::vector<Node> candidates = getNodes(cur_node,
                                                    surf_colors[0],
                                                    surf_colors[1]);
            
            if(candidates.empty()){
                throw GMDSException("No adjacent node on the curve");
            }
            if(candidates.size()==1){
                throw GMDSException("Only 1 adjacent node on the curve");
            }
            if(candidates.size()>2){
                std::cout<<"Colors "<< surf_colors[0]<<", "<< surf_colors[1]<<std::endl;
                std::cout<<"Candidates: ";
                for(auto c:candidates)
                    std::cout<<c.getID()<<" ";
                
                std::cout<<std::endl;
                throw GMDSException("Two many adjacent nodes on the curve (>2)");
            }
            //We get this edge and look towards which point to walk
            math::Vector v0(cur_pnt,candidates[0].getPoint());
            math::Vector v1(cur_pnt,candidates[1].getPoint());
            
            if(v0.dot(init_dir)>v1.dot(init_dir)){
                //walk toward the node 0
                if(v0.norm()>cur_dir.norm()){
                    //the end is already reached
                    finished = true;
                }
                next_node=candidates[0];
                
            }
            else{
                //walk toward the node 1
                if(v1.norm()>cur_dir.norm()){
                    //the end is already reached
                    finished = true;
                }
                next_node=candidates[1];
            }
        }
        std::vector<Face> cur_faces = getFaces(cur_node,next_node);
        if(finished){
            //before inserting the last node, we check that the current
            //point is in the list already
            if(cur_node!=ANewNodes.back()){
                std::cout<<"(7) PUSH NODE ON CURVE "<<cur_node.getID()<<std::endl;
                
                ANewNodes.push_back(cur_node);
                ANewFromPnt.push_back(false);
                for(auto f:cur_node.get<Face>()){
                    AData[f.getID()].push_back(ANewNodes.size()-1);
                }
            }
            
            math::Vector next_dir(next_node.getPoint(),ATo);
            
            if(next_dir.norm()<1e-5){
                //arrived on a node
                std::cout<<"(4) PUSH NODE ON CURVE "<<next_node.getID()<<std::endl;

                ATo = next_node.getPoint();
                ANewNodes.push_back(next_node);
                ANewFromPnt.push_back(true);
                //All the faces containing cur_node will be modified
                for(auto f:next_node.get<Face>()){
                    AData[f.getID()].push_back(ANewNodes.size()-1);
                }
                AOutFaces=next_node.get<Face>();
                return;
            }
            else {
                if(!AIsLast){
                    //As the node is supposed to be on the curve we move it
                    //to avoid wrong geometric computation afterward.
                    //We ensure to stay on the edge
                    std::vector<TCellID> f0_nodes, f1_nodes;
                    f0_nodes = cur_faces[0].getIDs<Node>();
                    f1_nodes = cur_faces[1].getIDs<Node>();
                    std::vector<TCellID> f01_nodes;
                    for(auto i0:f0_nodes){
                        for(auto i1:f1_nodes){
                            if(i0==i1){
                                f01_nodes.push_back(i0);
                            }
                        }
                    }
                    
                    math::Segment seg(m_mesh->get<Node>(f01_nodes[0]).getPoint(),
                                      m_mesh->get<Node>(f01_nodes[1]).getPoint());
                    ATo = seg.project(ATo);
                    Node to_insert = m_mesh->newNode(ATo);
                    std::cout<<"(2) CREATE NODE ON CURVE "<<to_insert.getID()<<std::endl;
                    ANewNodes.push_back(to_insert);
                    ANewFromPnt.push_back(true);
                    AData[cur_faces[0].getID()].push_back(ANewNodes.size()-1);
                    AData[cur_faces[1].getID()].push_back(ANewNodes.size()-1);
                    
                    AOutFaces=cur_faces;
                }
                return;
            }
        }
        else{
            std::cout<<"Add node "<<cur_node<<" in ANewNodes"<<std::endl;
            ANewNodes.push_back(cur_node);
            ANewFromPnt.push_back(false);
            //All the faces containing next_node will be modified
            for(auto f:cur_node.get<Face>()){
                AData[f.getID()].push_back(ANewNodes.size()-1);
            }
        }
        cur_node=next_node;
        
    }//While(!finished)
}
/*---------------------------------------------------------------------------*/
std::vector<Node> TriangularSurfaceManipulator::
getNodes(const Node& AFrom,
         const int ASurfColor0,
         const int ASurfColor1)
{
    std::vector<Node> to_return;
    Variable<int>* color=m_mesh->getVariable<int>(GMDS_FACE, "BND_SURFACE_COLOR");
    
    std::vector<Face> adj_faces = AFrom.get<Face>();
    std::set<TCellID> adj_nodes;
    for(auto f:adj_faces){
        std::vector<TCellID> nf_ids=f.getIDs<Node>();
        adj_nodes.insert(nf_ids.begin(),nf_ids.end());
    }
    for(auto n_id:adj_nodes){
        Node ni = m_mesh->get<Node>(n_id);
        if(ni==AFrom)
            continue;
        std::vector<Face> common = getFaces(AFrom,ni);
        int c0 = (*color)[common[0].getID()];
        int c1 = (*color)[common[1].getID()];
        if((c0==ASurfColor0 && c1==ASurfColor1) ||
           (c0==ASurfColor1 && c1==ASurfColor0)){
            to_return.push_back(ni);
        }
    }
    
    return to_return;
}
/*---------------------------------------------------------------------------*/
void TriangularSurfaceManipulator::
findOutFace(const Node& AFrom,
            const math::Point& ATo,
            const Face& ACurFace,
            Face& ANextFace)
{
    std::vector<Face> candidates = AFrom.get<Face>();
    math::Point from = AFrom.getPoint();
    
    math::Point intersection;
    double param_seg=0;
    double param_ray=0;
    std::cout<<"findOutFace from node: "<<AFrom.getID()<<" and face "<<ACurFace.getID()<<std::endl;
    for(auto c:candidates){
        if(c==ACurFace)
            continue;
        std::cout<<"\t face "<<c.getID();
        math::Vector normal = c.normal();
        std::vector<Node> nodes = c.get<Node>();
        math::Plane pl(nodes[0].getPoint(),
                       nodes[1].getPoint(),
                       nodes[2].getPoint());
        
        math::Vector dir(AFrom.getPoint(), pl.project(ATo));
        
        math::Ray r(from,dir);
        
        math::Segment opposite_seg;
        if(nodes[0]==AFrom){
            std::cout<<" try to intersect "
            <<nodes[1].getID()<<", "
            <<nodes[2].getID()<<std::endl;
            opposite_seg=math::Segment(nodes[1].getPoint(),
                                       nodes[2].getPoint());
        }
        else if(nodes[1]==AFrom){
            std::cout<<" try to intersect "
            <<nodes[0].getID()<<", "
            <<nodes[2].getID()<<std::endl;
            opposite_seg=math::Segment(nodes[0].getPoint(),
                                       nodes[2].getPoint());
        }
        else {// nodes[2]==AFrom
            std::cout<<" try to intersect "
            <<nodes[0].getID()<<", "
            <<nodes[1].getID()<<std::endl;
            opposite_seg=math::Segment(nodes[0].getPoint(),
                                       nodes[1].getPoint());
        }
        std::cout<<std::endl;
        std::cout<<"From "<<AFrom.getPoint()<<" to "<<pl.project(ATo)<<std::endl;
        if(r.intersect3D(opposite_seg, intersection,param_seg,param_ray)){
            std::cout<<" Intersect with "<<param_ray<<std::endl;
            if(param_ray>0){
                ANextFace = c;
                return;
            }
        }
    }
    
    throw GMDSException("No output face found");
}
/*---------------------------------------------------------------------------*/
std::vector<Face> TriangularSurfaceManipulator::
getFaces(const Node& AN1, const Node& AN2)
{
    std::vector<TCellID> f1 = AN1.getIDs<Face>();
    std::vector<TCellID> f2 = AN2.getIDs<Face>();
    std::vector<Face> f;
    for(auto i1:f1){
        for(auto i2:f2){
            if(i1==i2){
                f.push_back(m_mesh->get<Face>(i1));
            }
        }
    }
    return f;
}
/*---------------------------------------------------------------------------*/
void TriangularSurfaceManipulator::
triangulate(Face& AFace, std::vector<math::Point>& APoints)
{
    
}

/*---------------------------------------------------------------------------*/
std::vector<Face> TriangularSurfaceManipulator::
getFaces(const Face& ATri,
         const bool AOnEdge0,
         const bool AOnEdge1,
         const bool AOnEdge2)
{
    std::vector<Face> faces;
    std::vector<Node> n = ATri.get<Node>();
    if(AOnEdge0 && AOnEdge1){
        //means the point is located on node 2
        std::cout<<"GET ALL FACES OF "<<n[2].getID()<<std::endl;
        faces = n[2].get<Face>();
    }
    else if(AOnEdge0 && AOnEdge2){
        //means the point is located on node 1
        std::cout<<"GET ALL FACES OF "<<n[1].getID()<<std::endl;
        faces = n[1].get<Face>();
    }
    else if(AOnEdge1 && AOnEdge2){
        //means the point is located on node 0
        std::cout<<"GET ALL FACES OF "<<n[0].getID()<<std::endl;
        faces = n[0].get<Face>();
    }
    else if(AOnEdge0){
        //means the point is located on edge [n1,n2]
        faces = getFaces(n[1],n[2]);
    }
    else if(AOnEdge1){
        //means the point is located on edge [n0,n2]
        faces = getFaces(n[0],n[2]);
        
    }
    else if(AOnEdge2){
        //means the point is located on edge [n0,n1]
        faces = getFaces(n[0],n[1]);
        
    }
    else{
        //only in f
        faces.push_back(ATri);
    }
    return faces;
}
/*---------------------------------------------------------------------------*/
bool TriangularSurfaceManipulator::isIn(const math::Point& AP,
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
/*---------------------------------------------------------------------------*/
bool TriangularSurfaceManipulator::isIn(const math::Point& AP,
                                        const Face& ATri,
                                        std::vector<TCoord>& ACoord)
{
    
    std::vector<Node> n = ATri.get<Node>();
    std::vector<math::Point> p;
    p.resize(3);
    p[0] = n[0].getPoint();
    p[1] = n[1].getPoint();
    p[2] = n[2].getPoint();
    
    math::Plane pl(p[0],p[1],p[2]);
    math::Point proj = pl.project(AP);
    try{
        math::Vector V0(p[0],p[1]);
        math::Vector V1(p[0],p[2]);
        math::Vector V2(p[0],proj);
        
        float d00 = V0.dot(V0);
        float d01 = V0.dot(V1);
        float d11 = V1.dot(V1);
        float d20 = V2.dot(V0);
        float d21 = V2.dot(V1);
        if((d00 * d11 - d01 * d01)<1e-8)
            return false;
        
        float invDenom = 1.0 / (d00 * d11 - d01 * d01);
        ACoord[1] = (d11 * d20 - d01 * d21) * invDenom;
        ACoord[2] = (d00 * d21 - d01 * d20) * invDenom;
        ACoord[0] = 1.0 - ACoord[1] - ACoord[2];
        //math::Point::computeBarycentric(p, proj, ACoord);
    }
    catch(GMDSException& e){
        std::cout<<"Exception bary with triangle "<<ATri.getID()<<std::endl;
        std::cout<<"\t PO "<<n[0].getPoint()<<std::endl;
        std::cout<<"\t P1 "<<n[1].getPoint()<<std::endl;
        std::cout<<"\t P2 "<<n[2].getPoint()<<std::endl;
        std::cout<<"\t for point "<<AP<<std::endl;
        return false;
    }
    //==============================================================
    // If APntParam is "almost" into ATetParam, then we provide its
    // Barycentric coordinates
    //==============================================================
    double tolerance = 0.0001;
    if(ACoord[0]<-tolerance || ACoord[0]>1+tolerance ||
       ACoord[1]<-tolerance || ACoord[1]>1+tolerance ||
       ACoord[2]<-tolerance || ACoord[2]>1+tolerance )
        return false;
    
    return true;
}
/*---------------------------------------------------------------------------*/
void TriangularSurfaceManipulator::simplify(std::vector<Node>& AToKeep)
{
    //m is the triangular surface mesh we want to simplify
    
    //AToKeep are the points we want to preserver in this mesh. All the others
    //can be removed
    
    
}/*---------------------------------------------------------------------------*/
char TriangularSurfaceManipulator::orient3d(const gmds::math::Point& AP0,
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
void TriangularSurfaceManipulator::
writePoints(std::vector<math::Point >& ALoop)
{
    std::cout << "Write points (";
    std::cout<< ALoop.size()<<") ...";
    
    IGMesh mesh_loops(MeshModel(DIM3 | F | N | F2N));
    Variable<math::Vector>* v=mesh_loops.newVariable<math::Vector>(GMDS_NODE, "ORIENT");
    for(auto j=0; j<ALoop.size();j++){
        Node n1 = mesh_loops.newNode(ALoop[j]);
        Node n2 = mesh_loops.newNode(ALoop[(j+1)%ALoop.size()]);
        Face f = mesh_loops.newTriangle(n1,n1,n2);
        math::Vector v12(n1.getPoint(),n2.getPoint());
        (*v)[n1.getID()]=v12;
    }

    static int nb_file=0;
    std::stringstream file_name;
    file_name<<"INTERSECTION_LOOP" << nb_file++;
    VTKWriter<IGMesh> writer(mesh_loops);
    writer.write(file_name.str(), DIM3 | F| N);
    std::cout<<" done!"<<std::endl;
}

/*---------------------------------------------------------------------------*/


