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
#include <GMDS/Math/Numerics.h>
#include <GMDS/Math/Plane.h>
#include <GMDS/Math/Segment.h>

/*---------------------------------------------------------------------------*/
// STL File Headers
#include <set>
#include <algorithm>
/*---------------------------------------------------------------------------*/
// FRAME File Headers
#include "WhiskerWeaving.h"
/*---------------------------------------------------------------------------*/
using namespace gmds;
using namespace fhedo;
/*---------------------------------------------------------------------------*/
WhiskerWeaving::WhiskerWeaving():
m_mesh(MeshModel(DIM3|R|F|N|R2N|F2N|N2F)),
m_debug_file("ww_debug")
{
}
/*---------------------------------------------------------------------------*/
WhiskerWeaving::~WhiskerWeaving()
{
}
/*---------------------------------------------------------------------------*/
int WhiskerWeaving::execute(IGMesh* AM,
                            TCellID AF,
                            math::Vector3d& ADir)
{
    m_front_mark = m_mesh.getNewMark<Face>();

    prepareMesh(AM, AF, ADir);

    //==================================================================
    // We build hex and prism as much as possible
    //==================================================================
    std::vector<Loop> loops = buildLoops();
    for(auto& l:loops){
        computeData(l);
    }
    writeDebugMesh();

    std::vector<Loop> forbidden_loops;
    while(!loops.empty()){

        Loop l;
        bool dir;
        if(!select(loops,l,dir)){
            loops.clear();
            continue;

        }
        if(shrink(l,dir)){
            std::cout<<"Shrink done"<<std::endl;
            cleanFront();
            writeDebugMesh();
        }
        else{
            // l must be removed from
            forbidden_loops.push_back(l);
        }

        //all loops are built
        loops=buildLoops();
        //then forbidden loops are removed
        if(!forbidden_loops.empty()){
            for(auto i=0; i<loops.size(); i++){
                bool found=false;
                for(auto j=0; j<forbidden_loops.size() && !found;j++){
                    if(same(loops[i],forbidden_loops[j]))
                        found=true;
                }
                if(found){
                    //the current loop must be remove
                    if(i==loops.size()-1){
                        //last element
                        loops.pop_back();
                    }
                    else{
                        //it is switche by the last element
                        loops[i]=loops.back();
                        //and removed
                        loops.pop_back();
                        //and so the index is decremented
                        i--;
                    }

                }
            }
        }

        //data are computed for useful loops
        for(auto& l:loops){
            computeData(l);
        }

    }


    //==================================================================
    // If the mesh is not over, i.e. it remains front faces, we end
    // the process with tet and pyramids (if asked for)
    //==================================================================
    std::vector<TCellID> last_front_faces;
    for(IGMesh::face_iterator itf = m_mesh.faces_begin();
        !itf.isDone(); itf.next()){
        Face f = itf.value();
        if(m_mesh.isMarked(f, m_front_mark))
            last_front_faces.push_back(f.getID());
    }

    //Do we have to mesh some voids???
    if(!last_front_faces.empty()){
        std::cout<<"VOID MUST BE MESHED (size="
        <<last_front_faces.size()<<")"<<std::endl;
    }
    //clean marks
    m_mesh.unmarkAll<Face>(m_front_mark);
    m_mesh.freeMark<Face>(m_front_mark);

    return 0;
}

/*---------------------------------------------------------------------------*/
void WhiskerWeaving::cleanFront(){

    bool found = true;
    while(found){
        found=false;

        for(IGMesh::face_iterator it = m_mesh.faces_begin();
            !it.isDone() && !found; it.next()){
            Face f = it.value();

            //only front faces are considered
            if(!m_mesh.isMarked(f,m_front_mark))
                continue;

            std::vector<TCellID> adj_f = getAdjFaces(f);
            std::set<TCellID> set_adj_f;
            set_adj_f.insert(adj_f.begin(),adj_f.end());

            if(set_adj_f.size()!=adj_f.size()){
                found =true;
                //means one face share at least 3 nodes with
                //with f
                std::map<TCellID,int> nb_occ;
                for(auto i:set_adj_f){
                    nb_occ[i]=0;
                }
                for(auto i:adj_f){
                    nb_occ[i]++;
                }
                std::vector<TCellID> to_fuse;
                for(auto o:nb_occ){
                    if(o.second>1){
                        to_fuse.push_back(o.first);
                    }
                }
                if(to_fuse.size()>1 || to_fuse.empty()){
                    throw GMDSException("Unexpected configuration in front cleaning");
                }

                Face f2 = m_mesh.get<Face>(to_fuse[0]);
                std::vector<Node> n1 = f.get<Node>();
                std::vector<Node> n2 = f2.get<Node>();

                if(nb_occ[to_fuse[0]]>=3){
                    //means f and f2 share 3 edges so share all nodes
                    //n1==n2 N2F and F2N to handle
                    for(auto ni:n1){
                        ni.replace(f2,f);
                    }
                    m_mesh.deleteFace(f2);
                    //face f is no more in te front
                    m_mesh.unmark(f,m_front_mark);
                    writeDebugMesh();


                }
                else{
                    //share 3 nodes and 1 node is free. Found this one and merge.
                    Node alone1;
                    for(auto ni:n1){
                        bool found_ni=false;
                        for(auto nj:n2){
                            if(ni.getID()==nj.getID()){
                                found_ni = true;
                            }
                        }
                        if(!found_ni){
                            alone1=ni;
                        }
                    }
                    Node alone2;
                    for(auto ni:n2){
                        bool found_ni=false;
                        for(auto nj:n1){
                            if(ni.getID()==nj.getID()){
                                found_ni = true;
                            }
                        }
                        if(!found_ni){
                            alone2=ni;
                        }
                    }
                    std::vector<Node> common;
                    for(auto ni:n1){
                        if(ni!=alone1){
                            common.push_back(ni);
                        }
                    }
                    math::Plane pl(common[0].getPoint(),
                                   common[1].getPoint(),
                                   common[2].getPoint());
                    math::Point p1 = alone1.getPoint();
                    math::Point p2 = alone2.getPoint();
                    math::Point proj1 = pl.project(p1);
                    math::Point proj2 = pl.project(p2);
                    if(proj1.distance2(proj2)>0.5*p1.distance2(p2)){
                        found=false;
                        continue;
                    }

                    bool in_same_hex = false;
                    for(IGMesh::region_iterator itr = m_mesh.regions_begin();
                        !itr.isDone() && !in_same_hex; itr.next()){
                        Region r = itr.value();
                        bool found_alone2=false;
                        bool found_alone1=false;
                        for(auto i:r.getIDs<Node>()){
                            if(i==alone1.getID()){
                                found_alone1=true;
                            }
                            else if(i==alone2.getID()){
                                found_alone2=true;
                            }
                        }
                        if(found_alone1 && found_alone2){
                            in_same_hex=true;
                        }
                    }

                    if(in_same_hex){
                        found =false;
                        continue;
                    }

                    alone1.setPoint(0.5*(alone1.getPoint()+alone2.getPoint()));

                    for(auto ni:n2){
                        ni.remove(f2);
                    }
                    for(auto fi:alone2.get<Face>()){
                        fi.replace(alone2,alone1);
                        //by doing that, the face can be degenerated
                        std::set<TCellID> n_ids;
                        std::vector<TCellID> v_ids= fi.getIDs<Node>();
                        n_ids.insert(v_ids.begin(), v_ids.end());
                        if((fi.getType()==GMDS_TRIANGLE && n_ids.size()<3) ||
                           (fi.getType()==GMDS_QUAD && n_ids.size()==2)) {
                            for(auto x:n_ids){
                                m_mesh.get<Node>(x).remove(fi);
                            }
                            m_mesh.deleteFace(fi);
                        }
                        else if(fi.getType()==GMDS_QUAD && n_ids.size()==3) {
                            //the quad must be replace by a triangle
                            // or removed depending on the position of the duplicated
                            // point (consecutive-> triangle, not->remove)
                            std::vector<int> position;
                            for(auto i_p=0; i_p<v_ids.size();i_p++){
                                if(v_ids[i_p]==alone1.getID()){
                                    position.push_back(i_p);
                                }
                            }
                            if(position.size()!=2){
                                throw GMDSException("Topological issue during WW cleaning of quad");
                            }
                            if((position[1]==position[0]+1) ||
                               (position[0]==0 && position[1]==3)){
                                //becomes a triangle
                                TCellID i1 = *(n_ids.begin());
                                TCellID i2 = *(++n_ids.begin());
                                TCellID i3 = *(++(++(n_ids.begin())));
                                Face new_fi = m_mesh.newTriangle(i1,i2,i3);

                                for(auto x:n_ids){
                                    m_mesh.get<Node>(x).replace(fi,new_fi);
                                }
                                if(m_mesh.isMarked(fi,m_front_mark)){
                                    m_mesh.mark(new_fi, m_front_mark);
                                }
                                m_mesh.deleteFace(fi);

                            }
                            else{
                                //the face is removed
                                for(auto x:n_ids){
                                    m_mesh.get<Node>(x).remove(fi);
                                }
                                m_mesh.deleteFace(fi);
                            }
                        }


                        alone1.add(fi);
                    }

                    for(IGMesh::region_iterator itr = m_mesh.regions_begin();
                        !itr.isDone(); itr.next()){
                        Region r = itr.value();
                        bool found_alone2=false;
                        for(auto i:r.getIDs<Node>()){
                            if(i==alone2.getID()){
                                found_alone2=true;
                            }
                        }
                        if(found_alone2){
                            r.replace(alone2,alone1);
                        }
                    }
                    m_mesh.deleteNode(alone2);
                    m_mesh.deleteFace(f2);
                    //face f is no more in te front
                    m_mesh.unmark(f,m_front_mark);
                    writeDebugMesh();



                }
            }
        }
    }
}
/*---------------------------------------------------------------------------*/
void WhiskerWeaving::prepareMesh(IGMesh* AM,
                                 TCellID AF,
                                 math::Vector3d& ADir)
{
    math::Vector dir(ADir.X(),ADir.Y(),ADir.Z());

    //==================================================================
    // Node copy
    //==================================================================
    for(IGMesh::node_iterator it_n = AM->nodes_begin();
        !it_n.isDone(); it_n.next()){
        Node n = it_n.value();
        Node to = m_mesh.newNode(n.getPoint());

        m_from_bnd_nodes[n.getID() ] = to.getID();
        m_to_bnd_nodes  [to.getID()] = n.getID();

    }
    //==================================================================
    // Face copy
    //==================================================================
    for(IGMesh::face_iterator it_f = AM->faces_begin();
        !it_f.isDone(); it_f.next()){
        Face from = it_f.value();
        std::vector<TCellID> from_n = from.getIDs<Node>();
        std::vector<Node> to_n;
        to_n.reserve(from_n.size());
        for(auto i:from_n){
            to_n.push_back(m_mesh.get<Node>(m_from_bnd_nodes[i]));
        }
        Face to = m_mesh.newFace(to_n);
        //all faces are initially in the front
        m_mesh.mark(to,m_front_mark);
        for(auto n:to_n){
            n.add<Face>(to);
        }

        m_from_bnd_nodes[from.getID()] = to.getID();
        m_to_bnd_nodes  [to.getID()  ] = from.getID();

    }

    //==================================================================
    // Re-orient faces definition
    //==================================================================
    // Internal node numerotation of faces are not consistant, we
    // reorient them according to the initial face
    int done = m_mesh.getNewMark<Face>();
    TCellID seed_id = m_from_bnd_faces[AF];
    Face seed = m_mesh.get<Face>(seed_id);
    if(seed.normal().dot(dir)<0){
        //reverse the face numerotation
        std::vector<TCellID> n_ids = seed.getIDs<Node>();
        std::reverse(n_ids.begin(),n_ids.end());
        seed.set<Node>(n_ids);
    }
    m_mesh.mark<Face>(seed_id, done);

    std::vector<TCellID> to_propag;
    to_propag.push_back(seed_id);
    while(!to_propag.empty()){
        TCellID cur_id = to_propag.back();
        to_propag.pop_back();

        Face f = m_mesh.get<Face>(cur_id);
        //f is well-oriented
        std::vector<TCellID> fn = f.getIDs<Node>();

        for(auto i_n=0; i_n<fn.size(); i_n++){
            TCellID ni = fn[i_n];
            TCellID nj = fn[(i_n+1)%fn.size()];
            TCellID fij = getFace(f,ni,nj);
            if(!m_mesh.isMarked<Face>(fij, done)){
                //new face never reoriented
                m_mesh.mark<Face>(fij,done);
                to_propag.push_back(fij);

                //we reorient it if necessary, it must see edge i,j
                //in the opposite direction
                Face opp_ij = m_mesh.get<Face>(fij);
                reverse(opp_ij,nj,ni);
            }
        }

    }
    //We only negate the mask mark since all the faces have been traversed
    m_mesh.unmarkAll<Face>(done);
    m_mesh.freeMark<Face>(done);
}
/*---------------------------------------------------------------------------*/
void WhiskerWeaving::computeData(Loop& AL)
{
    //init all field to zero
    AL.init();

    //===================================================================
    // Check for self-intersection. Self-intersected loop are withdrawn
    //===================================================================
    std::set<TCellID> set_ids;
    set_ids.insert(AL.faces.begin(),AL.faces.end());
    if(set_ids.size()!=AL.faces.size()){
        AL.self_intersect=true;
        return;
    }

    //===================================================================
    int mark_loop = m_mesh.getNewMark<Face>();

    int nb_faces = AL.faces.size();

    AL.right_nodes[0].resize(nb_faces);
    AL.right_nodes[1].resize(nb_faces);
    AL.left_nodes[0].resize(nb_faces);
    AL.left_nodes[1].resize(nb_faces);
    AL.right_order_nodes.resize(nb_faces);
    AL.left_order_nodes.resize(nb_faces);
    for(auto i=0; i<nb_faces; i++){

        Face fi = m_mesh.get<Face>(AL.faces[i]);
        Face fj = m_mesh.get<Face>(AL.faces[(i+1)%nb_faces]);

        m_mesh.mark(fi, mark_loop);
        TCellID n0, n1;
        commonNodes(fi, fj, n0,n1);
        math::Point p0 = m_mesh.get<Node>(n0).getPoint();
        math::Point p1 = m_mesh.get<Node>(n1).getPoint();
        math::Point c = fi.center();
        math::Vector u(c,0.5*(p0+p1));
        math::Vector v(fi.normal());
        math::Vector w = u.cross(v);
        w.normalize();

        //w point to the right;
        std::vector<Node> ni = fi.get<Node>();
        int    best_aligned = -1;
        double best         = -10;
        for(auto i_n=0; i_n<ni.size(); i_n++){
            math::Point a = ni[i_n].getPoint();
            math::Point b = ni[(i_n+1)%ni.size()].getPoint();
            math::Vector ui(c,0.5*(a+b));
            ui.normalize();
            if(w.dot(ui)>best){
                best_aligned=i_n;
                best=w.dot(ui);
            }
        }

        TCellID right_nodes[2] = {
            ni[best_aligned].getID(),
            ni[(best_aligned+1)%ni.size()].getID()
        };

        TCellID left_nodes[2] = {NullID, NullID};

        getOtherNodesFace(fi, right_nodes[0], right_nodes[1],
                          left_nodes[0], left_nodes[1]);

        Face right_face = m_mesh.get<Face>(getFace(fi,right_nodes[0], right_nodes[1]));
        Face left_face  = m_mesh.get<Face>(getFace(fi,left_nodes[0] , left_nodes[1]));
        //===============================
        //We store the two nodes as being on the right of the loop
        AL.right_nodes[0][i]= right_nodes[0];
        AL.right_nodes[1][i]= right_nodes[1];
        AL.left_nodes[0][i] = left_nodes[0];
        AL.left_nodes[1][i] = left_nodes[1];
        if(right_nodes[0]==n0 || right_nodes[0]==n1 ){
            AL.right_order_nodes[i]=right_nodes[0];
        }
        else{
            AL.right_order_nodes[i]=right_nodes[1];
        }
        if(left_nodes[0]==n0 || left_nodes[0]==n1 ){
            AL.left_order_nodes[i]=left_nodes[0];
        }
        else{
            AL.left_order_nodes[i]=left_nodes[1];
        }
        //===============================
        // Flat index
        if(isFlatAngle(dihedralAngle(fi,fj)))
            AL.flat_index++;
        //===============================
        // Left sharp CONVEX
        if(isConvexAngle(dihedralAngle(fi,left_face)))
            AL.left_sharp++;
        //===============================
        // right sharp CONVEX
        if(isConvexAngle(dihedralAngle(fi,right_face)))
            AL.right_sharp++;
    }

    //Now we extract global side info (nb faces and triangles on each side)
    computeLoopSideInfo(AL,mark_loop, true);
    computeLoopSideInfo(AL,mark_loop, false);

//
//    for(auto i=0; i<nb_faces; i++){
//        m_mesh.unmark<Face>(AL.faces[i], mark_loop);
//    }
    m_mesh.unmarkAll<Face>(mark_loop);
    m_mesh.freeMark<Face>(mark_loop);

}

/*---------------------------------------------------------------------------*/
void WhiskerWeaving::
computeLoopSideInfo(Loop& AL, const int AM, const bool side)
{
    TCellID first_node_id = NullID;
    if(side==true){
        //to the right
        first_node_id = AL.right_nodes[0][0];
    }
    else{
        //to the left
        first_node_id = AL.left_nodes[0][0];
    }

    Node first_node = m_mesh.get<Node>(first_node_id);
    std::vector<TCellID> first_faces = first_node.getIDs<Face>();
    std::vector<TCellID> first_side_faces;
    for(auto f:first_faces){
        if(m_mesh.isMarked<Face>(f, m_front_mark) &&
           !m_mesh.isMarked<Face>(f, AM)){
            //means we have a face that is not in the loop
            first_side_faces.push_back(f);
        }
    }

    if(first_side_faces.empty()){
        throw GMDSException("Error, no starting faces in WW::computeLoopSideInfo");
    }

    int nb_faces = 0;
    int nb_triangles=0;
    int done = m_mesh.getNewMark<Face>();
    std::vector<TCellID> done_faces(1,first_side_faces[0]);
    std::vector<TCellID> to_proceed(1,first_side_faces[0]);
    m_mesh.mark<Face>(first_side_faces[0],done);

    while(!to_proceed.empty()){
        TCellID i = to_proceed.back();
        to_proceed.pop_back();
        Face fi = m_mesh.get<Face>(i);
        nb_faces++;
        if(fi.getType()==GMDS_TRIANGLE){
            nb_triangles++;
        }
        std::vector<TCellID> n = fi.getIDs<Node>();

        for(auto i_n=0; i_n<n.size(); i_n++){
            TCellID ni = n[i_n];
            TCellID nj = n[(i_n+1)%n.size()];
            TCellID fij = getFace(fi,ni,nj);
            if(!m_mesh.isMarked<Face>(fij, done) &&
               !m_mesh.isMarked<Face>(fij, AM)) {
                //new face never proceeded
                m_mesh.mark<Face>(fij,done);
                to_proceed.push_back(fij);
                done_faces.push_back(fij);
            }
        }


    }
    //unmark faces and free the mark
    for(auto i:done_faces){
        m_mesh.unmark<Face>(i,done);
    }
    m_mesh.freeMark<Face>(done);

    if(side==true){
        //to the right
        AL.right_nb_triangles = nb_triangles;
        AL.right_nb_faces     = nb_faces;
    }
    else{
        //to the left
        AL.left_nb_triangles = nb_triangles;
        AL.left_nb_faces     = nb_faces;
    }
}

/*---------------------------------------------------------------------------*/
void WhiskerWeaving::
getEnclosedFaces(const Loop& AL, const bool AToRight,
                 std::vector<TCellID>& AF)
{
    int done = m_mesh.getNewMark<Face>();

    TCellID first_node_id = NullID;
    if(AToRight==true){
        //to the right
        first_node_id = AL.right_nodes[0][0];
    }
    else{
        //to the left
        first_node_id = AL.left_nodes[0][0];
    }

    Node first_node = m_mesh.get<Node>(first_node_id);
    std::vector<TCellID> first_faces = first_node.getIDs<Face>();
    std::vector<TCellID> first_side_faces;
    for(auto f:first_faces){
        if(std::find(AL.faces.begin(),AL.faces.end(),f)==AL.faces.end()){
            //means we have a face that is not in the loop
            if(m_mesh.isMarked<Face>(f,m_front_mark))
               first_side_faces.push_back(f);
        }
    }

    if(first_side_faces.empty()){
        throw GMDSException("Error, no starting faces in WW::computeLoopSideInfo");
    }

    std::vector<TCellID> done_faces(1,first_side_faces[0]);
    std::vector<TCellID> to_proceed(1,first_side_faces[0]);
    m_mesh.mark<Face>(first_side_faces[0],done);

    while(!to_proceed.empty()){
        TCellID i = to_proceed.back();
        to_proceed.pop_back();
        Face fi = m_mesh.get<Face>(i);

        std::vector<TCellID> n = fi.getIDs<Node>();

        for(auto i_n=0; i_n<n.size(); i_n++){
            TCellID ni = n[i_n];
            TCellID nj = n[(i_n+1)%n.size()];
            TCellID fij = getFace(fi,ni,nj);
            if(!m_mesh.isMarked<Face>(fij, done) &&
               std::find(AL.faces.begin(),AL.faces.end(),fij)==AL.faces.end()) {
                //new face never proceeded
                m_mesh.mark<Face>(fij,done);
                to_proceed.push_back(fij);
                done_faces.push_back(fij);
            }
        }
    }
    AF = done_faces;

    for(auto f:done_faces){
        m_mesh.unmark<Face>(f,done);
    }
    m_mesh.freeMark<Face>(done);
}
/*---------------------------------------------------------------------------*/
bool WhiskerWeaving::isSharpAngle(const double& AAngle)
{
    return (!isFlatAngle(AAngle));
}
/*---------------------------------------------------------------------------*/
bool WhiskerWeaving::isFlatAngle(const double& AAngle)
{
    return (AAngle<45 && AAngle>-45);
}
/*---------------------------------------------------------------------------*/
bool WhiskerWeaving::isConvexAngle(const double& AAngle)
{
    return (AAngle>45);
}
/*---------------------------------------------------------------------------*/
bool WhiskerWeaving::isConcaveAngle(const double& AAngle)
{
    return (AAngle<-45);
}
/*---------------------------------------------------------------------------*/
std::vector<WhiskerWeaving::Loop> WhiskerWeaving::buildLoops()
{
    std::vector<Loop> loops;

    //==================================================================
    // We get all the quad faces in the front
    //==================================================================
    //only front quad faces are used
    std::vector<TCellID> quads;
    //Two colors are given to each quad:
    //  - one per direction
    //  - 1 means direction already done, (o nope)
    //  - color 1 is for nodes 0 and 1
    //  - color 2 is for nodes 1 and 2
    std::map<TCellID,int> c1;
    std::map<TCellID,int> c2;
    for(IGMesh::face_iterator it = m_mesh.faces_begin();
        !it.isDone(); it.next()){
        Face q = it.value();

        //only front faces are used
        if(!m_mesh.isMarked(q,m_front_mark))
            continue;

        //only quad faces are used
        if(q.getType()!=GMDS_QUAD)
            continue;

        //if a face share two edges with the same other we do not use it
        // to
        std::vector<TCellID> adj_q = getAdjFaces(q);
        std::set<TCellID> set_adj_q;
        set_adj_q.insert(adj_q.begin(),adj_q.end());
        if(set_adj_q.size()!=adj_q.size())
            continue;

        quads.push_back(q.getID());
        c1[q.getID()]=0;
        c2[q.getID()]=0;
    }

    //==================================================================
    // We build loops on this front
    //==================================================================
    for(auto q_id:quads){
        if(c1[q_id]==1 && c2[q_id]==1)
            continue; //nothing to do

        Face q = m_mesh.get<Face>(q_id);
        std::vector<TCellID> qn = q.getIDs<Node>();

        if(c1[q_id]==0){
            Loop l;
            if(buildLoop(q,qn[0],qn[1],c1,c2,l)){
                loops.push_back(l);
            }
        }
        if(c2[q_id]==0){
            Loop l;
            if(buildLoop(q,qn[1],qn[2],c1,c2,l)){
                loops.push_back(l);
            }
        }
    }
    return loops;
}

/*---------------------------------------------------------------------------*/
bool WhiskerWeaving::
select(const std::vector<Loop>& ALoops, Loop& AL, bool& AToRight)
{
    int nb_loops = ALoops.size();

    //============================================================
    // Step 1 - Prepare data associated to (loop, direction)
    //============================================================
    //one candidate = 1 loop id + 1 direction (true=right, false=left)
    std::vector<std::pair<int,bool> > candidates;
    //flatness : 0 (flat) - 1 only sharp edges.
    std::vector<double> flatness;
    // side_sharpness: all traversed edges are CONVEX (1) - 0 no convex
    std::vector<double> side_sharpness;
    //number of faces to shrink
    std::vector<int> nb_faces;
    //proportion of triangles among the faces to be shrunk
    std::vector<int> triangle_proportion;

    candidates.resize(2*nb_loops);
    flatness.resize(2*nb_loops);
    side_sharpness.resize(2*nb_loops);
    nb_faces.resize(2*nb_loops);
    triangle_proportion.resize(2*nb_loops);

    for(auto i=0; i<ALoops.size(); i++){
        Loop li = ALoops[i];
        int li_nb_faces = li.faces.size();


        //right
        candidates[2*i].first=i;
        candidates[2*i].second=true;
        flatness[2*i] = (double)li.flat_index/(double)li_nb_faces;
        side_sharpness[2*i] = (double)li.right_sharp/(double)li_nb_faces;
        nb_faces[2*i] = li.right_nb_faces;
        triangle_proportion[2*i] = (double)li.right_nb_triangles/(double)li.right_nb_faces;
        //left
        candidates[2*i+1].first=i;
        candidates[2*i+1].second=false;
        flatness[2*i+1] = (double)li.flat_index/(double)li_nb_faces;
        side_sharpness[2*i+1] = (double)li.left_sharp/(double)li_nb_faces;
        nb_faces[2*i+1] = li.left_nb_faces;
        triangle_proportion[2*i+1] = (double)li.left_nb_triangles/(double)li.left_nb_faces;

    }
    //============================================================
    // Step 2 - We shrink first the layers with a maximal side
    // sharpness
    //============================================================
    std::vector<int> sharp_filtered;
    int best_sharp = 0; //worst possible case
    for(auto i=0; i<candidates.size(); i++){
        if(ALoops[candidates[i].first].self_intersect)
            continue;

        if(side_sharpness[i]>best_sharp){
            sharp_filtered.clear();
            sharp_filtered.push_back(i);
            best_sharp=side_sharpness[i];
        }
        else if(side_sharpness[i]==best_sharp){
            sharp_filtered.push_back(i);
        }
        //otherwise we do nothing
    }
    if(best_sharp==0){
        //means no more valid candidate
        return false;
    }
    if(sharp_filtered.size()==1){
        //only one candidate!!!
        int c = sharp_filtered[0];
        AL = ALoops[candidates[c].first];
        AToRight = candidates[c].second;
        return true;
    }
    //============================================================
    // Step 3 - Now we select refine to the ones with the smallest
    // number of faces to shrink
    //============================================================
    std::vector<int> size_filter;
    int smallest_size = 10000000;
    for(auto i=0; i<sharp_filtered.size(); i++){
        int c_id = sharp_filtered[i];

        if(nb_faces[c_id]<smallest_size){
            size_filter.clear();
            size_filter.push_back(c_id);
            smallest_size=nb_faces[c_id];
        }
        else if(nb_faces[c_id]==smallest_size){
            size_filter.push_back(c_id);
        }
        //otherwise we do nothing
    }
    if(size_filter.size()>=1){
        //only one candidate!!!
        int c = size_filter[0];
        AL = ALoops[candidates[c].first];
        AToRight = candidates[c].second;
        return true;
    }

    for(auto l:ALoops){
        std::cout<<l<<std::endl;;
    }
    throw GMDSException("Unable to select a loop");

}
/*---------------------------------------------------------------------------*/
bool WhiskerWeaving::shrink(const Loop& AL, const bool AToRight)
{
    std::cout<<"=============== Shrink ==============="<<std::endl;
    std::cout<<AL<<std::endl;

    if(AToRight)
        std::cout<<"to right"<<std::endl;
    else
        std::cout<<"to left"<<std::endl;

    //==================================================================
    // 1) Define from and to bnd nodes
    //==================================================================
    std::map<TCellID, TCellID> node_mapping;
    std::vector<TCellID> ordered_node_from;
    //store the face going from the source to the target side for each
    //edge in the source mesh
    // Node with lowest id is always first.
    std::map<std::pair<TCellID,TCellID>, TCellID> n2face_mapping;

    if(AToRight){
        //from right to left
        ordered_node_from = AL.right_order_nodes;
        for(auto i=0; i<AL.faces.size(); i++){
            Face fi = m_mesh.get<Face>(AL.faces[i]);

            TCellID from_nodes[2] = {
                AL.right_nodes[0][i],
                AL.right_nodes[1][i]};
            TCellID to_nodes[2] = {
                AL.left_nodes[0][i],
                AL.left_nodes[1][i]};
            if(from_nodes[0]<from_nodes[1]){
                std::pair<TCellID,TCellID> e;
                e.first=from_nodes[0];
                e.second=from_nodes[1];
                n2face_mapping[e]=AL.faces[i];
            }
            else{
                std::pair<TCellID,TCellID> e;
                e.first=from_nodes[1];
                e.second=from_nodes[0];
                n2face_mapping[e]=AL.faces[i];
            }
            Node from_0 = m_mesh.get<Node>(from_nodes[0]);
            Node adj1, adj2;
            fi.getAdjacentNodes(from_0, adj1, adj2);
            if(adj1.getID()==from_nodes[1]){
                //now we look for the other mapping
                if(adj2.getID()==to_nodes[0]){
                    node_mapping[from_nodes[0]]=to_nodes[0];
                    node_mapping[from_nodes[1]]=to_nodes[1];
                }
                else{
                    node_mapping[from_nodes[0]]=to_nodes[1];
                    node_mapping[from_nodes[1]]=to_nodes[0];
                }

            }
            else if(adj2.getID()==from_nodes[1]){
                //now we look for the other mapping
                if(adj1.getID()==to_nodes[0]){
                    node_mapping[from_nodes[0]]=to_nodes[0];
                    node_mapping[from_nodes[1]]=to_nodes[1];
                }
                else{
                    node_mapping[from_nodes[0]]=to_nodes[1];
                    node_mapping[from_nodes[1]]=to_nodes[0];
                }
            }
            else{
                throw GMDSException("Impossible configuration for loop shrinking");
            }
        }

    }
    else{
        //from left to right
        ordered_node_from = AL.left_order_nodes;
        for(auto i=0; i<AL.faces.size(); i++){
            Face fi = m_mesh.get<Face>(AL.faces[i]);

            TCellID from_nodes[2] = {
                AL.left_nodes[0][i],
                AL.left_nodes[1][i]};
            TCellID to_nodes[2] = {
                AL.right_nodes[0][i],
                AL.right_nodes[1][i]};
            if(from_nodes[0]<from_nodes[1]){
                std::pair<TCellID,TCellID> e;
                e.first=from_nodes[0];
                e.second=from_nodes[1];
                n2face_mapping[e]=AL.faces[i];
            }
            else{
                std::pair<TCellID,TCellID> e;
                e.first=from_nodes[1];
                e.second=from_nodes[0];
                n2face_mapping[e]=AL.faces[i];
            }
            Node from_0 = m_mesh.get<Node>(from_nodes[0]);
            Node adj1, adj2;
            fi.getAdjacentNodes(from_0, adj1, adj2);
            if(adj1.getID()==from_nodes[1]){
                //now we look for the other mapping
                if(adj2.getID()==to_nodes[0]){
                    node_mapping[from_nodes[0]]=to_nodes[0];
                    node_mapping[from_nodes[1]]=to_nodes[1];
                }
                else{
                    node_mapping[from_nodes[0]]=to_nodes[1];
                    node_mapping[from_nodes[1]]=to_nodes[0];
                }

            }
            else if(adj2.getID()==from_nodes[1]){
                //now we look for the other mapping
                if(adj1.getID()==to_nodes[0]){
                    node_mapping[from_nodes[0]]=to_nodes[0];
                    node_mapping[from_nodes[1]]=to_nodes[1];
                }
                else{
                    node_mapping[from_nodes[0]]=to_nodes[1];
                    node_mapping[from_nodes[1]]=to_nodes[0];
                }
            }
            else{
                throw GMDSException("Impossible configuration for loop shrinking");
            }
        }

    }
    //we can deduced circle nodes from the mapping but with losing the cyclic order!!!
    int mark_source = m_mesh.getNewMark<Node>();
    double average_dist=0;
    std::vector<TCellID> from_circle_nodes, to_circle_nodes;
    from_circle_nodes.reserve(ordered_node_from.size());
    to_circle_nodes.reserve(ordered_node_from.size());
    for(auto id_from:ordered_node_from){
        //std::cout<<m.first<<" -> "<<m.second<<std::endl;
        from_circle_nodes.push_back(id_from);
        m_mesh.mark<Node>(id_from,mark_source);
        TCellID id_to = node_mapping[id_from];
        to_circle_nodes.push_back(id_to);
        math::Point p_s = m_mesh.get<Node>(id_from).getPoint();
        math::Point p_t = m_mesh.get<Node>(id_to).getPoint();
        average_dist +=p_s.distance(p_t);
    }
    average_dist /=node_mapping.size();

    //==================================================================
    // 2) Get source faces
    //==================================================================
    std::vector<TCellID> source_faces;
    getEnclosedFaces(AL, AToRight, source_faces);

    //==================================================================
    // 3) We compute the image of source nodes
    //==================================================================
    std::set<TCellID> source_nodes;
    //for each inner source node, target node position is precomputed
    //before creating a node mesh
    std::map<TCellID,math::Point> source_node_2_target_point;

    //two next structures are necesseray to detect impossible
    // shrinkings
    std::vector<math::Point> target_points;
    std::map<int,TCellID> target_point_2_source_node;

    //We collect all the nodes of the source surface
    for(auto f:source_faces){
        std::vector<TCellID> nf = m_mesh.get<Face>(f).getIDs<Node>();
        source_nodes.insert(nf.begin(),nf.end());
    }
    //source_nodes contains all the source nodes including those on the
    //circle loop
    for(auto s_id:source_nodes){
        if(m_mesh.isMarked<Node>(s_id,mark_source)){
            //we are on a loop node
            Node target_node = m_mesh.get<Node>(node_mapping[s_id]);
            //we store the existing point
            target_points.push_back(target_node.getPoint());
            target_point_2_source_node[target_points.size()-1]=s_id;
            continue;
        }
        m_mesh.mark<Node>(s_id,mark_source);
        //We compute an image of this point
        Node n = m_mesh.get<Node>(s_id);
        std::vector<Face> adj_f = n.get<Face>();
        //We keep only thos in source_faces;
        std::vector<Face> adj_source_f;
        for(auto f:adj_f){
            if(m_mesh.isMarked(f,m_front_mark))
                adj_source_f.push_back(f);
        }
        //image will be directed along the average direction
        math::Vector d(0,0,0);
        for(auto f:adj_source_f){
            d = d+f.normal();
        }
        d =d/adj_source_f.size();
        d.normalize();
        math::Point t = n.getPoint()+average_dist*d;
        source_node_2_target_point[s_id] = t;
        //we store the point
        target_points.push_back(t);
        target_point_2_source_node[target_points.size()-1]=s_id;
    }
    //==================================================================
    // 4) We check if we can shrink the loop
    //==================================================================
    bool can_shrink = true;


    //We get the quasi-planar opposite faces
    std::vector<TCellID> all_target_faces;
    getEnclosedFaces(AL, !AToRight, all_target_faces);

    //among all target faces, we keep those living close to the plane
    // defined by the target side nodes
    std::vector<math::Point> target_side_pnts;
    target_side_pnts.reserve(to_circle_nodes.size());
    for(auto id:to_circle_nodes){
        math::Point pi = m_mesh.get<Node>(id).getPoint();
        target_side_pnts.push_back(pi);
    }

    math::Point pl_pnt;
    math::Vector3d pl_dir;
    math::computeLeastSquarePlane(target_side_pnts, pl_pnt, pl_dir);
    math::Plane target_pl(pl_pnt,math::Vector(pl_dir.X(),
                                              pl_dir.Y(),
                                              pl_dir.Z()));

    double plane_tolerance = average_dist/10;
    double pnt_tolerance = average_dist/100;
    for(auto f_id=0; can_shrink && f_id<all_target_faces.size(); f_id++){
        Face f = m_mesh.get<Face>(all_target_faces[f_id]);
        std::vector<Node> n = f.get<Node>();
        bool on_plane = true;
        for(auto i_n=0; i_n<n.size(); i_n++){
            Node ni = n[i_n];
            if(target_pl.distance(ni.getPoint())>plane_tolerance){
                on_plane=false;
            }
        }
        if(on_plane){
//            std::cout<<"Potential target face on plane: "<<f.getID()<<std::endl;
            //means we have a potential target face, either it matches
            //and we reuse its node, or not and we can not shrink
            //We look for matching points;

            //but first we disclaimed the faces that do not lies in
            //target nodes polygon


            bool keep_going = isInPolygon(f.getID(), target_side_pnts, target_pl);


            if(keep_going){
                std::vector<TCellID> receed_source_nodes;
                receed_source_nodes.resize(n.size());
                for(auto i_n=0; i_n<n.size()&& can_shrink; i_n++){
                    Node ni = n[i_n];
                    math::Point pi = ni.getPoint();
                    bool not_found = true;
                    for(auto j=0; j<target_points.size() &&
                        not_found && can_shrink; j++){
                        if(target_points[j].distance(pi)<pnt_tolerance){
                            receed_source_nodes[i_n]=target_point_2_source_node[j];
                            not_found=false;
                        }
                    }
                    if(not_found){
                        std::cout<<"\t some point not found"<<std::endl;
                        can_shrink=false;
                    }
                }
                //We have found the preimage of all the points
                //Now we check if a face exists based on this points
                if(!existFace(source_faces, receed_source_nodes)){
                    can_shrink=false;
                    std::cout<<"\t no matching face"<<std::endl;

                }
//                else{
//                    std::cout<<"\t found a matching face"<<std::endl;
//                }
            }
        }//if onplane
    }

    //we use node_mapping

    //For each point we want to create, we look if is not too close
    //to faces which are on the other side of the shrinking loop
    //If
    for(auto n2p: source_node_2_target_point){
        math::Point t = n2p.second;
    }

    if(can_shrink==false){
        //We stop the process
        // Clean the source mark
       // std::cout<<"Chord "<<AL<<" is discarded"<<std::endl;
        for(auto s_id:source_nodes){
            m_mesh.unmark<Node>(s_id,mark_source);
        }
        m_mesh.freeMark<Node>(mark_source);
        return false;
    }

    for(auto n2p: source_node_2_target_point){
        math::Point t = n2p.second;
        Node target_node = m_mesh.newNode(t);
        node_mapping[n2p.first]=target_node.getID();
    }
    //==================================================================
    // 5) We build the new layer.
    //==================================================================
    for(auto id_s:source_faces){
        Face s_face = m_mesh.get<Face>(id_s);
        std::vector<Node> s_nodes = s_face.get<Node>();
        std::vector<Face> l_faces; //lateral faces
        // LATERAL FACES
        for(auto i_n=0; i_n<s_nodes.size(); i_n++){
            Node ni = s_nodes[i_n];
            Node nj = s_nodes[(i_n+1)%s_nodes.size()];
            std::pair<TCellID,TCellID> p_ij;
            if(ni.getID()<nj.getID()){
                p_ij.first  = ni.getID();
                p_ij.second = nj.getID();
            }
            else{
                p_ij.first  = nj.getID();
                p_ij.second = ni.getID();
            }
            Face fij;
            if(n2face_mapping.find(p_ij)==n2face_mapping.end()){
                //new face to create;
                Node n1 = m_mesh.get<Node>(p_ij.first);
                Node n2 = m_mesh.get<Node>(p_ij.second);
                Node n3 = m_mesh.get<Node>(node_mapping[p_ij.second]);
                Node n4 = m_mesh.get<Node>(node_mapping[p_ij.first]);
                fij =m_mesh.newQuad(n1,n2,n3,n4);
                n1.add(fij);
                n2.add(fij);
                n3.add(fij);
                n4.add(fij);

                n2face_mapping[p_ij]=fij.getID();
            }
            else{
                fij = m_mesh.get<Face>(n2face_mapping[p_ij]);
            }
            l_faces.push_back(fij);

        }

        // TARGET FACE
        Face t;
        if(s_face.getType()==GMDS_QUAD){
            TCellID i1 = node_mapping[s_nodes[0].getID()];
            TCellID i2 = node_mapping[s_nodes[1].getID()];
            TCellID i3 = node_mapping[s_nodes[2].getID()];
            TCellID i4 = node_mapping[s_nodes[3].getID()];

            t=m_mesh.newQuad(i1,i2,i3,i4);

            Node n[4] = {
                m_mesh.get<Node>(i1),
                m_mesh.get<Node>(i2),
                m_mesh.get<Node>(i3),
                m_mesh.get<Node>(i4)
            };
            for(auto ni:n){
                ni.add(t);
            }
            m_mesh.newHex(s_nodes[0].getID(),
                          s_nodes[1].getID(),
                          s_nodes[2].getID(),
                          s_nodes[3].getID(),
                          i1,i2,i3,i4);

        }
        else if(s_face.getType()==GMDS_TRIANGLE){
            TCellID i1 = node_mapping[s_nodes[0].getID()];
            TCellID i2 = node_mapping[s_nodes[1].getID()];
            TCellID i3 = node_mapping[s_nodes[2].getID()];

            t=m_mesh.newTriangle(i1,i2,i3);

            Node n[3] = {
                m_mesh.get<Node>(i1),
                m_mesh.get<Node>(i2),
                m_mesh.get<Node>(i3)
            };
            for(auto ni:n){
                ni.add(t);
            }
            m_mesh.newPrism3(s_nodes[0].getID(),
                             s_nodes[1].getID(),
                             s_nodes[2].getID(),
                             i1,i2,i3);
        }
        else{
            throw GMDSException("Only Q and T source faces are managed");
        }
        m_mesh.mark(t,m_front_mark);
        m_mesh.unmark(s_face,m_front_mark);

    }

    for(auto i=0; i<AL.faces.size(); i++){
        m_mesh.unmark<Face>(AL.faces[i],m_front_mark);
    }
    //==================================================================
    // Clean the source mark
    for(auto s_id:source_nodes){
        m_mesh.unmark<Node>(s_id,mark_source);
    }
    m_mesh.freeMark<Node>(mark_source);

    return true;
}

/*---------------------------------------------------------------------------*/
bool WhiskerWeaving:: isInPolygon(const gmds::TCellID AF,
                                  const std::vector<gmds::math::Point>& APolyPnts,
                                  const gmds::math::Plane& APolyPl)
{
    gmds::Face f = m_mesh.get<Face>(AF);
    std::vector<Node> n_f = f.get<Node>();
    double min_length = 100000000;
    for(auto i=0;i<n_f.size();i++){
        math::Point pi = n_f[i].getPoint();
        math::Point pj = n_f[(i+1)%n_f.size()].getPoint();
        double dij = pi.distance(pj);
        if(dij<min_length)
            min_length=dij;
    }
    double tol=0.05*min_length;

    //intersection require to have point in the right plane
    std::vector<gmds::math::Point> polygon;
    polygon.reserve(APolyPnts.size());
    for(auto p:APolyPnts){
        polygon.push_back(APolyPl.project(p));
    }
    math::Point center = math::Point::massCenter(APolyPnts);

    math::Vector v[2];
    v[0] =  math::Vector(center, 0.5*(APolyPnts[0]+ APolyPnts[1]));
    v[1] = v[0].cross(APolyPl.getNormal());

    math::Point inf_pnt[2] = {
        APolyPl.project(center+1000*v[0]),
        APolyPl.project(center+1000*v[1])
    };


    //for each point of f, we check if it is in or out of AF
    for(auto n:n_f){

        math::Point p = n.getPoint();
        bool ok = false;
        for(auto i=0; i<APolyPnts.size() && !ok; i++){
            math::Point pi = APolyPnts[i];
            if(p.distance(pi)<tol){
                ok=true;
            }
        }
        if(ok){
//            std::cout<<"\t close to polygon bnd"<<std::endl;
            continue;
        }
        //We check inside or outside

        math::Point proj_p = APolyPl.project(p);
        int nb_cut[2]={0,0};
        for(auto k=0;k<2;k++){
            math::Segment w(proj_p,inf_pnt[k]);
            for(auto i=0; i<polygon.size() ; i++){
                math::Point pi = polygon[i];
                math::Point pj = polygon[(i+1)%polygon.size()];
                math::Segment sij(pi,pj);

                math::Point intersection_pnt;
                double param_ij, param_w;
                if(sij.intersect3D(w, intersection_pnt, param_w,param_ij )&&
                   param_ij>0.01){
                    nb_cut[k]++;
                }
            }
        }

        if(nb_cut[0]%2==1 || nb_cut[1]%2==1)
            ok = true;

        if(!ok){
            return false;
        }


    }

    return true;
}
/*---------------------------------------------------------------------------*/
bool WhiskerWeaving::existFace(const std::vector<gmds::TCellID>& AF,
                               const std::vector<gmds::TCellID>& AN)
{
    for(auto i_f:AF){
        Face f = m_mesh.get<Face>(i_f);
        std::vector<TCellID> nf = f.getIDs<Node>();
        bool found = false;
        for(auto i:nf){
            for(auto j:AN){
                if(i==j){
                    found =true;
                }
            }
        }
        if (!found)
            return false;
    }
    return true;

}
/*---------------------------------------------------------------------------*/
bool WhiskerWeaving::buildLoop(const Face& AQ,
                               const TCellID& AN1,
                               const TCellID& AN2,
                               std::map<TCellID,int>& AC1,
                               std::map<TCellID,int>& AC2,
                               Loop& ALoop)
{
    TCellID end_edge[2]={AN1,AN2};
    //current face
    Face cur = AQ;
    //next edge
    TCellID cur_e[2]={AN1,AN2};
    bool done = false;
    do{
        //if the current face is not a quad, we stop the loop process
        if(cur.getType()!=GMDS_QUAD)
            return false;

        //if the current face share two edges with the same other
        //we do not use it
        std::vector<TCellID> adj_cur = getAdjFaces(cur);
        std::set<TCellID> set_adj_cur;
        set_adj_cur.insert(adj_cur.begin(),adj_cur.end());
        if(set_adj_cur.size()!=adj_cur.size())
            return false;

        //First we color cur according to cur_ed position into its local
        //node numbering (color 1 or 2)
        int direction=0;
        std::vector<TCellID> cur_nodes = cur.getIDs<Node>();

        if((cur_e[0]==cur_nodes[0] && cur_e[1]==cur_nodes[1])||
           (cur_e[0]==cur_nodes[1] && cur_e[1]==cur_nodes[0])){
            direction=1;
            cur_e[0]=cur_nodes[2];
            cur_e[1]=cur_nodes[3];
        }
        else if((cur_e[0]==cur_nodes[2] && cur_e[1]==cur_nodes[3])||
                (cur_e[0]==cur_nodes[3] && cur_e[1]==cur_nodes[2])){
            direction=1;
            cur_e[0]=cur_nodes[0];
            cur_e[1]=cur_nodes[1];
        }
        else if((cur_e[0]==cur_nodes[1] && cur_e[1]==cur_nodes[2])||
                (cur_e[0]==cur_nodes[2] && cur_e[1]==cur_nodes[1])){
            direction=2;
            cur_e[0]=cur_nodes[0];
            cur_e[1]=cur_nodes[3];
        }
        else {
            direction=2;
            cur_e[0]=cur_nodes[1];
            cur_e[1]=cur_nodes[2];
        }

        if(direction==1)
            AC1[cur.getID()]=1;
        else if(direction==2)
            AC2[cur.getID()]=1;
        else{
            throw GMDSException("Color issue for building a front loop");
        }

        ALoop.faces.push_back(cur.getID());

        TCellID next_face_id = getFace(cur,cur_e[0], cur_e[1]);
        cur = m_mesh.get<Face>(next_face_id);

        if(cur_e[0]==end_edge[0] && cur_e[1]==end_edge[1])
            done=true;
        else if(cur_e[0]==end_edge[1] && cur_e[1]==end_edge[0])
            done=true;
    }
    while(!done);

    return true;
}
/*---------------------------------------------------------------------------*/
void WhiskerWeaving::setDebugFile(const std::string& AName)
{
    m_debug_file = AName;
}
/*---------------------------------------------------------------------------*/
void WhiskerWeaving::reverse(Face& AF, TCellID AFrom, TCellID ATo)
{
    std::vector<TCellID> n = AF.getIDs<Node>();

    int i_from=-1;
    for(auto i=0; i<n.size();i++){
        if(n[i]==AFrom)
            i_from=i;
    }
    int i_next = (i_from+1)%n.size();

    if(n[i_next]!=ATo){
        std::reverse(n.begin(),n.end());
        AF.set<Node>(n);
    }
}
/*---------------------------------------------------------------------------*/
bool WhiskerWeaving::same(const Loop& AL1, const Loop& AL2)
{
    if(AL1.faces.size()!=AL2.faces.size()){
        return false;
    }

    std::list<TCellID> l1, l2;
    l1.insert(l1.end(), AL1.faces.begin(), AL1.faces.end());
    l2.insert(l2.end(), AL2.faces.begin(), AL2.faces.end());
    l1.sort();
    l2.sort();
    std::list<TCellID>::iterator it1=l1.begin(),it2=l2.begin();
    for(; it1!=l1.end() && it2!=l2.end(); it1++,it2++){
        if(*it1!=*it2)
            return false;
    }
    return true;
}
/*---------------------------------------------------------------------------*/
TCellID WhiskerWeaving::getFace(const Face& AF,
                                const TCellID AN1,
                                const TCellID AN2)
{
    Node n1 = m_mesh.get<Node>(AN1);
    Node n2 = m_mesh.get<Node>(AN2);
    std::vector<TCellID> f1 = n1.getIDs<Face>();
    std::vector<TCellID> f2 = n2.getIDs<Face>();

    for(auto i1:f1){
        for(auto i2:f2){
            if(i1==i2 && i1!=AF.getID() &&
               m_mesh.isMarked<Face>(i1,m_front_mark)){
                return i1;
            }
        }
    }

    return NullID;
}

/*---------------------------------------------------------------------------*/
std::vector<TCellID> WhiskerWeaving::getAdjFaces(const Face& AF)
{
    std::vector<TCellID> n = AF.getIDs<Node>();
    std::vector<TCellID> adj_f;

    for(auto i_n=0; i_n<n.size(); i_n++){
        TCellID ni = n[i_n];
        TCellID nj = n[(i_n+1)%n.size()];

        TCellID fij = getFace(AF,ni,nj);
        adj_f.push_back(fij);
    }
    return adj_f;
}

/*---------------------------------------------------------------------------*/
void WhiskerWeaving::getOtherNodesFace(const Face& AF,
                                       const TCellID AN1, const TCellID AN2,
                                       TCellID& AN3, TCellID& AN4)
{
    std::vector<TCellID> n = AF.getIDs<Node>();
    std::vector<TCellID> others;
    for(auto i:n){
        if(i!=AN1 && i!= AN2){
            others.push_back(i);
        }
    }
    if(others.size()!=2)
        throw GMDSException("Error in WhiskerWeaving::getOtherNodesFace");

    AN3 = others[0];
    AN4 = others[1];
}
/*---------------------------------------------------------------------------*/
double WhiskerWeaving::dihedralAngle(const Face& AF1,const Face& AF2)
{
    math::Vector a = AF1.normal();
    math::Vector b = AF2.normal();

    //compute the angle
    double ab =a.dot(b);
    math::Vector axb=a.cross(b);
    double angle = std::atan2(axb.norm(), ab);

    //check the right orientation
    TCellID n1=NullID,n2=NullID;
    commonNodes(AF1,AF2,n1,n2);
    math::Vector v2(0.5*(m_mesh.get<Node>(n1).getPoint()+
                         m_mesh.get<Node>(n2).getPoint()),
                    AF2.center());

    if(a.dot(v2)<0){
        angle = -angle;
    }

    double angle_deg =angle*math::Constants::INVPIDIV180;
    return angle_deg;
}
/*---------------------------------------------------------------------------*/
void WhiskerWeaving::commonNodes(const gmds::Face& AF1,
                                 const gmds::Face& AF2,
                                 TCellID& AN1, TCellID& AN2)
{
    std::vector<TCellID> common;
    std::vector<TCellID> n1= AF1.getIDs<Node>();
    std::vector<TCellID> n2= AF2.getIDs<Node>();
    for(auto i1:n1){
        for(auto i2:n2){
            if(i1==i2 ){
                common.push_back(i1);
            }
        }
    }

    if(common.size()!=2)
        throw GMDSException("Two adjacent faces should share 2 nodes!");

    AN1 = common[0];
    AN2 = common[1];
}

/*---------------------------------------------------------------------------*/
void WhiskerWeaving::writeDebugMesh()
{
    static int i=0;
    std::string file_name = m_debug_file+std::to_string(i);


    Variable<math::Vector>* inward = 0;
    try{
        inward = m_mesh.getVariable<math::Vector>(GMDS_FACE, "inward");
    }
    catch(GMDSException&){
        inward= m_mesh.newVariable<math::Vector>(GMDS_FACE, "inward");
    }

    Variable<int>* front = 0;
    try{
        front = m_mesh.getVariable<int>(GMDS_FACE, "front");
    }
    catch(GMDSException&){
        front= m_mesh.newVariable<int>(GMDS_FACE, "front");
    }

    for(IGMesh::face_iterator it_f = m_mesh.faces_begin();
        !it_f.isDone(); it_f.next()){
        Face f = it_f.value();
        math::Vector n = f.normal();
        (*inward)[f.getID()] = n;
        (*front )[f.getID()] = (m_mesh.isMarked(f,m_front_mark))?1:0;

    }

    VTKWriter<IGMesh> writer(m_mesh);
    writer.write(file_name, R|F);
    i++;
}
/*---------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& AS,
                         const WhiskerWeaving::Loop& AL)
{
    AS<<"Loop [";
    for(auto f:AL.faces){
        AS<<f<<" ";
    }
    AS<<"]\n";
    AS<<"\t Flat index: "<<AL.flat_index<<"\n";
    AS<<"\t Right (sharp, nb_faces, nb_tri): ("<<AL.right_sharp<<", "
    <<AL.right_nb_faces<<", "<<AL.right_nb_triangles<<")\n";
    AS<<"\t Left  (sharp, nb_faces, nb_tri): ("<<AL.left_sharp<<", "
    <<AL.left_nb_faces<<", "<<AL.left_nb_triangles<<")\n";
    return AS;
}
/*---------------------------------------------------------------------------*/
