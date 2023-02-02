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
/*---------------------------------------------------------------------------*/
// STL File Headers
#include <random>
#include <set>
/*---------------------------------------------------------------------------*/
// FRAME File Headers
#include "TetMeshManipulator.h"
/*---------------------------------------------------------------------------*/
using namespace gmds;
using namespace fhedo;
/*---------------------------------------------------------------------------*/
/** Local node numbering so that the tetrahedron formed with:
*   - local node  i
*   - m_local_tet_node2facet[i][0]
*   - m_local_tet_node2facet[i][1]
*   - m_local_tet_node2facet[i][2]
*   has the same orientation as the original tetrahedron for any vertex i.
*/
const int TetMeshManipulator::m_local_tet_node2facet[4][3] = {
    {1, 2, 3},
    {0, 3, 2},
    {3, 0, 1},
    {1, 0, 2}
};
/*---------------------------------------------------------------------------*/
TetMeshManipulator::TetMeshManipulator(IGMesh* AMesh)
{;}
/*---------------------------------------------------------------------------*/
TetMeshManipulator::~TetMeshManipulator()
{;}
/*---------------------------------------------------------------------------*/
void TetMeshManipulator::insertLoopOnSurface(IGMesh* ATriMesh,
                                             std::vector<math::Point>& ALoop)
{
    //m is a triangular surface mesh in which we insert the loop
    IGMesh* m = ATriMesh;
    
    
    

}
/*---------------------------------------------------------------------------*/
void TetMeshManipulator::insertPoints(IGMesh* AMesh,
                                      std::vector<math::Point>& APoints,
                                      const std::vector<int>& AClassification)
{
    bool use_tetgen = false;
    
    if(use_tetgen) {
        insertPointsTetgen(AMesh,
                           APoints,
                           AMesh);
    }
    else{
        GEO::PCK::initialize();
        m_mesh=AMesh;
        
        //=====================================================================
        // (1) Tet reordering
        //=====================================================================
        // First step is done to reorder all tetrahedral elements in such a way
        // that facet normals are all going outside
        correctTetOrientation();
        
        if(!isDelaunay())
            std::cout<<"WARNING: NOT DELAUNAY MESH"<<std::endl;
        
        //=====================================================================
        // (2) Outside tet are put outside to avoid boundary specific cases
        //=====================================================================
        buildOutTets();
        //=====================================================================
        // (0) Initializaiton of the cavity mark
        //=====================================================================
        m_cavity = m_mesh->newVariable<TCellID>(GMDS_REGION, "cavity");
        
        //=====================================================================
        // (3) Loop on the point to insert each of them
        //=====================================================================
        Region seed = AMesh->regions_begin().value();
        int nb_fails=0;
        int i_insert=0;
        Region next_r = this->seed();
        for(auto i=0; i<APoints.size(); i++){
            math::Point p= APoints[i];
            int c = AClassification[i];
            //        if(i_insert==5)
            //            exit(0);
            if(c==0 )
                continue;
            
            i_insert++;
            for(IGMesh::region_iterator itr = m_mesh->regions_begin(); !itr.isDone();
                itr.next()){
                (*m_cavity)[itr.value().getID()]=NullID;
            }
            std::cout<<"================================"<<std::endl;
            std::cout<<"Insert point "<<p<<" ("<<c<<")"<<std::endl;
            Region r = next_r;
            if(!insert(p,r,c, next_r))
                nb_fails++;
            
            //       correctTetOrientation();
            writeTetMesh();
            
        }
        std::cout<<"NB FAILS: "<<nb_fails<<std::endl;
        
        
        //=====================================================================
        // (4) Outside tet are now deleted
        //=====================================================================
        deleteOutTets();
        
        writeTetMesh();
        
        m_mesh->deleteVariable(GMDS_REGION, m_cavity);
        
    }
}

/*---------------------------------------------------------------------------*/
bool TetMeshManipulator::isDelaunay()
{
    for(IGMesh::region_iterator itr = m_mesh->regions_begin();
        !itr.isDone(); itr.next()){
        Region ri = itr.value();
        std::vector<Node> n = ri.get<Node>();
        for(int i=0;i<4;i++){
            Region r = getOppositeRegion(ri,
                                         n[(i+1)%4],
                                         n[(i+2)%4],
                                         n[(i+3)%4]);
            if(r.getID()==NullID || r.getID()==InfinityID )
                continue;
            std::vector<Node> nr = r.get<Node>();
            Node w;
            for(auto ni:nr){
                if(ni.getID()!=n[(i+1)%4].getID() &&
                   ni.getID()!=n[(i+2)%4].getID() &&
                   ni.getID()!=n[(i+3)%4].getID() )
                    w=ni;
            }
            if(isInConflict(w.getPoint(), ri))
                return false;
        }
    }
    
    return true;
}
/*---------------------------------------------------------------------------*/
void TetMeshManipulator::writeCavity(std::vector<Region>& ACav,
                                     std::vector<Region>& ACavBnd)
{
    static int i=0;
    std::string f_cav = "cavity_"+std::to_string(i);
    std::string f_bnd = "cavity_bnd_"+std::to_string(i);
    IGMesh mesh_cav(MeshModel(DIM3 | R | N | R2N));
    IGMesh mesh_bnd(MeshModel(DIM3 | R | N | R2N));
    std::map<TCellID,TCellID> n_map_cav;
    std::map<TCellID,TCellID> n_map_bnd;
    
    std::set<TCellID> n_in_cav, n_in_bnd;
    for(auto r:ACav){
        if(isInfinityTet(r))
            continue;
        
        std::vector<TCellID> nr = r.getIDs<Node>();
        for(auto i:nr)
            n_in_cav.insert(i);
    }
    
    for(auto r:ACavBnd){
        if(isInfinityTet(r))
            continue;
        
        std::vector<TCellID> nr = r.getIDs<Node>();
        for(auto i:nr)
            n_in_bnd.insert(i);
    }
    
    Variable<int>* cav = mesh_cav.newVariable<int>(GMDS_REGION, "CAV_ID");
    
    for(auto n:n_in_cav){
        Node c =mesh_cav.newNode(m_mesh->get<Node>(n).getPoint());
        n_map_cav[n]=c.getID();
    }
    for(auto n:n_in_bnd){
        Node c =mesh_bnd.newNode(m_mesh->get<Node>(n).getPoint());
        n_map_bnd[n]=c.getID();
    }
    
    for(auto r:ACav){
        std::vector<TCellID> nr = r.getIDs<Node>();
        if(isInfinityTet(r))
            continue;
        Region c = mesh_cav.newTet(n_map_cav[nr[0]],
                                   n_map_cav[nr[1]],
                                   n_map_cav[nr[2]],
                                   n_map_cav[nr[3]]);
        (*cav)[c.getID()]=r.getID();
    }
    
    for(auto r:ACavBnd){
        std::vector<TCellID> nr = r.getIDs<Node>();
        if(isInfinityTet(r))
            continue;
        mesh_bnd.newTet(n_map_bnd[nr[0]],
                        n_map_bnd[nr[1]],
                        n_map_bnd[nr[2]],
                        n_map_bnd[nr[3]]);
    }
    
    VTKWriter<IGMesh> w_cav(mesh_cav);
    w_cav.write(f_cav, DIM3 | R | N);
    VTKWriter<IGMesh> w_bnd(mesh_bnd);
    w_bnd.write(f_bnd, DIM3 | R | N);
    i++;
    

}

/*---------------------------------------------------------------------------*/
void TetMeshManipulator::writeTetMesh()
{
    static int i=0;
    std::string fileN = "Tet_improvement_"+std::to_string(i);
    IGMesh tmesh(MeshModel(DIM3 | R | N | R2N));
    std::map<TCellID,TCellID> n_map;
    Variable<int>* vn = tmesh.newVariable<int>(GMDS_NODE, "REF_ID");
    Variable<int>* vr = tmesh.newVariable<int>(GMDS_REGION, "REF_ID");
    for(IGMesh::node_iterator itn = m_mesh->nodes_begin();
        !itn.isDone();itn.next()){
        Node from = itn.value();
        if(from.getID()==InfinityID)
            continue;
        
//        std::vector<Region> rn = from.get<Region>();
//        std::cout<<"Node "<<from.getID()<<": ";
//        for(auto reg:rn){
//            std::cout<<reg.getID()<<" ";
//        }
//        std::cout<<std::endl;
        Node to = tmesh.newNode(from.getPoint());
        n_map[from.getID()] = to.getID();
        (*vn)[to.getID()] = from.getID();
    }
    int nb_regions = 0, nb_real_regions=0;
    for(IGMesh::region_iterator itr = m_mesh->regions_begin();
        !itr.isDone();itr.next()){
        nb_regions++;
        Region from = itr.value();
        std::vector<TCellID> from_ids = from.getIDs<Node>();
//        std::cout<<"Region "<<from.getID()<<": "
//        <<from_ids[0]<<" "
//        <<from_ids[1]<<" "
//        <<from_ids[2]<<" "
//        <<from_ids[3]<<std::endl;
        
        if(isInfinityTet(from))
            continue;
        nb_real_regions++;
        
        Region to =tmesh.newTet(n_map[from_ids[0]],
                                n_map[from_ids[1]],
                                n_map[from_ids[2]],
                                n_map[from_ids[3]]);
        (*vr)[to.getID()] = from.getID();
    }
    std::cout<<"Nb regions: "<<nb_real_regions<<"/"<<nb_regions<<std::endl;
    VTKWriter<IGMesh> w(tmesh);
    w.write(fileN, DIM3 | R | N);
    i++;

}
/*---------------------------------------------------------------------------*/
void TetMeshManipulator::correctTetOrientation()
{
    auto nb_orient=0;
    for(IGMesh::region_iterator itr = m_mesh->regions_begin();
        !itr.isDone(); itr.next()){
        Region r = itr.value();
        if(isInfinityTet(r))
            continue;
        std::vector<Node> n = r.get<Node>();
        
        math::Point p[4]={
            n[0].getPoint(),
            n[1].getPoint(),
            n[2].getPoint(),
            n[3].getPoint()
        };
        if(orient3d(p[0],p[1],p[2],p[3])==GEO::NEGATIVE){
            std::vector<Node> new_n;
            new_n.resize(4);
            new_n[0]= n[0];
            new_n[1]= n[1];
            new_n[2]= n[3];
            new_n[3]= n[2];
            r.set(new_n);
            nb_orient++;
        }
        
    }
    std::cout<<"Nb re-oriented tets: "<<nb_orient<<std::endl;
}

/*---------------------------------------------------------------------------*/
void TetMeshManipulator::buildOutTets()
{

    for(IGMesh::region_iterator itr = m_mesh->regions_begin();
        !itr.isDone(); itr.next()){

        Region r = itr.value();
        if(isInfinityTet(r))
            continue;
        
        std::vector<Node> n = r.get<Node>();
        
        for(auto i=0; i<4; i++){
            
            Node n0 = n[m_local_tet_node2facet[i][0]];
            Node n1 = n[m_local_tet_node2facet[i][1]];
            Node n2 = n[m_local_tet_node2facet[i][2]];
            
            Region opp_r = getOppositeRegion(r,n0,n1,n2);
            if(opp_r.getID()==NullID){
                //means we have a boundary face, so we create an out tet
                //ordered in the right way
                Region out_tet = m_mesh->newTet(InfinityID,
                                                n0.getID(),
                                                n2.getID(),
                                                n1.getID());
                //and built N->R
                n0.add(out_tet);
                n1.add(out_tet);
                n2.add(out_tet);
//                std::cout<<"Build out tet from "<<n0.getID()<<", "
//                <<n1.getID()<<" and "<<n2.getID()<<std::endl;
            }
        }
        
    }
}
/*---------------------------------------------------------------------------*/
void TetMeshManipulator::deleteOutTets()
{
    std::set<Region> slivers;
    for(IGMesh::region_iterator itr = m_mesh->regions_begin();
        !itr.isDone(); itr.next()){
        Region r = itr.value();
        std::vector<TCellID> n = r.getIDs<Node>();
        if(n.empty())
            m_mesh->deleteRegion(r);
        else if(n[0]==InfinityID ||
           n[1]==InfinityID ||
           n[2]==InfinityID ||
           n[3]==InfinityID ){
            //We look is the bnd real region is a sliver or not
            //if YES, it will be deleted at the end
            facet f = realFacet(r);
            Node n0 = m_mesh->get<Node>(f(0));
            Node n1 = m_mesh->get<Node>(f(1));
            Node n2 = m_mesh->get<Node>(f(2));
            Region real_r = getOppositeRegion(r, n0, n1, n2);
            Node n_opp;
            std::vector<Node> real_nodes = real_r.get<Node>();
            auto found_opp = false;
            for(int i=0; i<4 && !found_opp; i++){
                if(real_nodes[i].getID()!=n0.getID() &&
                   real_nodes[i].getID()!=n1.getID() &&
                   real_nodes[i].getID()!=n2.getID()){
                    n_opp = real_nodes[i];
                    found_opp=true;
                }
            }
            if(!found_opp)
                throw GMDSException("Topological error in out tet deletion");
            
            math::Point p =  n_opp.getPoint();
            math::Plane pl(n0.getPoint(), n1.getPoint(), n2.getPoint());
            math::Point pr = pl.project(p);
            math::Vector3d v01(n0.getPoint(), n1.getPoint());
            math::Vector3d v02(n0.getPoint(), n2.getPoint());
            math::Vector3d v12(n1.getPoint(), n2.getPoint());
            
            double min_dist = math::min3(v01.norm(),v02.norm(),v12.norm());

            if(p.distance(pr)<0.1*min_dist)
                slivers.insert(real_r);
            for(auto ni:n){
                if(ni!=InfinityID){
                     m_mesh->get<Node>(ni).remove(r);
                }
            }
            m_mesh->deleteRegion(r);
        }
    }
    std::cout<<"Nb bnd slivers: "<<slivers.size()<<std::endl;
    for(auto s:slivers){
        std::vector<Node> n = s.get<Node>();
        for(auto ni:n){
            if(ni.getID()!=InfinityID){
                ni.remove(s);
            }
        }
        m_mesh->deleteRegion(s);

    }
    
}
/*---------------------------------------------------------------------------*/
bool TetMeshManipulator::insert(math::Point &APnt,
                                Region& ATet,
                                const int AC,
                                Region& ANewTet)
{
    //if a wrong (or no tet) is provided in input, we select one randomly.
    if(ATet.getID()==NullID)
        ATet=seed();
    
    //    std::cout<<"Start from seed "<<ATet.getID()<<std::endl;
    
    //=====================================================================
    // (1) We find the tet in which APnt must be inserted
    //=====================================================================
    Region r = findTet(APnt,ATet, AC);
    //    std::cout<<"then from real region "<<r.getID()<<std::endl;
    math::Vector4d coord(0,0,0,0);
    if(isIn( APnt, ATet, coord)){
        //        std::cout<<"Coord: "<<coord[0]<<", "
        //        <<coord[1]<<", "
        //        <<coord[2]<<", "
//        <<coord[3]<<std::endl;
        for(auto i=0;i<4;i++){
            if(coord[i]>0.8){
                std::vector<TCellID> nids = ATet.getIDs<Node>();
                m_mesh->get<Node>(nids[i]).setPoint(APnt);
                return true;
            }
        }
        
    }
    
    //=====================================================================
    // (2) We check where the point is in r (inside, on a facet, on an
    //     edge, on a node?)
    //=====================================================================
    std::vector<Region> cavity = initFuzzyCavity(APnt, AC, r);
//    std::cout<<"Init cavity: ";
//    for(auto r:cavity) {
//        std::cout<<r.getID()<<" ";
//    }
    std::cout<<std::endl;
    if(cavity.empty()){
        return false;
    }
    //=====================================================================
    // In order to handle each region only once, we use the m_cavity
    // variable to mark the regions belonging to the same cavity
    // All the region will be marked with the id of the first cavity region
    // should be done with the thread id in concurrent programming
    ///====================================================================
    
    // The variable will be cleaned up by the calling method
    TCellID cavity_id = cavity[0].getID();
    
    for(auto in_cav:cavity){
        (*m_cavity)[in_cav.getID()]=cavity_id;
    }
    
    std::vector<Region> cavity_bnd=cavity;
//    //=====================================================================
//    // Cavity is expanded to include all the tet in conflict with APnt
//    
//    // No more expansion in our case of inserting element in a non-delaunay
//    // triangulation
//    //=====================================================================
//    
//    expandCavityND(APnt, AC, cavity, cavity_bnd, cavity_id);
//    std::cout<<"Cavity: ";
//    for(auto r:cavity)
//        std::cout<<" "<<r.getID();
//    std::cout<<std::endl;
//    std::cout<<"Cavity Bnd: ";
//    for(auto r:cavity_bnd)
//        std::cout<<" "<<r.getID();
//    std::cout<<std::endl;
//    //=====================================================================
//    // We try to stellate the cavity
//    //=====================================================================
//    writeCavity(cavity,cavity_bnd);
    bool ok = stellateCavity(APnt,cavity_bnd, cavity_id, ANewTet);
    


    //We clean the cavity mark
    for(auto r_cav:cavity){
        (*m_cavity)[r_cav.getID()]=NullID;
    }
    //=====================================================================
    // If the cavity replacement succeeds, we remove old tetrahedral
    // elements
    //=====================================================================
    if(ok){
        for(auto tet:cavity){
            std::vector<TCellID> n = tet.getIDs<Node>();
            for(auto ni:n){
                if(ni!=InfinityID)
                    m_mesh->get<Node>(ni).remove(tet);
            }
            
            //now the initial tet can be removed
            m_mesh->deleteRegion(tet);
                     std::cout<<"We delete "<<tet.getID()<<std::endl;
        }
    }
    return ok;
}
/*---------------------------------------------------------------------------*/

std::vector<gmds::Region> TetMeshManipulator::initCavity(math::Point& APnt,
                                                         const int AC,
                                                         const Region& AR)
{
    std::vector<Region> cavity;
    
    if(isInfinityTet(AR)){
        throw GMDSException("A cavity never starts from an infinity tet");
    }    //================================================================
    // WE COMPUTE EXACT FACE BELONGING FIRST
    //================================================================
    std::vector<Node> r_nodes = AR.get<Node>();
    //on_facet[i] = true means that the point relies on the plan define
    // by the face opposite to the ith node of r.
    int on_facet[4] = {false, false, false, false};
    int nb_on_facets=0;
    for(auto i=0;i<4;i++){
        //m_local_tet_node2facet is used to get a consistant local
        //numbering (and so well oriented)
        Node n[3]= {
            r_nodes[m_local_tet_node2facet[i][0]],
            r_nodes[m_local_tet_node2facet[i][1]],
            r_nodes[m_local_tet_node2facet[i][2]]
        };
        on_facet[i]=(orient3d(APnt,n[0],n[1],n[2])==GEO::ZERO);
        if(on_facet[i])
            nb_on_facets++;
    }
    //    std::cout<<on_facet[0]<<" ";
    //    std::cout<<on_facet[1]<<" ";
    //    std::cout<<on_facet[2]<<" ";
    //    std::cout<<on_facet[3]<<std::endl;
    // 3 facets can contain APnt
    if(nb_on_facets>2){
        //nothing to do the point is already in the mesh
        cavity.clear();
        return cavity;
    }
    //================================================================
    // CASE 1: INTO THE VOLUME
    //================================================================
    if(AC==3){
        if(on_facet[0] && on_facet[1] ){ // EDGE 1
            //We insert on edge [2,3]
            cavity = getRegions(r_nodes[2], r_nodes[3]);
        }
        else if(on_facet[0] && on_facet[2] ){// EDGE 2
            //We insert on edge [1,3]
            cavity = getRegions(r_nodes[1], r_nodes[3]);
        }
        else if(on_facet[0] && on_facet[3]){// EDGE 3
            //We insert on edge [1,2]
            cavity = getRegions(r_nodes[1], r_nodes[2]);
        }
        else if(on_facet[1] && on_facet[2]){// EDGE 4
            //We insert on edge [0,3]
            cavity = getRegions(r_nodes[0], r_nodes[3]);
        }
        else if(on_facet[1] && on_facet[3]){// EDGE 5
            //We insert on edge [0,2]
            cavity = getRegions(r_nodes[0], r_nodes[2]);
        }
        else if(on_facet[2] && on_facet[3]){// EDGE 6
            //We insert on edge [1,2]
            cavity = getRegions(r_nodes[0], r_nodes[1]);
        }
        else if(on_facet[0]){// FACE 1
            //We insert on face [1,2,3]
            cavity = getRegions(r_nodes[1], r_nodes[2], r_nodes[3]);
        }
        else if(on_facet[1]){// FACE 2
            //We insert on face [0,2,3]
            cavity = getRegions(r_nodes[0], r_nodes[2], r_nodes[3]);
        }
        else if(on_facet[2]){// FACE 3
            //We insert on face [0,1,3]
            cavity = getRegions(r_nodes[0], r_nodes[1], r_nodes[3]);
        }
        else if(on_facet[3]){// FACE 4
            //We insert on face [0,1,2]
            cavity = getRegions(r_nodes[0], r_nodes[1], r_nodes[2]);
        }
        else{
            //only in the tetrahedron r
            cavity.push_back(AR);
        }
        
    }
    //================================================================
    // CASE 2: ON A SURFACE OR A CURVE
    //================================================================
    else if (AC==2 || AC==1){
        std::vector<Region> adj_regions;
        std::cout<<"Exact facet nb= "<<nb_on_facets<<std::endl;
        if(nb_on_facets==0){
            //We encountered some numerical issues due to the wrong
            //location of the point to insert (out of the domain)
            
            std::vector<Region> adj = getRegions(AR);
            std::vector<Region> fuzzy_adj;
            
            for(auto adj_r:adj){
                if(isInfinityTet(adj_r)){
                    fuzzy_adj.push_back(adj_r);
                }
            }
            
            //geometrical correction must be applied for coarse
            //meshes
            for(auto adj_r:fuzzy_adj){
                double tol = 0.01;
                facet fa = realFacet(adj_r);
                math::Triangle t (m_mesh->get<Node>(fa(0)).getPoint(),
                                  m_mesh->get<Node>(fa(1)).getPoint(),
                                  m_mesh->get<Node>(fa(2)).getPoint());
                double dist =t.distance(APnt);
                std::cout<<"dist= "<<dist<<std::endl;
                //     APnt = t.project(APnt);
                if(dist<tol){
                    adj_regions.push_back(adj_r);
                    nb_on_facets++;
                }
            }
            
        }
        else{
            for(auto i=0;i<4;i++){
                if(on_facet[i]){
                    Node n[3]= {
                        r_nodes[m_local_tet_node2facet[i][0]],
                        r_nodes[m_local_tet_node2facet[i][1]],
                        r_nodes[m_local_tet_node2facet[i][2]]
                    };
                    adj_regions.push_back(getOppositeRegion(AR,
                                                            n[0],
                                                            n[1],
                                                            n[2]));
                }
            }
        }
        std::cout<<"After correction= "<<nb_on_facets<<std::endl;
        if(nb_on_facets==1){
            std::cout<<"here 1"<<std::endl;
            //we keep the current face in the cavity
            cavity.push_back(AR);
            cavity.push_back(adj_regions[0]);
        }
        else if(nb_on_facets==2){
            std::cout<<"here 2"<<std::endl;
            //We get the edge common to the two  cells and
            //we put all the regions incident to this edge into the
            //cavity
            std::vector<TCellID> n0 = adj_regions[0].getIDs<Node>();
            std::vector<TCellID> n1 = adj_regions[1].getIDs<Node>();
            TCellID edge_ids[2];
            auto i_e=0;
            for(auto id0:n0){
                if(id0!=InfinityID){
                    for(auto id1:n1){
                        if(id0==id1){
                            edge_ids[i_e++]=id0;
                        }
                    }
                }
            }
            cavity = getRegions(m_mesh->get<Node>(edge_ids[0]),
                                m_mesh->get<Node>(edge_ids[1]));
            
        }
        else {
            //We can not conclude, we must so geometrically
            // find the best location
            //We have 3 or 4 ghost facet
            bool find_cavity =false;
            double tol_dist = 1e-8;
            while (!find_cavity) {
                std::vector<facet> ghost_facet;
                for(auto a:adj_regions){
                    math::Triangle t;
                    facet fa;
                    if(isInfinityTet(a)){
                        fa = realFacet(a);
                        t =math::Triangle(m_mesh->get<Node>(fa(0)).getPoint(),
                                          m_mesh->get<Node>(fa(1)).getPoint(),
                                          m_mesh->get<Node>(fa(2)).getPoint());
                    }
                    else{
                        std::vector<Node> common;
                        common.reserve(3);
                        //real tet
                        std::vector<Node> a_nodes = a.get<Node>();
                        for(auto n1:r_nodes){
                            for(auto n2:a_nodes){
                                if(n1.getID()==n2.getID())
                                    common.push_back(n1);
                            }
                        }
                        t =math::Triangle(common[0].getPoint(),
                                          common[1].getPoint(),
                                          common[2].getPoint());
                        fa=facet(common[0].getID(),
                                 common[1].getID(),
                                 common[2].getID());
                        
                    }
                    double dist =t.distance(APnt);
                    //     APnt = t.project(APnt);
                    if(dist<tol_dist){
                        ghost_facet.push_back(fa);
                        
                    }
                }
                
                
                if(ghost_facet.size()==1){
                    find_cavity=true;
                    //We put the unique real tet adj to the facet
                    //which is AR
                    cavity.push_back(AR);
                }
                else if(ghost_facet.size()==2){
                    find_cavity=true;
                    //We get the edge common to the two infinity cells and
                    //we put all the regions incident to this edge into the
                    //cavity
                    TCellID n0[3] = {
                        ghost_facet[0](0),
                        ghost_facet[0](1),
                        ghost_facet[0](2)
                    };
                    TCellID n1[3] = {
                        ghost_facet[1](0),
                        ghost_facet[1](1),
                        ghost_facet[1](2)
                    };
                    TCellID edge_ids[2];
                    auto i_e=0;
                    for(auto id0:n0){
                        if(id0!=InfinityID){
                            for(auto id1:n1){
                                if(id0==id1){
                                    edge_ids[i_e++]=id0;
                                }
                            }
                        }
                    }
                    cavity = getRegions(m_mesh->get<Node>(edge_ids[0]),
                                        m_mesh->get<Node>(edge_ids[1]));
                }
                else if(ghost_facet.size()==3){
                    find_cavity=true;
                    cavity.clear();
                }
                
                if(!find_cavity){
                    tol_dist *=10;
                }
            }
            
        }
        
    }
    return cavity;
}
/*---------------------------------------------------------------------------*/
std::vector<gmds::Region>
TetMeshManipulator::initFuzzyCavity(math::Point& APnt,
                                    const int AC,
                                    const Region& AR)
{
    std::vector<Region> cavity;
    
    if(isInfinityTet(AR)){
        throw GMDSException("A cavity never starts from an infinity tet");
    }
    
    math::Vector4d coord(0,0,0,0);
    std::vector<Node> r_nodes = AR.get<Node>();
    isIn(APnt, AR, coord);
    
    
    //        std::cout<<"Coord: "<<coord[0]<<", "
    //        <<coord[1]<<", "
    //        <<coord[2]<<", "
    //        <<coord[3]<<std::endl;
    std::vector<int> sup045, sup033, sup08;
    for(auto i=0;i<4;i++){
        if(coord[i]>0.8){
            sup08.push_back(i);
        }
        if(coord[i]>0.4){
            sup045.push_back(i);
        }
        if(coord[i]>0.3){
            sup033.push_back(i);
        }
    }
    
    bool on_node =false;
    bool on_edge =false;
    bool on_face =false;
    std::vector<Node> cell_nodes;
    cell_nodes.reserve(3);
    if(sup08.size()==1){
        //THE POINT IS VIRTUALLY MOVED ONTO A NODE
        on_node=true;
        cell_nodes.push_back(r_nodes[sup08[0]]);
    }
    else if(sup045.size()==2){
        //THE POINT IS VIRTUALLY MOVED ONTO AN EDGE
        on_edge=true;
        cell_nodes.push_back(r_nodes[sup045[0]]);
        cell_nodes.push_back(r_nodes[sup045[1]]);
    }
    if(sup033.size()==3){
        //THE POINT IS VIRTUALLY MOVED ONTO AN FACE
        on_face=true;
        cell_nodes.push_back(r_nodes[sup033[0]]);
        cell_nodes.push_back(r_nodes[sup033[1]]);
        cell_nodes.push_back(r_nodes[sup033[2]]);

    }
    //================================================================
    // CASE 1: INTO THE VOLUME
    //================================================================
    
    if(on_node){
        cavity = cell_nodes[0].get<Region>();
    }
    else if(on_edge){
        cavity = getRegions(cell_nodes[0], cell_nodes[1]);
        
    }
    else if (on_face){
        cavity = getRegions(cell_nodes[0], cell_nodes[1],cell_nodes[0]);
    }
    else{
        //only in the tetrahedron r
        cavity.push_back(AR);
    }
    return cavity;
}
/*---------------------------------------------------------------------------*/
TetMeshManipulator::facet TetMeshManipulator::realFacet(const Region& ATet)
{
    std::vector<TCellID> n = ATet.getIDs<Node>();
    for(auto i=0; i<4; i++){
        if (n[i]==InfinityID){
            return facet(n[m_local_tet_node2facet[i][0]],
                         n[m_local_tet_node2facet[i][1]],
                         n[m_local_tet_node2facet[i][2]]);
        }
    }
    throw GMDSException("realFacet must be called for ghost tets only");
}

/*---------------------------------------------------------------------------*/
bool TetMeshManipulator::isInfinityTet(const Region& ATet) const
{
    std::vector<TCellID> n = ATet.getIDs<Node>();
    for(auto i:n){
        if (i==InfinityID){
            return true;
        }
    }
    return false;
}

/*---------------------------------------------------------------------------*/
bool TetMeshManipulator::contains(math::Point& AP,
                                  const Region& AR,
                                  const int AC)
{
    //================================================================
    // INFINITY TETRAHEDRON
    //================================================================
    if(isInfinityTet(AR)){
        facet bnd_f = realFacet(AR);
        if(AC==3){
            //should not be in this configuration, which is encountered
            //during the propagation process since we work on non-convex
            //meshes
            return false;
        }
        // Two cases:
        //================================================================
        // (1) the point to be inserted is very close since it
        //      was just outside of the domain
        if(AC==2){
            // (1.1) On a surface
            math::Point new_p;
            bool is_in = computeProj(AP,
                                     m_mesh->get<Node>(bnd_f(0)).getPoint(),
                                     m_mesh->get<Node>(bnd_f(1)).getPoint(),
                                     m_mesh->get<Node>(bnd_f(2)).getPoint(),
                                     new_p);
            if(is_in && new_p.distance(AP)<0.001){
                //We do not start from a current region but only from
                //the corresponding real one
                return  true;
//                getOppositeRegion(current_tet,
//                                          m_mesh->get<Node>(bnd_f(0)),
//                                          m_mesh->get<Node>(bnd_f(1)),
//                                          m_mesh->get<Node>(bnd_f(2)));
            }
            
        }
        else if(AC==1){
            // (1.1) On a curve
            math::Point new_p;
            math::Point p[3] ={
                m_mesh->get<Node>(bnd_f(0)).getPoint(),
                m_mesh->get<Node>(bnd_f(1)).getPoint(),
                m_mesh->get<Node>(bnd_f(2)).getPoint()
            };
            math::Segment s1(p[0],p[1]);
            math::Segment s2(p[1],p[2]);
            math::Segment s3(p[0],p[2]);
            //                std::cout<<"Point of triangle: "<<p[0]<<std::endl;
            //                std::cout<<"Point of triangle: "<<p[1]<<std::endl;
            //                std::cout<<"Point of triangle: "<<p[2]<<std::endl;
            
            math::Point proj[3] = {
                s1.project(AP),
                s2.project(AP),
                s3.project(AP)
            };
            //                std::cout<<"Projected point : "<<proj[0]<<std::endl;
            //                std::cout<<"Projected point : "<<proj[1]<<std::endl;
            //                std::cout<<"Projected point : "<<proj[2]<<std::endl;
            
            for(auto p:proj){
                std::cout<<p.distance(AP)<<std::endl;
                if(p.distance(AP)<0.01){
                    //We do not start from a current region but only from
                    //the corresponding real one
                    return  true;
//                    getOppositeRegion(current_tet,
//                                              m_mesh->get<Node>(bnd_f(0)),
//                                              m_mesh->get<Node>(bnd_f(1)),
//                                              m_mesh->get<Node>(bnd_f(2)));
                }
            }//for(auto p:proj)
            
        }
        
        return false;
    }
    else{
        std::vector<Node> nr = AR.get<Node>();
        //================================================================
        // REGULAR TETRAHEDRON
        //================================================================
        // Start from a random node
        auto i0 = random(4);
        bool is_out=false;
        //for each point starting from the i0th in current_tet
        //we chech the orientation of the opposite facet
        for(auto i=0;i<4 && !is_out;i++){
            auto c = (i0+i)%4;
            //m_local_tet_node2facet is used to get a consistant local
            //numbering (and so well oriented)
            Node n[3]= {
                nr[m_local_tet_node2facet[c][0]],
                nr[m_local_tet_node2facet[c][1]],
                nr[m_local_tet_node2facet[c][2]]
            };
            if(orient3d(AP,n[0],n[1],n[2])==GEO::NEGATIVE){
                
                is_out=true;
            }
        }//for(auto i=0;i<4 && !f....
     
        return (!is_out);
    }
    return false;
}
/*---------------------------------------------------------------------------*/
Region TetMeshManipulator::findTet(math::Point& APnt,
                                   const Region& ASeed,
                                   const int AC)
{
    gmds::Region current_tet = ASeed;
    TCellID previous_tet_id = NullID;
    std::cout<<"Current tet= "<<current_tet.getID()<<" and "<<APnt<<std::endl;
    int attempt = 0;
    while(attempt<5){
        
        bool found_next_tet=false;
        //================================================================
        // INFINITY TETRAHEDRON
        //================================================================
        if(isInfinityTet(current_tet)){
            facet bnd_f = realFacet(current_tet);
            
            if(AC==3){
                //should not be in this configuration, which is encountered
                //during the propagation process since we work on non-convex
                //meshes
                current_tet = seed();
                attempt++;
                found_next_tet=true;
            }
            // Two cases:
            //================================================================
            // (1) the point to be inserted is very close since it
            //      was just outside of the domain
            if(AC==2){
                // (1.1) On a surface
                math::Point new_p;
                bool is_in = computeProj(APnt,
                                         m_mesh->get<Node>(bnd_f(0)).getPoint(),
                                         m_mesh->get<Node>(bnd_f(1)).getPoint(),
                                         m_mesh->get<Node>(bnd_f(2)).getPoint(),
                                         new_p);
                if(is_in && new_p.distance(APnt)<0.001){
                    //We do not start from a current region but only from
                    //the corresponding real one
                    return  getOppositeRegion(current_tet,
                                              m_mesh->get<Node>(bnd_f(0)),
                                              m_mesh->get<Node>(bnd_f(1)),
                                              m_mesh->get<Node>(bnd_f(2)));
                }
                
            }
            else if(AC==1){
                // (1.1) On a curve
                math::Point new_p;
                math::Point p[3] ={
                    m_mesh->get<Node>(bnd_f(0)).getPoint(),
                    m_mesh->get<Node>(bnd_f(1)).getPoint(),
                    m_mesh->get<Node>(bnd_f(2)).getPoint()
                };
                math::Segment s1(p[0],p[1]);
                math::Segment s2(p[1],p[2]);
                math::Segment s3(p[0],p[2]);
                //                std::cout<<"Point of triangle: "<<p[0]<<std::endl;
                //                std::cout<<"Point of triangle: "<<p[1]<<std::endl;
                //                std::cout<<"Point of triangle: "<<p[2]<<std::endl;
                
                math::Point proj[3] = {
                    s1.project(APnt),
                    s2.project(APnt),
                    s3.project(APnt)
                };
                //                std::cout<<"Projected point : "<<proj[0]<<std::endl;
                //                std::cout<<"Projected point : "<<proj[1]<<std::endl;
                //                std::cout<<"Projected point : "<<proj[2]<<std::endl;
                
                for(auto p:proj){
                    std::cout<<p.distance(APnt)<<std::endl;
                    if(p.distance(APnt)<0.01){
                        //We do not start from a current region but only from
                        //the corresponding real one
                        return  getOppositeRegion(current_tet,
                                                  m_mesh->get<Node>(bnd_f(0)),
                                                  m_mesh->get<Node>(bnd_f(1)),
                                                  m_mesh->get<Node>(bnd_f(2)));
                    }
                }
                
            }
            //================================================================
            
            // (2) a new seed must be throw since the non-convexity
            //     of the domain must the reason why we are in this
            //     configuration
            current_tet = seed();
            attempt++;
            found_next_tet=true;
            

//            if(AC==3){
//                //should not be in this configuration, which is encountered
//                //during the propagation process since we work on non-convex
//                //meshes
//                current_tet = seed();
//                attempt++;
//                found_next_tet=true;
//            }
//            // Two cases:
//            //================================================================
//            // (1) the point to be inserted is very close since it
//            //      was just outside of the domain
//            if(AC==2){
//                // (1.1) On a surface
//                math::Point new_p;
//                bool is_in = computeProj(APnt,
//                                         m_mesh->get<Node>(bnd_f(0)).getPoint(),
//                                         m_mesh->get<Node>(bnd_f(1)).getPoint(),
//                                         m_mesh->get<Node>(bnd_f(2)).getPoint(),
//                                         new_p);
//                if(is_in && new_p.distance(APnt)<0.001){
//                    //We do not start from a current region but only from
//                    //the corresponding real one
//                    return  getOppositeRegion(current_tet,
//                                              m_mesh->get<Node>(bnd_f(0)),
//                                              m_mesh->get<Node>(bnd_f(1)),
//                                              m_mesh->get<Node>(bnd_f(2)));
//                }
//                
//            }
//            else if(AC==1){
//                // (1.1) On a curve
//                math::Point new_p;
//                math::Point p[3] ={
//                    m_mesh->get<Node>(bnd_f(0)).getPoint(),
//                    m_mesh->get<Node>(bnd_f(1)).getPoint(),
//                    m_mesh->get<Node>(bnd_f(2)).getPoint()
//                };
//                math::Segment s1(p[0],p[1]);
//                math::Segment s2(p[1],p[2]);
//                math::Segment s3(p[0],p[2]);
////                std::cout<<"Point of triangle: "<<p[0]<<std::endl;
////                std::cout<<"Point of triangle: "<<p[1]<<std::endl;
////                std::cout<<"Point of triangle: "<<p[2]<<std::endl;
//                
//                math::Point proj[3] = {
//                    s1.project(APnt),
//                    s2.project(APnt),
//                    s3.project(APnt)
//                };
////                std::cout<<"Projected point : "<<proj[0]<<std::endl;
////                std::cout<<"Projected point : "<<proj[1]<<std::endl;
////                std::cout<<"Projected point : "<<proj[2]<<std::endl;
//
//                for(auto p:proj){
//                    std::cout<<p.distance(APnt)<<std::endl;
//                    if(p.distance(APnt)<0.01){
//                        //We do not start from a current region but only from
//                        //the corresponding real one
//                        return  getOppositeRegion(current_tet,
//                                                  m_mesh->get<Node>(bnd_f(0)),
//                                                  m_mesh->get<Node>(bnd_f(1)),
//                                                  m_mesh->get<Node>(bnd_f(2)));
//                    }
//                }
//
//            }
//            //================================================================
//
//            // (2) a new seed must be throw since the non-convexity
//            //     of the domain must the reason why we are in this
//            //     configuration
//            current_tet = seed();
//            attempt++;
//            found_next_tet=true;
//
            
        } //if(isInfinityTet(current_tet))
        else{
            std::vector<Node> nr = current_tet.get<Node>();
            //================================================================
            // REGULAR TETRAHEDRON
            //================================================================
            // Start from a random node
            auto i0 = random(4);
            
            //for each point starting from the i0th in current_tet
            //we chech the orientation of the opposite facet
            for(auto i=0;i<4 && !found_next_tet;i++){
                auto c = (i0+i)%4;
                //m_local_tet_node2facet is used to get a consistant local
                //numbering (and so well oriented)
                Node n[3]= {
                    nr[m_local_tet_node2facet[c][0]],
                    nr[m_local_tet_node2facet[c][1]],
                    nr[m_local_tet_node2facet[c][2]]
                };
                if(orient3d(APnt,n[0],n[1],n[2])==GEO::NEGATIVE){
                    
                    Region next_r = getOppositeRegion(current_tet,
                                                      n[0],n[1],n[2]);
                    std::cout<<"- from "<<current_tet<<" to "<<next_r<<std::endl;;
                    if(next_r.getID()==previous_tet_id)
                        return current_tet;
                    
                    previous_tet_id = current_tet.getID();
                    current_tet = next_r;
                    //jump to next tet
                    found_next_tet=true;
                    
                }
            }//for(auto i=0;i<4 && !f....
            
        }
        if(!found_next_tet){
            //We have found the tet
            return current_tet;
        }
    }
    
    //So we go through all the mesh
    for(IGMesh::region_iterator itr= m_mesh->regions_begin();
        !itr.isDone(); itr.next()){
        Region r = itr.value();
        if(contains(APnt,r,AC)){
            if(isInfinityTet(r)){
                facet bnd_f = realFacet(r);
                return  getOppositeRegion(r,
                                          m_mesh->get<Node>(bnd_f(0)),
                                          m_mesh->get<Node>(bnd_f(1)),
                                          m_mesh->get<Node>(bnd_f(2)));
            }
            else{
                return r;
            }
        }
        
    }
    throw GMDSException("Impossible to find a tet containing a point");
    
}
/*---------------------------------------------------------------------------*/
void TetMeshManipulator::expandCavity(const math::Point &APnt,
                                      std::vector<Region>& ACav,
                                      std::vector<Region>& ACavBnd,
                                      const int ACavID)
{
    std::set<TCellID> cavity_bnd_set;
    std::set<TCellID> cavity_set;
    //We list all the regions adjacent to a region of ACav and which is not
    // in ACav right yet
    std::vector<Region> grewing_cavity;
    for(auto r:ACav){
        cavity_set.insert(r.getID());
        if(isInfinityTet(r)){
            cavity_bnd_set.insert(r.getID());
        }
        else{
            std::vector<Region> adj_r = getRegions(r);
            
            for(auto a:adj_r){
                std::cout<<"Test "<<a.getID()<<std::endl;
                
                if(isInfinityTet(a))
                    std::cout<<"\t "<<a.getID()<<": infinity"<<std::endl;
                else
                    std::cout<<"\t "<<a.getID()<<": regular"<<std::endl;
                
                std::vector<Node> an = a.get<Node>();
                for(auto node: an){
                    std::cout<<"\t "<<node<<std::endl;
                }
                if((*m_cavity)[a.getID()]!=ACavID){
                    std::cout<<"\t not already in cav"<<std::endl;
                    if(isInConflict(APnt, a)){
                        std::cout<<"\t and in conflict"<<std::endl;
                        grewing_cavity.push_back(a);
                        (*m_cavity)[a.getID()]=ACavID;
                    }
                    else
                        cavity_bnd_set.insert(r.getID());
                }
            }
        }
    }
    std::cout<<"----- grewing now------"<<std::endl;
    while(!grewing_cavity.empty()){
        //We pick a region which is conflict by definition
        Region current = grewing_cavity.back();
        grewing_cavity.pop_back();
        cavity_set.insert(current.getID());
        if(isInfinityTet(current)){
            cavity_bnd_set.insert(current.getID());
        }
        else{
            std::vector<Region> adj_r = getRegions(current);
            for(auto a:adj_r){
                std::cout<<"Test "<<a.getID()<<std::endl;
                if(isInfinityTet(a))
                    std::cout<<"\t "<<a.getID()<<": infinity"<<std::endl;
                else
                    std::cout<<"\t "<<a.getID()<<": regular"<<std::endl;
                
                
                
                if((*m_cavity)[a.getID()]!=ACavID){
                    std::cout<<"\t not already in cav"<<std::endl;
                    
                    if(isInConflict(APnt, a)){
                        std::cout<<"\t and in conflict"<<std::endl;
                        grewing_cavity.push_back(a);
                        (*m_cavity)[a.getID()]=ACavID;
                    }
                    else
                        cavity_bnd_set.insert(current.getID());
                }
            }
        }
    }
    
    //We fill in the Cavity Boundary parameter
    ACavBnd.resize(cavity_bnd_set.size());
    auto i=0;
    for(auto c:cavity_bnd_set){
        ACavBnd[i++]=m_mesh->get<Region>(c);
    }
    
    //and we also update the cavity itself //to be improved;
    ACav.resize(cavity_set.size());
    i=0;
    for(auto c:cavity_set){
        ACav[i++]=m_mesh->get<Region>(c);
    }
}

/*---------------------------------------------------------------------------*/
void TetMeshManipulator::expandCavityND(const math::Point &APnt,
                                        const int AClass,
                                        std::vector<Region>& ACav,
                                        std::vector<Region>& ACavBnd,
                                        const int ACavID)
{
    std::set<TCellID> cavity_bnd_set;
    std::set<TCellID> cavity_set;
    //We list all the regions adjacent to a region of ACav and which is not
    // in ACav right yet
    for(auto r:ACav){
        cavity_set.insert(r.getID());
        std::vector<Region> adj_r = getRegions(r);
        std::cout<<"Current region "<<r.getID()<<std::endl;
        if(isInfinityTet(r)){
            std::cout<<"\t infinite, we do nothing"<<std::endl;
            cavity_bnd_set.insert(r.getID());
            
        }
        else{
            std::vector<Node> r_n = r.get<Node>();
            for(auto a:adj_r){
                std::cout<<"\t adj region "<<a.getID()<<std::endl;
                
                if(isInfinityTet(a))
                    std::cout<<"\t "<<a.getID()<<": infinity"<<std::endl;
                else
                    std::cout<<"\t "<<a.getID()<<": regular"<<std::endl;
                
                if((*m_cavity)[a.getID()]!=ACavID){
                    std::cout<<"\t not already in cav"<<std::endl;
                    bool is_inf = isInfinityTet(a);
                    bool added = false;
                    //We had it in cav depending on some tet quality
                    if(is_inf && isInConflict(APnt, a)){
                        std::cout<<"\t and in conflict"<<std::endl;
                        cavity_bnd_set.insert(a.getID());
                        cavity_set.insert(a.getID());
                        (*m_cavity)[a.getID()]=ACavID;
                        added=true;
                    }
                    else if(!is_inf){
                        //We compare the quality of tets before and after insertion
                        // 2 tets vs. 3 tets
                        std::vector<Node> a_n = a.get<Node>();
                        math::Point common[3];
                        Node common_node[3];
                        math::Point a_node;
                        int i_c=0;
                        for(auto ni:a_n){
                            if(ni!=r_n[0] && ni!=r_n[1] &&
                               ni!=r_n[2] && ni!=r_n[3]){
                                a_node=ni.getPoint();
                            }
                            else{
                                common[i_c]=ni.getPoint();
                                common_node[i_c++]=ni;

                            }
                        }
                        
                        double cur_qual=0;
                        double swap_qual=0;
                        
                        double q1 = quality(APnt, common[0], common[1], common[2]);
                        double q2 = quality(a_node, common[0], common[2], common[1]);
                        cur_qual= math::min2(q1,q2);
                        
                        
                        double q3 = quality(a_node, APnt, common[0], common[1]);
                        double q4 = quality(a_node, APnt, common[1], common[2]);
                        double q5 = quality(a_node, APnt, common[2], common[0]);
                        if(AClass==3){
                            swap_qual= math::min3(q3,q4,q5);
                        }
                        else if(AClass==1){
                           
                            // Indeed, 2 must be bnd slivers that we do not look at, so only
                            //one is relevant (the max of three of them)
                            swap_qual= math::max3(q3,q4,q5);
                        }
                        else if(AClass==2){
                            //One is flat, so we remove it
                            if(q5<q3 && q5<q4)
                                swap_qual= math::min2(q3,q4);
                            else if(q4<q3 && q4<q5)
                                swap_qual= math::min2(q3,q5);
                            else
                                swap_qual= math::min2(q4,q5);
                        }
                        //best quality with swap, we do it
                        if(swap_qual>cur_qual){
                            std::cout<<"\t and in conflict"<<std::endl;
                            cavity_bnd_set.insert(a.getID());
                            cavity_set.insert(a.getID());
                            
                            (*m_cavity)[a.getID()]=ACavID;
                            //moreover we add in the cavity the infinity region of a
                            //sharing and edge with r
                            std::vector<Region> adj_a = getRegions(a);
                            for(auto cur_adj:adj_a){
                                if(isInfinityTet(cur_adj)){
                                    facet f_real = realFacet(cur_adj);
                                    int nb_common=0;
                                    for(auto n: common_node){
                                        if(n.getID()==f_real(0) || n.getID()==f_real(1) ||
                                           n.getID()==f_real(2))
                                            nb_common++;
                                    }
                                    if(nb_common>1){
                                        cavity_bnd_set.insert(cur_adj.getID());
                                        cavity_set.insert(cur_adj.getID());
                                         (*m_cavity)[cur_adj.getID()]=ACavID;
                                    }
                                        
                                }
                            }
                            added=true;
                        }
                        
                    }
                    if(!added){
                        std::cout<<"\t so current region IN"<<std::endl;
                        cavity_bnd_set.insert(r.getID());
                    }
                }
            }
        }
    }
    
    //We fill in the Cavity Boundary parameter
    ACavBnd.resize(cavity_bnd_set.size());
    auto i=0;
    for(auto c:cavity_bnd_set){
        ACavBnd[i++]=m_mesh->get<Region>(c);
    }
    
    //and we also update the cavity itself //to be improved;
    ACav.resize(cavity_set.size());
    i=0;
    for(auto c:cavity_set){
        ACav[i++]=m_mesh->get<Region>(c);
    }
}
/*---------------------------------------------------------------------------*/
void TetMeshManipulator::checkAndCorrectCavity(const math::Point &APnt,
                                               std::vector<Region>& ACav)
{
    //  throw GMDSException("computeCavity not yet implemented");
}
/*---------------------------------------------------------------------------*/
bool TetMeshManipulator::isInConflict(const math::Point& APnt,
                                      Region& AR)
{
    std::vector<Node> n = AR.get<Node>();
    
    auto infinity_node = -1;
    for(auto i=0; i<4; i++){
        if(n[i].getID()==InfinityID){
            infinity_node=i;
        }
    }
    if(infinity_node!=-1){
//        //The circumcenter has barycentric coordinates
//        //a^2( b^2 + c^2 -a^2):\;b^2(c^2 + a^2-b^2):\;c^2(a^2 + b^2 - c^2)
//        //where a, b, c are edge lengths (BC, CA, AB respectively) of the triangle.
//        //see https://en.wikipedia.org/wiki/Circumscribed_circle#Barycentric_coordinates
//        
//        facet f = realFacet(AR);
//        math::Point A = m_mesh->get<Node>(f(0)).getPoint();
//        math::Point B = m_mesh->get<Node>(f(1)).getPoint();
//        math::Point C = m_mesh->get<Node>(f(2)).getPoint();
//        math::Vector AB(A,B);
//        math::Vector AC(A,C);
//        math::Vector BC(B,C);
//        double a = BC.norm();
//        double b = AC.norm();
//        double c = AB.norm();
//        double a2 = a*a;
//        double b2 = b*b;
//        double c2 = c*c;
//        math::Point center(a2*(b2+c2-a2)*A+
//                           b2*(c2+a2-b2)*B+
//                           c2*(a2+b2-c2)*C);
//        math::Vector normal = AB.cross(AC);
//        normal.normalize();
//        
//        double radius  = center.distance(A);
//        
//        math::Point D = center + radius*normal;
//        double co[4][3] = {
//            {A.X(), A.Y(), A.Z()},
//            {B.X(), B.Y(), B.Z()},
//            {C.X(), C.Y(), C.Z()},
//            {D.X(), D.Y(), D.Z()}
//        };
//        
//        double pnt[3] = {APnt.X(), APnt.Y(), APnt.Z()};
//        
//        return (GEO::PCK::in_sphere_3d_SOS(co[0], co[1], co[2], co[3], pnt) > 0);
//    }
        //We are in an out tet, where the infinity point is always first
        math::Point p[4] = {
            n[0].getPoint(),
            n[1].getPoint(),
            n[2].getPoint(),
            n[3].getPoint()
        };
        p[infinity_node] = APnt;
//        math::Point p[4] = {
//            APnt,
//            n[m_local_tet_node2facet[infinity_node][0]].getPoint(),
//            n[m_local_tet_node2facet[infinity_node][1]].getPoint(),
//            n[m_local_tet_node2facet[infinity_node][2]].getPoint()
//        };
        
        char s = orient3d(p[0],p[1],p[2],p[3]);
        
        return(s == GEO::POSITIVE || s== GEO::ZERO);
    }
    else
    {
        std::vector<char> signs;
        return isIn(APnt, AR, signs);
//                    
//        //We are in a regular true tet
//        math::Point p[4] = {
//            n[0].getPoint(),
//            n[1].getPoint(),
//            n[2].getPoint(),
//            n[3].getPoint()
//        };
//        
//        double c[4][3] = {
//            {p[0].X(), p[0].Y(), p[0].Z()},
//            {p[1].X(), p[1].Y(), p[1].Z()},
//            {p[2].X(), p[2].Y(), p[2].Z()},
//            {p[3].X(), p[3].Y(), p[3].Z()}
//        };
//        
//        double pnt[3] = {APnt.X(), APnt.Y(), APnt.Z()};
//        
//        return (GEO::PCK::in_sphere_3d_SOS(c[0], c[1], c[2], c[3], pnt) > 0);
    }
    

    return false;
}
/*---------------------------------------------------------------------------*/
bool TetMeshManipulator::stellateCavity(const math::Point& APnt,
                                        std::vector<Region>& ACavBnd,
                                        const int ACavID,
                                        Region& ANewTet)
{
//    std::cout<<"CAV ID "<<ACavID<<std::endl;
    Node c = m_mesh->newNode(APnt);
    TCellID c_id = c.getID();
    //=====================================================================
    //Each tet of ACavBnd is in the cavity and has at least one facet
    // on the cavity boundary
    //=====================================================================
    std::vector<Region> added;

    for(auto tet:ACavBnd){
        std::vector<TCellID> n_ids = tet.getIDs<Node>();
        auto infinity_node = -1;
        for(auto i=0; i<4; i++){
            if(n_ids[i] == InfinityID){
                infinity_node=i;
            }
        }
        std::cout<<"Tet "<<tet.getID()<<" in the cavity"<<std::endl;
        if(isInfinityTet(tet)){
            std::cout<<"> INFINITY REGION"<<std::endl;
        }
        else{
            std::cout<<"> REGULAR REGION"<<std::endl;
        }
        //================================================================
        // Case 1 - Out tet
        //================================================================
        if(infinity_node!=-1){
            //We are in an out tet
            //only one real face which is in the cavity necessary
            // and 3 adjacent ghost regions that can be in the cavity too
            Node bnd[3] = {
                m_mesh->get<Node>(n_ids[m_local_tet_node2facet[infinity_node][0]]),
                m_mesh->get<Node>(n_ids[m_local_tet_node2facet[infinity_node][1]]),
                m_mesh->get<Node>(n_ids[m_local_tet_node2facet[infinity_node][2]])
            };
            bool put_tetra[3]={true, true, true};
            for(auto i=0; i<3; i++){
                std::vector<Region> regs_i = getRegions(bnd[i],bnd[(i+1)%3]);
                //only one ghost region share this edge with tet
                Region ghost_i;
                for(auto j=0; j<regs_i.size() && put_tetra[i];j++){
                    Region rj = regs_i[j];
                    if(rj.getID()!= tet.getID() &&
                       isInfinityTet(rj) &&
                       (*m_cavity)[rj.getID()]==ACavID){
                        put_tetra[i]=false;
                    }
                }
            }

            for(auto i=0; i<3; i++){
                if(put_tetra[i]){
                    Region ti = m_mesh->newTet(InfinityID,
                                               bnd[i].getID(),
                                               bnd[(i+1)%3].getID(),
                                               c_id);
                    
                    added.push_back(ti);
                }
            }
        }
        //================================================================
        // Case 2 - Regulat tet
        //================================================================
        else{
            // Each of the cavity bnd facet will generate a new tet
            for(auto i=0; i<4; i++){
                Node bnd[3] = {
                    m_mesh->get<Node>(n_ids[m_local_tet_node2facet[i][0]]),
                    m_mesh->get<Node>(n_ids[m_local_tet_node2facet[i][1]]),
                    m_mesh->get<Node>(n_ids[m_local_tet_node2facet[i][2]])
                };
                Region ri = getOppositeRegion(tet, bnd[0], bnd[1], bnd[2]);
                std::cout<<"BND "<<bnd[0].getID()
                <<" "<<bnd[1].getID()
                <<" "<<bnd[2].getID();
                std::cout<<" -> r opp: "<<ri.getID()<<std::endl;


                //ri is out of the cavity if it is not marked with ACavID
                // A NullID region should not happen
                if(ri.getID()==NullID){
                    throw GMDSException("Null region should not happen");
                }
                if((*m_cavity)[ri.getID()]==ACavID){
                    std::cout<<"IN CAVITY REGION"<<std::endl;
                    continue;
                }
//                if(isInfinityTet(ri)){
//                    std::cout<<"INFINITY REGION"<<std::endl;
//                    math::Point pc = c.getPoint();
//                    math::Point p0 =bnd[0].getPoint();
//                    math::Point p1 =bnd[1].getPoint();
//                    math::Point p2 =bnd[2].getPoint();
//                    if(pc.areCoplanar(p0,p1,p2))
//                        continue;
//                }
                //Otherwise, bnd[0], bnd[1] and bnd[2] define a bnd face
                Region ni = m_mesh->newTet(c_id,
                                           bnd[0].getID(),
                                           bnd[1].getID(),
                                           bnd[2].getID());
                std::cout<<"Create region "<<c_id<<" "
                <<bnd[0].getID()<<" "
                <<bnd[1].getID()<<" "
                <<bnd[2].getID()<<std::endl;
                //this region must be added at the end of the cavity process
                added.push_back(ni);
            }

        }//else
        
    }
    for(auto r: added){
        std::vector<TCellID> n = r.getIDs<Node>();
        for(auto ni:n){
            if(ni!=InfinityID)
                m_mesh->get<Node>(ni).add(r);
        }
    }
    ANewTet = added[0];
    return true;
}
/*---------------------------------------------------------------------------*/
std::vector<Region> TetMeshManipulator::
getRegions(const Node& AN1, const Node& AN2, const Node& AN3)
{
    std::vector<Region> common;
    std::vector<TCellID> r[3] = {
        AN1.getIDs<Region>(),
        AN2.getIDs<Region>(),
        AN3.getIDs<Region>()
    };
    for(auto r0:r[0]){
        for(auto r1:r[1]){
            for(auto r2:r[2]){
                if(r0==r1 && r0==r2){
                    common.push_back(m_mesh->get<Region>(r0));
                }
            }
        }
    }
    return common;
}
/*---------------------------------------------------------------------------*/
std::vector<Region> TetMeshManipulator::
getRegions(const Region& AR)
{
    std::vector<Region> adj;
    if(isInfinityTet(AR)){
        //one real region
        adj.resize(4);
        facet f = realFacet(AR);
        Node n[3] = {
            m_mesh->get<Node>(f(0)),
            m_mesh->get<Node>(f(1)),
            m_mesh->get<Node>(f(2))
        };
        adj[0] = getOppositeRegion(AR, n[0], n[1], n[2]);
//        return adj;
        //other ghost regions are only sharing an edge with AR.
        for(auto i=0; i<3; i++){
            std::vector<Region> regs_i = getRegions(n[i],n[(i+1)%3]);
            //only one ghost region share this edge with AR
            Region ghost_i;
            bool found_gi = false;
            for(auto j=0; j<regs_i.size() && !found_gi;j++){
                if(regs_i[j].getID()!= AR.getID() &&
                   isInfinityTet(regs_i[j])){
                    found_gi=true;
                    ghost_i = regs_i[j];
                }
            }
            if(!found_gi)
                throw GMDSException("Not found an adj. ghost cell");
            
            adj[i+1]=ghost_i;
        }
    }
    else {
        adj.resize(4);

        //max number of adj region, region with NullID are returned for boundary
        // facets
        
        std::vector<Node> n = AR.get<Node>();
        
        adj[0] = getOppositeRegion(AR, n[0], n[1],n[2]);
        adj[1] = getOppositeRegion(AR, n[0], n[1],n[3]);
        adj[2] = getOppositeRegion(AR, n[0], n[2],n[3]);
        adj[3] = getOppositeRegion(AR, n[1], n[2],n[3]);
    }
    return adj;
}
/*---------------------------------------------------------------------------*/
std::vector<Region> TetMeshManipulator::
getRegions(const Node& AN1, const Node& AN2)
{
    std::vector<Region> common;
    std::vector<TCellID> r[2] = {
        AN1.getIDs<Region>(),
        AN2.getIDs<Region>()
    };
    for(auto r0:r[0]){
        for(auto r1:r[1]){
            if(r0==r1 ){
                common.push_back(m_mesh->get<Region>(r0));
            }
        }
    }
    return common;
}

/*---------------------------------------------------------------------------*/
void TetMeshManipulator::removeNodeViaEdgeContraction(Node& AN1, Node& AN2)
{
    std::vector<Region> n1_reg = AN1.get<Region>();
    
    // Regions common to AN1 and AN2 are going to vanish due to the edge
    // contraction
    std::vector<Region> to_remove_reg = getRegions(AN1,AN2);
    
    //Some regions of AN1 must be connected to AN2 afterward
    std::vector<Region> n1_reg_to_keep;
    
    n1_reg_to_keep.reserve(n1_reg.size()-to_remove_reg.size());
    
    for(auto r:n1_reg){
        bool found_r=false;
        for(auto i=0; i<to_remove_reg.size() && !found_r;i++){
            if(to_remove_reg[i]==r){
                found_r=true;
            }
        }
        if(!found_r){
            n1_reg_to_keep.push_back(r);
        }
    }
    
    // We remove the regions adjacent to [AN1, AN2]
    for(auto r:to_remove_reg){
        std::vector<Node> nodes_r = r.get<Node>();
        for(auto ni:nodes_r){
            ni.remove(r);
        }
        m_mesh->deleteRegion(r);
    }
    
    for(auto r: n1_reg_to_keep){
        //r is connected to AN2 now and vice versa
        AN2.add(r);
        r.replace(AN1, AN2);
    }
    
    m_mesh->deleteNode(AN1);
    
}
/*---------------------------------------------------------------------------*/
bool TetMeshManipulator::isIn(const math::Point APnt,
                              const Region& ATet,
                              math::Vector4d& ACoord) const
{
//    const double EPST = -1.e-14;
//    const double EPSR =  1.e+14;
//    
//    
//    std::vector<Node> n = ATet.get<Node>();
//
//    /* compute barycentrics */
//    math::Point p0 = n[0].getPoint();
//    math::Point p1 = n[1].getPoint();
//    math::Point p2 = n[2].getPoint();
//    math::Point p3 = n[3].getPoint();
//    double	bx,by,bz,cx,cy,cz,dx,dy,dz,vx,vy,vz,apx,apy,apz;
//    double	epsra,vol1,vol2,vol3,vol4,dd;
//    
//    /* barycentric */
//    bx  = p1.X() - p0.X();
//    by  = p1.Y() - p0.Y();
//    bz  = p1.Z() - p0.Z();
//    cx  = p2.X() - p0.X();
//    cy  = p2.Y() - p0.Y();
//    cz  = p2.Z() - p0.Z();
//    dx  = p3.X() - p0.X();
//    dy  = p3.Y() - p0.Y();
//    dz  = p3.Z() - p0.Z();
//    
//    /* test volume */
//    vx  = cy*dz - cz*dy;
//    vy  = cz*dx - cx*dz;
//    vz  = cx*dy - cy*dx;
//    
//    epsra = EPST*(bx*vx + by*vy + bz*vz);
//    apx = APnt.X() - p0.X();
//    apy = APnt.Y() - p0.Y();
//    apz = APnt.Z() - p0.Z();
//    
//    /* p in 2 */
//    vol2  = apx*vx + apy*vy + apz*vz;
//    if ( epsra > vol2 )  return(false);
//    
//    /* p in 3 */
//    vx  = by*apz - bz*apy;
//    vy  = bz*apx - bx*apz;
//    vz  = bx*apy - by*apx;
//    vol3 = dx*vx + dy*vy + dz*vz;
//    if ( epsra > vol3 )  return(false);
//    
//    /* p in 4 */
//    vol4 = -cx*vx - cy*vy - cz*vz;
//    if ( epsra > vol4 )  return(false);
//    
//    /* p in 1 */
//    vol1 = -epsra * EPSR - vol2 - vol3 - vol4;
//    if ( epsra > vol1 )  return(false);
//    
//    dd = vol1+vol2+vol3+vol4;
//    if ( dd != 0.0 )  dd = 1.0 / dd;
//    ACoord[0] = vol1 * dd;
//    ACoord[1] = vol2 * dd;
//    ACoord[2] = vol3 * dd;
//    ACoord[3] = vol4 * dd;
//    
//    return(true);
//
    //==============================================================
    // First we compute the location of APntParam ito ATetParam
    //==============================================================
    double coeff[4]={0, 0, 0, 0};
    
    std::vector<Node> n = ATet.get<Node>();
    try{
        math::Point::computeBarycentric(n[0].getPoint(),
                                        n[1].getPoint(),
                                        n[2].getPoint(),
                                        n[3].getPoint(),
                                        APnt,
                                        coeff[0],coeff[1],
                                        coeff[2],coeff[3]);
    }
    catch(GMDSException& e){
        std::cout<<"Exception bary with region "<<ATet.getID()<<std::endl;
        std::cout<<"PO "<<n[0].getPoint()<<std::endl;
        std::cout<<"P1 "<<n[1].getPoint()<<std::endl;
        std::cout<<"P2 "<<n[2].getPoint()<<std::endl;
        std::cout<<"P3 "<<n[3].getPoint()<<std::endl;
        return false;
    }
    ACoord[0] = coeff[0];
    ACoord[1] = coeff[1];
    ACoord[2] = coeff[2];
    ACoord[3] = coeff[3];

    //==============================================================
    // If APntParam is "almost" into ATetParam, then we provide its
    // Barycentric coordinates
    //==============================================================
    double tolerance = 0;
    if(coeff[0]<tolerance || coeff[1]<tolerance ||
       coeff[2]<tolerance || coeff[3]<tolerance)
        return false;
    
    return true;
}
/*---------------------------------------------------------------------------*/
Region TetMeshManipulator::getOppositeRegion(const Region& AR,
                                             const Node&   AN0,
                                             const Node&   AN1,
                                             const Node&   AN2)
{
    //==============================================================
    // We look for regions adjacent to all opp_nodes;
    //==============================================================
    std::vector<Region> common = getRegions(AN0, AN1, AN2);
    
    if(common.size()==2){
        if(common[0].getID()==AR.getID()){
            return common[1];
        }
        else if(common[1].getID()==AR.getID()){
            return common[0];
        }
        else {
            throw GMDSException("TetMeshManipulator::getOppositeRegion (1)");
        }
    }
    else if(common.size()==1){
        if(common[0]!=AR){
            throw GMDSException("TetMeshManipulator::getOppositeRegion (2)");
        }
        return Region();
    }
    else{
        std::cout<<AN0.getPoint()<<std::endl;
        std::cout<<AN1.getPoint()<<std::endl;
        std::cout<<AN2.getPoint()<<std::endl;
        for(auto c:common){
            std::cout<<c.getID()<<" ";
            if(isInfinityTet(c))
                std::cout<<"infinity";
            else
                std::cout<<"real";
            std::cout<<std::endl;
        }
        throw GMDSException("TetMeshManipulator::getOppositeRegion (3)");
    }
    
    throw GMDSException("TetMeshManipulator::getOppositeRegion (4)");
}
/*---------------------------------------------------------------------------*/
Region TetMeshManipulator::seed(){
    /* The implementation of seed was initially recursive but it can lead
     * to memory issues due to the great number of infinite regions*/
    // MUst be optimized by looking only into real tet.
    std::cout<<"New seed"<<std::endl;
    int seed_id;
    do{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, m_mesh->getNbRegions()-1);
        seed_id=dis(gen);
    }
    while(!m_mesh->has<Region>(seed_id) ||
          isInfinityTet(m_mesh->get<Region>(seed_id)));

    //The seed is a real cell now
    return m_mesh->get<Region>(seed_id);

}
/*---------------------------------------------------------------------------*/
int TetMeshManipulator::random(const int& AMax) const {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, AMax-1);
    
    return  dis(gen);
}

/*---------------------------------------------------------------------------*/
void TetMeshManipulator::insertPointsTetgen(IGMesh* AInMesh,
                                            std::vector<math::Point>& APoints,
                                            IGMesh* AOutMesh)
{
    
    tetgenio in_mesh, between_mesh, out_mesh, points;
    //    in_mesh.initialize();
    //    out_mesh.initialize();
    //
    //    points.initialize();
    
    /** map of gmds node ids to tetgen point ids */
    std::map<TCellID,int> gmds_to_tetgen;
    /** map of tetgen point ids to gmds node ids */
    std::map<int,TCellID> tetgen_to_gmds;
    //====================================================================
    // 1) AMesh is copied into a tetgen mesh.
    //====================================================================
    
    //the new mesh is loaded
    loadTetMesh(AInMesh, in_mesh, gmds_to_tetgen, tetgen_to_gmds);
    
    int nb_points_init= AInMesh->getNbNodes();
    //the point to be inserted too
    loadNodes(APoints, points);
    std::cout<<"IN: "<<AInMesh->getNbNodes()<<" + "<<APoints.size()<<std::endl;
    //====================================================================
    // 2) Tetgen mesh is now refined
    //====================================================================
    std::string param_str = "riq";
    char *param = new char[param_str.length()+1];
    strcpy(param, param_str.c_str());
    // do stuff
    tetrahedralize(param, &in_mesh, &between_mesh, &points);
    std::cout<<" INTERMEDIATE MESH: "<<between_mesh.numberofpoints<<std::endl;
    
    
    auto nb_points = between_mesh.numberofpoints;
    between_mesh.pointmarkerlist= new int[nb_points ];
    
    for (auto i=0; i<nb_points; i++){
        if(i<nb_points_init)
            between_mesh.pointmarkerlist[i]=-1;
        else
            between_mesh.pointmarkerlist[i]=0;
    }
    
    
    
    std::string param_str2 = "rRV";
    char *param2 = new char[param_str.length()+1];
    strcpy(param2, param_str2.c_str());
    // do stuff
    tetrahedralize(param2, &between_mesh, &out_mesh, &points);
    std::cout<<" OUT MESH: "<<out_mesh.numberofpoints<<std::endl;
    
    //====================================================================
    // 3) Now we fill in the new mesh AOutMesh
    //====================================================================
    /** map of tetgen point ids to gmds node ids */
    std::map<int,Node> node_id;
    
    for(auto i=0; i<out_mesh.numberofpoints; i++){
        
        double *coord = &out_mesh.pointlist[i * 3];
        math::Point pi(coord[0],coord[1],coord[2]);
        Node ni = AOutMesh->newNode(pi);
        node_id[i]=ni;
        
    }
    for (auto i=0; i<out_mesh.numberoftetrahedra; i++) {
        int *plist = &(out_mesh.tetrahedronlist[i * 4]);
        AOutMesh->newTet(node_id[plist[0]],
                         node_id[plist[1]],
                         node_id[plist[2]],
                         node_id[plist[3]]);
        
    }
    
    VTKWriter<IGMesh> writer(*AOutMesh);
    writer.write("AFTER_TETGEN", DIM3 | R | N);
    
    
    delete [] param;
}

/*---------------------------------------------------------------------------*/
void TetMeshManipulator::tetrahedralizeTetgen(gmds::IGMesh* ABndMesh,
                                              gmds::IGMesh* ATetMesh)
{
    
    tetgenio in_mesh, out_mesh, points;
    //    in_mesh.initialize();
    //    out_mesh.initialize();
    //
    //    points.initialize();
    
    /** map of gmds node ids to tetgen point ids */
    std::map<TCellID,int> gmds_to_tetgen;
    /** map of tetgen point ids to gmds node ids */
    std::map<int,TCellID> tetgen_to_gmds;
    //====================================================================
    // 1) AMesh is copied into a tetgen mesh.
    //====================================================================
    
    //the new mesh is loaded
    loadTriangularMesh(ABndMesh, in_mesh, gmds_to_tetgen, tetgen_to_gmds);
    
    int nb_points_init= ABndMesh->getNbNodes();
    //====================================================================
    // 2) Tetgen mesh is now refined
    //====================================================================
    std::string param_str = "pYq";
    char *param = new char[param_str.length()+1];
    strcpy(param, param_str.c_str());
    // do stuff
    tetrahedralize(param, &in_mesh, &out_mesh);
    std::cout<<" TETGEN MESH NB POINTS: "<<out_mesh.numberofpoints<<std::endl;
    
    
    auto nb_points = out_mesh.numberofpoints;
    out_mesh.pointmarkerlist= new int[nb_points ];
    
    for (auto i=0; i<nb_points; i++){
        if(i<nb_points_init)
            out_mesh.pointmarkerlist[i]=-1;
        else
            out_mesh.pointmarkerlist[i]=0;
    }
    
    //====================================================================
    // 3) Now we fill in the new mesh AOutMesh
    //====================================================================
    /** map of tetgen point ids to gmds node ids */
    std::map<int,Node> node_id;
    
    for(auto i=0; i<out_mesh.numberofpoints; i++){
        
        double *coord = &out_mesh.pointlist[i * 3];
        math::Point pi(coord[0],coord[1],coord[2]);
        Node ni = ATetMesh->newNode(pi);
        node_id[i]=ni;
        
    }
    for (auto i=0; i<out_mesh.numberoftetrahedra; i++) {
        int *plist = &(out_mesh.tetrahedronlist[i * 4]);
        ATetMesh->newTet(node_id[plist[0]],
                         node_id[plist[1]],
                         node_id[plist[2]],
                         node_id[plist[3]]);
        
    }
    
    static int cav_nb = 0;
    VTKWriter<IGMesh> writer(*ATetMesh);
    writer.write("AFTER_TETGEN_"+std::to_string(cav_nb++), DIM3 | R | N);
    
    
    delete [] param;
}
/*---------------------------------------------------------------------------*/
void TetMeshManipulator::loadTetMesh(IGMesh* AMesh,
                                     tetgenio& ATetMesh,
                                     std::map<TCellID,int>& AGMDS2Tet,
                                     std::map<int,TCellID>& ATet2GMDS)
{
    int nverts = AMesh->getNbNodes();
    // Allocate memory for 'tetgenio'
    if (nverts > 0) {
        ATetMesh.numberofpoints = nverts;
        ATetMesh.pointlist = new REAL[nverts * 3];
    }
    
    // Read the node list.
    int point_index =0;
    
    for (IGMesh::node_iterator itn=AMesh->nodes_begin();
         !itn.isDone(); itn.next()) {
        Node n = itn.value();
        math::Point p = n.getPoint();
        
        double *coord = &ATetMesh.pointlist[point_index * 3];
        coord[0] = p.X();
        coord[1] = p.Y();
        coord[2] = p.Z();
        
        //update map between tetgen and gmds meshes
        AGMDS2Tet[n.getID()  ] = point_index;
        ATet2GMDS[point_index] = n.getID();
        
        //go to the next index in m_tet_mesh.pointlist
        point_index++;
    }
    
    
    int ntets = AMesh->getNbRegions();
    
    if (ntets > 0) {
        // It is a tetrahedral mesh.
        ATetMesh.numberoftetrahedra = ntets;
        ATetMesh.numberofcorners = 4;
        ATetMesh.numberoftetrahedronattributes = 1;
        ATetMesh.tetrahedronlist = new int[ntets * 4];
        ATetMesh.tetrahedronattributelist = new REAL[ntets];
    }
    int tet_index =0;
    for (IGMesh::region_iterator itr=AMesh->regions_begin();
         !itr.isDone(); itr.next()) {
        Region r = itr.value();
        std::vector<TCellID> n = r.getIDs<Node>();
        int *plist = &(ATetMesh.tetrahedronlist[tet_index * 4]);
        for (auto i=0; i<4; i++) {
            plist[i] = AGMDS2Tet[n[i]];
        }
        //go to the next index in m_tet_mesh.tetrahedronlist
        tet_index++;
    }
    // the firstnumber of the index.
    ATetMesh.firstnumber = 0;
}

/*---------------------------------------------------------------------------*/
void TetMeshManipulator::loadTriangularMesh(IGMesh* AMesh,
                                            tetgenio& ATriMesh,
                                            std::map<TCellID,int>& AGMDS2Tet,
                                            std::map<int,TCellID>& ATet2GMDS)
{
    int nverts = AMesh->getNbNodes();
    
    tetgenio::facet *facet;
    tetgenio::polygon *p;
    
    // All indices start from 1.
    ATriMesh.firstnumber = 0;
    
    ATriMesh.numberofpoints = nverts;
    
    
    
    // Allocate memory for 'tetgenio'
    ATriMesh.numberofpoints = nverts;
    ATriMesh.pointlist = new REAL[nverts * 3];
    
    // Read the node list.
    int point_index =0;
    
    for (IGMesh::node_iterator itn=AMesh->nodes_begin();
         !itn.isDone(); itn.next()) {
        Node n = itn.value();
        math::Point p = n.getPoint();
        
        double *coord = &ATriMesh.pointlist[point_index * 3];
        coord[0] = p.X();
        coord[1] = p.Y();
        coord[2] = p.Z();
        
        //update map between tetgen and gmds meshes
        AGMDS2Tet[n.getID()  ] = point_index;
        ATet2GMDS[point_index] = n.getID();
        
        //go to the next index in m_tet_mesh.pointlist
        point_index++;
    }

    //Must be modified using the surface-color to define a facet made
    //of many triangles
    int ntris = AMesh->getNbFaces();
    ATriMesh.numberoffacets = ntris;
    ATriMesh.facetlist = new tetgenio::facet[ATriMesh.numberoffacets];
//    ATriMesh.facetmarkerlist = new int[ATriMesh.numberoffacets];
    
    int tri_index =0;
    for (IGMesh::face_iterator itf=AMesh->faces_begin();
         !itf.isDone(); itf.next()) {
        Face f = itf.value();
        
        facet = &ATriMesh.facetlist[tri_index];
        facet->numberofpolygons = 1;
        facet->polygonlist = new tetgenio::polygon[facet->numberofpolygons];
        facet->numberofholes = 0;
        facet->holelist = NULL;
        p = &facet->polygonlist[0];
        p->numberofvertices = 3;
        p->vertexlist = new int[p->numberofvertices];

        std::vector<TCellID> n = f.getIDs<Node>();
        for (auto i=0; i<3; i++) {
            p->vertexlist[i] = AGMDS2Tet[n[i]];
        }
        
//        ATriMesh.facetmarkerlist[tri_index]=0;
        //go to the next index in m_tet_mesh.tetrahedronlist
        tri_index++;
    }
    
}
/*---------------------------------------------------------------------------*/
void TetMeshManipulator::loadNodes(std::vector<math::Point>& APnts,
                                   tetgenio& ATetPnts)
{
    ATetPnts.mesh_dim = 3;
    ATetPnts.numberofpointattributes = 0;  // no point attribute.

    int nb_points = APnts.size();
    // Allocate memory for 'tetgenio'
    if (nb_points > 0) {
        ATetPnts.numberofpoints = nb_points;
        ATetPnts.pointlist = new REAL[nb_points * 3];
        ATetPnts.pointattributelist=NULL;
        ATetPnts.pointparamlist=NULL;
    }
    
    for (auto i=0; i<APnts.size(); i++){
        math::Point p = APnts[i];
        
        ATetPnts.pointlist[(i*3)  ] = p.X();
        ATetPnts.pointlist[(i*3)+1] = p.Y();
        ATetPnts.pointlist[(i*3)+2] = p.Z();
    }
    
}
/*---------------------------------------------------------------------------*/
bool TetMeshManipulator::isIn(const math::Point& AP,
                              const Region& AR,
                              std::vector<char>& AS) const
{
    AS.resize(4);
    std::vector<Node> n = AR.get<Node>();
    math::Point p[4] ={
        n[0].getPoint(),
        n[1].getPoint(),
        n[2].getPoint(),
        n[3].getPoint()
    };
    double pnt[3] = {AP.X(), AP.Y(), AP.Z()};
    double p0[3] = {p[0].X(), p[0].Y(), p[0].Z()};
    double p1[3] = {p[1].X(), p[1].Y(), p[1].Z()};
    double p2[3] = {p[2].X(), p[2].Y(), p[2].Z()};
    double p3[3] = {p[3].X(), p[3].Y(), p[3].Z()};

    // Check for flat tetrahedra
    math::Vector3d v01(p[0],p[1]);
    math::Vector3d v02(p[0],p[2]);
    math::Vector3d v03(p[0],p[3]);
    if (std::abs(v03.dot(v02.cross(v03)))<1e-15) {
        return false;
    }
    
    AS[0] = GEO::PCK::orient_3d(pnt, p1 , p2 , p3 );
    AS[1] = GEO::PCK::orient_3d(p0 , pnt, p2 , p3 );
    AS[2] = GEO::PCK::orient_3d(p0 , p1 , pnt, p3 );
    AS[3] = GEO::PCK::orient_3d(p0 , p1 , p2 , pnt);
    
    if ((AS[0] >= 0 && AS[1] >= 0 && AS[2] >= 0 && AS[3] >= 0) ||
        (AS[0] <= 0 && AS[1] <= 0 && AS[2] <= 0 && AS[4] <= 0)) {
        return true;
    }
    
    return false; 
}

/*---------------------------------------------------------------------------*/
bool TetMeshManipulator::same(math::Point& AP1, math::Point& AP2)
{
    return (AP1.X()==AP2.X()) && (AP1.Y()==AP2.Y()) && (AP1.Z()==AP2.Z());
}

/*---------------------------------------------------------------------------*/
bool TetMeshManipulator::colinear(math::Point& AP1,
                                  math::Point& AP2,
                                  math::Point& AP3)
{
    // Colinearity is tested by using four coplanarity
    // tests with four points that are not coplanar.
    static const double q000[3] = {0.0, 0.0, 0.0};
    static const double q001[3] = {0.0, 0.0, 1.0};
    static const double q010[3] = {0.0, 1.0, 0.0};
    static const double q100[3] = {1.0, 0.0, 0.0};
    
    double p1[3]={AP1.X(),AP1.Y(),AP1.Z()};
    double p2[3]={AP2.X(),AP2.Y(),AP2.Z()};
    double p3[3]={AP3.X(),AP3.Y(),AP3.Z()};
    
    return (GEO::PCK::orient_3d(p1, p2, p3, q000) == GEO::ZERO &&
            GEO::PCK::orient_3d(p1, p2, p3, q001) == GEO::ZERO &&
            GEO::PCK::orient_3d(p1, p2, p3, q010) == GEO::ZERO &&
            GEO::PCK::orient_3d(p1, p2, p3, q100) == GEO::ZERO);
}

/*---------------------------------------------------------------------------*/
char TetMeshManipulator::orient3d(const gmds::math::Point& AP0,
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
char TetMeshManipulator::orient3d(const gmds::math::Point& AP,
                                  const gmds::Node& AN1,
                                  const gmds::Node& AN2,
                                  const gmds::Node& AN3)
{
    double p0[3]={AP.X(),AP.Y(),AP.Z()};
    double p1[3]={AN1.getPoint().X(),AN1.getPoint().Y(),AN1.getPoint().Z()};
    double p2[3]={AN2.getPoint().X(),AN2.getPoint().Y(),AN2.getPoint().Z()};
    double p3[3]={AN3.getPoint().X(),AN3.getPoint().Y(),AN3.getPoint().Z()};
    
    return GEO::PCK::orient_3d(p0, p1, p2, p3);
}

/*---------------------------------------------------------------------------*/
bool TetMeshManipulator::computeProj(const math::Point& AP,
                                     const math::Point& AT1,
                                     const math::Point& AT2,
                                     const math::Point& AT3,
                                     math::Point& ANewP)
{
    //We project AP onto the plane defined by AT1, AT2 and AT3
    math::Plane plane(AT1,AT2,AT3);
    math::Point p = plane.project(AP);
    
    // We look if p is insie or outside of the triangle defined by AT1,
    // AT2 and AT3
    math::Vector normal = plane.getNormal();
    
    char ori[3] ={
        orient3d(p, AT1 , AT2 , AT1+normal),
        orient3d(p, AT2 , AT3 , AT2+normal),
        orient3d(p, AT3 , AT1 , AT3+normal)
    };
    if ((ori[0] >= 0 && ori[1] >= 0 && ori[2] >= 0 ) ||
        (ori[0] <= 0 && ori[1] <= 0 && ori[2] <= 0 ) ) {
        ANewP=p;
        return true;
    }
    return false;
}

/*---------------------------------------------------------------------------*/
double TetMeshManipulator::quality(const math::Point& AN1, const math::Point& AN2,
                                   const math::Point& AN3, const math::Point& AN4)
{
    //=====================================================================
    // See V. N. Parthasarathy, C. M. Graichen, and A. F. Hathaway. A
    // Comparison of Tetrahedron Quality Measures. Finite Elements in
    // Analysis and Design 15(3):255–261, January 1994.
    //=====================================================================
    math::Point p[4]={AN1, AN2, AN3, AN4};
    
    double  vol = math::Tetrahedron(p[0], p[2], p[1], p[3]).getVolume();
    
    double h[6] ={
        math::Vector3d(p[0],p[1]).norm(),
        math::Vector3d(p[0],p[2]).norm(),
        math::Vector3d(p[0],p[3]).norm(),
        math::Vector3d(p[1],p[2]).norm(),
        math::Vector3d(p[1],p[3]).norm(),
        math::Vector3d(p[2],p[3]).norm()
    };
    double hs2 =0;
    for(auto hi:h){
        hs2 += hi*hi;
    }
    double hs = sqrt(hs2);
    return abs((6*sqrt(2.0)*vol)/(hs*hs*hs));
    //*6*sqrt(2) to get highest quality at 1 (equilateral tet)
    
}
/*---------------------------------------------------------------------------*/


