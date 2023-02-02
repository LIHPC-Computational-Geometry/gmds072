/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * The GMDS library is a computer program whose purpose is to provide a set of
 * functionnalities to represent and handle any type of meshes (2D, 3D,
 * triangles, tetrahedra, quad, hexa, polygons, polyhedra, etc.) and write
 * meshing algorithms. So it gathers many mathematical objects like points,
 * segment, quaternions, etc. and basic algorithms useful to build more evolved
 * ones.
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
/** \file    DistanceFieldBuilder2D.cpp
 *  \author  F. LEDOUX
 *  \date    08/31/2010
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/Algo/DistanceFieldBuilder2D.h>
#include <GMDS/Math/Numerics.h>
#include <set>
#include <map>
#include <list>
#include <climits>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
DistanceFieldBuilder2D::DistanceFieldBuilder2D(IGMesh* AMesh)
:m_mesh(AMesh), m_distance(0), m_cross_field(0)
{
    m_adv_field[0]=0;
    m_adv_field[1]=0;
    m_adv_dist[0]=0;
    m_adv_dist[1]=0;
}
/*----------------------------------------------------------------------------*/
DistanceFieldBuilder2D::~DistanceFieldBuilder2D()
{}
/*----------------------------------------------------------------------------*/
bool DistanceFieldBuilder2D::isValid()
{
    //F2N and N2F connectivites are required
    MeshModel model = m_mesh->getModel();
    if (!model.has(N2F))
    {
        std::cout << "Error in the 2D distance field computation, N2F connectivity must be available" << std::endl;
        return false;
    }
    
    //The mesh faces must be triangles only
    IGMesh::face_iterator it_faces = m_mesh->faces_begin();
    for (; !it_faces.isDone(); it_faces.next()){
        Face f = it_faces.value();
        if (f.getType() != GMDS_TRIANGLE)
        {
            std::cout << "Error in the 2D distance field computation, ";
            std::cout << "only triangular faces are supported" << std::endl;
            return false;
        }
    }
    
    
    return true;
    
}

/*----------------------------------------------------------------------------*/
std::vector<Node> DistanceFieldBuilder2D::
getAdjacentNodes(Node& ANode, const int AMarkNode, const int AMarkFace)
{
    std::set<Node> adj_nodes;
    std::vector<Face>  adj_faces = ANode.get<Face>();
    
    for (unsigned int i= 0; i < adj_faces.size(); i++)
    {
        Face f_i = adj_faces[i];
        if (m_mesh->isMarked(f_i, AMarkFace))
        {
            std::vector<Node> nodes_fi = f_i.get<Node>();
            for (unsigned int j = 0; j < nodes_fi.size(); j++)
            {
                Node n_j = nodes_fi[j];
                if (m_mesh->isMarked(n_j, AMarkNode) && n_j.getID() != ANode.getID())
                {
                    adj_nodes.insert(n_j);
                }
                
            }//for (unsigned int j = 0; j < nodes_fi.size(); j++)
            
        }//if (m_mesh->isMarked(f_i, AMarkFace))
        
    }//for (unsigned int i= 0; i < adj_faces.size(); i++)
    
    std::vector<Node> adj;
    adj.insert(adj.end(), adj_nodes.begin(), adj_nodes.end());
    
    return adj;
}
/*----------------------------------------------------------------------------*/
std::vector<TCellID> DistanceFieldBuilder2D::getInsertionOrder() const
{
    return m_insertion_order;
}
/*----------------------------------------------------------------------------*/
double DistanceFieldBuilder2D::
extrapolateDistance(Node& AN, Node& AN1, Node& AN2,
                    Variable<double>*& ADist,
                     const bool AOnlyOne)
{
    std::cout<<"# Compute distance in "<<AN.getID()<<" from "
    <<AN1.getID()<<" and "<<AN2.getID()<<std::endl;
    
    double d1 = (*ADist)[AN1.getID()];
    double d2 = (*ADist)[AN2.getID()];
    //  std::cout<<"# with "<<d1<<", "<<d2<<std::endl;
    
    math::Vector v12(AN1.getPoint(), AN2.getPoint());
    math::Vector v13(AN1.getPoint(), AN.getPoint());
    
    math::Vector normal = v12.cross(v13);
    
    math::Vector x = v12.getNormalize();
    math::Vector y = x.cross(normal);
    y.normalize();
    
    if(y.dot(v13)<0)
        y= y.opp();
    
    TCoord x2 = v12.norm();
    TCoord x3 = x.dot(v13);
    TCoord y3 = y.dot(v13);
    
    TCoord c = (d2-d1) / x2;
    TCoord d3 = d1 + c*x3 + sqrt(1 - c*c)*y3;
    //std::cout<<"# (x2, x3, y3, c, d3)=("<<x2<<", "<<x3<<", "<<y3<<", "<<c<<", "<<d3<<std::endl;
    //std::cout<<"#     --> "<<d3<<std::endl;
    
    return d3;
    
}

/*----------------------------------------------------------------------------*/
double DistanceFieldBuilder2D::
extrapolateLInfDistance(Node& AN,
                        Node&AN1,
                        Node& AN2,
                        Variable<double>*& ADist,
                        const bool AOnlyOne)
{
    TCoord d3 = 0;
    math::Cross2D ref_cross = (*m_cross_field)[AN.getID()];
    
    std::vector<math::Vector> ref_vectors = ref_cross.componentVectors();

    if(AOnlyOne){//only AN1 is available
     
        
        std::cout<<"compute in "<<AN.getID()<<" from ONLY "<<AN1.getID()<<std::endl;
        exit(0);
        
    }
    else {
        //  std::cout<<"# Compute distance in "<<AN.getID()<<" from "
        //<<AN1.getID()<<" and "<<AN2.getID()<<std::endl;
        
        std::cout<<"compute in "<<AN.getID()<<" from "<<AN1.getID()
            <<" and "<<AN2.getID()<<std::endl;
        //  std::cout<<"# with "<<d1<<", "<<d2<<std::endl;
        if(AN.getID()==91 || AN.getID()==805 || AN.getID()==789){
            std::cout<<"check"<<std::endl;
        }
        
        math::Vector adv1_u = (*m_adv_field[0])[AN1.getID()];
        math::Vector adv1_v = (*m_adv_field[1])[AN1.getID()];
        double dist1_u =(*m_adv_dist[0])[AN1.getID()];
        double dist1_v =(*m_adv_dist[1])[AN1.getID()];
        std::cout<<"adv1_u: "<<adv1_u<<std::endl;
        std::cout<<"adv1_v: "<<adv1_v<<std::endl;
        math::Vector adv2_u = (*m_adv_field[0])[AN2.getID()];
        math::Vector adv2_v = (*m_adv_field[1])[AN2.getID()];
        double dist2_u =(*m_adv_dist[0])[AN2.getID()];
        double dist2_v =(*m_adv_dist[1])[AN2.getID()];
        std::cout<<"adv2_u: "<<adv2_u<<std::endl;
        std::cout<<"adv2_v: "<<adv2_v<<std::endl;

        math::Vector u[2], v[2];
        double du[2], dv[2];
        if(fabs(adv1_u.dot(adv2_u)) > fabs(adv1_u.dot(adv2_v))){
            //match u1 and u2
            u[0]=adv1_u;
            u[1]=adv2_u;
            
            du[0]=dist1_u;
            du[1]=dist2_u;

            v[0]=adv1_v;
            v[1]=adv2_v;
            
            dv[0]=dist1_v;
            dv[1]=dist2_v;

        }
        else{
            u[0]=adv1_u;
            u[1]=adv2_v;
            
            du[0]=dist1_u;
            du[1]=dist2_v;

            v[0]=adv1_v;
            v[1]=adv2_u;
            
            dv[0]=dist1_v;
            dv[1]=dist2_u;
        }
        if(u[0].dot(u[1])<0)
            u[1] = u[1].opp();
        
        if(v[0].dot(v[1])<0)
            v[1] = v[1].opp();
        
        std::cout<<"u[0]: "<<u[0]<<std::endl;
        std::cout<<"u[1]: "<<u[1]<<std::endl;
        std::cout<<"v[0]: "<<v[0]<<std::endl;
        std::cout<<"v[1]: "<<v[1]<<std::endl;

        //now I look for the ref vector that maximize the dot product with
        //v1
        math::Vector u_N =ref_vectors[0];
        math::Vector v_N =ref_vectors[0];
        TCoord dot_u = u_N.dot(u[0]);
        TCoord dot_v = v_N.dot(v[0]);
        for(int i=1; i<4; i++){
            math::Vector ref_i =ref_vectors[i];
            TCoord dot_ui = ref_i.dot(u[0]);
            TCoord dot_vi = ref_i.dot(v[0]);
            if(dot_ui>dot_u){
                dot_u=dot_ui;
                u_N = ref_i;
            }
            if(dot_vi>dot_v){
                dot_v=dot_vi;
                v_N = ref_i;
            }
        }
        std::cout<<"u_N: "<<u_N<<std::endl;
        std::cout<<"v_N: "<<v_N<<std::endl;

        //compute the relative distance
        math::Vector v1(AN1.getPoint(), AN.getPoint());
        math::Vector v2(AN2.getPoint(), AN.getPoint());
        double norm0 = v1.norm();
        double norm1 = v2.norm();
        double du_0 =(v1.dot(u_N));
        double du_1 =(v2.dot(u_N));
        
        std::cout<<"du[0]= "<<du[0]<<std::endl;
        std::cout<<"du[1]= "<<du[1]<<std::endl;
        std::cout<<"du_0 = "<<du_0<<std::endl;
        std::cout<<"du_1 = "<<du_1<<std::endl;
        
        double dv_0 =(v1.dot(v_N));
        double dv_1 =(v2.dot(v_N));
        std::cout<<"dv[0]= "<<dv[0]<<std::endl;
        std::cout<<"dv[1]= "<<dv[1]<<std::endl;
        std::cout<<"dv_0 = "<<dv_0<<std::endl;
        std::cout<<"dv_1 = "<<dv_1<<std::endl;

        double du01;
        if (du_0>=0 && du_1>=0) {
            du01 = (norm1*(du_0+du[0]) + norm0*(du_1+du[1]))/(norm0+norm1);
        }
        else if (du_0<=0 && du_1<=0){
            du01 = (norm1*(-du_0+du[0]) + norm0*(-du_1+du[1]))/(norm0+norm1);
        }
        else if (du_0>=0){
            du01 = du_0+du[0];
        }
        else {
            du01 = du_1+du[1];
        }
        
        double dv01;
        if (dv_0>=0 && dv_1>=0) {
            dv01 = (norm1*(dv_0+dv[0]) + norm0*(dv_1+dv[1]))/(norm0+norm1);
        }
        else if (dv_0<=0 && dv_1<=0){
            dv01 = (norm1*(-dv_0+dv[0]) + norm0*(-dv_1+dv[1]))/(norm0+norm1);
        }
        else if (du_0>=0){
            dv01 = dv_0+dv[0];
        }
        else {
            dv01 = dv_1+dv[1];
        }
        
        
        double du_N = du01;
        double dv_N = dv01;

        
        (*m_adv_field[0])[AN.getID()] = u_N;
        (*m_adv_field[1])[AN.getID()] = v_N;
        (*m_adv_dist[0])[AN.getID()]  = du_N;
        (*m_adv_dist[1])[AN.getID()]  = dv_N;
        std::cout<<"  --> ("<<du_N<<", "<<dv_N<<")"<<std::endl;
        d3=math::max2(du_N, dv_N);
        
        
        //std::cout<<"# (x2, x3, y3, c, d3)=("<<x2<<", "<<x3<<", "<<y3<<", "<<c<<", "<<d3<<std::endl;
        std::cout<<"#     --> "<<d3<<std::endl;
    }
    return d3;
    
}
/*---------------------------------------------------------------------------*/
bool DistanceFieldBuilder2D::
belongToTheSameFace(Node& AN1, Node& AN2, Node& AN3, Face& AF)
{
    std::vector<TCellID> common_faces_12;
    std::vector<TCellID> common_faces_123;
    std::vector<TCellID> fs1 = AN1.getIDs<Face>();
    std::vector<TCellID> fs2 = AN2.getIDs<Face>();
    for (unsigned int i = 0; i < fs1.size(); i++)
    {
        TCellID id1 = fs1[i];
        bool found = false;
        for (unsigned int j = 0; !found && j < fs2.size(); j++)
        {
            TCellID id2 = fs2[j];
            if (id1 == id2)
            {
                found = true;
                common_faces_12.push_back(id1);
            }
        }
    }
    fs2 = AN3.getIDs<Face>();
    
    for (unsigned int i = 0; i < common_faces_12.size(); i++)
    {
        TCellID id1 = common_faces_12[i];
        bool found = false;
        for (unsigned int j = 0; !found && j < fs2.size(); j++)
        {
            TCellID id2 = fs2[j];
            if (id1 == id2)
            {
                found = true;
                common_faces_123.push_back(id1);
            }
        }
    }
    if (common_faces_123.size() == 1){
        AF = m_mesh->get<Face>(common_faces_123[0]);
        return true;
    }
    return false;
}
/*----------------------------------------------------------------------------*/
void DistanceFieldBuilder2D::
getNodesInTheSameFace(const std::vector<Node>& AIN,
                      std::vector<std::pair<TCellID, TCellID> >& AOUT)
{
    for (unsigned int i = 0; i < AIN.size(); i++)
    {
        Node ni = AIN[i];
        std::vector<TCellID> fsi = ni.getIDs<Face>();
        for (unsigned int j = i+1; j < AIN.size(); j++)
        {
            Node nj = AIN[j];
            std::vector<TCellID> fsj = nj.getIDs<Face>();
            
            //we look for a face that would be both in fsi and fsj
            bool found = false;
            for (unsigned int ki = 0; !found && ki < fsi.size(); ki++)
            {
                TCellID id1 = fsi[ki];
                for (unsigned int kj = 0; !found && kj < fsj.size(); kj++)
                {
                    TCellID id2 = fsj[kj];
                    if (id1 == id2)
                    {
                        found = true;
                    }
                }
            }
            if (found)
            {
                AOUT.push_back(std::pair<TCellID, TCellID>(ni.getID(), nj.getID()));
            }
        }
    }
}
/*----------------------------------------------------------------------------*/
void DistanceFieldBuilder2D::
initNarrowBand(std::vector<Node>& AFrom,
               std::multimap<double, Node>& ANarrowBand,
               std::map<TCellID, double>& AExtrapolateDistance,
               const int AMarkNodeToWorkOn,
               const int AMarkNodeAlive,
               const int AMarkNodeNarrow,
               const int AMarkFace)
{
    std::vector<Node> source_nodes;
    for (unsigned int i = 0; i < AFrom.size(); i++) {
        Node n_i = AFrom[i];
        std::vector<Node> adj_nodes =
            getAdjacentNodes(n_i, AMarkNodeToWorkOn, AMarkFace);
        for (unsigned int j = 0; j < adj_nodes.size(); j++) {
            Node n_j = adj_nodes[j];
            // we look for nodes that are not yet in the narrow band or alive
            if (m_mesh->isMarked(n_j, AMarkNodeNarrow) ||
                m_mesh->isMarked(n_j, AMarkNodeAlive))
                continue;
            
            std::vector<Node> adj_nodes_j = getAdjacentNodes(n_j, AMarkNodeToWorkOn, AMarkFace);
            std::vector<Node> ext_nodes;
            int nb_alive_adj_nodes = 0;
            for (unsigned int k = 0; k < adj_nodes_j.size(); k++) {
                Node n_k = adj_nodes_j[k];
                if (m_mesh->isMarked(n_k, AMarkNodeAlive))
                {
                    nb_alive_adj_nodes++;
                    ext_nodes.push_back(n_k);
                }
            } //for (unsigned int k = 0; k < adj_nodes_j.size(); k++)
            
            if (nb_alive_adj_nodes == 1) {
                source_nodes.push_back(ext_nodes[0]);
                           }//if (nb_alive_adj_nodes == 1)
            else if (nb_alive_adj_nodes == 2) {
                Face commonFace;
                if (belongToTheSameFace(n_j, ext_nodes[0], ext_nodes[1], commonFace))
                {
                    m_mesh->mark(n_j, AMarkNodeNarrow);
                    
                    //we extrapolate the distance in n_j
                    double d_j = (this->*m_extrapolate_func_ptr)(n_j, ext_nodes[0], ext_nodes[1], m_distance,false);
                    if (d_j < 0)
                    {
                        d_j = 0;
                        std::cout << "Warning, we get a negative distance (2D)"
                        << std::endl;
                        AExtrapolateDistance[n_j.getID()] = 10000;
                        ANarrowBand.insert(std::pair<double, Node>(d_j, n_j));
                    }
                    else
                    {
                        AExtrapolateDistance[n_j.getID()] = d_j;
                        ANarrowBand.insert(std::pair<double, Node>(d_j, n_j));
                    }
                }
            }//if (nb_alive_adj_nodes == 2)
            else if (nb_alive_adj_nodes > 2) {
                
                //we look for 2 nodes of ext_nodes that belongs to the same face
                //we extrapolate the distance in n_j
                std::vector<std::pair<TCellID, TCellID> > node_pairs;
                getNodesInTheSameFace(ext_nodes, node_pairs);
                if (!node_pairs.empty())
                {
                    m_mesh->mark(n_j, AMarkNodeNarrow);
                    Node from1 = m_mesh->get<Node>(node_pairs[0].first);
                    Node from2 = m_mesh->get<Node>(node_pairs[0].second);
                    double d_j = (this->*m_extrapolate_func_ptr)(n_j, from1, from2, m_distance,false);
                    if (d_j < 0)
                    {
                        d_j = 0;
                        std::cout << "Warning, we get a negative distance (2D)"
                        << std::endl;
                        AExtrapolateDistance[n_j.getID()] = 10000;
                        ANarrowBand.insert(std::pair<double, Node>(d_j, n_j));
                    }
                    else
                    {
                        AExtrapolateDistance[n_j.getID()] = d_j;
                        ANarrowBand.insert(std::pair<double, Node>(d_j, n_j));
                    }
                    
                }
            } //else if (nb_alive_adj_nodes > 2)
        }
    }
    
    std::cout<<"We deal with the source nodes now"<<std::endl;
    for(unsigned int i=0; i<source_nodes.size(); i++){
        Node source = source_nodes[i];
        std::vector<Edge> adj_edges;
        source.get<Edge>(adj_edges);
        
        //=======================================================
        //We get  all the adjacent nodes
        //=======================================================
        std::vector<Node> adj_nodes;
        for(unsigned int i_e=0; i_e<adj_edges.size();i_e++){
            std::vector<Node> current_nodes = adj_edges[i_e].get<Node>();
            if(current_nodes[0].getID()==source.getID())
                adj_nodes.push_back(current_nodes[1]);
            else
                adj_nodes.push_back(current_nodes[0]);

        }
        //=======================================================
        // We pair directions
        //=======================================================
        math::Cross2D ref_cross = (*m_cross_field)[adj_nodes[0].getID()];
        std::vector<math::Vector> ref_vectors = ref_cross.componentVectors();
        math::Vector ref_u = ref_vectors[0];
        math::Vector ref_v = ref_vectors[1];
        (*m_adv_field[0])[adj_nodes[0].getID()]=ref_u;
        (*m_adv_field[1])[adj_nodes[0].getID()]=ref_v;
    
        for(unsigned int j=1; j<adj_nodes.size();j++)
        {
            Node n_j = adj_nodes[j];
            
            math::Cross2D c_j = (*m_cross_field)[n_j.getID()];
            std::vector<math::Vector> vec_j = c_j.componentVectors();

            math::Vector u_j = vec_j[0];
            math::Vector v_j = vec_j[0];
            double du = u_j.dot(ref_u);
            double dv = v_j.dot(ref_v);
            for(int k=1;k<4;k++){
                
                math::Vector u_tmp = vec_j[k];
                double duk = u_tmp.dot(ref_u);
                double dvk = u_tmp.dot(ref_v);
                if(duk>du){
                    du=duk;
                    u_j=u_tmp;
                }
                if(dvk>dv){
                    dv=dvk;
                    v_j=u_tmp;
                }
            }
            
            (*m_adv_field[0])[n_j.getID()]=u_j;
            (*m_adv_field[1])[n_j.getID()]=v_j;
        }
        
        //=======================================================
        // Directions are paired, initial distance can be
        // computed
        //=======================================================
        for(unsigned int j=0; j<adj_nodes.size();j++)
        {
            Node n_j = adj_nodes[j];
            m_mesh->mark(n_j, AMarkNodeNarrow);
            
            math::Vector dir_j(source.getPoint(), n_j.getPoint());
            // we keep the vector in mind
            math::Vector uj =(*m_adv_field[0])[n_j.getID()];
            math::Vector vj =(*m_adv_field[1])[n_j.getID()];

            
            TCoord du_inf = fabs(uj.dot(dir_j));
            TCoord dv_inf = fabs(vj.dot(dir_j));
            
            
            (*m_adv_dist[0])[n_j.getID()]=du_inf;
            (*m_adv_dist[1])[n_j.getID()]=dv_inf;
            
            double d_j = math::max2(du_inf, dv_inf);
            
            
            if (d_j < 0) {
                d_j = 0;
                std::cout << "Warning, we get a negative distance (2D)"
                << std::endl;
                AExtrapolateDistance[n_j.getID()] = 10000;
                ANarrowBand.insert(std::pair<double, Node>(d_j, n_j));
            }
            else
            {
                AExtrapolateDistance[n_j.getID()] = d_j;
                ANarrowBand.insert(std::pair<double, Node>(d_j, n_j));
            }
        }

    }
    //	std::cout << "Narrow Band Size: " << ANarrowBand.size() << std::endl;
}

/*----------------------------------------------------------------------------*/
void DistanceFieldBuilder2D::advanceDistanceFront(
                                                  std::multimap<double,Node>& ANarrowBand,
                                                  std::map<TCellID, double>& AExtrapolateDistance,
                                                  const int AMarkNodeToWorkOn,
                                                  const int AMarkNodeAlive,
                                                  const int AMarkNodeNarrow,
                                                  const int AMarkFace)
{
    while (!ANarrowBand.empty())
    {
        //we take the node with the smallest distance computed in
        std::multimap<double, Node>::iterator it = ANarrowBand.begin();
        Node trial_node = it->second;
        double trial_dist = it->first;
        (*m_distance)[trial_node.getID()] = trial_dist;
        m_insertion_order.push_back(trial_node.getID());
        ANarrowBand.erase(it);
        
        //the trial node becomes alive ...
        m_mesh->mark(trial_node, AMarkNodeAlive);
        //... and is removed from the narrow band
        m_mesh->unmark(trial_node, AMarkNodeNarrow);
        
        //adjacent nodes that are not alive or in tha narrow band are added
        //in the narrow band if they are now adjacent to 2 alive nodes
        std::vector<Node> trial_adj_nodes =
        getAdjacentNodes(trial_node, AMarkNodeToWorkOn,AMarkFace);
        int nb_neigbor_alive_or_narrow = 0;
        for (unsigned int j = 0; j < trial_adj_nodes.size(); j++)
        {
            Node n_j = trial_adj_nodes[j];
            // we look for nodes that are not yet in the narrow band or alive
            if (m_mesh->isMarked(n_j, AMarkNodeNarrow) || m_mesh->isMarked(n_j, AMarkNodeAlive))
            {
                nb_neigbor_alive_or_narrow++;
                continue;
            }
            std::vector<Node> adj_nodes_j = getAdjacentNodes(n_j, AMarkNodeToWorkOn,AMarkFace);
            std::vector<Node> ext_nodes;
            int nb_alive_adj_nodes = 0;
            for (unsigned int k = 0; k < adj_nodes_j.size(); k++)
            {
                Node n_k = adj_nodes_j[k];
                if (m_mesh->isMarked(n_k, AMarkNodeAlive))
                {
                    nb_alive_adj_nodes++;
                    ext_nodes.push_back(n_k);
                }
            }
            if (nb_alive_adj_nodes == 2)
            {
                // we have to be sure that n_j, ext_nodes[0], ext_nodes[1]
                // belong to the same face. Iy yes we can put n_j in the
                // narrow band.
                Face commonFace;
                if (belongToTheSameFace(n_j, ext_nodes[0], ext_nodes[1], commonFace))
                {
                    
                    m_mesh->mark(n_j, AMarkNodeNarrow);
                    
                    //we extrapolate the distance in n_j
                    double d_j = (this->*m_extrapolate_func_ptr)(n_j, ext_nodes[0], ext_nodes[1], m_distance,false);
                    if (d_j < 0)
                    {
                        d_j = 0;
                        std::cout << "Warning, we get a negative distance (2D)"
                        << std::endl;
                        AExtrapolateDistance[n_j.getID()] = 10000;
                        ANarrowBand.insert(std::pair<double, Node>(d_j, n_j));
                    }
                    else
                    {
                        AExtrapolateDistance[n_j.getID()] = d_j;
                        ANarrowBand.insert(std::pair<double, Node>(d_j, n_j));
                    }
                }
            }
            else if (nb_alive_adj_nodes > 2)
            {
                
                //we look for 2 nodes of ext_nodes that belongs to the same face
                //we extrapolate the distance in n_j
                std::vector<std::pair<TCellID, TCellID> > node_pairs;
                getNodesInTheSameFace(ext_nodes, node_pairs);
                if (!node_pairs.empty())
                {
                    m_mesh->mark(n_j, AMarkNodeNarrow);
                    Node from1 = m_mesh->get<Node>(node_pairs[0].first);
                    Node from2 = m_mesh->get<Node>(node_pairs[0].second);
                    double d_j = (this->*m_extrapolate_func_ptr)(n_j, from1, from2, m_distance,false);
                    if (d_j < 0)
                    {
                        d_j = 0;
                        std::cout << "Warning, we get a negative distance (2D)"
                        << std::endl;
                        AExtrapolateDistance[n_j.getID()] = 10000;
                        ANarrowBand.insert(std::pair<double, Node>(d_j, n_j));
                    }
                    else
                    {
                        AExtrapolateDistance[n_j.getID()] = d_j;
                        ANarrowBand.insert(std::pair<double, Node>(d_j, n_j));
                    }
                    
                }
            }
        }
        
    }
}

/*----------------------------------------------------------------------------*/
void DistanceFieldBuilder2D::
internalComputeDistance(std::vector<Node>& AFrom,
                        std::vector<Node>& AToCompute,
                        const int AMarkFaceOnBoundary)
{
    //======================================================================
    // Boolean marks definition
    //======================================================================
    // a mark is given to the node to be computed
    int node_to_use = m_mesh->getNewMark<Node>();
    //all the nodes a value is computed are marked node_alive
    int node_alive  = m_mesh->getNewMark<Node>();
    //all the nodes that are in the narrow band
    int node_narrow  = m_mesh->getNewMark<Node>();
    
    //======================================================================
    // Initialization - ALIVE NODES
    //======================================================================
    //we mark all the nodes that we want to compute a distance on
    for (unsigned int i = 0; i < AToCompute.size(); i++)
    {
        Node n_i = AToCompute[i];
        m_mesh->mark(n_i, node_to_use);
    }
    
    for (unsigned int i = 0; i < AFrom.size(); i++)
    {
        Node n_i = AFrom[i];
        //initial nodes are evaluated to zero
        (*m_distance)[n_i.getID()] = 0;
        // the node is also marked as done
        m_mesh->mark(n_i, node_alive);
        m_mesh->mark(n_i, node_to_use);
    }
    // now, all the nodes where a distance is assigned are marked alive,
    // and all the others are not marked alive.
    
    //======================================================================
    // Initialization - NARROW BAND NODES
    //======================================================================
    // We build the narrow band list, which contains all the nodes we are going
    // to work on. They are adjacent to at 2 alive nodes, while being
    // not alive.
    std::multimap<double,Node> narrow_band;
    // we keep in mind the extrapolate distance for any node in the narrow
    // band
    std::map<TCellID,double> extrapolate_distance;
    
    
    initNarrowBand(AFrom, narrow_band, extrapolate_distance,
                   node_to_use, node_alive, node_narrow, AMarkFaceOnBoundary);
    
    
    //======================================================================
    // Advancing front loop
    //======================================================================
    advanceDistanceFront(narrow_band, extrapolate_distance,
                         node_to_use, node_alive, node_narrow, AMarkFaceOnBoundary);
    
    //======================================================================
    // Cleaning of the Boolean marks
    //======================================================================
    m_mesh->unmarkAll<Node>(node_alive);
    m_mesh->unmarkAll<Node>(node_narrow);
    m_mesh->unmarkAll<Node>(node_to_use);
    m_mesh->freeMark<Node>(node_alive);
    m_mesh->freeMark<Node>(node_narrow);
    m_mesh->freeMark<Node>(node_to_use);
}

/*----------------------------------------------------------------------------*/
Variable<double>* DistanceFieldBuilder2D::
computeDistance(std::vector<Node>& AFrom, std::vector<Node>& AToCompute,
                const int AMarkFaceOnBoundary)
{
    m_insertion_order.clear();
    m_distance = m_mesh->newVariable<double>(GMDS_NODE, "distance");
    m_extrapolate_func_ptr = &DistanceFieldBuilder2D::extrapolateDistance;
    internalComputeDistance(AFrom,AToCompute,AMarkFaceOnBoundary);
    return m_distance;
}
/*----------------------------------------------------------------------------*/
Variable<double>* DistanceFieldBuilder2D::
computeLInfDistanceFromCrossField(std::vector<Node>& AFrom,
                                  std::vector<Node>& AToCompute,
                                  Variable<math::Cross2D>* AField,
                                  const int AMarkFaceOnBoundary)
{
    m_insertion_order.clear();
    m_distance = m_mesh->newVariable<double>(GMDS_NODE, "distance");
    m_cross_field = AField;
    m_adv_field[0]= m_mesh->newVariable<math::Vector>(GMDS_NODE, "adv_dir0");
    m_adv_field[1]= m_mesh->newVariable<math::Vector>(GMDS_NODE, "adv_dir1");
    m_adv_dist[0]= m_mesh->newVariable<double>(GMDS_NODE, "adv_dist0");
    m_adv_dist[1]= m_mesh->newVariable<double>(GMDS_NODE, "adv_dist1");
    m_adv_order= m_mesh->newVariable<int>(GMDS_NODE, "advancing_order");
    m_extrapolate_func_ptr = &DistanceFieldBuilder2D::extrapolateLInfDistance;
    internalComputeDistance(AFrom,AToCompute,AMarkFaceOnBoundary);
    
    for(unsigned int i=0;i<m_insertion_order.size();i++){
        (*m_adv_order)[m_insertion_order[i]]=i;
    }
    return m_distance;
}

/*----------------------------------------------------------------------------*/
