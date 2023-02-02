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
#include <GMDS/Utils/Timer.h>
#include <GMDS/Algo/BoundaryOperator.h>
#include <GMDS/Utils/Log.h>
/*---------------------------------------------------------------------------*/
// FHeDo File Headers
#include "FieldGenerator.h"

/*---------------------------------------------------------------------------*/
using namespace gmds;
using namespace fhedo;
///*----------------------------------------------------------------------------*/
//// for hlbgs optimization
//Eigen::VectorXd feasible_objective_vector;
///*---------------------------------------------------------------------------*/
//// for alglib optimization
//math::SHVector ALGLIB_SH_Target;
///*----------------------------------------------------------------------------*/
///*
// * IN N -> Nb variable to handle. Here 3(alpha, beta, gamma)
// * IN x -> 3-sized vectors containing the Euler angles
// * IN prev_x -> Not used - contains previous x value
// * OUT func-> function to minimize
// * OUT grad -> grad(func)
// */
//
//void evalFuncAndGrad(const alglib::real_1d_array &AX,
//                     double &AFunc,
//                     alglib::real_1d_array &AGrad,
//                     void *APtr) {
//    double alpha = AX[0];
//    double beta  = AX[1];
//    double gamma = AX[2];
//    
//    math::SHFrame f;
//    f.rotateXYZ(alpha, beta, gamma);
//    
//    math::SHVector sh_diff = f.a() - ALGLIB_SH_Target;
//    
//    AFunc = sh_diff.squaredNorm();
//    
//    Eigen::Matrix<double, 3, 9> diff;
//    
//    const math::SHMatrix RX  = math::SH::RX (alpha);
//    const math::SHMatrix RY  = math::SH::RY (beta );
//    const math::SHMatrix RZ  = math::SH::RZ (gamma);
//    const math::SHMatrix DRX = math::SH::DRX(alpha);
//    const math::SHMatrix DRY = math::SH::DRY(beta );
//    const math::SHMatrix DRZ = math::SH::DRZ(gamma);
//    const math::SHVector I   = math::SH::Identity();
//    
//    math::SHVector DSH_alpha = DRX*  RY *  RZ * I;
//    math::SHVector DSH_beta  = RX * DRY *  RZ * I;
//    math::SHVector DSH_gamma = RX * RY  * DRZ * I;
//    
//    diff << DSH_alpha.transpose(), DSH_beta .transpose(), DSH_gamma.transpose();
//    
//    Eigen::Vector3d deriv = 2 * diff * sh_diff;
//    AGrad[0] = deriv[0];
//    AGrad[1] = deriv[1];
//    AGrad[2] = deriv[2];
//    
//}
///*----------------------------------------------------------------------------*/
///*
// * IN N -> Nb variable to handle. Here 3(alpha, beta, gamma)* nb free nodes
// * IN x -> 3-sized vectors containing the Euler angles * nb free nodes
// * IN prev_x -> Not used - contains previous x value
// * OUT func-> function to minimize
// * OUT grad -> grad(func)
// */
//
//void evalGlobalFuncAndGrad(int N, double* x, double *prev_x,
//                           double* func, double* grad)
//{
//    int nb_free_nodes = N/3;
//    for(int i =0; i< nb_free_nodes; i++){
//        
//        double alpha = x[3*i];
//        double beta  = x[3*i+1];
//        double gamma = x[3*i+2];
//        
//        math::SHFrame f;
//        f.rotateXYZ(alpha, beta, gamma);
//        
//        math::SHVector ref_i;
//        for(int d=0;d<9;d++)
//            ref_i[d]= feasible_objective_vector[9*i+d];
//        math::SHVector sh_diff = f.a() - ref_i;
//        
//        *func = sh_diff.squaredNorm();
//        
//        Eigen::Matrix<double, 3, 9> diff;
//        
//        const math::SHMatrix RX  = math::SH::RX (alpha);
//        const math::SHMatrix RY  = math::SH::RY (beta );
//        const math::SHMatrix RZ  = math::SH::RZ (gamma);
//        const math::SHMatrix DRX = math::SH::DRX(alpha);
//        const math::SHMatrix DRY = math::SH::DRY(beta );
//        const math::SHMatrix DRZ = math::SH::DRZ(gamma);
//        const math::SHVector I   = math::SH::Identity();
//        
//        math::SHVector DSH_alpha = DRX*  RY *  RZ * I;
//        math::SHVector DSH_beta  = RX * DRY *  RZ * I;
//        math::SHVector DSH_gamma = RX * RY  * DRZ * I;
//        
//        diff << DSH_alpha.transpose(), DSH_beta .transpose(), DSH_gamma.transpose();
//        
//        Eigen::Vector3d deriv = 2 * diff * sh_diff;
//        grad[3*i] = deriv[0];
//        grad[3*i+1] = deriv[1];
//        grad[3*i+2] = deriv[2];
//    }
//    
//}
//
///*----------------------------------------------------------------------------*/
//void displayHLBGFSIteration(int iter, int call_iter, double *x, double* f, double *g, double* gnorm)
//{
//    std::cout << iter << ": " << call_iter << " " << *f << " " << *gnorm << std::endl;
//}
/*---------------------------------------------------------------------------*/
FieldGenerator::
FieldGenerator(FieldSolverStrategyItf* ASolver,
               IGMesh* AMesh,
               const ParamsGlobal& AGParam,
               const ParamsFrameField& AParam,
               const ParamsMark& ABM):
m_solver(ASolver),
m_mesh(AMesh),
m_global_params(AGParam),
m_params(AParam),
m_bm(ABM),
m_nb_inner_nodes(0),
m_nb_fixed_nodes(0),
m_nb_surface_nodes(0),
m_nb_free_nodes(0),
m_nb_smooth_edges(0)
{
    /** Definition of variable on the mesh nodes */
    m_harmonic_field = m_mesh->newVariable<math::SHarmonicL4>(GMDS_NODE, "SHL4"    );
    m_chart_field    = m_mesh->newVariable<math::Chart      >(GMDS_NODE, "SHChart" );
    m_ordering       = m_mesh->newVariable<int              >(GMDS_NODE, "Ordering");
    
}
/*---------------------------------------------------------------------------*/
double FieldGenerator::computeFieldEnergy(){

    double e=0;
    for(IGMesh::edge_iterator ite=m_mesh->edges_begin();
        !ite.isDone();
        ite.next()) {
        Edge ei = ite.value();
        std::vector<TCellID> n_ids=ei.getIDs<Node>();
        
        math::SHarmonicL4 sh_e = ((*m_harmonic_field)[n_ids[0]]-
                                  (*m_harmonic_field)[n_ids[1]]);
        e+=(sh_e).dot(sh_e);
        
    }

    return e;
}
/*---------------------------------------------------------------------------*/
void FieldGenerator::execute(){
    
    //==================================================================
    // INITIALIZATION, boundary cells shoud be marked outside of this
    // algorithm (beforehand)
    //==================================================================
    // Detection of the boundary cells
    
    initMarks();
    Log::mng()<< "FF GEN > Marks boundary \n";
    Log::mng().flush();
    markBoundaryCells();
    
    //==================================================================
    // NODES ARE REORDERED TO ASSEMBLE THE SYSTEM MATRIX
    //==================================================================
    // Node ordering
    Log::mng()<< "FF GEN > Node sorting \n";
    Log::mng().flush();
    sortNodes();
    
    // Parse edges
    //    std::cout << "> Edge parsing " << std::endl;
    //    m_nb_smooth_edges = 0;
    //    m_nb_strict_smooth_edges=0;
    //    for(IGMesh::edge_iterator ite = m_mesh->edges_begin();
    //        !ite.isDone();ite.next()){
    //        Edge e = ite.value();
    //        std::vector<Node> e_nodes = e.get<Node>();
    //        if((!m_mesh->isMarked(e_nodes[0],m_markNodeFixed)) ||
    //           (!m_mesh->isMarked(e_nodes[1],m_markNodeFixed))){
    //            m_mesh->mark(e,m_markEdgeSmooth);
    //            m_nb_smooth_edges++;
    //        }
    //        if(!m_mesh->isMarked(e_nodes[0],m_markNodeFixed) &&
    //           !m_mesh->isMarked(e_nodes[1],m_markNodeFixed)){
    //            m_nb_strict_smooth_edges++;
    //        }
    //    }
    
    m_nb_smooth_edges = m_mesh->getNbEdges();

    //==================================================================
    // BOUNDARY DATA, FIXED SH ON CURVES AND POINTS AND FIXED NORMALS
    // ON SURFACES
    //==================================================================
    Log::mng()<< "FF GEN > Boundary constraint assignment \n";
    Log::mng().flush();

    computeBoundaryData();
    
    if(m_global_params.with_debug_files){
        writeSolution();
    }
    
    int nb_nodes     = m_mesh->getNbNodes();
    int nb_bnd_nodes = m_nb_fixed_nodes+m_nb_surface_nodes;
    m_nb_inner_nodes = nb_nodes-nb_bnd_nodes;
    m_nb_free_nodes  = nb_nodes;//-m_nb_fixed_nodes;
    
    //==================================================================
    // Initialization of most of the mandatory field for the solver
    m_solver->setNbFreeNodes(m_nb_free_nodes);
    m_solver->setMesh(m_mesh);
    m_solver->setMarkNodeLocked(m_markNodeFixed);
    m_solver->setHarmonicVariable(m_harmonic_field);
    m_solver->setChartVariable(m_chart_field);
    m_solver->setOrderingVariable(m_ordering);
    
    m_solver->setH0(&m_H0);
    m_solver->setH4(&m_H4);
    m_solver->setH8(&m_H8);
    m_solver->setSurfaceNodes(&m_surface_nodes);
    //==================================================================
    // INITIALIZATION STEP
    //==================================================================
    std::cout<<"========================================" << std::endl;
    std::cout<<"Field Initialization"<<std::endl;
  

    
    //  9D representation coef for each node
    //+ 2D alignment vector for bound. nodes

    int nb_unknowns =9*m_nb_free_nodes + 2*m_nb_surface_nodes;
    
    m_solver->setNbUnknowns(nb_unknowns);

   

    buildSystem(0);
    std::cout << "> Solving system with "<<nb_unknowns<<" unknowns"
    << std::endl;
    m_solver->solve();
    m_solver->getFeasibleSolution();
    writeSolution();

    /* Cleanup */
    std::cout << "> Cleaning system"<< std::endl;
    m_solver->clean();

    //==================================================================
    // SMOOTHING
    //==================================================================
    int nb_iterations = m_params.smoothing_nb_iter;
    
    //  9D representation coef for each node
    //+ 2D alignment vector for bound. nodes
    //+ 3D feasibility constraints, only for optimization
    nb_unknowns =9*m_nb_free_nodes + 2*m_nb_surface_nodes
    + 3*m_nb_free_nodes;

    m_solver->setNbUnknowns(nb_unknowns);
    
    
    bool stop =false;
    int iteration_i=0;
    for(; !stop && iteration_i<nb_iterations; iteration_i++){


        buildSystem(iteration_i);
        std::cout << "> Solving system with "<<nb_unknowns<<" unknowns"
        << std::endl;
        m_solver->solve();
        m_solver->getFeasibleSolution();


//        if(iteration_i<10 || iteration_i%10==0){
//            std::cout<<">>>>>>>>>> Write for iteration "
//            <<iteration_i<<std::endl;
//            writeSolution();
//        }
        std::cout << "> Cleaning system"<< std::endl;
        m_solver->clean();
    }
    
    std::cout<<"Converge after solving "
    <<iteration_i<<" least-squared systems"<<std::endl;
    writeSolution();
}
/*---------------------------------------------------------------------------*/
void FieldGenerator::
buildSystem(const int AI){
    std::cout << "> Building system"<< std::endl;

    m_solver->initializeAssembly();
    
    if(AI>0){
        m_solver->setX();
    }
    
    m_solver->addLockedTerms();
    
    //We fill in the first 9*nb_nodes lines of A and b
    std::cout << "\t Add smoothing terms to A and b" << std::endl;
    m_solver->addSmoothingTerms();

    //We fill in the next 2*nb_bnd_nodes lines of A and b
    std::cout << "\t Add boundary constraint terms to A and b" << std::endl;
    m_solver->addBoundaryConstraints();
    
    if(AI>0){
        m_solver->addLocalOptimConstraints();
    }
    
    m_solver->finalizeAssembly();
    
}
/*---------------------------------------------------------------------------*/
void FieldGenerator::initMarks(){
    //create node marks
    m_markNodeOnBnd = m_mesh->getNewMark<Node>();
    m_markNodeFixed = m_mesh->getNewMark<Node>();
    
    //create edge marks
    m_markEdgeSmooth = m_mesh->getNewMark<Edge>();
    
 }
/*---------------------------------------------------------------------------*/
void FieldGenerator::cleanMarks()
{
    //clean node marks
    m_mesh->unmarkAll<Node>(m_markNodeOnBnd);
    m_mesh->unmarkAll<Node>(m_markNodeFixed);
    //clean edge marks
    m_mesh->unmarkAll<Edge>(m_markEdgeSmooth);
    
    //free node marks
    m_mesh->freeMark<Node>(m_markNodeOnBnd);
    m_mesh->freeMark<Node>(m_markNodeFixed);
    //free edge marks
    m_mesh->freeMark<Edge>(m_markEdgeSmooth);    
}
/*---------------------------------------------------------------------------*/
void FieldGenerator::markBoundaryCells(){
   
    
    IGMesh::node_iterator it_node = m_mesh->nodes_begin();
    for (; !it_node.isDone(); it_node.next())
    {
        Node current_node = it_node.value();
        if (m_mesh->isMarked(current_node, m_bm.mark_edge_on_surf) ||
            m_mesh->isMarked(current_node, m_bm.mark_node_on_curv) ||
            m_mesh->isMarked(current_node, m_bm.mark_node_on_pnt)) {
            m_mesh->mark(current_node, m_markNodeOnBnd);
            
        }
        if (m_mesh->isMarked(current_node, m_bm.mark_node_on_surf) &&
            (!m_mesh->isMarked(current_node, m_bm.mark_node_on_curv)) &&
            (!m_mesh->isMarked(current_node, m_bm.mark_node_on_pnt))){
            m_surface_nodes.push_back(current_node);
        }
        
        if (m_mesh->isMarked(current_node, m_bm.mark_node_on_curv) ||
            m_mesh->isMarked(current_node, m_bm.mark_node_on_pnt)) {
            m_mesh->mark(current_node, m_markNodeFixed);
        }
    }
}

/*---------------------------------------------------------------------------*/
void FieldGenerator::sortNodes()
{
    int i=0;
    
    for (IGMesh::node_iterator it_node = m_mesh->nodes_begin();
         !it_node.isDone(); it_node.next()) {
        Node current_node = it_node.value();
        if(m_mesh->isMarked(current_node, m_bm.mark_node_on_surf) &&
           (!m_mesh->isMarked(current_node, m_markNodeFixed)))
        {
            (*m_ordering)[current_node.getID()] = i++;
        }
    }

    //==================================================================
    //STEP 1- Traversal to list ALL boundary nodes
    //==================================================================
    
    // We add the fixed nodes in the variables
    for (IGMesh::node_iterator it_node = m_mesh->nodes_begin();
         !it_node.isDone(); it_node.next()) {
        Node current_node = it_node.value();
        if(m_mesh->isMarked(current_node, m_markNodeFixed))
        {
            (*m_ordering)[current_node.getID()] = i++;
        }
    }
    
    //==================================================================
    //STEP 2- Traversal to list inner nodes
    //==================================================================
    for (IGMesh::node_iterator it_node = m_mesh->nodes_begin();
         !it_node.isDone(); it_node.next()) {
        Node current_node = it_node.value();
        if(!m_mesh->isMarked(current_node, m_markNodeOnBnd)){
            (*m_ordering)[current_node.getID()] = i++;
        }
    }
    
}

/*---------------------------------------------------------------------------*/
void FieldGenerator::computeBoundaryData()
{
    m_nb_fixed_nodes=0;
    m_nb_surface_nodes=0;
    
    
    std::vector<Node> nodes_on_points;
    //====================================================================
    //STEP 1- First traversal to initialize frames on curve and surfaces,
    // nodes classsified on vertices are also kept in mind since frames
    // will be assigned to them in a second step
    //====================================================================
    for (IGMesh::node_iterator it_node = m_mesh->nodes_begin();
         !it_node.isDone(); it_node.next()) {
        Node current_node = it_node.value();
        
        //======================================================
        // Node on a geometrical point
        if(m_mesh->isMarked(current_node, m_bm.mark_node_on_pnt)){
            nodes_on_points.push_back(current_node);
            m_nb_fixed_nodes++;
            
        }
        //======================================================
        // Node on a geometrical curve
        else if( m_mesh->isMarked(current_node, m_bm.mark_node_on_curv)&&
                !m_mesh->isMarked(current_node, m_bm.mark_node_on_pnt)){
            initHarmonicOnCurve(current_node);
            m_nb_fixed_nodes++;
        }
        //======================================================
        // Node on a geometrical surface
        else if(m_mesh->isMarked(current_node, m_bm.mark_node_on_surf)){
            //We are on a surface node
            initHarmonicOnSurface(current_node);
            m_nb_surface_nodes++;
        }
        //Otherwise nothing to do.
    }
    for(unsigned int i=0; i<nodes_on_points.size();i++){
        Node current_node = nodes_on_points[i];
        initHarmonicOnPoint(current_node);
    }
    
}
/*---------------------------------------------------------------------------*/
void  FieldGenerator::initHarmonicOnPoint(Node &ANode)
{
    std::vector<Node> curve_nodes;
    std::vector<Edge> curve_edges = getEdgesOnCurve(ANode);
    
    for (unsigned int j = 0; j < curve_edges.size(); j++) {
        Edge ej = curve_edges[j];
        curve_nodes.push_back(getNeighboorOn(ANode, ej));
    }
    
    std::vector<math::Quaternion> adj_quat;
    std::vector<TCoord> adj_coef;
    adj_quat.resize(curve_nodes.size());
    adj_coef.resize(curve_nodes.size());
    for (unsigned int j = 0; j < curve_nodes.size(); j++){
        math::Chart cj = (*m_chart_field)[curve_nodes[j].getID()];
        math::Quaternion qj(cj);
        adj_quat[j]=qj;
        adj_coef[j]=1;
    }
    
    math::Quaternion q= math::Quaternion::mean(adj_quat, adj_coef);
    math::Chart c(q);
    
    (*m_harmonic_field)[ANode.getID()] = math::SHarmonicL4(c);
    (*m_chart_field)   [ANode.getID()] = c;

}
/*---------------------------------------------------------------------------*/
void  FieldGenerator::initHarmonicOnSurface(Node &ANode)
{
    TCellID node_id = ANode.getID();
    BoundaryOperator boundaryOp(m_mesh);
    math::Vector n = boundaryOp.getOutputNormalOfABoundaryNode(ANode);
    
    
    //find Euler angle to rotate Oz onto n
    //OZ -> n right rotation, checked several times and debugged
    // NOT n-> OZ
    //    vec3 nr(n.X(),n.Y(),n.Z());
    //    AxisAngleRot rot = bring_z_to_vec(nr);
    //
    //    SphericalHarmonicL4 h0(1, 0, 0 , 0, 0, 0 ,0, 0, 0);
    //    SphericalHarmonicL4 h4(0, 0, 0 , 0, 1, 0 ,0, 0, 0);
    //    SphericalHarmonicL4 h8(0, 0, 0 , 0, 0, 0 ,0, 0, 1);
    //    h0.Rot(rot.v());
    //    h4.Rot(rot.v());
    //    h8.Rot(rot.v());
    //    m_H0[node_id]= math::SHarmonicL4( h0[0],h0[1],h0[2],h0[3],h0[4],h0[5],h0[6],h0[7],h0[8]);
    //    m_H4[node_id]= math::SHarmonicL4(h4[0],h4[1],h4[2],h4[3],h4[4],h4[5],h4[6],h4[7],h4[8]);
    //    m_H8[node_id]= math::SHarmonicL4( h8[0],h8[1],h8[2],h8[3],h8[4],h8[5],h8[6],h8[7],h8[8]);
    
    math::AxisAngleRotation rot =
    math::AxisAngleRotation::alignZ(n);
    
    math::Matrix<3,3,double> m_rot =rot.toRotationMatrix();
    math::Vector3d xyz = m_rot.eulerAngles(0,1,2);
    
    double alpha = xyz.X();
    double beta  = xyz.Y();
    double gamma = xyz.Z();
    //With new classes
    math::SHL4Matrix m_sh = math::fromXYZEulerAngle(alpha,beta,gamma);
    
    math::Vector9d v0=m_sh*math::SHarmonicL4(1,0,0,0,0,0,0,0,0);
    math::Vector9d v4=m_sh*math::SHarmonicL4(0,0,0,0,1,0,0,0,0);
    math::Vector9d v8=m_sh*math::SHarmonicL4(0,0,0,0,0,0,0,0,1);
    
    m_H0[node_id] = math::SHarmonicL4(v0[0],v0[1],v0[2],v0[3],v0[4],
                                      v0[5],v0[6],v0[7],v0[8]);
    
    m_H4[node_id] = math::SHarmonicL4(v4[0],v4[1],v4[2],v4[3],v4[4],
                                      v4[5],v4[6],v4[7],v4[8]);
    
    m_H8[node_id] = math::SHarmonicL4(v8[0],v8[1],v8[2],v8[3],v8[4],
                                      v8[5],v8[6],v8[7],v8[8]);
    
}
/*---------------------------------------------------------------------------*/
void FieldGenerator::initHarmonicOnCurve(Node &ANode)
{
    //for each node on a geometrical curve, we compute
    //its associated quaternion
    
    Node current_node = ANode;
    
    math::Chart current_chart;
    
    //current_node is on a geometric curve
    //We get adjacent edges that are classified onto
    //a geometric curve. We must get one or two edges
    //at most.
    std::vector<Edge> ridges = getEdgesOnCurve(current_node);
    math::Vector v1, newN1, newN2;
    if (ridges.size() == 1){
        Edge current_edge = ridges[0];
        Node node1 = getNeighboorOn(current_node, current_edge);
        Node node2 = current_node;
        
        //we build the direction vector of the current edge
        math::Point p1 = node1.getPoint();
        math::Point p2 = node2.getPoint();
        v1 = math::Vector(p1, p2);
        v1.normalize();
        //We retrieve adjacent boundary faces and we get their normal vectors
        std::vector<Face> boundary_faces = getFacesOnSurface(current_edge);
        
        math::Vector n1 = boundary_faces[0].normal();
        math::Vector n2 = boundary_faces[1].normal();
        n1.normalize();
        n2.normalize();
        
        //They are projected onto the tangent vector
        n1 = n1 - (v1.dot(n1) * v1);
        n2 = n2 - (v1.dot(n2) * v1);
        n1.normalize();
        n2.normalize();
        
        newN1 = n1;
        newN2 = n2;
    }
    else if (ridges.size() == 2){
        //With 2 adajcent edges on the curve, we compute average values
        
        Edge edge1 = ridges[0];
        Edge edge2 = ridges[1];
        
        Node node1 = getNeighboorOn(current_node, edge1);
        Node node2 = getNeighboorOn(current_node, edge2);
        //we build the direction vector linking the two adjacent nodes
        math::Point p1 = node1.getPoint();
        math::Point p2 = node2.getPoint();
        v1 = math::Vector(p1, p2);
        v1.normalize();
        
        //Normal are computed in the same way
        std::vector<Face> faces1 = getFacesOnSurface(edge1);
        std::vector<Face> faces2 = getFacesOnSurface(edge2);
        Face SelectedFace[4];
        SelectedFace[0] = faces1[0];
        SelectedFace[1] = faces1[1];
        SelectedFace[2] = faces2[0];
        SelectedFace[3] = faces2[1];
        
        math::Vector n1 = SelectedFace[0].normal();
        math::Vector n2 = SelectedFace[1].normal();
        math::Vector n3 = SelectedFace[2].normal();
        math::Vector n4 = SelectedFace[3].normal();
        
        n1.normalize();
        n2.normalize();
        n3.normalize();
        n4.normalize();
        
        
        //Calcul des produits scalaires
        const TCoord prodScal1 = abs(n1.dot(n3));
        const TCoord prodScal2 = abs(n1.dot(n4));
        
        //Correction
        if (prodScal1 < prodScal2)
        {
            Face fTmp = SelectedFace[2];
            SelectedFace[2] = SelectedFace[3];
            SelectedFace[3] = fTmp;
        }
        
        //Recuperation des vecteurs normaux.
        
        math::Vector normal1 = SelectedFace[0].normal();
        math::Vector normal2 = SelectedFace[1].normal();
        {
            math::Vector normal3 = SelectedFace[2].normal();
            math::Vector normal4 = SelectedFace[3].normal();
            //Mise dans le mememe sens des normales
            TCoord scal;
            scal = normal1.dot(normal3);
            if (scal < 0.0)
                normal3 = normal3.opp();
            scal = normal2.dot(normal4);
            if (scal < 0.0)
                normal4 = normal4.opp();
            
            //Moyenne des normales
            normal1 = normal1 + normal3;
            normal2 = normal2 + normal4;
            
            normal1.normalize();
            normal2.normalize();
        }
        //Projection des normales sur la tangente
        normal1 = normal1 - (v1.dot(normal1)* v1);
        normal2 = normal2 - (v1.dot(normal2)* v1);
        //Normalisation des vecteurs projetes
        normal1.normalize();
        normal2.normalize();
        
        newN1 = normal1;
        newN2 = normal2;
    }
    else{
        std::cout << "Nb ridges for an edge adjacent to node "
        << current_node.getID() << ": " << ridges.size() << std::endl;
        throw GMDSException("A ridge node has an illegal number of edges.");
    }
    //Creation of N1 prime
    math::Vector N1Prime = v1.cross(newN1);
    if (N1Prime.dot(newN2) < 0.0)
        N1Prime = N1Prime.opp();
    
    //Computation of v2 from N1' and N2
    math::Vector v2 = N1Prime + newN2;
    v2.normalize();
    
    //Computation of v3 from v1 and v2.
    const math::Vector v3 = v1.cross(v2);
    
    current_chart = math::Chart(v1, v2, v3);
    
    (*m_harmonic_field)[ANode.getID()] = math::SHarmonicL4(current_chart);
    (*m_chart_field)[ANode.getID()] = current_chart;
    
}
/*---------------------------------------------------------------------------*/
std::vector<Edge>  FieldGenerator::
getEdgesOnCurve(Node& ANode) const {
    std::vector<Edge> edges_on_curve;
    std::vector<Edge> adj_edges = ANode.get<Edge>();
    for (unsigned int i = 0; i < adj_edges.size(); i++)
    {
        Edge ei = adj_edges[i];
        if (m_mesh->isMarked(ei, m_bm.mark_edge_on_curv))
            edges_on_curve.push_back(ei);
    }
    return edges_on_curve;
}
/*---------------------------------------------------------------------------*/
std::vector<Face>  FieldGenerator::
getFacesOnSurface(Edge& AEdge) const {
    std::vector<Face> faces_on_surf;
    std::vector<Face> adj_faces = AEdge.get<Face>();
    for (unsigned int i = 0; i < adj_faces.size(); i++)
    {
        Face fi = adj_faces[i];
        if (m_mesh->isMarked(fi, m_bm.mark_face_on_surf))
            faces_on_surf.push_back(fi);
    }
    return faces_on_surf;
}
/*---------------------------------------------------------------------------*/
Node  FieldGenerator::getNeighboorOn(Node& AN, Edge& AE) const {
    std::vector<Node> nodes = AE.get<Node>();
    if (nodes[0].getID() == AN.getID())
        return nodes[1];
    
    return nodes[0];
}
/*---------------------------------------------------------------------------*/
void FieldGenerator::writeSolution(){
    static int nb_file = 0;
    
    IGMesh::node_iterator it = m_mesh->nodes_begin();
    double bound = 100000;
    double x_min = bound;
    double y_min = bound;
    double z_min = bound;
    double x_max = -bound;
    double y_max = -bound;
    double z_max = -bound;
    for (; !it.isDone(); it.next())
    {
        Node n = it.value();
        math::Point p = n.getPoint();
        if (p.X() < x_min)
            x_min = p.X();
        if (p.X() > x_max)
            x_max = p.X();
        
        if (p.Y() < y_min)
            y_min = p.Y();
        if (p.Y() > y_max)
            y_max = p.Y();
        
        if (p.Z() < z_min)
            z_min = p.Z();
        if (p.Z() > z_max)
            z_max = p.Z();
    }
    double dist_x = x_max - x_min;
    double dist_y = y_max - y_min;
    double dist_z = z_max - z_min;
    
    double cube_size = 0;
    if (dist_x <= dist_y && dist_x <= dist_z){
        cube_size = dist_x;
    }
    else if (dist_y <= dist_x && dist_y <= dist_z){
        cube_size = dist_y;
    }
    else
        cube_size = dist_z;
    
    cube_size /= 20;
    
    MeshModel model_cube(DIM3 | R | N | R2N);
    IGMesh mesh_cube(model_cube);
    
    Variable<int>* v = mesh_cube.newVariable<int>(GMDS_REGION, "classification");
    
    
    for (it = m_mesh->nodes_begin(); !it.isDone(); it.next())
    {
        Node n = it.value();
        math::Point center = n.getPoint();
        
        math::Chart c = (*m_chart_field)[n.getID()];

        math::Vector3d evx(c.X()[0],c.X()[1],c.X()[2]);
        math::Vector3d evy(c.Y()[0],c.Y()[1],c.Y()[2]);
        math::Vector3d evz(c.Z()[0],c.Z()[1],c.Z()[2]);
        
        evx.normalize();
        evy.normalize();
        evz.normalize();
        
        if(evx.X()>100){
            std::cout<<"Error Nan for node "<<n.getID()<<std::endl;
            continue;
        }
        math::Vector vx(evx.X(),evx.Y(),evx.Z());
        math::Vector vy(evy.X(),evy.Y(),evy.Z());
        math::Vector vz(evz.X(),evz.Y(),evz.Z());
        
        math::Point p1 = center + (vx + vy - vz)*cube_size;
        
        Node n1 = mesh_cube.newNode(p1);
        math::Point p2 = center + (vx - vy - vz)*cube_size;
        Node n2 = mesh_cube.newNode(p2);
        math::Point p3 = center + (vx + vy + vz).opp()*cube_size;
        Node n3 = mesh_cube.newNode(p3);
        math::Point p4 = center + (vy - vx - vz)*cube_size;
        Node n4 = mesh_cube.newNode(p4);
        
        math::Point p5 = center + (vx + vy + vz)*cube_size;
        Node n5 = mesh_cube.newNode(p5);
        math::Point p6 = center + (vx - vy + vz)*cube_size;
        Node n6 = mesh_cube.newNode(p6);
        math::Point p7 = center + (vx + vy - vz).opp()*cube_size;
        Node n7 = mesh_cube.newNode(p7);
        math::Point p8 = center + (vy - vx + vz)*cube_size;
        Node n8 = mesh_cube.newNode(p8);
        //        std::cout<<"Cube from "<<std::endl;
        //        std::cout<<"\t "<<p1<<std::endl;
        //        std::cout<<"\t "<<p2<<std::endl;
        //        std::cout<<"\t "<<p3<<std::endl;
        //        std::cout<<"\t "<<p4<<std::endl;
        //        std::cout<<"\t "<<p5<<std::endl;
        //        std::cout<<"\t "<<p6<<std::endl;
        //        std::cout<<"\t "<<p7<<std::endl;
        //        std::cout<<"\t "<<p8<<std::endl;
        Region r = mesh_cube.newHex(n1, n2, n3, n4, n5, n6, n7, n8);
        if(m_mesh->isMarked(n, m_bm.mark_node_on_pnt))
            (*v)[r.getID()]=0;
        else if(m_mesh->isMarked(n, m_bm.mark_node_on_curv))
            (*v)[r.getID()]=1;
        else if(m_mesh->isMarked(n, m_bm.mark_node_on_surf))
            (*v)[r.getID()]=2;
        else
            (*v)[r.getID()]=3;
        
    }
    VTKWriter<IGMesh> writer_cube(mesh_cube);
    
    std::stringstream file_name_cube;
    file_name_cube<<m_global_params.output_dir+"/SH_DEBUG_CUBE_" << nb_file;
    writer_cube.write(file_name_cube.str(), DIM3 | R | N);
    
    
    
    Variable<int>* var_sing = 0;
    try{
        var_sing = m_mesh->getVariable<int>(GMDS_REGION, "sing_tet");
    }
    catch (GMDSException& e){
        var_sing = m_mesh->newVariable<int>(GMDS_REGION, "sing_tet");
    }
    
    IGMesh::region_iterator itr = m_mesh->regions_begin();
    
    int nbColoredTet = 0;
    //=========================================================================
    // INTERN SKELETON CREATION
    //=========================================================================
    // FIRST LOOP ON REGIONS TO GET ALL THE 3-SING. TETS
    //=========================================================================
    for (; !itr.isDone(); itr.next()){
        Region current_region = itr.value();
        std::vector<Node> nodes = current_region.get<Node>();
        bool onPnt = false;
        for (unsigned int i_node = 0; i_node < nodes.size(); i_node++){
            Node ni = nodes[i_node];
            if (m_mesh->isMarked(ni, m_bm.mark_node_on_pnt))
                onPnt = true;
        }
        if (onPnt){
            (*var_sing)[current_region.getID()] = 0;
        }
        else
        {
            std::vector<TCellID> node_ids = current_region.getIDs<Node>();
            math::Quaternion q[4];
            for(int i_n=0; i_n<4; i_n++){
                math::Chart ci= (*m_chart_field)[node_ids[i_n]];
                q[i_n]= math::Quaternion(ci);
            }
            
            int sing_type = math::Quaternion::testSingularity(q[0],q[1],q[2],q[3]);
            if (sing_type != 0)
                nbColoredTet++;
            
            (*var_sing)[current_region.getID()] = sing_type;
        }
    }
    std::cout << "Nb. colored tetrahedra: " << nbColoredTet << std::endl;
    
    
    
    VTKWriter<IGMesh> writer(*m_mesh);
    
    std::string file_name=m_global_params.output_dir+"/SH_DEBUG_"+to_string(nb_file);
    writer.write(file_name, DIM3 | R| F | N);
    
    
    nb_file++;
}

/*---------------------------------------------------------------------------*/


