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
// STL Header files
#include <iostream>
#include <map>
#include <set>
/*---------------------------------------------------------------------------*/
#include "CrossFieldSmoother2D.h"
/*---------------------------------------------------------------------------*/
// HLBFGS Header file
#include "HLBFGS.h"
#include <GMDS/IO/VTKWriter.h>
/*---------------------------------------------------------------------------*/
using namespace gmds;
/*---------------------------------------------------------------------------*/
namespace hlbfgs_for_cross2D {
  /*---------------------------------------------------------------------------*/
  // GLOBAL VARIABLES USED FOR THE CALL TO HLBFGS
  /*---------------------------------------------------------------------------*/
  gmds::IGMesh* HLBFGS_mesh = 0;
  // mark to identify the nodes that must be steady
  int HLBFGS_mark_steady = 0;
  // mark to identify the nodes that can be modified
  int HLBFGS_mark_smooth = 0;
  // mark to identify the nodes that can be modified
  double HLBFGS_global_area = 0;
  //the ordered list of candidate nodes
  std::vector<gmds::Node> HLBFGS_candidates;
  //the index of a candidate node in m_candidates
  std::map<gmds::TCellID, int> HLBFGS_candidates_index;
  //the list of edges that connect candidates these edges will be used to compute
  // smoothing values between end nodes
  std::vector<gmds::Edge> HLBFGS_edges;
  //for each candidate node, we precompute an ordered list of neigbours
  std::map<gmds::TCellID, std::vector<gmds::TCellID> > HLBFGS_vicinity_index;
  //... and some associated weights
  std::map<gmds::TCellID, std::vector<double> > HLBFGS_vicinity_weight;
  std::map<gmds::TCellID, double > HLBFGS_boundary_x;

  void writeMesh(const int i, double* x){
    double cube_size = 0.4;

    gmds::MeshModel model_cube(DIM3 | F | N | F2N);
    gmds::IGMesh mesh_cube(model_cube);
    for (gmds::IGMesh::node_iterator it = HLBFGS_mesh->nodes_begin(); !it.isDone(); it.next()) {
      gmds::Node n = it.value();
      gmds::math::Point center = n.getPoint();
      gmds::math::Cross2D current_cross;
      if(HLBFGS_mesh->isMarked(n, HLBFGS_mark_steady)) {
	current_cross =  HLBFGS_boundary_x[n.getID()];
      }
      else if (HLBFGS_mesh->isMarked(n, HLBFGS_mark_smooth)){ 
	current_cross = gmds::math::Cross2D(x[HLBFGS_candidates_index[n.getID()]]);
      }
      std::vector<gmds::math::Vector> current_vectors = current_cross.componentVectors();
      gmds::math::Vector vx = current_vectors[0];
      gmds::math::Vector vy = current_vectors[1];
      gmds::math::Point p1 = center + (vx + vy )*cube_size;
      gmds::Node n1 = mesh_cube.newNode(p1);
      gmds::math::Point p2 = center + (vx - vy)*cube_size;
      gmds::Node n2 = mesh_cube.newNode(p2);
      gmds::math::Point p3 = center + (vx + vy).opp()*cube_size;
      gmds::Node n3 = mesh_cube.newNode(p3);
      gmds::math::Point p4 = center + (vy - vx)*cube_size;
      gmds::Node n4 = mesh_cube.newNode(p4);

      mesh_cube.newQuad(n1, n2, n3, n4);
    }
    gmds::VTKWriter<gmds::IGMesh> writer_cube(mesh_cube);
    std::string filename = "smoothing."+std::to_string((long long int)(i));

    writer_cube.write(filename, DIM3 | F | N);
    
  }
  /*----------------------------------------------------------------------------*/
  /*
   * IN N -> nb variables  (nb vertices in our case)
   * IN x -> unknowns (size=N)
   * IN prev_x -> NOT USED previous value of x
   * OUT func-> pointer on the computed function
   * OUT grad -> grad(func)
   */
  void evalF(int N, 
	     double* x, 
	     double* prev_x, 
	     double* func, 
	     double* grad)
  {
  //  size_t nvars = N;
//
    //Size of "grad" is already equal to N, since it points on a vector
    //initialized with the right size in HLBFGS

    *func = 0.0;

    for (unsigned int i = 0; i < N; i++)
      grad[i] = 0.0;

    // Nodes are traversed to build the system
    for (unsigned int i_edge = 0; i_edge < HLBFGS_edges.size(); i_edge++) {
      Edge e_i =  HLBFGS_edges[i_edge];
      std::vector<Node> edge_nodes = e_i.get<Node>();
      Node ni = edge_nodes[0];
      Node nj = edge_nodes[1];
    
      int indexi = HLBFGS_candidates_index[ni.getID()];
      int indexj = HLBFGS_candidates_index[nj.getID()];
    

      double alphai=0;
      bool boundary_i = false;
      bool boundary_j = false;
      if( HLBFGS_mesh->isMarked(ni,HLBFGS_mark_smooth))
	alphai = x[indexi];
      else {//marked steady
	alphai = HLBFGS_boundary_x[ni.getID()];  
	boundary_i=true;
      }

      double alphaj=0;
      if( HLBFGS_mesh->isMarked(nj,HLBFGS_mark_smooth))
	alphaj = x[indexj];
      else {//marked steady
	alphaj = HLBFGS_boundary_x[nj.getID()];
	boundary_j=true;
      }

      std::vector<Face> adj_faces = e_i.get<Face>();
      double w=0;
      for(unsigned int i_face=0; i_face<adj_faces.size();i_face++){
	Face f_i = adj_faces[i_face];
	w +=0.33333333*f_i.area();
      }
      double cos_ij = cos( alphai - alphaj );
      double sin_ij = sin( alphai - alphaj );

      double cos_ji = cos( alphaj - alphai );
      double sin_ji = sin( alphaj - alphai );

      *func += -0.5*w*(cos_ij+cos_ji);
      if(!boundary_i)
	grad[indexi] += 2.0*w*sin_ij;
      if(!boundary_j)
	grad[indexj] += 2.0*w*sin_ji;
    }

    // for (unsigned int a = 0; a < N; a++) {
    //   Node    n_a =  HLBFGS_candidates[a];
    //   double alpha_a = x[a];
    //   TCellID n_a_id =  n_a.getID();
   
    //   std::vector<gmds::TCellID> ring_a = HLBFGS_vicinity_index [n_a_id];
    //   std::vector<double> weight_a      = HLBFGS_vicinity_weight[n_a_id];

    //   //   std::cout<<"Ref node "<<n_a_id
    //   //	     <<", alpha: "<<alpha_a<<std::endl; 
    //   //================================================================
    //   // 1) We compute the local contribution to f and grad(f)
    //   //================================================================  
    //   double func_i = 0.0;
    //   double grad_i = 0.0;
    //   for(unsigned int i=0;i<ring_a.size();i++) {
    //     Node n_i = HLBFGS_mesh->get<Node>(ring_a[i]);
    //     double alpha_i=0.0;
    //     /*     std::cout<<"Node "<<n_i.getID();
    //     if(HLBFGS_mesh->isMarked(n_i, HLBFGS_mark_steady)) {
    // 	alpha_i =  HLBFGS_boundary_x[ring_a[i]];
    // 	std::cout<<" boundary";
    //     }
    //     else {
    // 	alpha_i = x[ring_a[i]];
    // 	std::cout<<" inner";
    //     }
    //     std::cout<<", with alpha "<<alpha_i;*/
    //     double w_i     = weight_a[i];
    //     //      std::cout<<", and w: "<<w_i<<std::endl;

    //     //      std::cout<<"cos(a-i), sin (a-i): "<<cos(alpha_a-alpha_i)<<", "
    //     //	       <<sin(alpha_a-alpha_i)<<std::endl;
    //     double cos_ai = -cos(alpha_a-alpha_i);
    //     double sin_ai = -sin(alpha_a-alpha_i);

    //     func_i += cos_ai;
    //     grad_i += 2.0*sin_ai;
    //     //  func_i += 0.5*w_i*(1-cos(alpha_a-alpha_i));
    //     //      grad_i += 2.0*w_i*sin(alpha_a-alpha_i);
    //   }//for(unsigned int k=0;k<cyclic_faces.size();k++) {


    //   //================================================================
    //   // 2) We contribute to global f and grad(f)
    //   //================================================================
    //   *func += func_i;
    //   grad[a] += grad_i;

    // }//for (unsigned int i = 0; i < HLBFGS_candidates.size(); i++)

  }
  /*----------------------------------------------------------------------------*/
  void newiteration(int iter, int call_iter, double *x, double* f, double *g,
		    double* gnorm)
  {
    std::cout << iter << ": " 
	      << call_iter << " " 
	      << *f << " " 
	      << *gnorm << std::endl;

    long long int i = iter;
    //   writeMesh(iter,x);
  }
} //namespace hlbfgs_for_cross2D 
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
CrossFieldSmoother2D::
CrossFieldSmoother2D(gmds::IGMesh* AMesh,
		     gmds::Variable<gmds::math::Cross2D>*& AField)
  : m_mesh(AMesh), m_field(AField), m_marks_initialized(false),
    m_mark_smooth(0), m_mark_steady(0)
{ }
/*---------------------------------------------------------------------------*/
void CrossFieldSmoother2D::initMarks(const int AToSmooth, const int ASteady)
{
  std::cout<<"marques IN :"<<AToSmooth
	   <<", "<<ASteady <<std::endl;
  m_mark_smooth = AToSmooth;
  m_mark_steady = ASteady;
  hlbfgs_for_cross2D::HLBFGS_mark_smooth = AToSmooth;
  hlbfgs_for_cross2D::HLBFGS_mark_steady = ASteady;
  std::cout<<"marques IN :"<<m_mark_smooth
	   <<", "<<m_mark_steady <<std::endl;
  std::cout<<"marques :"<<hlbfgs_for_cross2D::HLBFGS_mark_smooth
	   <<", "<<hlbfgs_for_cross2D::HLBFGS_mark_steady <<std::endl;
  m_marks_initialized = true;
}
/*---------------------------------------------------------------------------*/
void CrossFieldSmoother2D::execute( const int AMaxIteration,
				    const double AMaxDelta)
{
  if (!m_marks_initialized)
    throw GMDSException("Boolean marks were not initialized");

  if (AMaxIteration<1)
    throw GMDSException("The maximum number of iterations must be strictly positive");

  if (AMaxDelta<=0)
    throw GMDSException("The maximum delta must be strictly positive");

  hlbfgs_for_cross2D::HLBFGS_mesh = m_mesh;
  //======================================================================
  // 1) We get the total area of the mesh we work on   
  //======================================================================
  hlbfgs_for_cross2D::HLBFGS_global_area = computeTotalArea();

  //======================================================================
  // 2) We initialize the set of nodes to work on with HLBFGS
  //======================================================================
  std::cout<<" Initialization of data "<<std::endl;
  initOptimizationData();
  std::cout<<"\t DONE"<<std::endl;
  //======================================================================
  // 3) We initialize the constraint system
  //======================================================================
  int reference_angles_size =  hlbfgs_for_cross2D::HLBFGS_candidates.size();
  double* ref_angles = new double[reference_angles_size];

  // Euler angles are computed for every defined vertex (including ridges)
  for (unsigned int i = 0; i <  hlbfgs_for_cross2D::HLBFGS_candidates.size(); i++)
    {
      Node n =  hlbfgs_for_cross2D::HLBFGS_candidates[i];
      TCellID n_id = n.getID();
      //we store the candidate index of n_i
      hlbfgs_for_cross2D::HLBFGS_candidates_index[n_id] = i;

      math::Cross2D c2d = (*m_field)[n_id];
      ref_angles[i] = c2d.referenceAngle();
    }

  //Definition of the support Edges
  hlbfgs_for_cross2D::HLBFGS_edges.clear();
  for (unsigned int i = 0; i < hlbfgs_for_cross2D::HLBFGS_candidates.size(); i++){
    Node n = hlbfgs_for_cross2D::HLBFGS_candidates[i];

    std::vector<Edge> adj_edges = n.get<Edge>();

    for (unsigned int j = 0; j < adj_edges.size(); j++) {
	Edge e_j = adj_edges[j];

	std::vector<Node> e_j_nodes = e_j.get<Node>();
	Node other_node =
	  (e_j_nodes[0].getID() == n.getID()) ? e_j_nodes[1] : e_j_nodes[0];

	if (m_mesh->isMarked(other_node, m_mark_smooth) || 
	    m_mesh->isMarked(other_node, m_mark_steady) ) {
	  if (other_node.getID() < n.getID()) //to store the edge only once
	    hlbfgs_for_cross2D::HLBFGS_edges.push_back(e_j);
	}
      }
  }
  //======================================================================
  // 4) System solving
  //======================================================================
  //Initialization of HLBFGS
  double parameter[20];
  int info[20];
  INIT_HLBFGS(parameter, info);
  int N =  hlbfgs_for_cross2D::HLBFGS_candidates.size(); // Nb variable, one angle per node
  int M = 7;//Default choice for the optimization algorithm
  int T = 0;
  int num_iter = 10000;
  info[3] = 1; //0 is default value, but 1 is recommended is the HLBGS documentation
  info[4] = num_iter;
  info[5] = 1; //print messages
  info[6] = T;
  info[7] = 0;
  info[10] = 0;
  info[11] = 1;
  /*  parameter[0] = 0.1;
  parameter[1] = 0.1;
  parameter[3] = 0.1;
  parameter[4] = 0.1;
  */

  std::cout << "=== Before HLBFGS Smoothing" << std::endl;

  /*  int nannb = 0;
  for (unsigned int i = 0; i<N; i++){
    double d = ref_angles[i];
    if (!(d>-5) && !(d < 5)){
      nannb++;
    }

  }
  std::cout << "Nb NAN in reference angles: " << nannb << std::endl;
  */
  HLBFGS(N, M, &ref_angles[0],  
	 hlbfgs_for_cross2D::evalF, 0, 
	 HLBFGS_UPDATE_Hessian,
	 hlbfgs_for_cross2D::newiteration, 
	 parameter, info);
  std::cout << "=== After smoothing, nb iter : " << info[2]
	    << " for " << info[1] << " evaluations" <<std::endl;

  std::cout << "2D crosses rebuild" << std::endl;
  rebuild2DCrosses(ref_angles);
  std::cout << "2D crosses rebuild done" << std::endl;


  delete[] ref_angles;
}

/*----------------------------------------------------------------------------*/
void CrossFieldSmoother2D::rebuild2DCrosses(double*& ARefAngles)
{

  for (unsigned int i = 0; i < hlbfgs_for_cross2D::HLBFGS_candidates.size(); i++) {
    Node n = hlbfgs_for_cross2D::HLBFGS_candidates[i];
    TCoord a = ARefAngles[i];
    (*m_field)[n.getID()] =math::Cross2D(a);
  }
}

/*---------------------------------------------------------------------------*/
void CrossFieldSmoother2D::initOptimizationData()
{ 
  //HLBFSG global data structure are cleared first
  hlbfgs_for_cross2D::HLBFGS_candidates.clear();
  hlbfgs_for_cross2D::HLBFGS_vicinity_index.clear();
  hlbfgs_for_cross2D::HLBFGS_vicinity_weight.clear();  

  //======================================================================
  // 1) We gather the node we will work on
  //======================================================================
  IGMesh::node_iterator it_nodes = m_mesh->nodes_begin();
  for (; !it_nodes.isDone(); it_nodes.next()) {
    Node n = it_nodes.value();
    if (m_mesh->isMarked(n, m_mark_smooth))
      hlbfgs_for_cross2D::HLBFGS_candidates.push_back(it_nodes.value());
    else if (m_mesh->isMarked(n, m_mark_steady)){
      hlbfgs_for_cross2D::HLBFGS_boundary_x[n.getID()] 
	= (*m_field)[n.getID()].referenceAngle();
    }
  }
  
  std::cout << std::endl << "Nb candidates (" 
	    << hlbfgs_for_cross2D::HLBFGS_candidates.size()
	    << " / " << m_mesh->getNbNodes() << ")" << std::endl;

  //======================================================================
  // 2) some data useful to solve the HLBFGS system are stored
  //======================================================================
  for (unsigned int a = 0; a<hlbfgs_for_cross2D::HLBFGS_candidates.size(); a++) {

    Node    n_a =   hlbfgs_for_cross2D::HLBFGS_candidates[a];
    std::vector<TCellID> ring_a;
    std::vector<double> weight_a;
   
    // We build the oriented ring of nodes around n_a
    std::vector<Face> adj_faces = n_a.get<Face>();
    Face first_face   = adj_faces[0];
    std::vector<Node> current_nodes = first_face.get<Node>();

    Node n_i;
    if(current_nodes[0].getID()==n_a.getID())
      n_i = current_nodes[1];
    else
      n_i = current_nodes[0];

    Face current_face = adj_faces[0];
    Face next_face;
    do{
      //get the next face
      bool found_next_face=false;
      for(unsigned int i=0; !found_next_face && i<adj_faces.size();i++){
	Face f_i = adj_faces[i];
	//we do not consider the current face, we need the next one
	if(f_i.getID()==current_face.getID())
	  continue;

	std::vector<TCellID> node_ids = f_i.getIDs<Node>();
	bool found_na = false;
	bool found_ni = false;
	for(unsigned int j=0; j<node_ids.size(); j++){
	  if(node_ids[j]== n_a.getID())
	    found_na=true;
	  else if (node_ids[j] == n_i.getID())
	    found_ni=true;
	}
	if(found_ni && found_na){
	  found_next_face = true;
	  next_face = f_i;
	}
      }// for(unsigned int i=0;!found_next_face && i<adj_faces.size();i++){
      //We've got the next face      
      double f_cur_area = current_face.area();
      double f_nex_area = next_face.area();
      double w_i = (f_cur_area+f_nex_area)/3;//(3*hlbfgs_for_cross2D::HLBFGS_global_area);

      //data relative to n_a, n_i are stored
      ring_a.push_back(n_i.getID());
      weight_a.push_back(w_i);
      
      //we compute the next n_i and current_face;
      current_face = next_face;
      Node prev_n_i = n_i;
      std::vector<Node> current_nodes = current_face.get<Node>();
      for(unsigned int i=0; i<current_nodes.size();i++){
	TCellID current_id = current_nodes[i].getID();
	if(current_id!=prev_n_i.getID() && current_id != n_a.getID())
	  n_i = current_nodes[i];
      }
    }
    while(next_face.getID()!=first_face.getID());

    //now we update the HLBFGS maps
    hlbfgs_for_cross2D::HLBFGS_vicinity_index[n_a.getID()] = ring_a;
    hlbfgs_for_cross2D::HLBFGS_vicinity_weight[n_a.getID()] = weight_a;
  }

}
/*---------------------------------------------------------------------------*/
double CrossFieldSmoother2D::computeTotalArea()
{
  double area = 0.0;
  IGMesh::face_iterator it_faces = m_mesh->faces_begin();
  for(;!it_faces.isDone(); it_faces.next()) {
    Face f  = it_faces.value();
    std::vector<Node> f_nodes = f.get<Node>();
    bool valid = true;
    for(unsigned int i=0; valid && i<f_nodes.size();i++){
      if(!m_mesh->isMarked(f_nodes[i],m_mark_smooth) && 
	 !m_mesh->isMarked(f_nodes[i],m_mark_steady) )
	valid = false;
    }
    if(valid)
      area += f.area();
  }
  return area;
}
/*---------------------------------------------------------------------------*/
