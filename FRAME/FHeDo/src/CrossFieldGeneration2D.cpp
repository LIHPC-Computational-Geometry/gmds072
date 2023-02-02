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
#include <iostream>
/*---------------------------------------------------------------------------*/
#include <GMDS/Algo/BoundaryOperator.h>
#include <GMDS/Algo/DistanceFieldBuilder2D.h>
/*---------------------------------------------------------------------------*/
#include <GMDS/IO/VTKWriter.h>
#include <sstream>
/*---------------------------------------------------------------------------*/
#include "CrossFieldGeneration2D.h"
#include "LevelSetCross2D.h"
#include "StabilityBallCross2D.h"
#include "LaplaceCross2D.h"
#include "FunctionalApproachCross2D.h"
/*---------------------------------------------------------------------------*/
using namespace gmds;
/*---------------------------------------------------------------------------*/
CrossFieldGeneration2D::CrossFieldGeneration2D(IGMesh* AMesh)
  : m_mesh(AMesh)
{
  try{
    m_cross_field_2D = m_mesh->getVariable<math::Cross2D>(GMDS_NODE, "cross");
  }
  catch (GMDSException& e){
    std::cout << e.what() << std::endl;
    std::cout << "\t -> generation of the cross variable on nodes"<<std::endl;
    m_cross_field_2D = m_mesh->newVariable<math::Cross2D>(GMDS_NODE, "cross");
  }

  m_debug_output = ".";
}
/*---------------------------------------------------------------------------*/
void CrossFieldGeneration2D::setDebugPrefix(const std::string& AName)
{
  m_debug_output = AName;
}
/*---------------------------------------------------------------------------*/
void CrossFieldGeneration2D::execute(const Strategy AStrategy)
{
  std::cout << "======================================"<< std::endl;
  std::cout<<" Starting 2D cross field generation " <<std::endl;
  std::cout << "======================================"<< std::endl;

  std::cout << "Boolean Marks' initialization " << std::endl; 
  initMarks(); 
  std::cout << "  -->  DONE" << std::endl;

  //==================================================================
  //STEP 1 - We mark all the nodes, edges, faces classified on
  // boundary entities (points, curves, surfaces)
  //==================================================================
  std::cout << "======================================"<< std::endl;
  std::cout << "Mark boundary cells" << std::endl;
  markBoundaryCells();
  std::cout << "  -->  DONE" << std::endl;

  //==================================================================
  //STEP 2 - We store all the nodes classified on curves and surfaces
  // in STL vectors
  //=================================================================
  std::cout << "======================================"<< std::endl;
  std::cout << " Storage of nodes classified on curves " << std::endl;
  IGMesh::node_iterator it_nodes = m_mesh->nodes_begin();
  m_curve_nodes.clear();
  m_surf_nodes.clear();
  int c=0;
  for (; !it_nodes.isDone(); it_nodes.next()) {
    Node n = it_nodes.value();
    c++;
    if (m_mesh->isMarked(n, m_markNodeOnCurv)) 
      m_curve_nodes.push_back(n);
    else if (!m_mesh->isMarked(n, m_markIsolated) &&
	     !m_mesh->isMarked(n, m_markNodeOnPnt) )
      m_surf_nodes.push_back(n);
  }
  std::cout<<"nodes on curves  : "<<m_curve_nodes.size()<<std::endl;
  std::cout<<"nodes on surfaces: "<<m_surf_nodes.size()<<std::endl;
  std::cout<<"total nb. nodes  : "<<m_mesh->getNbNodes()<<" / "<<c<<std::endl;

  std::cout << "    DONE" << std::endl;

  //==================================================================
  //STEP 2 - For nodes on curves, we compute crosses from the 
  // geometric information we have
  //==================================================================
  std::cout << "======================================"<< std::endl;

  std::cout << "Initialization of crosses along curves" << std::flush;
  initCrossesOnCurves(); 

  std::cout << "    DONE" << std::endl;
  
  
  //==================================================================
  //STEP 3 - A cross is associated to each node classified on 
  // a geometric point
  //==================================================================
  std::cout << "======================================"<< std::endl;

  std::cout << "Initialization of crosses on points" << std::flush;
  initCrossesOnPoints(); 
  std::cout << "    DONE" << std::endl;
  writeForDebug();  
  //==================================================================
  // Now, depending on the strategy, the algorithm changes
  //==================================================================
  if(AStrategy==laplace_solve)
    buildFieldViaLaplaceEDP();
  else if(AStrategy== level_sets)
      buildFieldViaLevelSets();
  else if(AStrategy== functional_approach)
          buildFieldViaFunctionalApproach();
  else
    throw GMDSException("Wrong strategy value for 2D cross field generation");
 //==================================================================
  // A final smoothing step can be performed if required
  //==================================================================

  // if(AStrategy!= laplace_solve) {
  //   std::cout << "======================================"<< std::endl;
  //   std::cout << " smoothing" << std::flush;
  //   smoothAll();
  //   std::cout << "    DONE" << std::endl;
  // }
    
  // //==================================================================
  // // DEBUGGING STEP - Tet coloring
  // //==================================================================
  // std::cout << "======================================"<< std::endl;
  // std::cout << " Coloring simplices" << std::endl;
  // colorSimplices();
  // std::cout << "    DONE" << std::endl;

  //==================================================================
  // CLEANING - Boolean marks are cleaned
  //==================================================================
  cleanMarks();
  writeForDebug();  
}

/*---------------------------------------------------------------------------*/
void CrossFieldGeneration2D::buildFieldViaLaplaceEDP()
{
    std::vector<Face> mesh_faces;
    mesh_faces.resize(m_mesh->getNbFaces());
    IGMesh::face_iterator it_f = m_mesh->faces_begin();
    int f_index=0;
    for(;!it_f.isDone();it_f.next()) {
        mesh_faces[f_index++]= it_f.value();
    }
    LaplaceCross2D algo(m_mesh, m_cross_field_2D,
                        m_curve_nodes, m_surf_nodes,
                        mesh_faces);
    algo.execute();
}
/*---------------------------------------------------------------------------*/
void CrossFieldGeneration2D::buildFieldViaFunctionalApproach()
{
    std::vector<Face> mesh_faces;
    mesh_faces.resize(m_mesh->getNbFaces());
    IGMesh::face_iterator it_f = m_mesh->faces_begin();
    int f_index=0;
    for(;!it_f.isDone();it_f.next()) {
        mesh_faces[f_index++]= it_f.value();
    }
    FunctionalApproachCross2D algo(this, m_mesh, m_cross_field_2D,
                                   m_curve_nodes, m_surf_nodes,
                                   mesh_faces);
    algo.execute();
}
/*---------------------------------------------------------------------------*/
void CrossFieldGeneration2D::buildFieldViaLevelSets()
{
  // only nodes on curves are given at the initialization
  // nodes on points have crosses without any meaning
  LevelSetCross2D algo(m_mesh, m_cross_field_2D,
		       m_curve_nodes, m_surf_nodes,
		       m_markFace);
  algo.execute();
}
/*---------------------------------------------------------------------------*/
void CrossFieldGeneration2D::cleanMarks()
{
  m_mesh->unmarkAll<Node>(m_markNodeOnCurv);
  m_mesh->unmarkAll<Node>(m_markNodeOnPnt);
  m_mesh->unmarkAll<Node>(m_markIsolated);
  m_mesh->unmarkAll<Edge>(m_markEdgeOnCurv);
  m_mesh->unmarkAll<Face>(m_markFace);  

  m_mesh->freeMark<Node>(m_markNodeOnCurv);
  m_mesh->freeMark<Node>(m_markNodeOnPnt);
  m_mesh->freeMark<Node>(m_markIsolated);
  m_mesh->freeMark<Edge>(m_markEdgeOnCurv);
  m_mesh->freeMark<Face>(m_markFace); 
}
/*---------------------------------------------------------------------------*/
void CrossFieldGeneration2D::initMarks()
{
  m_markNodeOnCurv = m_mesh->getNewMark<Node>();
  m_markNodeOnPnt  = m_mesh->getNewMark<Node>();
  m_markIsolated   = m_mesh->getNewMark<Node>();
  m_markEdgeOnCurv = m_mesh->getNewMark<Edge>();
  m_markFace       = m_mesh->getNewMark<Face>();

}
/*---------------------------------------------------------------------------*/
void CrossFieldGeneration2D::markBoundaryCells()
{
  BoundaryOperator boundaryOp(m_mesh);
  if (!boundaryOp.isValid())
    {
      std::cout << "Invalid model for boundary operations" << std::endl;
      throw GMDSException("Invalid model for boundary operations");
    }

  int mark_edge_on_surf = m_mesh->getNewMark<Edge>();
  int mark_node_on_surf = m_mesh->getNewMark<Node>();

  boundaryOp.markCellOnGeometry(m_markFace,
				mark_edge_on_surf,
				mark_node_on_surf,
				m_markEdgeOnCurv, 
				m_markNodeOnCurv, 
				m_markNodeOnPnt,
				m_markIsolated);


  m_mesh->unmarkAll<Node>(mark_node_on_surf);
  m_mesh->unmarkAll<Edge>(mark_edge_on_surf);

  m_mesh->freeMark<Node>(mark_node_on_surf);
  m_mesh->freeMark<Edge>(mark_edge_on_surf);
}
/*---------------------------------------------------------------------------*/
std::vector<Edge> CrossFieldGeneration2D::
getEdgesOnCurve(const Node& ANode) const
{
  std::vector<Edge> edges_on_curve;
  std::vector<Edge> adj_edges = ANode.get<Edge>();
  for (unsigned int i = 0; i < adj_edges.size(); i++)
    {
      Edge ei = adj_edges[i];
      if (m_mesh->isMarked(ei, m_markEdgeOnCurv))
	edges_on_curve.push_back(ei);
    }
  return edges_on_curve;
}
/*---------------------------------------------------------------------------*/
Node CrossFieldGeneration2D::
getNeighboorOn(const Node& ANode, const Edge& AEdge) const
{
  std::vector<Node> nodes = AEdge.get<Node>();
  if (nodes[0].getID() == ANode.getID())
    return nodes[1];

  return nodes[0];
}
/*---------------------------------------------------------------------------*/
void  CrossFieldGeneration2D::initCrossesOnCurves()
{
  //for each node on a geometrical curve, we compute
  //its associated 2D cross
  IGMesh::node_iterator it_node = m_mesh->nodes_begin();
  for (; !it_node.isDone(); it_node.next())
    {
      Node current_node = it_node.value();
      
      if (m_mesh->isMarked(current_node, m_markNodeOnCurv) &&
	  !m_mesh->isMarked(current_node, m_markNodeOnPnt)){

	//current_node is on a geometric curve
	//We get adjacent edges that are classified onto
	//a geometric curve. We must get one or two edges
	//at most.
	std::vector<Edge> ridges = getEdgesOnCurve(current_node);
	math::Vector newN;
	if (ridges.size() == 1) {
	  Edge current_edge = ridges[0];
	  Node node1 = getNeighboorOn(current_node, current_edge);
	  Node node2 = current_node;

	  //we build the direction vector of the current edge
	  math::Point p1 = node1.getPoint();
	  math::Point p2 = node2.getPoint();
	  math::Vector v1 = math::Vector(p1, p2);
	  v1.normalize();

	  newN = v1;
	}
	else if (ridges.size() == 2) {
	  //With 2 adajcent edges on the curve, we compute average values

	  Edge edge1 = ridges[0];
	  Edge edge2 = ridges[1];

	  Node node1 = getNeighboorOn(current_node, edge1);
	  Node node2 = getNeighboorOn(current_node, edge2);

	  //we build the average between adjacent edges
	  math::Point  p1 = node1.getPoint();
	  math::Point  p  = current_node.getPoint();
	  math::Point  p2 = node2.getPoint();
	  math::Vector v1 = math::Vector(p1, p);
	  math::Vector v2 = math::Vector(p, p2);
	  v1.normalize();
	  v2.normalize();

	  newN = v1+v2;
	  newN.normalize();
	}
	else{
	  std::cout << "Nb ridges for an edge adjacent to node "
		    << current_node.getID() << ": " << ridges.size() << std::endl;
	  throw GMDSException("A ridge node has an illegal number of edges.");
	}
	std::vector<Face> adj_faces;
	current_node.get<Face>(adj_faces);
	Face adj_f = adj_faces[0];
	math::Vector normal  = adj_f.normal();
	
	math::Vector v1 = newN;

	//Computation of v3 from v1 and v2.
	math::Vector v2 = v1.cross(normal);
	math::Cross2D c(v1,v2);

	(*m_cross_field_2D)[current_node.getID()] = c;
      }//if (m_mesh->isMarked(current_node, m_markNodeOnCurv))

    }//for (; !it_node.isDone(); it_node.next())
}
/*----------------------------------------------------------------------------*/
void CrossFieldGeneration2D::initCrossesOnPoints()
{
  IGMesh::node_iterator it_node = m_mesh->nodes_begin();
  for (; !it_node.isDone(); it_node.next())
    {
      Node current_node = it_node.value();

      if (m_mesh->isMarked(current_node, m_markNodeOnPnt))
	{
	  //we initialize the value at points, but we do not use it after.
	  (*m_cross_field_2D)[current_node.getID()] = math::Cross2D();
	}

    }
}
/*----------------------------------------------------------------------------*/
void CrossFieldGeneration2D::smoothAll()
{
  int mark_smooth = m_mesh->getNewMark<Node>();
  IGMesh::node_iterator it_nodes = m_mesh->nodes_begin();
  for (; !it_nodes.isDone(); it_nodes.next())
    {
      Node n = it_nodes.value();
      if (!m_mesh->isMarked(n, m_markNodeOnCurv) &&
	  !m_mesh->isMarked(n, m_markNodeOnPnt))
	m_mesh->mark(n, mark_smooth);

    }
  smooth(mark_smooth);

  m_mesh->unmarkAll<Node>(mark_smooth);
  m_mesh->freeMark<Node>(mark_smooth);



} 
/*----------------------------------------------------------------------------*/
void CrossFieldGeneration2D::smooth(const int AMark)
{

  //	SmoothingHLBFGS smoother(m_mesh, m_cross_field_2D, m_surf_normal);
  /*	FrameFieldLaplacianSmoothing smoother(m_mesh, m_cross_field_2D, m_surf_normal);
	smoother.initBoundaryMarks(m_markNodeOnPnt, m_markNodeOnCurv, m_markNodeOnSurf);
	smoother.selectNodes(AMark);
	smoother.execute();
  */
}
/*----------------------------------------------------------------------------*/
void CrossFieldGeneration2D::colorSimplices()
{
  // Variable<int>* var_sing = 0;
  // try{
  //   var_sing = m_mesh->getVariable<int>(GMDS_REGION, "sing_tet");
  // }
  // catch (GMDSException& e){
  //   var_sing = m_mesh->newVariable<int>(GMDS_REGION, "sing_tet");
  // }
		    
  // IGMesh::region_iterator it = m_mesh->regions_begin();
		      
  // int nbOfCluster = 0;
  // int nbColoredTet = 0;
  // //=========================================================================
  // // INTERN SKELETON CREATION
  // //=========================================================================
  // // FIRST LOOP ON REGIONS TO GET ALL THE 3-SING. TETS
  // //=========================================================================
  // for (; !it.isDone(); it.next()){
  //   Region current_region = it.value();
  //   std::vector<Node> nodes = current_region.get<Node>();
  //   bool onPnt = false;
  //   for (unsigned int i_node = 0; i_node < nodes.size(); i_node++)
  //     {
  // 	Node ni = nodes[i_node];
  // 	if (m_mesh->isMarked(ni, m_markNodeOnPnt))
  // 	  onPnt = true;

  //     }
  //   if (onPnt)
  //     {
  // 	(*var_sing)[current_region.getID()] = 0;
  //     }
  //   else
  //     {
  // 	std::vector<TCellID> nodeIDs = current_region.getIDs<Node>();
  // 	int ID1 = nodeIDs[0];
  // 	int ID2 = nodeIDs[1];
  // 	int ID3 = nodeIDs[2];
  // 	int ID4 = nodeIDs[3];
						
  // 	int singTypeTmp = math::Quaternion::testSingularity((*m_cross_field)[ID1],
  // 							    (*m_cross_field)[ID2],
  // 							    (*m_cross_field)[ID3],
  // 							    (*m_cross_field)[ID4]);
  // 	if (singTypeTmp != 0)
  // 	  nbColoredTet++;
						    
  // 	(*var_sing)[current_region.getID()] = singTypeTmp;
  //     }
  // }
  // std::cout << "Nb. colored tetrahedra: " << nbColoredTet << std::endl;
}

/*---------------------------------------------------------------------------*/
void CrossFieldGeneration2D::writeForDebug(const std::string AFileName)
{
  static int nb_file = 0;
  VTKWriter<IGMesh> writer(*m_mesh);
  if (AFileName != "")
    { 
      std::stringstream file_name;
      file_name <<m_debug_output<<"_"<<AFileName;
      writer.write(file_name.str(), DIM3 | F | N);
    }
  else
    {
      std::stringstream file_name;
      file_name <<m_debug_output<<"_FFG_2D_Debug_" << nb_file;
      writer.write(file_name.str(), DIM3 | F | N);
      nb_file++;
    }

  IGMesh::node_iterator it = m_mesh->nodes_begin();
  double x_min = 100000;
  double y_min = 100000;
  double x_max = -100000;
  double y_max = -100000;
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


    }
  double dist_x = x_max - x_min;
  double dist_y = y_max - y_min;

  double cube_size = 0;
  if (dist_x <= dist_y ){
    cube_size = dist_x;
  }
  else {
    cube_size = dist_y;
  }

  cube_size /= 20;

  MeshModel model_cube(DIM3 | F | N | F2N);
  IGMesh mesh_cube(model_cube);
  for (it = m_mesh->nodes_begin(); !it.isDone(); it.next())
    {
      Node n = it.value();
      math::Point center = n.getPoint();
      //      if (m_mesh->isMarked(n, m_mark_alive))
{
	math::Cross2D current_cross = (*m_cross_field_2D)[n.getID()];

	std::vector<math::Vector> current_vectors = current_cross.componentVectors();
	math::Vector vx = current_vectors[0];
	math::Vector vy = current_vectors[1];
	math::Point p1 = center + (vx + vy )*cube_size;
	Node n1 = mesh_cube.newNode(p1);
	math::Point p2 = center + (vx - vy)*cube_size;
	Node n2 = mesh_cube.newNode(p2);
	math::Point p3 = center + (vx + vy).opp()*cube_size;
	Node n3 = mesh_cube.newNode(p3);
	math::Point p4 = center + (vy - vx)*cube_size;
	Node n4 = mesh_cube.newNode(p4);

	mesh_cube.newQuad(n1, n2, n3, n4);
      }
    }
  VTKWriter<IGMesh> writer_cube(mesh_cube);

  std::stringstream file_name_cube;
  file_name_cube<<m_debug_output <<"FFG_2D_Debug_Cube_" << nb_file;
  writer_cube.write(file_name_cube.str(), DIM3 | F | N);

}
