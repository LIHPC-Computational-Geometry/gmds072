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
 * SmoothingField2DTest.h
 *
 *  Created on: 01 march 2015
 *      Author: ledouxf
 */

/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IGMesh.h>
#include <GMDS/IG/IGMeshDoctor.h>
#include <GMDS/IO/MeditReader.h>
#include <GMDS/IO/VTKReader.h>
#include <GMDS/IO/VTKWriter.h>
#include <GMDS/Algo/BoundaryOperator.h>
/*----------------------------------------------------------------------------*/
#include "CrossFieldGeneration2D.h"
#include "CrossFieldSmoother2D.h"
/*----------------------------------------------------------------------------*/
#include <iostream> 
#include <vector>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
class SmoothingField2DTest: public ::testing::Test 
{
 protected:
  SmoothingField2DTest(){;}
  virtual ~SmoothingField2DTest(){;}
};
/*----------------------------------------------------------------------------*/
TEST_F(SmoothingField2DTest, DISABLED_square_holes_level_sets)
{
  MeshModel model(DIM3 | F | E | N | F2N | E2N | F2E | E2F | N2F | N2E);
  IGMesh mesh(model);
  
  //==================================================================
  // MESH READING
  //==================================================================
  MeditReader<IGMesh> reader(mesh);
  reader.read("Samples/square_holes.mesh", DIM3|F|N);
  //==================================================================
  // MESH TOPOLOGY PREPARATION
  //==================================================================3
  IGMeshDoctor doctor(&mesh);
  doctor.buildEdgesAndX2E();
  doctor.updateUpwardConnectivity();
  doctor.orient2DFaces();

  EXPECT_TRUE(mesh.getNbRegions()==0);  
  EXPECT_TRUE(mesh.getNbFaces()!=0);  
  EXPECT_TRUE(mesh.getNbEdges()!=0);
  EXPECT_TRUE(mesh.getNbNodes()!=0);
 //==================================================================
  // FRAME FIELD GENERATION
  //=================================================================
  CrossFieldGeneration2D algo(&mesh);
  algo.setDebugPrefix("SquareHoleLevelSetSmooth");
  algo.execute(CrossFieldGeneration2D::level_sets);

 //==================================================================
  // FRAME FIELD SMOOTHING
  //================================================================= 
  std::cout<<"===== SMOOTHING ====="<<std::endl;
  gmds::BoundaryOperator boundaryOp(&mesh); 

  int edge_curve_mark = mesh.getNewMark<gmds::Edge>();
  int node_curve_mark = mesh.getNewMark<gmds::Node>();
  int node_point_mark = mesh.getNewMark<gmds::Node>();
  int i_mark = mesh.getNewMark<gmds::Node>();
  int a_mark = mesh.getNewMark<gmds::Node>();
  int b_mark = mesh.getNewMark<gmds::Node>();
  std::vector<TCellID> b_ids;
  boundaryOp.markCellsOnCurves(edge_curve_mark, node_curve_mark);
  boundaryOp.markNodesOnPoint(edge_curve_mark, node_curve_mark, node_point_mark);
  boundaryOp.markAloneNodes(a_mark);
  
  gmds::IGMesh::node_iterator it_node = mesh.nodes_begin();
  for(;!it_node.isDone();it_node.next()){
    Node n = it_node.value();
    if(mesh.isMarked(n,node_curve_mark) && !mesh.isMarked(n,node_point_mark))
      mesh.mark(n, b_mark);
    else if (!mesh.isMarked(n,node_curve_mark) && !mesh.isMarked(n,node_point_mark))
      mesh.mark(n, i_mark);
  }
  
  gmds::Variable<gmds::math::Cross2D>* field  = 
    mesh.getVariable<math::Cross2D>(GMDS_NODE, "cross");
  // std::cout.precision(16);
  //  std::cout<<std::scientific;
  CrossFieldSmoother2D smooth(&mesh, field);
  smooth.initMarks(i_mark, b_mark);
  smooth.execute();

  mesh.unmarkAll<Node>(a_mark);
  mesh.unmarkAll<Node>(b_mark);
  mesh.unmarkAll<Node>(i_mark);
  mesh.unmarkAll<Node>(node_curve_mark);
  mesh.unmarkAll<Node>(node_point_mark);
  mesh.unmarkAll<Edge>(edge_curve_mark);

  mesh.freeMark <Node>(a_mark);
  mesh.freeMark <Node>(b_mark);
  mesh.freeMark <Node>(i_mark); 
  mesh.freeMark<Node>(node_curve_mark);
  mesh.freeMark<Node>(node_point_mark);
  mesh.freeMark<Edge>(edge_curve_mark);

  gmds::IGMesh::node_iterator it = mesh.nodes_begin();
  double x_min = 100000;
  double y_min = 100000;
  double x_max = -100000;
  double y_max = -100000;
  for (; !it.isDone(); it.next())
    {
      Node n = it.value();
      gmds::math::Point p = n.getPoint();
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

  gmds::MeshModel model_cube(DIM3 | F | N | F2N);
  gmds::IGMesh mesh_cube(model_cube);
  for (it = mesh.nodes_begin(); !it.isDone(); it.next())
    {
      gmds::Node n = it.value();
      gmds::math::Point center = n.getPoint();
      //      if (m_mesh->isMarked(n, m_mark_alive))
      {
	gmds::math::Cross2D current_cross = (*field)[n.getID()];
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
    }
  //  gmds::VTKWriter<gmds::IGMesh> writer_cube(mesh_cube);
  //  writer_cube.write("final_level_sets", DIM3 | F | N);
}
/*----------------------------------------------------------------------------*/
TEST_F(SmoothingField2DTest, square_holes_laplace)
{
  MeshModel model(DIM3 | F | E | N | F2N | E2N | F2E | E2F | N2F | N2E);
  IGMesh mesh(model);
  
  //==================================================================
  // MESH READING
  //==================================================================
  MeditReader<IGMesh> reader(mesh);
  reader.read("Samples/square_holes.mesh", DIM3|F|N);
  //==================================================================
  // MESH TOPOLOGY PREPARATION
  //==================================================================3
  IGMeshDoctor doctor(&mesh);
  doctor.buildEdgesAndX2E();
  doctor.updateUpwardConnectivity();
  doctor.orient2DFaces();

  EXPECT_TRUE(mesh.getNbRegions()==0);  
  EXPECT_TRUE(mesh.getNbFaces()!=0);  
  EXPECT_TRUE(mesh.getNbEdges()!=0);
  EXPECT_TRUE(mesh.getNbNodes()!=0);
 //==================================================================
  // FRAME FIELD GENERATION
  //=================================================================
  CrossFieldGeneration2D algo(&mesh);
  algo.setDebugPrefix("SquareHoleLaplaceSmooth");
  algo.execute(CrossFieldGeneration2D::laplace_solve);

 //==================================================================
  // FRAME FIELD SMOOTHING
  //================================================================= 
  std::cout<<"===== SMOOTHING ====="<<std::endl;
  gmds::BoundaryOperator boundaryOp(&mesh); 

  int edge_curve_mark = mesh.getNewMark<gmds::Edge>();
  int node_curve_mark = mesh.getNewMark<gmds::Node>();
  int node_point_mark = mesh.getNewMark<gmds::Node>();
  int i_mark = mesh.getNewMark<gmds::Node>();
  int a_mark = mesh.getNewMark<gmds::Node>();
  int b_mark = mesh.getNewMark<gmds::Node>();
  std::vector<TCellID> b_ids;
  boundaryOp.markCellsOnCurves(edge_curve_mark, node_curve_mark);
  boundaryOp.markNodesOnPoint(edge_curve_mark, node_curve_mark, node_point_mark);
  boundaryOp.markAloneNodes(a_mark);
  
  gmds::IGMesh::node_iterator it_node = mesh.nodes_begin();
  for(;!it_node.isDone();it_node.next()){
    Node n = it_node.value();
    if(mesh.isMarked(n,node_curve_mark) && !mesh.isMarked(n,node_point_mark))
      mesh.mark(n, b_mark);
    else if (!mesh.isMarked(n,node_curve_mark) && !mesh.isMarked(n,node_point_mark))
      mesh.mark(n, i_mark);
  }
  
  gmds::Variable<gmds::math::Cross2D>* field  = 
    mesh.getVariable<math::Cross2D>(GMDS_NODE, "cross");
  // std::cout.precision(16);
  //  std::cout<<std::scientific;
  CrossFieldSmoother2D smooth(&mesh, field);
  smooth.initMarks(i_mark, b_mark);
  smooth.execute();

  mesh.unmarkAll<Node>(a_mark);
  mesh.unmarkAll<Node>(b_mark);
  mesh.unmarkAll<Node>(i_mark);
  mesh.unmarkAll<Node>(node_curve_mark);
  mesh.unmarkAll<Node>(node_point_mark);
  mesh.unmarkAll<Edge>(edge_curve_mark);

  mesh.freeMark <Node>(a_mark);
  mesh.freeMark <Node>(b_mark);
  mesh.freeMark <Node>(i_mark); 
  mesh.freeMark<Node>(node_curve_mark);
  mesh.freeMark<Node>(node_point_mark);
  mesh.freeMark<Edge>(edge_curve_mark);


  gmds::IGMesh::node_iterator it = mesh.nodes_begin();
  double x_min = 100000;
  double y_min = 100000;
  double x_max = -100000;
  double y_max = -100000;
  for (; !it.isDone(); it.next())
    {
      Node n = it.value();
      gmds::math::Point p = n.getPoint();
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

  gmds::MeshModel model_cube(DIM3 | F | N | F2N);
  gmds::IGMesh mesh_cube(model_cube);
  for (it = mesh.nodes_begin(); !it.isDone(); it.next())
    {
      gmds::Node n = it.value();
      gmds::math::Point center = n.getPoint();
      //      if (m_mesh->isMarked(n, m_mark_alive))
      {
	gmds::math::Cross2D current_cross = (*field)[n.getID()];
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
    }
  //  gmds::VTKWriter<gmds::IGMesh> writer_cube(mesh_cube);
  //  writer_cube.write("final_laplace", DIM3 | F | N);
}

/*----------------------------------------------------------------------------*/
TEST_F(SmoothingField2DTest, square_holes_perturb_laplace)
{
  MeshModel model(DIM3 | F | E | N | F2N | E2N | F2E | E2F | N2F | N2E);
  IGMesh mesh(model);
  
  //==================================================================
  // MESH READING
  //==================================================================
  MeditReader<IGMesh> reader(mesh);
  reader.read("Samples/square_holes_perturb.mesh", DIM3|F|N);
  //==================================================================
  // MESH TOPOLOGY PREPARATION
  //==================================================================3
  IGMeshDoctor doctor(&mesh);
  doctor.buildEdgesAndX2E();
  doctor.updateUpwardConnectivity();
  doctor.orient2DFaces();

  EXPECT_TRUE(mesh.getNbRegions()==0);  
  EXPECT_TRUE(mesh.getNbFaces()!=0);  
  EXPECT_TRUE(mesh.getNbEdges()!=0);
  EXPECT_TRUE(mesh.getNbNodes()!=0);
 //==================================================================
  // FRAME FIELD GENERATION
  //=================================================================
  CrossFieldGeneration2D algo(&mesh);
  algo.setDebugPrefix("SquareHolePerturbSmoothLaplace");
  algo.execute(CrossFieldGeneration2D::laplace_solve);

 //==================================================================
  // FRAME FIELD SMOOTHING
  //================================================================= 
  std::cout<<"===== SMOOTHING ====="<<std::endl;
  gmds::BoundaryOperator boundaryOp(&mesh); 

  int edge_curve_mark = mesh.getNewMark<gmds::Edge>();
  int node_curve_mark = mesh.getNewMark<gmds::Node>();
  int node_point_mark = mesh.getNewMark<gmds::Node>();
  int i_mark = mesh.getNewMark<gmds::Node>();
  int a_mark = mesh.getNewMark<gmds::Node>();
  int b_mark = mesh.getNewMark<gmds::Node>();
  std::vector<TCellID> b_ids;
  boundaryOp.markCellsOnCurves(edge_curve_mark, node_curve_mark);
  boundaryOp.markNodesOnPoint(edge_curve_mark, node_curve_mark, node_point_mark);
  boundaryOp.markAloneNodes(a_mark);
  
  gmds::IGMesh::node_iterator it_node = mesh.nodes_begin();
  for(;!it_node.isDone();it_node.next()){
    Node n = it_node.value();
    if(mesh.isMarked(n,node_curve_mark) && !mesh.isMarked(n,node_point_mark))
      mesh.mark(n, b_mark);
    else if (!mesh.isMarked(n,node_curve_mark) && !mesh.isMarked(n,node_point_mark))
      mesh.mark(n, i_mark);
  }
  
  gmds::Variable<gmds::math::Cross2D>* field  = 
    mesh.getVariable<math::Cross2D>(GMDS_NODE, "cross");
  // std::cout.precision(16);
  //  std::cout<<std::scientific;
  CrossFieldSmoother2D smooth(&mesh, field);
  smooth.initMarks(i_mark, b_mark);
  smooth.execute();

  mesh.unmarkAll<Node>(a_mark);
  mesh.unmarkAll<Node>(b_mark);
  mesh.unmarkAll<Node>(i_mark);
  mesh.unmarkAll<Node>(node_curve_mark);
  mesh.unmarkAll<Node>(node_point_mark);
  mesh.unmarkAll<Edge>(edge_curve_mark);

  mesh.freeMark <Node>(a_mark);
  mesh.freeMark <Node>(b_mark);
  mesh.freeMark <Node>(i_mark); 
  mesh.freeMark<Node>(node_curve_mark);
  mesh.freeMark<Node>(node_point_mark);
  mesh.freeMark<Edge>(edge_curve_mark);


  gmds::IGMesh::node_iterator it = mesh.nodes_begin();
  double x_min = 100000;
  double y_min = 100000;
  double x_max = -100000;
  double y_max = -100000;
  for (; !it.isDone(); it.next())
    {
      Node n = it.value();
      gmds::math::Point p = n.getPoint();
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

  gmds::MeshModel model_cube(DIM3 | F | N | F2N);
  gmds::IGMesh mesh_cube(model_cube);
  for (it = mesh.nodes_begin(); !it.isDone(); it.next())
    {
      gmds::Node n = it.value();
      gmds::math::Point center = n.getPoint();
      //      if (m_mesh->isMarked(n, m_mark_alive))
      {
	gmds::math::Cross2D current_cross = (*field)[n.getID()];
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
    }
  //  gmds::VTKWriter<gmds::IGMesh> writer_cube(mesh_cube);
  //  writer_cube.write("final_laplace", DIM3 | F | N);
}
/*----------------------------------------------------------------------------*/
TEST_F(SmoothingField2DTest, hecht_laplace)
{
    MeshModel model(DIM3 | F | E | N | F2N | E2N | F2E | E2F | N2F | N2E);
    IGMesh mesh(model);
    
    //==================================================================
    // MESH READING
    //==================================================================
    MeditReader<IGMesh> reader(mesh);
    reader.read("Samples/Hecht_example.mesh", DIM3|F|N);
    //==================================================================
    // MESH TOPOLOGY PREPARATION
    //==================================================================3
    IGMeshDoctor doctor(&mesh);
    doctor.buildEdgesAndX2E();
    doctor.updateUpwardConnectivity();
    doctor.orient2DFaces();
    
    EXPECT_TRUE(mesh.getNbRegions()==0);
    EXPECT_TRUE(mesh.getNbFaces()!=0);
    EXPECT_TRUE(mesh.getNbEdges()!=0);
    EXPECT_TRUE(mesh.getNbNodes()!=0);
    //==================================================================
    // FRAME FIELD GENERATION
    //=================================================================
    CrossFieldGeneration2D algo(&mesh);
    algo.setDebugPrefix("HechtExampleSmoothLaplace");
    algo.execute(CrossFieldGeneration2D::laplace_solve);
    
    //==================================================================
    // FRAME FIELD SMOOTHING
    //=================================================================
    std::cout<<"===== SMOOTHING ====="<<std::endl;
    gmds::BoundaryOperator boundaryOp(&mesh);
    
    int edge_curve_mark = mesh.getNewMark<gmds::Edge>();
    int node_curve_mark = mesh.getNewMark<gmds::Node>();
    int node_point_mark = mesh.getNewMark<gmds::Node>();
    int i_mark = mesh.getNewMark<gmds::Node>();
    int a_mark = mesh.getNewMark<gmds::Node>();
    int b_mark = mesh.getNewMark<gmds::Node>();
    std::vector<TCellID> b_ids;
    boundaryOp.markCellsOnCurves(edge_curve_mark, node_curve_mark);
    boundaryOp.markNodesOnPoint(edge_curve_mark, node_curve_mark, node_point_mark);
    boundaryOp.markAloneNodes(a_mark);
    
    gmds::IGMesh::node_iterator it_node = mesh.nodes_begin();
    for(;!it_node.isDone();it_node.next()){
        Node n = it_node.value();
        if(mesh.isMarked(n,node_curve_mark) && !mesh.isMarked(n,node_point_mark))
            mesh.mark(n, b_mark);
        else if (!mesh.isMarked(n,node_curve_mark) && !mesh.isMarked(n,node_point_mark) &&
                 !mesh.isMarked(n,a_mark))
            mesh.mark(n, i_mark);
    }
    
    gmds::Variable<gmds::math::Cross2D>* field  =
    mesh.getVariable<math::Cross2D>(GMDS_NODE, "cross");
    // std::cout.precision(16);
    //  std::cout<<std::scientific;
    CrossFieldSmoother2D smooth(&mesh, field);
    smooth.initMarks(i_mark, b_mark);
    smooth.execute();
    
    mesh.unmarkAll<Node>(a_mark);
    mesh.unmarkAll<Node>(b_mark);
    mesh.unmarkAll<Node>(i_mark);
    mesh.unmarkAll<Node>(node_curve_mark);
    mesh.unmarkAll<Node>(node_point_mark);
    mesh.unmarkAll<Edge>(edge_curve_mark);
    
    mesh.freeMark <Node>(a_mark);
    mesh.freeMark <Node>(b_mark);
    mesh.freeMark <Node>(i_mark);
    mesh.freeMark<Node>(node_curve_mark);
    mesh.freeMark<Node>(node_point_mark);
    mesh.freeMark<Edge>(edge_curve_mark);
    
    
    gmds::IGMesh::node_iterator it = mesh.nodes_begin();
    double x_min = 100000;
    double y_min = 100000;
    double x_max = -100000;
    double y_max = -100000;
    for (; !it.isDone(); it.next())
    {
        Node n = it.value();
        gmds::math::Point p = n.getPoint();
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
    
    gmds::MeshModel model_cube(DIM3 | F | N | F2N);
    gmds::IGMesh mesh_cube(model_cube);
    for (it = mesh.nodes_begin(); !it.isDone(); it.next())
    {
        gmds::Node n = it.value();
        gmds::math::Point center = n.getPoint();
        //      if (m_mesh->isMarked(n, m_mark_alive))
        {
            gmds::math::Cross2D current_cross = (*field)[n.getID()];
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
    }
    //  gmds::VTKWriter<gmds::IGMesh> writer_cube(mesh_cube);
    //  writer_cube.write("final_hecht", DIM3 | F | N);
}

/*----------------------------------------------------------------------------*/
TEST_F(SmoothingField2DTest, square_45degree_laplace)
{
    MeshModel model(DIM3 | F | E | N | F2N | E2N | F2E | E2F | N2F | N2E);
    IGMesh mesh(model);
    
    //==================================================================
    // MESH READING
    //==================================================================
    MeditReader<IGMesh> reader(mesh);
    reader.read("Samples/square_45degree.mesh", DIM3|F|N);
    //==================================================================
    // MESH TOPOLOGY PREPARATION
    //==================================================================3
    IGMeshDoctor doctor(&mesh);
    doctor.buildEdgesAndX2E();
    doctor.updateUpwardConnectivity();
    doctor.orient2DFaces();
    
    EXPECT_TRUE(mesh.getNbRegions()==0);
    EXPECT_TRUE(mesh.getNbFaces()!=0);
    EXPECT_TRUE(mesh.getNbEdges()!=0);
    EXPECT_TRUE(mesh.getNbNodes()!=0);
    //==================================================================
    // FRAME FIELD GENERATION
    //=================================================================
    CrossFieldGeneration2D algo(&mesh);
    algo.setDebugPrefix("Square45Laplace");
    algo.execute(CrossFieldGeneration2D::laplace_solve);
    
    //==================================================================
    // FRAME FIELD SMOOTHING
    //=================================================================
    std::cout<<"===== SMOOTHING ====="<<std::endl;
    gmds::BoundaryOperator boundaryOp(&mesh);
    
    int edge_curve_mark = mesh.getNewMark<gmds::Edge>();
    int node_curve_mark = mesh.getNewMark<gmds::Node>();
    int node_point_mark = mesh.getNewMark<gmds::Node>();
    int i_mark = mesh.getNewMark<gmds::Node>();
    int a_mark = mesh.getNewMark<gmds::Node>();
    int b_mark = mesh.getNewMark<gmds::Node>();
    std::vector<TCellID> b_ids;
    boundaryOp.markCellsOnCurves(edge_curve_mark, node_curve_mark);
    boundaryOp.markNodesOnPoint(edge_curve_mark, node_curve_mark, node_point_mark);
    boundaryOp.markAloneNodes(a_mark);
    
    gmds::IGMesh::node_iterator it_node = mesh.nodes_begin();
    for(;!it_node.isDone();it_node.next()){
        Node n = it_node.value();
        if(mesh.isMarked(n,node_curve_mark) && !mesh.isMarked(n,node_point_mark))
            mesh.mark(n, b_mark);
        else if (!mesh.isMarked(n,node_curve_mark) && !mesh.isMarked(n,node_point_mark) &&
                 !mesh.isMarked(n,a_mark))
            mesh.mark(n, i_mark);
    }
    
    gmds::Variable<gmds::math::Cross2D>* field  =
    mesh.getVariable<math::Cross2D>(GMDS_NODE, "cross");
    // std::cout.precision(16);
    //  std::cout<<std::scientific;
    CrossFieldSmoother2D smooth(&mesh, field);
    smooth.initMarks(i_mark, b_mark);
    smooth.execute();
    
    mesh.unmarkAll<Node>(a_mark);
    mesh.unmarkAll<Node>(b_mark);
    mesh.unmarkAll<Node>(i_mark);
    mesh.unmarkAll<Node>(node_curve_mark);
    mesh.unmarkAll<Node>(node_point_mark);
    mesh.unmarkAll<Edge>(edge_curve_mark);
    
    mesh.freeMark <Node>(a_mark);
    mesh.freeMark <Node>(b_mark);
    mesh.freeMark <Node>(i_mark);
    mesh.freeMark<Node>(node_curve_mark);
    mesh.freeMark<Node>(node_point_mark);
    mesh.freeMark<Edge>(edge_curve_mark);
    
    
    gmds::IGMesh::node_iterator it = mesh.nodes_begin();
    double x_min = 100000;
    double y_min = 100000;
    double x_max = -100000;
    double y_max = -100000;
    for (; !it.isDone(); it.next())
    {
        Node n = it.value();
        gmds::math::Point p = n.getPoint();
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
    
    gmds::MeshModel model_cube(DIM3 | F | N | F2N);
    gmds::IGMesh mesh_cube(model_cube);
    for (it = mesh.nodes_begin(); !it.isDone(); it.next())
    {
        gmds::Node n = it.value();
        gmds::math::Point center = n.getPoint();
        //      if (m_mesh->isMarked(n, m_mark_alive))
        {
            gmds::math::Cross2D current_cross = (*field)[n.getID()];
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
    }
    //  gmds::VTKWriter<gmds::IGMesh> writer_cube(mesh_cube);
    //  writer_cube.write("final_hecht", DIM3 | F | N);
}

/*----------------------------------------------------------------------------*/
TEST_F(SmoothingField2DTest, half_donut_laplace)
{
    MeshModel model(DIM3 | F | E | N | F2N | E2N | F2E | E2F | N2F | N2E);
    IGMesh mesh(model);
    
    //==================================================================
    // MESH READING
    //==================================================================
    MeditReader<IGMesh> reader(mesh);
    reader.read("Samples/half_donut.mesh", DIM3|F|N);
    //==================================================================
    // MESH TOPOLOGY PREPARATION
    //==================================================================3
    IGMeshDoctor doctor(&mesh);
    doctor.buildEdgesAndX2E();
    doctor.updateUpwardConnectivity();
    doctor.orient2DFaces();
    
    EXPECT_TRUE(mesh.getNbRegions()==0);
    EXPECT_TRUE(mesh.getNbFaces()!=0);
    EXPECT_TRUE(mesh.getNbEdges()!=0);
    EXPECT_TRUE(mesh.getNbNodes()!=0);
    //==================================================================
    // FRAME FIELD GENERATION
    //=================================================================
    CrossFieldGeneration2D algo(&mesh);
    algo.setDebugPrefix("HalfDonut");
    algo.execute(CrossFieldGeneration2D::laplace_solve);
    
    //==================================================================
    // FRAME FIELD SMOOTHING
    //=================================================================
    std::cout<<"===== SMOOTHING ====="<<std::endl;
    gmds::BoundaryOperator boundaryOp(&mesh);
    
    int edge_curve_mark = mesh.getNewMark<gmds::Edge>();
    int node_curve_mark = mesh.getNewMark<gmds::Node>();
    int node_point_mark = mesh.getNewMark<gmds::Node>();
    int i_mark = mesh.getNewMark<gmds::Node>();
    int a_mark = mesh.getNewMark<gmds::Node>();
    int b_mark = mesh.getNewMark<gmds::Node>();
    std::vector<TCellID> b_ids;
    boundaryOp.markCellsOnCurves(edge_curve_mark, node_curve_mark);
    boundaryOp.markNodesOnPoint(edge_curve_mark, node_curve_mark, node_point_mark);
    boundaryOp.markAloneNodes(a_mark);
    
    gmds::IGMesh::node_iterator it_node = mesh.nodes_begin();
    for(;!it_node.isDone();it_node.next()){
        Node n = it_node.value();
        if(mesh.isMarked(n,node_curve_mark) && !mesh.isMarked(n,node_point_mark))
            mesh.mark(n, b_mark);
        else if (!mesh.isMarked(n,node_curve_mark) && !mesh.isMarked(n,node_point_mark) &&
                 !mesh.isMarked(n,a_mark))
            mesh.mark(n, i_mark);
    }
    
    gmds::Variable<gmds::math::Cross2D>* field  =
    mesh.getVariable<math::Cross2D>(GMDS_NODE, "cross");
    // std::cout.precision(16);
    //  std::cout<<std::scientific;
    CrossFieldSmoother2D smooth(&mesh, field);
    smooth.initMarks(i_mark, b_mark);
    smooth.execute();
    
    mesh.unmarkAll<Node>(a_mark);
    mesh.unmarkAll<Node>(b_mark);
    mesh.unmarkAll<Node>(i_mark);
    mesh.unmarkAll<Node>(node_curve_mark);
    mesh.unmarkAll<Node>(node_point_mark);
    mesh.unmarkAll<Edge>(edge_curve_mark);
    
    mesh.freeMark <Node>(a_mark);
    mesh.freeMark <Node>(b_mark);
    mesh.freeMark <Node>(i_mark);
    mesh.freeMark<Node>(node_curve_mark);
    mesh.freeMark<Node>(node_point_mark);
    mesh.freeMark<Edge>(edge_curve_mark);
    
   
}
/*----------------------------------------------------------------------------*/
TEST_F(SmoothingField2DTest, nautilus_laplace)
{
  MeshModel model(DIM3 | F | E | N | F2N | E2N | F2E | E2F | N2F | N2E);
  IGMesh mesh(model);
  
  //==================================================================
  // MESH READING
  //==================================================================
  MeditReader<IGMesh> reader(mesh);
  reader.read("Samples/nautilus.mesh", DIM3|F|N);
  //==================================================================
  // MESH TOPOLOGY PREPARATION
  //==================================================================3
  IGMeshDoctor doctor(&mesh);
  doctor.buildEdgesAndX2E();
  doctor.updateUpwardConnectivity();
  doctor.orient2DFaces();

  EXPECT_TRUE(mesh.getNbRegions()==0);  
  EXPECT_TRUE(mesh.getNbFaces()!=0);  
  EXPECT_TRUE(mesh.getNbEdges()!=0);
  EXPECT_TRUE(mesh.getNbNodes()!=0);
 //==================================================================
  // FRAME FIELD GENERATION
  //=================================================================
  CrossFieldGeneration2D algo(&mesh);
  algo.setDebugPrefix("NautilusSmoothLaplace");
  algo.execute(CrossFieldGeneration2D::laplace_solve);

 //==================================================================
  // FRAME FIELD SMOOTHING
  //================================================================= 
  std::cout<<"===== SMOOTHING ====="<<std::endl;
  gmds::BoundaryOperator boundaryOp(&mesh); 

  int edge_curve_mark = mesh.getNewMark<gmds::Edge>();
  int node_curve_mark = mesh.getNewMark<gmds::Node>();
  int node_point_mark = mesh.getNewMark<gmds::Node>();
  int i_mark = mesh.getNewMark<gmds::Node>();
  int a_mark = mesh.getNewMark<gmds::Node>();
  int b_mark = mesh.getNewMark<gmds::Node>();
  std::vector<TCellID> b_ids;
  boundaryOp.markCellsOnCurves(edge_curve_mark, node_curve_mark);
  boundaryOp.markNodesOnPoint(edge_curve_mark, node_curve_mark, node_point_mark);
  boundaryOp.markAloneNodes(a_mark);
  
  gmds::IGMesh::node_iterator it_node = mesh.nodes_begin();
  for(;!it_node.isDone();it_node.next()){
    Node n = it_node.value();
    if(mesh.isMarked(n,node_curve_mark) && !mesh.isMarked(n,node_point_mark))
      mesh.mark(n, b_mark);
    else if (!mesh.isMarked(n,node_curve_mark) && !mesh.isMarked(n,node_point_mark) &&
	     !mesh.isMarked(n,a_mark))
      mesh.mark(n, i_mark);
  }
  
  gmds::Variable<gmds::math::Cross2D>* field  = 
    mesh.getVariable<math::Cross2D>(GMDS_NODE, "cross");
  // std::cout.precision(16);
  //  std::cout<<std::scientific;
  CrossFieldSmoother2D smooth(&mesh, field);
  smooth.initMarks(i_mark, b_mark);
  smooth.execute();

  mesh.unmarkAll<Node>(a_mark);
  mesh.unmarkAll<Node>(b_mark);
  mesh.unmarkAll<Node>(i_mark);
  mesh.unmarkAll<Node>(node_curve_mark);
  mesh.unmarkAll<Node>(node_point_mark);
  mesh.unmarkAll<Edge>(edge_curve_mark);

  mesh.freeMark <Node>(a_mark);
  mesh.freeMark <Node>(b_mark);
  mesh.freeMark <Node>(i_mark); 
  mesh.freeMark<Node>(node_curve_mark);
  mesh.freeMark<Node>(node_point_mark);
  mesh.freeMark<Edge>(edge_curve_mark);


  gmds::IGMesh::node_iterator it = mesh.nodes_begin();
  double x_min = 100000;
  double y_min = 100000;
  double x_max = -100000;
  double y_max = -100000;
  for (; !it.isDone(); it.next())
    {
      Node n = it.value();
      gmds::math::Point p = n.getPoint();
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

  gmds::MeshModel model_cube(DIM3 | F | N | F2N);
  gmds::IGMesh mesh_cube(model_cube);
  for (it = mesh.nodes_begin(); !it.isDone(); it.next())
    {
      gmds::Node n = it.value();
      gmds::math::Point center = n.getPoint();
      //      if (m_mesh->isMarked(n, m_mark_alive))
      {
	gmds::math::Cross2D current_cross = (*field)[n.getID()];
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
    }
  //  gmds::VTKWriter<gmds::IGMesh> writer_cube(mesh_cube);
  //  writer_cube.write("final_nautilus", DIM3 | F | N);
}
/*----------------------------------------------------------------------------*/
TEST_F(SmoothingField2DTest, half_disk_laplace)
{
  MeshModel model(DIM3 | F | E | N | F2N | E2N | F2E | E2F | N2F | N2E);
  IGMesh mesh(model);
  
  //==================================================================
  // MESH READING
  //==================================================================
  MeditReader<IGMesh> reader(mesh);
  reader.read("Samples/half_disk.mesh", DIM3|F|N);
  //==================================================================
  // MESH TOPOLOGY PREPARATION
  //==================================================================3
  IGMeshDoctor doctor(&mesh);
  doctor.buildEdgesAndX2E();
  doctor.updateUpwardConnectivity();
  doctor.orient2DFaces();

  EXPECT_TRUE(mesh.getNbRegions()==0);  
  EXPECT_TRUE(mesh.getNbFaces()!=0);  
  EXPECT_TRUE(mesh.getNbEdges()!=0);
  EXPECT_TRUE(mesh.getNbNodes()!=0);
 //==================================================================
  // FRAME FIELD GENERATION
  //=================================================================
  CrossFieldGeneration2D algo(&mesh);
  algo.setDebugPrefix("HalfDiskSmoothLaplace");
  algo.execute(CrossFieldGeneration2D::laplace_solve);

 //==================================================================
  // FRAME FIELD SMOOTHING
  //================================================================= 
  std::cout<<"===== SMOOTHING ====="<<std::endl;
  gmds::BoundaryOperator boundaryOp(&mesh); 

  int edge_curve_mark = mesh.getNewMark<gmds::Edge>();
  int node_curve_mark = mesh.getNewMark<gmds::Node>();
  int node_point_mark = mesh.getNewMark<gmds::Node>();
  int i_mark = mesh.getNewMark<gmds::Node>();
  int a_mark = mesh.getNewMark<gmds::Node>();
  int b_mark = mesh.getNewMark<gmds::Node>();
  std::vector<TCellID> b_ids;
  boundaryOp.markCellsOnCurves(edge_curve_mark, node_curve_mark);
  boundaryOp.markNodesOnPoint(edge_curve_mark, node_curve_mark, node_point_mark);
  boundaryOp.markAloneNodes(a_mark);
  
  gmds::IGMesh::node_iterator it_node = mesh.nodes_begin();
  for(;!it_node.isDone();it_node.next()){
    Node n = it_node.value();
    if(mesh.isMarked(n,node_curve_mark) && !mesh.isMarked(n,node_point_mark))
      mesh.mark(n, b_mark);
    else if (!mesh.isMarked(n,node_curve_mark) && !mesh.isMarked(n,node_point_mark) &&
	     !mesh.isMarked(n,a_mark))
      mesh.mark(n, i_mark);
  }
  
  gmds::Variable<gmds::math::Cross2D>* field  = 
    mesh.getVariable<math::Cross2D>(GMDS_NODE, "cross");
  // std::cout.precision(16);
  //  std::cout<<std::scientific;
  CrossFieldSmoother2D smooth(&mesh, field);
  smooth.initMarks(i_mark, b_mark);
  smooth.execute();

  mesh.unmarkAll<Node>(a_mark);
  mesh.unmarkAll<Node>(b_mark);
  mesh.unmarkAll<Node>(i_mark);
  mesh.unmarkAll<Node>(node_curve_mark);
  mesh.unmarkAll<Node>(node_point_mark);
  mesh.unmarkAll<Edge>(edge_curve_mark);

  mesh.freeMark <Node>(a_mark);
  mesh.freeMark <Node>(b_mark);
  mesh.freeMark <Node>(i_mark); 
  mesh.freeMark<Node>(node_curve_mark);
  mesh.freeMark<Node>(node_point_mark);
  mesh.freeMark<Edge>(edge_curve_mark);


  gmds::IGMesh::node_iterator it = mesh.nodes_begin();
  double x_min = 100000;
  double y_min = 100000;
  double x_max = -100000;
  double y_max = -100000;
  for (; !it.isDone(); it.next())
    {
      Node n = it.value();
      gmds::math::Point p = n.getPoint();
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

  gmds::MeshModel model_cube(DIM3 | F | N | F2N);
  gmds::IGMesh mesh_cube(model_cube);
  for (it = mesh.nodes_begin(); !it.isDone(); it.next())
    {
      gmds::Node n = it.value();
      gmds::math::Point center = n.getPoint();
      //      if (m_mesh->isMarked(n, m_mark_alive))
      {
	gmds::math::Cross2D current_cross = (*field)[n.getID()];
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
    }
  //  gmds::VTKWriter<gmds::IGMesh> writer_cube(mesh_cube);
  //  writer_cube.write("final_nautilus", DIM3 | F | N);
}
/*----------------------------------------------------------------------------*/
