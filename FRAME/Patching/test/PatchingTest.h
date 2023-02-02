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
 * PatchingTest.h
 *
 *  Created on: June 17, 2015
 *      Author: ledouxf
 */

/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IGMesh.h>
#include <GMDS/IG/IGMeshDoctor.h>
#include <GMDS/IO/VTKReader.h>
#include <GMDS/IO/VTKWriter.h>
#include <GMDS/Algo/BoundaryOperator.h>
/*----------------------------------------------------------------------------*/
#include "Patching2D.h"
/*----------------------------------------------------------------------------*/
#include <iostream> 
#include <vector>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
class PatchingTest: public ::testing::Test 
{
 protected:
  PatchingTest(){;}
  virtual ~PatchingTest(){;}
};
/*----------------------------------------------------------------------------*/
TEST_F(PatchingTest, DISABLED_quad)
{
    MeshModel model(DIM3 | F | E | N | F2N | E2N | F2E | E2F | N2F | N2E);
    IGMesh mesh(model);
    
    //==================================================================
    // MESH READING
    //==================================================================
    VTKReader<IGMesh> reader(mesh);
    reader.read("Samples/quad_field.vtk");
    //==================================================================
    // MESH TOPOLOGY PREPARATION
    //==================================================================
    IGMeshDoctor doctor(&mesh);
    doctor.buildEdgesAndX2E();
    doctor.updateUpwardConnectivity();
    doctor.orient2DFaces();
    
    EXPECT_TRUE(mesh.getNbRegions()==0);
    EXPECT_TRUE(mesh.getNbFaces()!=0);
    EXPECT_TRUE(mesh.getNbEdges()!=0);
    EXPECT_TRUE(mesh.getNbNodes()!=0);
    
    //==================================================================
    // INIT MARKs FOR BOUNDARY NODES
    //==================================================================
    int edge_curve_mark = mesh.getNewMark<gmds::Edge>();
    int node_curve_mark = mesh.getNewMark<gmds::Node>();
    int node_point_mark = mesh.getNewMark<gmds::Node>();
    
    BoundaryOperator boundaryOp(&mesh);
    boundaryOp.markCellsOnCurves(edge_curve_mark, node_curve_mark);
    boundaryOp.markNodesOnPoint(edge_curve_mark, node_curve_mark, node_point_mark);
    
    
    //==================================================================
    // CROSS FIELD EXTRACTION FROM THE IMPORTED FILE
    //==================================================================
    Variable<math::Vector>* field_X =
    mesh.getVariable<math::Vector>(gmds::GMDS_NODE,"cross_X");
    
    Variable<math::Cross2D>* field =
    mesh.newVariable<math::Cross2D>(gmds::GMDS_NODE,"c");
    
    IGMesh::node_iterator it_nodes = mesh.nodes_begin();
    math::Vector OX(1,0,0);
    while(!it_nodes.isDone()){
        Node n = it_nodes.value();
        math::Vector vx = (*field_X)[n.getID()];
        (*field)[n.getID()]= math::Cross2D(4*vx.angle(OX));
        it_nodes.next();
    }
    
    //==================================================================
    // PATCHING ALGO USED TO GENERATE FULL-QUAD MESHES
    //==================================================================
    Patching2D algo(&mesh, field);
    algo.setDebugPrefix("PatchingQuad");
    algo.initMarks(node_point_mark, node_curve_mark, edge_curve_mark);
    algo.execute();
    //==================================================================
    // MARKS CLEANING
    //==================================================================
    mesh.unmarkAll<Node>(node_curve_mark);
    mesh.unmarkAll<Node>(node_point_mark);
    mesh.unmarkAll<Edge>(edge_curve_mark);
    mesh.freeMark<Node>(node_curve_mark);
    mesh.freeMark<Node>(node_point_mark);
    mesh.freeMark<Edge>(edge_curve_mark);
    
}
/*----------------------------------------------------------------------------*/
TEST_F(PatchingTest, DISABLED_half_disk)
{
    MeshModel model(DIM3 | F | E | N | F2N | E2N | F2E | E2F | N2F | N2E);
    IGMesh mesh(model);
    
    //==================================================================
    // MESH READING
    //==================================================================
    VTKReader<IGMesh> reader(mesh);
    reader.read("Samples/half_disk.vtk");
    //==================================================================
    // MESH TOPOLOGY PREPARATION
    //==================================================================
    IGMeshDoctor doctor(&mesh);
    doctor.buildEdgesAndX2E();
    doctor.updateUpwardConnectivity();
    doctor.orient2DFaces();
    
    EXPECT_TRUE(mesh.getNbRegions()==0);
    EXPECT_TRUE(mesh.getNbFaces()!=0);
    EXPECT_TRUE(mesh.getNbEdges()!=0);
    EXPECT_TRUE(mesh.getNbNodes()!=0);
    
    //==================================================================
    // INIT MARKs FOR BOUNDARY NODES
    //==================================================================
    int edge_curve_mark = mesh.getNewMark<gmds::Edge>();
    int node_curve_mark = mesh.getNewMark<gmds::Node>();
    int node_point_mark = mesh.getNewMark<gmds::Node>();
    
    BoundaryOperator boundaryOp(&mesh);
    boundaryOp.markCellsOnCurves(edge_curve_mark, node_curve_mark);
    boundaryOp.markNodesOnPoint(edge_curve_mark, node_curve_mark, node_point_mark);
    
    
    //==================================================================
    // CROSS FIELD EXTRACTION FROM THE IMPORTED FILE
    //==================================================================
    Variable<math::Vector>* field_X =
    mesh.getVariable<math::Vector>(gmds::GMDS_NODE,"cross_X");
    
    Variable<math::Cross2D>* field =
    mesh.newVariable<math::Cross2D>(gmds::GMDS_NODE,"c");
    
    IGMesh::node_iterator it_nodes = mesh.nodes_begin();
    math::Vector OX(1,0,0);
    while(!it_nodes.isDone()){
        Node n = it_nodes.value();
        math::Vector vx = (*field_X)[n.getID()];
        (*field)[n.getID()]= math::Cross2D(4*vx.angle(OX));
        it_nodes.next();
    }
    
    //==================================================================
    // PATCHING ALGO USED TO GENERATE FULL-QUAD MESHES
    //==================================================================
    Patching2D algo(&mesh, field);
    algo.setDebugPrefix("PatchingHalfCircle");
    algo.initMarks(node_point_mark, node_curve_mark, edge_curve_mark);
    algo.execute();
    //==================================================================
    // MARKS CLEANING
    //==================================================================
    mesh.unmarkAll<Node>(node_curve_mark);
    mesh.unmarkAll<Node>(node_point_mark);
    mesh.unmarkAll<Edge>(edge_curve_mark);
    mesh.freeMark<Node>(node_curve_mark);
    mesh.freeMark<Node>(node_point_mark);
    mesh.freeMark<Edge>(edge_curve_mark);
    
}
/*----------------------------------------------------------------------------*/
