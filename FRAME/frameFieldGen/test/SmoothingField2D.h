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
/*----------------------------------------------------------------------------*/
#include "CrossFieldGeneration2D.h"
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
TEST_F(SmoothingField2DTest, square_holes_level_sets)
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
  // FRAME FIELD GENERATION
  //==================================================================
  CrossFieldGeneration2D algo(&mesh);
  algo.setDebugPrefix("SquareHolesLevelSets");
  algo.execute(CrossFieldGeneration2D::level_sets);
}
/*----------------------------------------------------------------------------*/
TEST_F(SmoothingField2DTest, DISABLED_square_holes_coarse_laplace)
{
  MeshModel model(DIM3 | F | E | N | F2N | E2N | F2E | E2F | N2F | N2E);
  IGMesh mesh(model);
  
  //==================================================================
  // MESH READING
  //==================================================================
  MeditReader<IGMesh> reader(mesh);
  reader.read("Samples/square_holes_coarse.mesh", DIM3|F|N);
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
  // FRAME FIELD GENERATION
  //==================================================================
  CrossFieldGeneration2D algo(&mesh);
  algo.setDebugPrefix("SquareHolesCoarseLaplace");
  algo.execute(CrossFieldGeneration2D::laplace_solve);
}
/*----------------------------------------------------------------------------*/
