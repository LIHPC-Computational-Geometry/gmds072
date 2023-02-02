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
/*----------------------------------------------------------------------------*/
#include<string>
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IG.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
class GroupTest: public ::testing::Test {

  protected:
	GroupTest(){;}
    virtual ~GroupTest(){;}
};
/*----------------------------------------------------------------------------*/
TEST_F(GroupTest,surfaceCreation) {
	MeshModel mod = DIM3|N|F|F2N;
	IGMesh mesh(mod);

	Node t[3];
	Face f[1];

	t[0] = mesh.newNode(0.,0.,0.);
	t[1] = mesh.newNode(0.,1.,0.);
	t[2] = mesh.newNode(1.,1.,0.);

	f[0] = mesh.newTriangle(t[0],t[1],t[2]);

	std::string surfaceName = "bob";
	IGMesh::surface& surf = mesh.newSurface(surfaceName.c_str());
	surf.add(f[0]);

	EXPECT_EQ(1,mesh.getNbSurfaces());
}
/*----------------------------------------------------------------------------*/
TEST_F(GroupTest,surfaceDelete) {
	MeshModel mod = DIM3|N|F|F2N;
	IGMesh mesh(mod);

	Node t[4];
	Face f[2];
	t[0] = mesh.newNode(0.,0.,0.);
	t[1] = mesh.newNode(0.,1.,0.);
	t[2] = mesh.newNode(1.,1.,0.);
	t[3] = mesh.newNode(1.,1.,2.);

	f[0] = mesh.newTriangle(t[0],t[1],t[2]);
	f[1] = mesh.newTriangle(t[0],t[1],t[3]);

	std::string surfaceName = "bob";
	IGMesh::surface& surf = mesh.newSurface(surfaceName.c_str());
	surf.add(f[0]);

	mesh.deleteSurface(surf);

	EXPECT_EQ(0,mesh.getNbSurfaces());

	surfaceName = "bill";
	IGMesh::surface& surf1 = mesh.newSurface(surfaceName.c_str());
	surf1.add(f[0]);
	EXPECT_EQ(1, mesh.getNbSurfaces());
	surfaceName = "bab";
	IGMesh::surface& surf2 = mesh.newSurface(surfaceName.c_str());
	surf2.add(f[1]);

	EXPECT_EQ(2, mesh.getNbSurfaces());

	IGMesh::surfaces_iterator its   = mesh.surfaces_begin();

	while(its!=mesh.surfaces_end())
	{
		mesh.deleteSurface(*its);
		its   = mesh.surfaces_begin();
	}

	EXPECT_EQ(0,mesh.getNbSurfaces());

}
/*----------------------------------------------------------------------------*/
TEST_F(GroupTest,surfaceBoundingBox) {
	MeshModel mod = DIM3|N|F|F2N;
	IGMesh mesh(mod);

	Node t[3];
	Face f[1];

	t[0] = mesh.newNode(0.,0.,0.);
	t[1] = mesh.newNode(0.,1.,0.);
	t[2] = mesh.newNode(1.,1.,0.);

	f[0] = mesh.newTriangle(t[0],t[1],t[2]);

	std::string surfaceName = "bob";
	IGMesh::surface& surf = mesh.newSurface(surfaceName.c_str());
	surf.add(f[0]);

	TCoord minCoords[3];
	TCoord maxCoords[3];
	surf.computeBoundingBox(minCoords,maxCoords);

	EXPECT_TRUE(minCoords[0]==0.);
	EXPECT_TRUE(minCoords[1]==0.);
	EXPECT_TRUE(minCoords[2]==0.);
	EXPECT_TRUE(maxCoords[0]==1.);
	EXPECT_TRUE(maxCoords[1]==1.);
	EXPECT_TRUE(maxCoords[2]==0.);
}
/*----------------------------------------------------------------------------*/
