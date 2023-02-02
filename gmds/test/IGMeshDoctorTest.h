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
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IG.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
class IGMeshDoctorTest: public ::testing::Test {

  protected:
	IGMeshDoctorTest(){;}
    virtual ~IGMeshDoctorTest(){;}
};
/*----------------------------------------------------------------------------*/
TEST_F(IGMeshDoctorTest,testTET) {
	MeshModel mod = DIM3|N|F|R|F2N|R2F|R2N;
	IGMesh mesh(mod);

	Node 	t[5];
	Region 	r[2];

	t[0] = mesh.newNode(0.,0.,0.);
	t[1] = mesh.newNode(0.,1.,0.);
	t[2] = mesh.newNode(1.,1.,0.);

	t[3] = mesh.newNode(0.5,0.5,-1.);
	t[4] = mesh.newNode(0.5,0.5,1.);

	r[0] = mesh.newTet(t[0],t[1],t[2],t[3]);
	r[1] = mesh.newTet(t[0],t[1],t[2],t[4]);

	EXPECT_EQ(0,mesh.getNbFaces());

	IGMeshDoctor doc(&mesh);
	doc.buildFacesAndR2F();

	EXPECT_EQ(7, mesh.getNbFaces());
}
/*----------------------------------------------------------------------------*/
TEST_F(IGMeshDoctorTest,testHEX) {
	MeshModel mod = DIM3|N|F|R|F2N|R2F|R2N;
	IGMesh mesh(mod);
	const int si = 3;
	const int sj = 2;
	const int sk = 2;
	Node 	t[si][sj][sk];
	Region 	r[si-1][sj-1][sk-1];

	/* nodes creation */
	for(int i=0;i<si;i++)
		for(int j=0;j<sj;j++)
			for(int k=0;k<sk;k++){
				t[i][j][k] = mesh.newNode(i,j,k);
			}

	for(int i=1;i<si;i++)
		for(int j=1;j<sj;j++)
			for(int k=1;k<sk;k++){
				r[i-1][j-1][k-1] = mesh.newHex(t[i-1][j-1][k],t[i-1][j][k],t[i][j][k],t[i][j-1][k],
						t[i-1][j-1][k-1],t[i-1][j][k-1],t[i][j][k-1],t[i][j-1][k-1]);
			}

	EXPECT_EQ(0,mesh.getNbFaces());

	IGMeshDoctor doc(&mesh);
	doc.buildFacesAndR2F();

	EXPECT_EQ(11,mesh.getNbFaces());
}
/*----------------------------------------------------------------------------*/
TEST_F(IGMeshDoctorTest,testPYRAMID_1) {
	MeshModel mod = DIM3|N|F|R|F2N|R2F|R2N;
	IGMesh mesh(mod);

	Node	t[5];
	Region	r;

	t[0] = mesh.newNode(0.,0.,0.);
	t[1] = mesh.newNode(0.,1.,0.);
	t[2] = mesh.newNode(1.,1.,0.);
	t[3] = mesh.newNode(1.,0.,0.);

	t[4] = mesh.newNode(0.5,0.5,-1.);


	r = mesh.newPyramid(t[3],t[2],t[1],t[0],t[4]);

		EXPECT_EQ(0, mesh.getNbFaces());

	IGMeshDoctor doc(&mesh);
	doc.buildFacesAndR2F();

	EXPECT_EQ(5, mesh.getNbFaces());
	EXPECT_EQ(1, mesh.getNbQuadrilaterals());
	EXPECT_EQ(4, mesh.getNbTriangles());

}
/*----------------------------------------------------------------------------*/
TEST_F(IGMeshDoctorTest,testPYRAMID_2) {
	MeshModel mod = DIM3|N|F|R|R2F|R2N|F2N;
	IGMesh mesh(mod);

	Node	t[6];
	Region	r[2];

	t[0] = mesh.newNode(0.,0.,0.);
	t[1] = mesh.newNode(0.,1.,0.);
	t[2] = mesh.newNode(1.,1.,0.);
	t[3] = mesh.newNode(1.,0.,0.);

	t[4] = mesh.newNode(0.5,0.5,-1.);
	t[5] = mesh.newNode(0.5,0.5,1.);

	r[0] = mesh.newPyramid(t[3],t[2],t[1],t[0],t[4]);
	r[1] = mesh.newPyramid(t[0],t[1],t[2],t[3],t[5]);

	EXPECT_EQ(0, mesh.getNbFaces());

	IGMeshDoctor doc(&mesh);
	doc.buildFacesAndR2F();

	EXPECT_EQ(9, mesh.getNbFaces());
	EXPECT_EQ(1, mesh.getNbQuadrilaterals());
	EXPECT_EQ(8, mesh.getNbTriangles());

}
/*----------------------------------------------------------------------------*/
