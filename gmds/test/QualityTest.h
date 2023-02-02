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
class QualityTest: public ::testing::Test {

  protected:
	QualityTest(){;}
    virtual ~QualityTest(){;}
};
/*----------------------------------------------------------------------------*/
TEST_F(QualityTest,hexJacobian) {
	MeshModel mod = DIM2|N|R|R2N;
	IGMesh mesh(mod);
	IGMeshQualityEvaluation qualityEval;

	Node n1 = mesh.newNode(0,0,0);
	Node n2 = mesh.newNode(1,0,0);
	Node n3 = mesh.newNode(1,1,0);
	Node n4 = mesh.newNode(0,1,0);
	Node n5 = mesh.newNode(0,0,1);
	Node n6 = mesh.newNode(1,0,1);
	Node n7 = mesh.newNode(1,1,1);
	Node n8 = mesh.newNode(0,1,1);

	Region r = mesh.newHex(n1,n2,n3,n4,n5,n6,n7,n8);
	EXPECT_EQ(1,qualityEval.jacobian(r));
	EXPECT_EQ(1,qualityEval.scaledJacobian(r));
}
/*----------------------------------------------------------------------------*/
TEST_F(QualityTest,hexJacobian2) {
	MeshModel mod = DIM2|N|R|R2N;
	IGMesh mesh(mod);
	IGMeshQualityEvaluation qualityEval;

	Node n1 = mesh.newNode(0,0,0);
	Node n2 = mesh.newNode(1,0,0);
	Node n3 = mesh.newNode(1,1,0);
	Node n4 = mesh.newNode(0,1,0);
	Node n5 = mesh.newNode(0,0,2);
	Node n6 = mesh.newNode(1,0,2);
	Node n7 = mesh.newNode(1,1,2);
	Node n8 = mesh.newNode(0,1,2);

	Region r = mesh.newHex(n1,n2,n3,n4,n5,n6,n7,n8);
	EXPECT_EQ(2,qualityEval.jacobian(r));
	EXPECT_EQ(1,qualityEval.scaledJacobian(r));
}
/*----------------------------------------------------------------------------*/
TEST_F(QualityTest,hexJacobian3) {
	MeshModel mod = DIM2|N|R|R2N;
	IGMesh mesh(mod);
	IGMeshQualityEvaluation qualityEval;

	Node n1 = mesh.newNode(0,0,0);
	Node n2 = mesh.newNode(4,0,0);
	Node n3 = mesh.newNode(4,3,0);
	Node n4 = mesh.newNode(0,3,0);
	Node n5 = mesh.newNode(0,0,2);
	Node n6 = mesh.newNode(4,0,2);
	Node n7 = mesh.newNode(4,3,2);
	Node n8 = mesh.newNode(0,3,2);

	Region r = mesh.newHex(n1,n2,n3,n4,n5,n6,n7,n8);
	EXPECT_EQ(1,qualityEval.scaledJacobian(r));
}
/*----------------------------------------------------------------------------*/
TEST_F(QualityTest,hexJacobian4) {
	MeshModel mod = DIM2|N|R|R2N;
	IGMesh mesh(mod);
	IGMeshQualityEvaluation qualityEval;

	TCoord perturb = 0.8;
	Node n1 = mesh.newNode(0,0,0);
	Node n2 = mesh.newNode(1,0,0);
	Node n3 = mesh.newNode(perturb,perturb,0);
	Node n4 = mesh.newNode(0,1,0);
	Node n5 = mesh.newNode(0,0,1);
	Node n6 = mesh.newNode(1,0,1);
	Node n7 = mesh.newNode(perturb,perturb,1);
	Node n8 = mesh.newNode(0,1,1);

	Region r = mesh.newHex(n1,n2,n3,n4,n5,n6,n7,n8);

	EXPECT_TRUE(qualityEval.scaledJacobian(r)<1 && qualityEval.scaledJacobian(r)>0);
	perturb = 0.6;
	n3.X()=perturb;
	n3.Y()=perturb;
	n7.X()=perturb;
	n7.Y()=perturb;
	EXPECT_TRUE(qualityEval.scaledJacobian(r)<1 && qualityEval.scaledJacobian(r)>0);

	perturb = 0.2;
	n3.X()=perturb;
	n3.Y()=perturb;
	n7.X()=perturb;
	n7.Y()=perturb;

	EXPECT_TRUE(qualityEval.scaledJacobian(r)<0);

	perturb = -0.2;
	n3.X()=perturb;
	n3.Y()=perturb;
	n7.X()=perturb;
	n7.Y()=perturb;

	EXPECT_TRUE(qualityEval.scaledJacobian(r)<0);

	perturb = 0.5;
	n3.X()=perturb;
	n3.Y()=perturb;
	n7.X()=perturb;
	n7.Y()=perturb;

	EXPECT_TRUE(isZero(qualityEval.scaledJacobian(r)));
}
/*----------------------------------------------------------------------------*/
TEST_F(QualityTest,quadMinAngle) {
	MeshModel mod = DIM2|N|F|F2N;
	IGMesh mesh(mod);
	IGMeshQualityEvaluation qualityEval;

	Node n1 = mesh.newNode(0,0,0);
	Node n2 = mesh.newNode(1,0,0);
	Node n3 = mesh.newNode(1,1,0);
	Node n4 = mesh.newNode(0,1,0);

	Face f = mesh.newQuad(n1,n2,n3,n4);
	EXPECT_TRUE(areEquals(math::Constants::PIDIV2,qualityEval.minAngle(f)));

	n1.Y()=n1.Y()-1;
	EXPECT_TRUE(areEquals(math::Constants::PIDIV4,qualityEval.minAngle(f)));

}
/*----------------------------------------------------------------------------*/
TEST_F(QualityTest,quadMinAngle2) {
	MeshModel mod = DIM2|N|F|F2N;
	IGMesh mesh(mod);
	IGMeshQualityEvaluation qualityEval;

	Node n1 = mesh.newNode(0,0,0);
	Node n2 = mesh.newNode(1,0,0);
	Node n3 = mesh.newNode(1,1,0);
	Node n4 = mesh.newNode(0,1,0);

	Face f = mesh.newQuad(n4,n3,n2,n1);
	EXPECT_TRUE(areEquals(math::Constants::PIDIV2,qualityEval.minAngle(f)));

	n1.Y()=n1.Y()-1;
	EXPECT_TRUE(areEquals(math::Constants::PIDIV4,qualityEval.minAngle(f)));

}
