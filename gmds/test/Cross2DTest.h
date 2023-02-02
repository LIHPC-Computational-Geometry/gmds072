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
#include <GMDS/Math/Cross2D.h>
#include <GMDS/Math/Numerics.h>
/*----------------------------------------------------------------------------*/
class Cross2DTest: public ::testing::Test {
protected:
	Cross2DTest(){;}
	virtual ~Cross2DTest(){;}
};

using namespace gmds;

/*----------------------------------------------------------------------------*/
TEST_F(Cross2DTest,initVec1) {
  math::Vector v1(1,0,0);
  math::Vector v2(0,1,0);
  math::Cross2D c(v1,v2);

  EXPECT_EQ(0,c.referenceAngle());
}

/*----------------------------------------------------------------------------*/
TEST_F(Cross2DTest,initVec2) {
  math::Vector v1(0,1,0);
  math::Vector v2(-1,0,0);
  math::Cross2D c(v1,v2);

  EXPECT_EQ(0,c.referenceAngle());
}

/*----------------------------------------------------------------------------*/
TEST_F(Cross2DTest,initVec3) {
  math::Vector v1(-1,0,0);
  math::Vector v2(0,-1,0);
  math::Cross2D c(v1,v2);

  EXPECT_EQ(0,c.referenceAngle());
}

/*----------------------------------------------------------------------------*/
TEST_F(Cross2DTest,initVec4) {
  math::Vector v1(1,0,0);
  math::Vector v2(0,-1,0);
  math::Cross2D c(v1,v2);

  EXPECT_EQ(0,c.referenceAngle());
}

/*----------------------------------------------------------------------------*/
TEST_F(Cross2DTest,initVecFail) {
  math::Vector v1(1,0,0);
  math::Vector v2(0.1,1,0);
  EXPECT_THROW({
      math::Cross2D c(v1,v2);

    }, GMDSException);

}

/*----------------------------------------------------------------------------*/
TEST_F(Cross2DTest,initVectorsVsRef) {
  math::Vector v_ref(1,0,0);
  math::Cross2D c(v_ref);

  EXPECT_EQ(0,c.referenceAngle());

  math::Vector ref = c.referenceVector();
  
  EXPECT_EQ(1,ref.X());
  EXPECT_EQ(0,ref.Y());
  EXPECT_EQ(0,ref.Z());

  std::vector<math::Vector> c_vecs = c.componentVectors();

  EXPECT_EQ(4,c_vecs.size());

  for(unsigned int i=0; i<c_vecs.size();i++){
    math::Vector current_vec = c_vecs[i];
    TCoord current_angle = current_vec.angle(v_ref);
    EXPECT_NEAR(0.0,math::modulo(current_angle,math::Constants::PIDIV2), 1e-12);
  }
}

/*----------------------------------------------------------------------------*/
TEST_F(Cross2DTest,Cross2DIndex1) {
  math::Cross2D c1(0);
  math::Cross2D c2(0);
  math::Cross2D c3(0);
 
  int index = math::Cross2D::index(c1,c2,c3);
  EXPECT_EQ(0, index);
}

/*----------------------------------------------------------------------------*/
TEST_F(Cross2DTest,Cross2DIndex2) {
  math::Cross2D c1(15);
  math::Cross2D c2(15);
  math::Cross2D c3(15);
 
  int index = math::Cross2D::index(c1,c2,c3);
  EXPECT_EQ(0, index);
}

/*----------------------------------------------------------------------------*/
TEST_F(Cross2DTest,Cross2DIndex3) {
  math::Cross2D c1(math::Constants::PI);
  math::Cross2D c2(math::Constants::PI*1.1);
  math::Cross2D c3(math::Constants::PI);
 
  int index = math::Cross2D::index(c1,c2,c3);
  EXPECT_EQ(0, index);
}

/*----------------------------------------------------------------------------*/
TEST_F(Cross2DTest,Cross2DIndex4) {
  math::Cross2D c1(1);
  math::Cross2D c2(-2.1);
  math::Cross2D c3(-1.3);
 
  int index = math::Cross2D::index(c1,c2,c3);
  EXPECT_EQ(0, index);
}

/*----------------------------------------------------------------------------*/
TEST_F(Cross2DTest,Cross2DIndex5) {
  math::Cross2D c1(math::Constants::PIDIV2);
  math::Cross2D c2(0);
  math::Cross2D c3(math::Constants::PIDIV2*3);
 
  int index = math::Cross2D::index(c3,c2,c1);
  EXPECT_EQ(1, index);
}

/*----------------------------------------------------------------------------*/
TEST_F(Cross2DTest,Cross2DIndex6) {
  math::Cross2D c1(math::Constants::PIDIV2);
  math::Cross2D c2(0);
  math::Cross2D c3(math::Constants::PIDIV2*3);
 
  int index = math::Cross2D::index(c1,c2,c3);
  EXPECT_EQ(-1, index);
}
