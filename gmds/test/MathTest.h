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
#include <GMDS/Math/Hexahedron.h>
#include <GMDS/Math/Point.h>
#include <GMDS/Math/Segment.h>
#include <GMDS/Math/Ray.h>
#include <GMDS/Math/Triangle.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
class MathTest: public ::testing::Test {

  protected:
	MathTest(){;}
    virtual ~MathTest(){;}
};
/*----------------------------------------------------------------------------*/
TEST_F(MathTest,hexahedronTriangleIntersection) {
	
	math::Point pnts[8];

	pnts[0] = math::Point(1, 2, 4.3);
	pnts[1] = math::Point(1.1, 2, 4.3);
	pnts[2] = math::Point(1.1, 2.25, 4.3);
	pnts[3] = math::Point(1, 2.25, 4.3);
	pnts[4] = math::Point(1, 2, 4.6);
	pnts[5] = math::Point(1.1, 2, 4.6);
	pnts[6] = math::Point(1.1, 2.25, 4.6);
	pnts[7] = math::Point(1, 2.25, 4.6);
	
	math::Hexahedron hex(pnts[0],pnts[1],pnts[2],pnts[3],
			pnts[4],pnts[5],pnts[6],pnts[7]);

	math::Point pntsbis[3];
	pntsbis[0] = math::Point(1, 1.73205, 4);
	pntsbis[1] = math::Point(1, 1.73205, 4.5);
	pntsbis[2] = math::Point(0.72555, 1.86375, 4.16573);

	math::Triangle triangle(pntsbis[0],pntsbis[1],pntsbis[2]);

	bool res = hex.intersect(triangle);

	EXPECT_EQ(false,res);
}
/*----------------------------------------------------------------------------*/
TEST_F(MathTest,TriangleTriangleIntersection) {
	math::Point pnts[3];
	pnts[0] = math::Point(1, 2.25, 4.3);
	pnts[1] = math::Point(1, 2, 4.3);
	pnts[2] = math::Point(1, 2, 4.6);

	math::Triangle triangle(pnts[0],pnts[1],pnts[2]);

	math::Point pntsbis[3];
        pntsbis[0] = math::Point(1, 1.73205, 4);
        pntsbis[1] = math::Point(1, 1.73205, 4.5);
        pntsbis[2] = math::Point(0.72555, 1.86375, 4.16573);

        math::Triangle triangle2(pntsbis[0],pntsbis[1],pntsbis[2]);

	bool res = triangle.intersect(triangle2);

        EXPECT_EQ(false,res);
}
/*----------------------------------------------------------------------------*/
TEST_F(MathTest,TriangleTriangleIntersection_2) {
        math::Point pnt0(-0.492124, -0.712611, -0.5);
        math::Point pnt1(-0.620854, -0.603772, -0.5);
        math::Point pnt2(-0.596846, -0.627515, -0.5);
        math::Triangle triangle(pnt0,pnt1,pnt2);

	math::Point pnt3(-0.7, -0.8, -0.5);
        math::Point pnt4(-0.6, -0.8, -0.5);
        math::Point pnt5(-0.6, -0.7, -0.5);
        math::Triangle triangle2(pnt3,pnt4,pnt5);

        bool res = triangle.intersect(triangle2);

        EXPECT_EQ(false,res);
}
/*----------------------------------------------------------------------------*/
TEST_F(MathTest,TriangleTriangleIntersection_3) {
	math::Point pnt0(0, 0.5, 0.5);
        math::Point pnt1(0, 0, 0);
        math::Point pnt2(0, 1, 0);
        math::Triangle triangle(pnt0,pnt1,pnt2);

        math::Point pnt3(0, 1.1, 0.1);
        math::Point pnt4(0, 1, 0.1);
        math::Point pnt5(0, 1, 0.2);
        math::Triangle triangle2(pnt3,pnt4,pnt5);

        bool res = triangle.intersect(triangle2);

        EXPECT_EQ(false,res);
}
/*----------------------------------------------------------------------------*/
TEST_F(MathTest,TriangleSegmentIntersection) {
        math::Point pnts[3];
        pnts[0] = math::Point(1, 2.25, 4.3);
        pnts[1] = math::Point(1, 2, 4.3);
        pnts[2] = math::Point(1, 2, 4.6);

        math::Triangle triangle(pnts[0],pnts[1],pnts[2]);

        math::Point pntsbis[3];
        pntsbis[0] = math::Point(1, 1.73205, 4);
        pntsbis[1] = math::Point(1, 1.73205, 4.5);
//        pntsbis[2] = math::Point(0.72555, 1.86375, 4.16573);

        math::Segment segment(pntsbis[0],pntsbis[1]);

        bool res = triangle.intersect(segment);

        EXPECT_EQ(false,res);
}
/*----------------------------------------------------------------------------*/
TEST_F(MathTest,TriangleSegmentIntersection2D) {
        math::Point pnt3(-0.7, -0.8, 0.);
        math::Point pnt4(-0.6, -0.8, 0.);
        math::Point pnt5(-0.6, -0.7, 0.);
        math::Triangle triangle2(pnt3,pnt4,pnt5);

	math::Point pnt0(-0.492124, -0.712611, 0.);
	math::Point pnt1(-0.620854, -0.603772, 0.);
        math::Segment segment(pnt0,pnt1);

        bool res = triangle2.intersect2D(segment);

        EXPECT_EQ(false,res);
}
/*----------------------------------------------------------------------------*/
TEST_F(MathTest,TriangleSegmentIntersection2D_2) {
        math::Point pnt3(1.1, 0.1, 0);
        math::Point pnt4(1., 0.1, 0);
        math::Point pnt5(1., 0.2, 0);
        math::Triangle triangle2(pnt3,pnt4,pnt5);

        math::Point pnt0(1, 0, 0.);
        math::Point pnt1(0.5, 0.5, 0.);
        math::Segment segment(pnt0,pnt1);

        bool res = triangle2.intersect2D(segment);

        EXPECT_EQ(false,res);
}
/*----------------------------------------------------------------------------*/
TEST_F(MathTest,TriangleRayIntersection) {
	math::Point pnt0(0, -0.4, 0);
	math::Point pnt1(0.47629, -0.479422, 0);
	math::Point pnt2(0.312907, -0.766642, 0);
	math::Triangle triangle(pnt0,pnt1,pnt2);

	math::Point pnt(0.263162, -0.546133, 5);
	math::Vector dir(0, 0, -0.15228);
	math::Ray ray(pnt,dir);	

	bool res = triangle.intersect(ray);
	
	EXPECT_EQ(true,res);
}
/*----------------------------------------------------------------------------*/
TEST_F(MathTest,HexahedronTriangleIntersection_2) {
        math::Point pnt0(-0.492124, -0.712611, -0.5);
        math::Point pnt1(-0.620854, -0.603772, -0.5);
        math::Point pnt2(-0.596846, -0.627515, -0.5);
        math::Triangle triangle(pnt0,pnt1,pnt2);

        math::Point pnta(-0.7, -0.8, -0.6);
        math::Point pntb(-0.6, -0.8, -0.6);
        math::Point pntc(-0.6, -0.7, -0.6);
        math::Point pntd(-0.7, -0.7, -0.6);
        math::Point pnte(-0.7, -0.8, -0.5);
        math::Point pntf(-0.6, -0.8, -0.5);
        math::Point pntg(-0.6, -0.7, -0.5);
        math::Point pnth(-0.7, -0.7, -0.5);

        math::Hexahedron hex(pnta,pntb,pntc,pntd,pnte,pntf,pntg,pnth);

        bool res = hex.intersect(triangle);

        EXPECT_EQ(false,res);
}
/*----------------------------------------------------------------------------*/
TEST_F(MathTest,HexahedronTriangleIntersection_3) {
        math::Point pnt0(0, 0.5, 0.5);
        math::Point pnt1(0, 0, 0);
        math::Point pnt2(0, 1, 0);
        math::Triangle triangle(pnt0,pnt1,pnt2);

        math::Point pnta(0, 1, 0.1);
        math::Point pntb(0.1, 1, 0.1);
        math::Point pntc(0.1, 1.1, 0.1);
        math::Point pntd(0, 1.1, 0.1);
        math::Point pnte(0, 1, 0.2);
        math::Point pntf(0.1, 1, 0.2);
        math::Point pntg(0.1, 1.1, 0.2);
        math::Point pnth(0, 1.1, 0.2);

        math::Hexahedron hex(pnta,pntb,pntc,pntd,pnte,pntf,pntg,pnth);

        bool res = hex.intersect(triangle);

        EXPECT_EQ(false,res);
}
/*----------------------------------------------------------------------------*/
TEST_F(MathTest,HexahedronScaledJacobian) {
        math::Point pnta(0, 0, 0);
        math::Point pntb(1, 0, 0);
        math::Point pntc(1, 1, 0);
        math::Point pntd(0, 1, 0);
        math::Point pnte(0, 0, 1);
        math::Point pntf(1, 0, 1);
        math::Point pntg(1, 1, 1);
        math::Point pnth(0, 1, 1);

        math::Hexahedron hex(pnta,pntb,pntc,pntd,pnte,pntf,pntg,pnth);

        double res = hex.computeScaledJacobian();

        EXPECT_EQ(1,res);
}
/*----------------------------------------------------------------------------*/
TEST_F(MathTest,VectorAngle) {
  math::Vector v1(1,0,0);
  math::Vector v2(0,-1,0);
 
  double angle1 = v1.angle(v2);
  double angle2 = v2.angle(v1);
   
  EXPECT_NEAR(angle1, angle2, 1e-12);
}
/*----------------------------------------------------------------------------*/
TEST_F(MathTest,VectorOrientedAngle) {
  math::Vector v1(1,0,0);
  math::Vector v2(0,-1,0);
 
  double angle1 = v1.orientedAngle(v2);
  double angle2 = v2.orientedAngle(v1);
   
  EXPECT_TRUE(angle1<0);
  EXPECT_TRUE(angle2>0);
  EXPECT_NEAR(angle1, -angle2, 1e-12);
}
/*----------------------------------------------------------------------------*/
TEST_F(MathTest,VectorSelfAngle) {
  math::Vector v(1,1,0);
 
  double angle1 = v.angle(v);
  double angle2 = v.orientedAngle(v);
   
  EXPECT_NEAR(0.0, angle1, 1e-6);
  EXPECT_NEAR(0.0, angle2, 1e-6);
}
/*----------------------------------------------------------------------------*/
TEST_F(MathTest,VectorOrderingAngles) {
  math::Vector v(1,0,0);
  math::Vector v1(1,1,0);
  math::Vector v2(1,-1,0);
 
  double angle1   = v.angle(v1);
  double angle2   = v.angle(v2);
  double angle1_o = v.orientedAngle(v1);
  double angle2_o = v.orientedAngle(v2);

  EXPECT_NEAR(angle2, angle1, 1e-6);

  EXPECT_TRUE(angle1_o>angle2_o);

  //  EXPECT_TRUE(angle2>(3*math::Constants:PIDIV3));

}
/*----------------------------------------------------------------------------*/
