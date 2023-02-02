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
/*
 * Ray.cpp
 *
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/Math/Ray.h>
/*----------------------------------------------------------------------------*/
#include <GMDS/Math/Segment.h>
#include <GMDS/Math/Line.h>
#include <GMDS/Math/Plane.h>
#include <GMDS/Math/Triangle.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
  /*--------------------------------------------------------------------------*/
  namespace math{
    /*------------------------------------------------------------------------*/
    const Point& Ray::getPoint() const
    {
      return m_pnt;
    }
    /*------------------------------------------------------------------------*/
    const Vector& Ray::getDir() const {
      return m_dir;
    }
    /*------------------------------------------------------------------------*/
    const Vector& Ray::getDirUnit() const {
      if(m_isDirUnit)
	return m_dirUnit;

      m_isDirUnit = true;
      m_dirUnit = m_dir.getNormalize();
      return m_dirUnit;
    }
    /*------------------------------------------------------------------------*/
    Ray& Ray::operator=(const Ray& ARay) {
      m_pnt = ARay.m_pnt;
        m_dir = ARay.m_dir;
        m_isDirUnit = ARay.m_isDirUnit;
        m_dirUnit = ARay.m_dirUnit;
        return *this;
    }
      /*----------------------------------------------------------------------------*/
      bool Ray::intersect2D(const Segment& AS, Point& AP, double& AParam) const {
          
          Point  p1= AS.getPoint(0);
          Point  p2= AS.getPoint(1);
          
          Point src_pnt = m_pnt;
          Point dir_pnt = m_pnt + m_dirUnit;
          
          Vector v12(p1,p2);
          if(v12.isColinear(m_dirUnit)){
              return false; // No intersection
          }
          
          double x1 = p1.X();
          double y1 = p1.Y();
          double x2 = p2.X();
          double y2 = p2.Y();
          double xs = src_pnt.X();
          double ys = src_pnt.Y();
          double xd = dir_pnt.X();
          double yd = dir_pnt.Y();
          
          double D = x1*(yd-ys) + x2*(ys-yd) + xd*(y2-y1) + xs* (y1-y2);
          double N = x1*(yd-ys) + xs*(y1-yd) + xd*(ys-y1);
          
          //param in [p1,p2]
          double s = N/D;
          
          N =  -(x1*(ys- y2) + x2*(y1-ys) + xs*(y2-y1));
          
          //param in [src_pnt,dir_pnt]
          double t = N/D;
          
          if ((0.0 <= s) && (s <= 1.0) && (0.0 < t)){
              AP = p1 + s*(p2-p1);
              AParam = 1-s;
              return true;
          }
          return false;
          
      }
      
      /*----------------------------------------------------------------------------*/
      bool Ray::intersect3D(const Segment& AS, Point& AP, double& AParamSeg,
                            double& AParamRay) const {
          // see details in http://www.lucidarme.me/?p=1872
          
          Point  p1= AS.getPoint(0);
          Point  p2= AS.getPoint(1);
          
          Point p3 = m_pnt;
          Point p4 = m_pnt + m_dirUnit;

          Vector v13(p1,p3);
          Vector v12(p1,p2);
          Vector v34(p3,p4);
          
          double tol = 0.01;
          
          if(v12.isColinear(m_dirUnit)){
//              std::cout<<"\t colinear"<<std::endl;
              return false; // No intersection
          }

          double coplanar = std::abs(v13.dot(v12.cross(v34)));

          if(coplanar<tol){//Coplanar
              Vector v12_c_v34 = v12.cross(v34);
              double t = (v13.cross(v34)).dot(v12_c_v34)/(v12_c_v34.dot(v12_c_v34));
              if(0<=t && t<=1){
                  AParamSeg = t;
                  AP = p1 +t*v12;
                  AParamRay = (v13.cross(v12)).dot(v12_c_v34)/(v12_c_v34.dot(v12_c_v34));

                  return true;
              }
          }
          else
              std::cout<<"\t not coplanar: "<<std::abs(v13.dot(v12.cross(v34)))<<std::endl;
          return false;
          
      }
      /*----------------------------------------------------------------------------*/
      bool Ray::intersect2D(const Ray& AR, Point& AP) const {
          
          Point  p1= AR.getPoint();
          Point  p2= p1 + AR.getDirUnit();
          
          Point src_pnt = m_pnt;
          Point dir_pnt = m_pnt + m_dirUnit;
          
      Vector v12(p1,p2);
      if(v12.isColinear(m_dirUnit)){
	return false; // No intersection
      }

      double x1 = p1.X();
      double y1 = p1.Y();
      double x2 = p2.X();
      double y2 = p2.Y();
      double xs = src_pnt.X();
      double ys = src_pnt.Y();
      double xd = dir_pnt.X();
      double yd = dir_pnt.Y();

      double D = x1*(yd-ys) + x2*(ys-yd) + xd*(y2-y1) + xs* (y1-y2);
      double N = x1*(yd-ys) + xs*(y1-yd) + xd*(ys-y1);

      //param in [p1,p2]
      double s = N/D;

      N =  -(x1*(ys- y2) + x2*(y1-ys) + xs*(y2-y1));

      //param in [src_pnt,dir_pnt]
      double t = N/D;

      if ((0.0 <= s)  && (0.0 <= t)){
	AP = p1 + s*(p2-p1);
	return true;
      }
      return false;

    }
/*----------------------------------------------------------------------------*/
bool Ray::intersect3D(const Plane& APlane, Point& AP) const
{
	// check whether the line intersects the plane
	Line ln(m_pnt,m_dirUnit);
	double param(0.);
	if(!ln.intersect3D(APlane,AP,param)) {
		return false;
	}
	
	// check whether the intersection point is on the ray or on the other side
	if(param < 0.) {
		return false;
	}
	
	return true;
}
/*----------------------------------------------------------------------------*/
bool Ray::intersect3D(const Triangle& ATri, Point& AP) const
{
        // check whether the line intersects the triangle
        Line ln(m_pnt,m_dirUnit);
        double param(0.);
        if(!ln.intersect3D(ATri,AP,param)) {
                return false;
        }

        // check whether the intersection point is on the ray or on the other side
        if(param < 0.) {
                return false;
        }

        return true;
}
/*----------------------------------------------------------------------------*/
Point
Ray::project(const Point& AP) const
{
	gmds::math::Line line(m_pnt,m_dirUnit);
	gmds::math::Point projectedPoint = line.project(AP);

	if(AP.distance2(projectedPoint) < AP.distance2(m_pnt)) {
		return projectedPoint;
	} else {
		return m_pnt;
	}
}
/*---------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& AStr, const Ray& ARay){
        AStr<<"("<<ARay.m_pnt<<" | "<<ARay.m_dirUnit<<")";
                return AStr;
}
    /*----------------------------------------------------------------------------*/
  } // namespace math
  /*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
