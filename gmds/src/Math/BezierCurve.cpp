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
 * BezieCurve.cpp
 *
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/Math/BezierCurve.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
  /*--------------------------------------------------------------------------*/
  namespace math{
    /*------------------------------------------------------------------------*/
    BezierCurve::BezierCurve(const Point& AP1, const Point& AP2, 
			     const Point& AP3) {
      m_control_points.resize(3);
      m_control_points[0] = AP1;
      m_control_points[1] = AP2;
      m_control_points[2] = AP3;
    }
    /*------------------------------------------------------------------------*/
    BezierCurve::BezierCurve(const Point& AP1, const Point& AP2, 
			     const Point& AP3, const Point& AP4) {
      m_control_points.resize(4);
      m_control_points[0] = AP1;
      m_control_points[1] = AP2;
      m_control_points[2] = AP3;
      m_control_points[3] = AP4;
    }
    /*------------------------------------------------------------------------*/
    BezierCurve:: BezierCurve(const Point& AP1, const Vector& AV1, 
			      const Point& AP2, const Vector& AV2){ 
      m_control_points.resize(4);
      m_control_points[0] = AP1;
      m_control_points[1] = AP1 + (1.0/3.0)*AV1;
      m_control_points[2] = AP2 - (1.0/3.0)*AV2;
      m_control_points[3] = AP2;
    }
    
    /*------------------------------------------------------------------------*/
    BezierCurve::BezierCurve(const std::vector<Point>& APts) {
      m_control_points = APts;
    }
    /*------------------------------------------------------------------------*/
    Point BezierCurve::operator()(const double& AT) const {
      if(AT>1.0 || AT<0.0)
	throw GMDSException("BezierCurve Query out of the range [0,1]");
      
     
      std::vector<Point> current_control = m_control_points;
      while(current_control.size()>1){

	std::vector<Point> next_control;
	const int deg = current_control.size();

	next_control.resize(deg-1);
	for(int i=0; i<deg-1; i++){
	  next_control[i] = (1-AT)*current_control[i]+ AT*current_control[i+1];
	}
       
	current_control = next_control;
      }
      
      return current_control[0];
    }
     /*------------------------------------------------------------------------*/
    std::vector<Point> BezierCurve:: getDiscretization(const int ANb) const {
      if(ANb<1)
	throw GMDSException("BezierCurve discretization impossible with this parameter");
            
     
      std::vector<Point> points;
      points.resize(ANb+1);
      double step = 1.0/ANb;
      double val =0;
      for(int i=0; i<=ANb; i++){
	points[i] = this->operator()(val);
	val+=step;
      }
      return points;
    }
  
    /*----------------------------------------------------------------------------*/
  } // namespace math
  /*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
