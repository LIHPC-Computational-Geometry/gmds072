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
 * BezierCurve.h
 *
 *  Created on: 07/02/2015
 *      Author: F. Ledoux
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_BEZIERCURVE_H_
#define GMDS_MATH_BEZIERCURVE_H_
/*----------------------------------------------------------------------------*/
#include <cmath>
#include <vector>
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/CommonTypes.h>
#include <GMDS/Math/Point.h>
#include <GMDS/Math/Vector.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
  /*--------------------------------------------------------------------------*/
  namespace math {
    /*------------------------------------------------------------------------*/
    /** \class BezierCurve
     *  \brief Defines a Bezier curve in 3D. Underlying computations are based
     *         on the simple de Casteljau algorithm
     */
    class EXPORT_GMDS BezierCurve {
      
    
    public:
      
      /*------------------------------------------------------------------------*/
      /** \brief Constructor of a quadratic bezier curve from 3 control points.
       * 
       * \param AP1 first control point
       * \param AP2 second control point
       * \param AP3 third control point
       */
      BezierCurve(const Point& AP1, const Point& AP2, const Point& AP3);

      /*------------------------------------------------------------------------*/
      /** \brief Constructor of a cubic bezier curve from 4 control points.
       * 
       * \param AP1 first control point
       * \param AP2 second control point
       * \param AP3 third control point
       * \param AP4 fourth control point
       */
      BezierCurve(const Point& AP1, const Point& AP2, const Point& AP3,
		  const Point& AP4);


      /*------------------------------------------------------------------------*/
      /** \brief Constructor of a cubic bezier curve from 2 end points and two
       *         tangent vectors
       * 
       * \param AP1 first end point
       * \param AV1 first derivative vector at AP1
       * \param AP2 second end point
       * \param AV2 second derivative vector at AP1
       */
      BezierCurve(const Point& AP1, const Vector& AV1, 
		  const Point& AP2, const Vector& AV2);


      /*------------------------------------------------------------------------*/
      /** \brief Constructor of a bezier curve from an ordered set of control 
       *         points.
       * 
       * \param APts the set of control points to define the curve
       */
      BezierCurve(const std::vector<Point>& APts);

      /*------------------------------------------------------------------------*/
      /** \brief Returns the point located on the parametric curve wit
       * 
       * \param AT the parameter in [0..1]
       *
       * \return The point located at (*this)(AT)
       */
      Point operator()(const double& AT) const;

      /*------------------------------------------------------------------------*/
      /** \brief Returns a set of point that discretize the curve in ANb segment
       *         in the parametric space
       * 
       * \param ANb the number of segment we want to have
       *
       * \return The set of point discretizing (*this);
       */
      std::vector<Point> getDiscretization(const int ANb) const;
      
    private:
      /** ordered list of control points of the curve */
      std::vector<Point> m_control_points;
      
    };
    /*----------------------------------------------------------------------------*/
  }
  /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_BEZIERCURVE_H_ */
/*----------------------------------------------------------------------------*/
