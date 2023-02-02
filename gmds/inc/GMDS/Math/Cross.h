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
/*-----------------------------------------------------------------*/
/*
 * Cross.h
 *
 *  Created on: 03/01/2015
 */
/*-----------------------------------------------------------------*/
#ifndef GMDS_MATH_CROSS_H_
#define GMDS_MATH_CROSS_H_
/*-----------------------------------------------------------------*/
#include "GMDS/Math/Vector.h"
#include "GMDS/Math/Quaternion.h"
#include "GMDS/Math/Chart.h"
#include <cmath>
/*-----------------------------------------------------------------*/
#include <iostream>
/*-----------------------------------------------------------------*/
using namespace std;
namespace gmds {
  namespace math {
    /*-----------------------------------------------------------------*/
    class EXPORT_GMDS Cross
    {
    private:
      /// A cross is represented by two orthogonal vectors. 
      /* Using only one representation vector is only possible if we
       * can define a reference vector int the whole world (not done
       * in 3D, try to use transport equation along surface?)
       */
      math::Vector m_x;
      math::Vector m_y;

    public:

      /*------------------------------------------------------------*/
      /* Constructor
       *  Defines a cross with (1,0,0) and (0,1,0)
       */
      Cross();
      /*------------------------------------------------------------*/
      /* Constructor
       *  Defines a cross from two orthogonal vectors.
       */
      Cross(math::Vector& AV1, math::Vector& AV2);
      /*------------------------------------------------------------*/
      /* Constructor
       * Defines a cross from a quaternion and a normal vector to be
       * aligned with
       */
      Cross(const Quaternion& AQ, const math::Vector& AN);

      /*------------------------------------------------------------*/
      /* returns the right-hand chart (X,Y, XxY)
       */
      math::Chart chart() const;

      /*------------------------------------------------------------*/
      /*
       *  return the vector, which is the closest of AN between m_x,
       * m_y, -m_x, -m_y.
       */
      math::Vector closestVector(const math::Vector& AN);

      /*------------------------------------------------------------*/
      /*
       *  return the vector x
       */
      math::Vector X() const { return m_x; }
      /*------------------------------------------------------------*/
      /*
       *  return the vector y
       */
      math::Vector Y() const { return m_y; }
    };


    EXPORT_GMDS ostream & operator << (ostream & op_g, const Cross & op_d);
  }
  /*-----------------------------------------------------------------*/
}
/*-----------------------------------------------------------------*/
#endif /* GMDS_MATH_CROSS_H_ */
/*-----------------------------------------------------------------*/
