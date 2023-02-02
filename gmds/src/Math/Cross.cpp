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
 * cross.cpp
 *
 *  Created on: Sep 05, 2014
 *      Author: franck Ledoux
 */
/*-----------------------------------------------------------------*/
#include "GMDS/Math/Cross.h"
#include "GMDS/Math/Quaternion.h"
#include "GMDS/Math/Chart.h"
#include "GMDS/Math/Vector.h"
/*-----------------------------------------------------------------*/
#include <iostream>
#include <math.h>
namespace gmds {
  /*-----------------------------------------------------------------*/
  namespace math{
    /*-----------------------------------------------------------------*/
    Cross::Cross()
    {
      m_x = math::Vector(1,0,0);
      m_y = math::Vector(0,1,0);
    }

    /*-----------------------------------------------------------------*/
    Cross::Cross(math::Vector& AV1, math::Vector& AV2)
    {
      if (AV1.dot(AV2) != 0.0)
	throw GMDSException("A CROSS object can only be built from 2 orthogonal vectors");

      m_x = AV1;
      m_y = AV2;
      m_x.normalize();
      m_y.normalize();

    }

    /*-----------------------------------------------------------------*/
    Cross::
    Cross(const  Quaternion& AQ, const math::Vector& AN)
    {
      Quaternion q = AQ;
      // we rotate q to be aligned with AN
      q.alignWith(AN);
      Chart tq(q);
      if (fabs(AN.dot(tq.X())) >0.5)
	{
	  m_x = tq.Y();
	  m_y = tq.Z();
	}
      else if (fabs(AN.dot(tq.Y())) >0.5)
	{
	  m_x = tq.X();
	  m_y = tq.Z();
	}
      else
	{
	  m_x = tq.X();
	  m_y = tq.Y();
	}
      m_x.normalize();
      m_y.normalize();
      //std::cout << "normal: " << AN << std::endl;
      //std::cout << " -> Cross X: " << m_x << std::endl;
      //std::cout << " -> Cross Y: " << m_y << std::endl;
    }
    /*-----------------------------------------------------------------*/
    math::Vector Cross::
    closestVector(const math::Vector& AN)
    {
      math::Vector v = AN;
      v.normalize();
      math::Vector x_opp(-m_x.X(), -m_x.Y(), -m_x.Z());
      math::Vector y_opp(-m_y.X(), -m_y.Y(), -m_y.Z());

      math::Vector result = m_x;
      double valX = v.dot(m_x);
      double valXOpp = v.dot(x_opp);
      double valY = v.dot(m_y);
      double valYOpp = v.dot(y_opp);

      double val = valX;
      if (valXOpp > val)
	{
	  val = valXOpp;
	  result = x_opp;
	}
      if (valY > val)
	{
	  val = valY;
	  result = m_y;
	}
      if (valYOpp > val)
	{
	  result = y_opp;
	}

      return result;
    }
    /*-----------------------------------------------------------------*/

    ostream & operator << (ostream & str, const Cross & c)
    {
      str << "Cross (" << c.X() << ", " << c.Y() << ")";
      return str;
    }
  }
  /*-----------------------------------------------------------------*/
}
/*-----------------------------------------------------------------*/
