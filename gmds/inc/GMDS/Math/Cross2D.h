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
 * Cross2D.h
 *
 *  Created on: 03/19/2015
 *      Author: ledouxf
 */
/*-----------------------------------------------------------------*/
#ifndef GMDS_MATH_CROSS_2D_H_
#define GMDS_MATH_CROSS_2D_H_
/*-----------------------------------------------------------------*/
#include "GMDS/Math/Vector.h"
#include "GMDS/Math/Quaternion.h"
#include "GMDS/Math/Chart.h"
#include <cmath>
/*-----------------------------------------------------------------*/
#include <iostream>
/*-----------------------------------------------------------------*/
using namespace std;
/*-----------------------------------------------------------------*/
namespace gmds {
  /*-----------------------------------------------------------------*/
  namespace math {
    /*-----------------------------------------------------------------*/
    class EXPORT_GMDS Cross2D
    {

    public:

      static int index(const Cross2D& AC1,
		       const Cross2D& AC2,
		       const Cross2D& AC3);
      /*------------------------------------------------------------*/
      /* Constructor
       *  Defines a cross aligned with (1,0,0)
       */
      Cross2D();
      /*------------------------------------------------------------*/
      /* Constructor
       *  Defines a cross from two orthogonal vectors
       */
      Cross2D(Vector& AV1, Vector& AV2);
      /*------------------------------------------------------------*/
      /* Constructor
       * Defines from a reference vector
       */
      Cross2D(const Vector& ARefVec);
      /*------------------------------------------------------------*/
      /* Constructor
       * Defines from a reference angle
       */
      Cross2D(const TCoord& ARefAngle);

      /*------------------------------------------------------------*/
      /* Copy  Constructor
       */
      Cross2D(const Cross2D& AC);

      /*------------------------------------------------------------*/
      /* Return the vector of this, which is the closest of AN
       */
      Vector closestComponentVector(const Vector& AN) const;

      /*------------------------------------------------------------*/
      /* Return the vector representation of this with 4 orthogonal
       * vectors
       */
      std::vector<Vector> componentVectors() const;

      /*------------------------------------------------------------*/
      /* Compute and stores the cross vectors
       */
      void computeComponentVectors();

      /*------------------------------------------------------------*/
      /* Indicates if cross vectors are stored or not
       */
      bool storeComponentVectors() const;

      /*------------------------------------------------------------*/
      /* Return the reference vector
       */
      Vector referenceVector() const;

      /*------------------------------------------------------------*/
      /* Return the reference angle
       */
      TCoord referenceAngle() const {return m_angle;}; 
      /*------------------------------------------------------------*/
      /* Return the reference angle
       */
      TCoord orientedReferenceAngle() const 
      {
	if(m_angle>Constants::PI)
	  return Constants::PI2-m_angle;
	else
	  return m_angle;
      }; 

      /*------------------------------------------------------------*/
      /* Return the angle between this and AC
       */
      TCoord angle(const Cross2D& AC) const;

      /*------------------------------------------------------------------------*/
      /** \brief  Overloaded operator+ to create a new point from 2 points
       */
      friend EXPORT_GMDS Cross2D operator+(const Cross2D&, const Cross2D&);

      /*------------------------------------------------------------------------*/
      /** \brief A pondered mean
       * This methods implement a k-pass direct mean pondered by the weights
       *
       * \param ACrosses the crosses we want to compute the mean
       * \param Aweights weights associated to each cross in the mean computation
       * \param ANbSteps the max number of step to converge (default = 5)
       */
      static EXPORT_GMDS Cross2D mean(const vector<Cross2D> & ACrosses, 
				      const vector<TCoord> & AWeights,
				      const TInt ANbSteps=5);
      static EXPORT_GMDS Cross2D mean(const Cross2D & AC1, 
				      const TCoord& AW1,
				      const Cross2D & AC2, 
				      const TCoord& AW2);
    private:
      // A 2D cross is only represented by an angle, which is the
      // angle between its representation vector and the (OX) axis
      TCoord m_angle;

      // the four vector that defines the 2D cross can be kept in
      // mind or computed on the fly. By default they are not 
      // stored
      bool m_known_vectors;
      // storage container of cross vectors
      std::vector<Vector> m_vectors;
    };


    EXPORT_GMDS ostream & operator << (ostream & AStr, 
				       const Cross2D & AC);
  }//namespace math
  /*-----------------------------------------------------------------*/
} //namespace gmds
/*-----------------------------------------------------------------*/
#endif /* GMDS_MATH_CROSS_2D_H_ */
/*-----------------------------------------------------------------*/
