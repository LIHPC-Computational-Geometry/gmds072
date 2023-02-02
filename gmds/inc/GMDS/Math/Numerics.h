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
 * Numerics.h
 *
 *  Created on:March 01, 2015
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_NUMERICS_H_
#define GMDS_MATH_NUMERICS_H_
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/CommonTypes.h>
#include <GMDS/Math/Constants.h>
#include <GMDS/Math/Point.h>
#include <GMDS/Math/Matrix.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
  /*----------------------------------------------------------------------------*/
  namespace math{
    
    /*------------------------------------------------------------------------*/
    /** \brief Compute AVal modulo AMod
     */
    TCoord modulo(const TCoord AVal, const TCoord AMod);
    /*------------------------------------------------------------------------*/
    /** \brief Compute AVal modulo 2xPI
     */
    TCoord modulo2PI(const TCoord AVal); 
    /*------------------------------------------------------------------------*/
    /** \brief Returns the min value between AV1, AV2 and AV3
     */
    TCoord min3(const TCoord AV1, const TCoord AV2, const TCoord AV3);
    /*------------------------------------------------------------------------*/
    /** \brief Returns the max value between AV1, AV2 and AV3
     */
    TCoord max3(const TCoord AV1, const TCoord AV2, const TCoord AV3);
    /*------------------------------------------------------------------------*/
    /** \brief Returns the min value between AV1 and AV2
     */
    TCoord min2(const TCoord AV1, const TCoord AV2);
    /*------------------------------------------------------------------------*/
    /** \brief Returns the max value between AV1 and AV2
     */
    TCoord max2(const TCoord AV1, const TCoord AV2);
    /*------------------------------------------------------------------------*/
    /** \brief Returns true if AV1==AV2 [AEpsilon]
     */
    bool near(const TCoord AV1, TCoord AV2, TCoord AEpsilon = Constants::EPSILON);

      /*------------------------------------------------------------------------*/
      /** \brief Returns the stiffness matrix for the triangle elment defined by
       *         AP1, AP2 and AP3.
       *
       *    WARNING ONLY 2D!!!!! (x,y)
       *
       * \param AP1 a first point
       * \param AP2 a second point
       * \param AP3 a third point
       *
       * \return A 3x3 the stiffness matrix of (AP1,AP2,AP3)
       */
      Matrix<3,3,double> stiffnessMatrix2D(const Point& AP1,
                                           const Point& AP2,
                                           const Point& AP3);
      
      /*------------------------------------------------------------------------*/
      /** \brief Compute an average plane going through the set of points \p AP
       *
       * \param[in]  AP           the set of points the plane must fit to
       * \param[out] APlanePnt    one point of the plane
       * \param[out] APlaneNormal normal vector to the plane
       *
       */
      void computeLeastSquarePlane(const std::vector<Point>& AP,
                                   Point& APlanePnt,
                                   math::Vector3d& APlaneNormal);
      
      
      /*------------------------------------------------------------------------*/
      /** \brief Solve a 2nd degree polynom ax2+bx+c=0
       *
       * \param[in]  AA coeff of X2
       * \param[in]  AB coeff of X
       * \param[in]  AC constant coeff
       * \param[out] AX solutions
       *
       * \return the number of real solutions (0, 1 or 2)
       */
      int solve2ndDegreePolynomial( const double& AA,
                                   const double& AB,
                                   const double& AC,
                                   std::vector<double>& AX);
      
      /*------------------------------------------------------------------------*/
      /** \brief Considering tetrahedron T defined by points \p AP0, \p AP1,
       *         \p AP2 and \p AP3, this function returs the cotangent weight
       *         for edge[\p AP0, \p AP1] in T.
       *
       * \param[in]  AP0 a first tet point
       * \param[in]  AP1 a second tet point
       * \param[in]  AP2 a third tet point
       * \param[in]  AP3 a fourth tet point
       *
       * \return the cotangent weight of edge [\p AP0, \p AP1] in this tet.
       */
      
      double cotangentWeight(const Point& AP1,
                             const Point& AP2,
                             const Point& AP3,
                             const Point& AP4);
      
      /*------------------------------------------------------------------------*/
      /** \brief Considering tetrahedron T defined by points \p AP0, \p AP1,
       *         \p AP2 and \p AP3, this function returs the dihedral angle
       *         for edge[\p AP0, \p AP1] in T.
       *
       * \param[in]  AP0 a first tet point
       * \param[in]  AP1 a second tet point
       * \param[in]  AP2 a third tet point
       * \param[in]  AP3 a fourth tet point
       *
       * \return the dihedral angle of edge [\p AP0, \p AP1] in this tet.
       */
      
      double dihedralAngle(const Point& AP1,
                             const Point& AP2,
                             const Point& AP3,
                             const Point& AP4);
      
      
      
      /*----------------------------------------------------------------------------*/
  } // namespace math
  /*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_CONSTANTS_H_ */
/*----------------------------------------------------------------------------*/
