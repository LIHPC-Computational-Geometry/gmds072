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
 * Constants.h
 *
 *  Created on: 6 févr. 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_CONSTANTS_H_
#define GMDS_MATH_CONSTANTS_H_
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/CommonTypes.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
  /*----------------------------------------------------------------------------*/
  namespace math{
    /*----------------------------------------------------------------------------*/
    namespace Constants{
      const TCoord EPSILON = (TCoord)1.e-8;
      const TCoord PI = (TCoord)3.1415926535897931;
      const TCoord PI2 = (TCoord)6.2831853071795862;
      const TCoord PIDIV2 = (TCoord)1.5707963267948966;
      const TCoord PIDIV4 = (TCoord)0.78539816339744828;
      const TCoord PIDIV6 = (TCoord)0.52359877559829882;
      const TCoord PIDIV8 = (TCoord)0.39269908169872414;
      const TCoord INVPIDIV180 = (TCoord)57.295779513082323;
      const TCoord PIDIV180 = (TCoord)0.017453292519943295;

    } // namespace Constants
      /*----------------------------------------------------------------------------*/
  } // namespace math
  /*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_CONSTANTS_H_ */
/*----------------------------------------------------------------------------*/
