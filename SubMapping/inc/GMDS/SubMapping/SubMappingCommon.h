/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * This software is a computer program whose purpose is to provide a set of
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
 * data to be ensured and, more generally, to use and operate it in the 
 * same conditions as regards security. 
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    SubMappingCommon.h
 *  \author  legoff
 *  \date    01/02/2016
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_SUBMAPPING_SUBMAPPINGCOMMON_H_
#define GMDS_SUBMAPPING_SUBMAPPINGCOMMON_H_
/*----------------------------------------------------------------------------*/
// stl File Headers
/*----------------------------------------------------------------------------*/
#include <cmath>
/*----------------------------------------------------------------------------*/
// External File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
// CaGe File Headers
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace submapping{
/*----------------------------------------------------------------------------*/
const double PI = 4*atan(1.);
/*----------------------------------------------------------------------------*/
/** Values for ideal angle classification*/
/*----------------------------------------------------------------------------*/
const double END_VALUE = PI/2;
const double SIDE_VALUE = PI;
const double CORNER_VALUE = 3*PI/2;
const double REVERSAL_VALUE = 2*PI;
/*----------------------------------------------------------------------------*/
/** Values of angle classification*/
/*----------------------------------------------------------------------------*/
const int END = 1;
const int SIDE = 0;
const int CORNER = -1;
const int REVERSAL = -2;
/*----------------------------------------------------------------------------*/
/** Print vertices classification*/
/*----------------------------------------------------------------------------*/
std::map<int, std::string> verticesClassificationStr = {
	std::pair<int,std::string>(END,std::string("END")),
	std::pair<int,std::string>(SIDE,std::string("SIDE")),
	std::pair<int,std::string>(CORNER,std::string("CORNER")),
	std::pair<int,std::string>(REVERSAL,std::string("REVERSAL"))
};
/*----------------------------------------------------------------------------*/
/** Values for curves classification*/
/*----------------------------------------------------------------------------*/
typedef enum {
   PLUS_I,
   PLUS_J,
   MINUS_I,
   MINUS_J
 } CurveClassification;
/*----------------------------------------------------------------------------*/
/** Print curves classification*/
/*----------------------------------------------------------------------------*/
std::map<int, std::string> curvesClassificationStr = {
        std::pair<int,std::string>(PLUS_I,std::string("i+")),
        std::pair<int,std::string>(PLUS_J,std::string("j+")),
        std::pair<int,std::string>(MINUS_I,std::string("i-")),
        std::pair<int,std::string>(MINUS_J,std::string("j-"))
};
/*----------------------------------------------------------------------------*/
} // end namespace submapping
/*----------------------------------------------------------------------------*/
} // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_SUBMAPPING_SUBMAPPINGCOMMON_H_ */
