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
 * FACFacetedGeomReadAndWrite.h
 *
 *  Created on: 16 oct 2014
 *      Author: legoff
 */

/*----------------------------------------------------------------------------*/
#ifndef GMDS_GEOM_FACFACETEDGEOMREADANDWRITE_H_
#define GMDS_GEOM_FACFACETEDGEOMREADANDWRITE_H_
/*----------------------------------------------------------------------------*/
#include <string>
/*----------------------------------------------------------------------------*/
#include "GMDS/CAD/FacetedGeomManager.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace geom{
/*----------------------------------------------------------------------------*/
/** \class  FACFacetedGeomReadAndWrite
 *  \brief  This class initialize a faceted geometric model from a FAC file and
 *  		export in FAC files.
 */
	class EXPORT_GMDS FACFacetedGeomReadAndWrite {
public:

	/*------------------------------------------------------------------------*/
	/** \brief  Default constructor
	 */
	FACFacetedGeomReadAndWrite();

	/*------------------------------------------------------------------------*/
	/** \brief  Default destructor
	 */virtual ~FACFacetedGeomReadAndWrite();


	 /*------------------------------------------------------------------------*/
	 /** \brief  Importation of a FAC File into a faceted geometric manager
	  *
	  *  \param AGeomMng the geometric manager that will get the geometric
	  *  				 faceted entities stored in the file "AFileName"
	  *
	  *  \param AFile   the name of the file to be imported
	  */
	 void import( geom::FacetedGeomManager& AGeomMng,
			 const std::string& 	   AFile);



	 /*------------------------------------------------------------------------*/
	 /** \brief  export the model stored in AGeomMng into a FAC file for
	  * 		 vizualisation.
	  *
	  *  \param AGeomMng the geometric manager that will get the geometric
	  *  				 faceted entities stored in the file "AFileName"
	  *
	  *
	  *  \param AFile the output file
	  */
	 void exportFAC( geom::FacetedGeomManager& AGeomMng,
			 const std::string& AFile = "faceted_model");


};
/*----------------------------------------------------------------------------*/
} // namespace geom
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_GEOM_FACFACETEDGEOMREADANDWRITE_H_ */
/*----------------------------------------------------------------------------*/
