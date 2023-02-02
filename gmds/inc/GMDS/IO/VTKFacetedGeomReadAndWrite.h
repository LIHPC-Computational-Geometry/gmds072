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
 * VTKFacetedGeomReadAndWrite.h
 *
 *  Created on: 29 août 2014
 *      Author: ledouxf
 */

/*----------------------------------------------------------------------------*/
#ifndef GMDS_VTKFACETEDGEOMREADANDWRITE_H_
#define GMDS_VTKFACETEDGEOMREADANDWRITE_H_
/*----------------------------------------------------------------------------*/
#include <string>
/*----------------------------------------------------------------------------*/
#include "GMDS/CAD/FacetedGeomManager.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
/** \class  VTKFacetedGeomReadAndWrite
 *  \brief  This class initialize a faceted geometric model from a VTK file and
 *  		export in VTK files.
 */
class VTKFacetedGeomReadAndWrite {
public:

	/*------------------------------------------------------------------------*/
	/** \brief  Default constructor
	 */
	VTKFacetedGeomReadAndWrite();

	/*------------------------------------------------------------------------*/
	/** \brief  Default destructor
	 */virtual ~VTKFacetedGeomReadAndWrite();


	 /*------------------------------------------------------------------------*/
	 /** \brief  Importation of a VTK File into a faceted geometric manager
	  *
	  *  \param AGeomMng the geometric manager that will get the geometric
	  *  				 faceted entities stored in the file "AFileName"
	  *
	  *  \param AFile   the name of the file to be imported
	  */
	 void import( geom::FacetedGeomManager& AGeomMng,
			 const std::string& 	   AFile);

	 /*------------------------------------------------------------------------*/
         /** \brief  Build a faceted geometric manager from a VTK file
          *
          *  \param AGeomMng the geometric manager that will be built from
	  *				from the mesh in the file				
          *
          *  \param AFile   the name of the file to be imported
	  *  \param ASingleSurface build the model with only one surface
          */
         void build( geom::FacetedGeomManager& AGeomMng,
                         const std::string&        AFile,
			 const bool ASingleSurface=false);

	 /*------------------------------------------------------------------------*/
	 /** \brief  export the model stored in AGeomMng into a VTK file for
	  * 		 vizualisation.
	  *
	  *  \param AGeomMng the geometric manager that will get the geometric
	  *  				 faceted entities stored in the file "AFileName"
	  *
	  *
	  *  \param AFile the output file
	  */
	 void exportVTK( geom::FacetedGeomManager& AGeomMng,
			 const std::string& AFile = "faceted_model");


};
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_VTKFACETEDGEOMREADANDWRITE_H_ */
/*----------------------------------------------------------------------------*/
