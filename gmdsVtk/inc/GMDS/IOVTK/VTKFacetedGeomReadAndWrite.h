/*----------------------------------------------------------------------------*/
/*
 * VTKFacetedGeomReadAndWrite.h
 *
 *  Created on: 29 août 2014
 *      Author: ledouxf
 */

/*----------------------------------------------------------------------------*/
#ifndef VTKFACETEDGEOMREADANDWRITE_H_
#define VTKFACETEDGEOMREADANDWRITE_H_
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
#endif /* VTKFACETEDGEOMREADANDWRITE_H_ */
/*----------------------------------------------------------------------------*/
