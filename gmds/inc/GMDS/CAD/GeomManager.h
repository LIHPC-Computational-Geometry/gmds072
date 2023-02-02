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
/** \file    GeomManager.h
 *  \author  F. LEDOUX
 *  \date    30/06/11
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_GEOM_GEOMMANAGER_H_
#define GMDS_GEOM_GEOMMANAGER_H_
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include "GMDS/Utils/Exception.h"
#include "GMDS/Utils/CommonTypes.h"
#include "GMDS/CAD/GeomPoint.h"
#include "GMDS/CAD/GeomCurve.h"
#include "GMDS/CAD/GeomSurface.h"
#include "GMDS/CAD/GeomVolume.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace geom{
/*----------------------------------------------------------------------------*/
/** \class GeomManager
 *  \brief This interface gathers the factory methods required for every
 *  		geometric model, the access to all the geom entities stored in the
 *  		model as the responsability to delete geometric entities.
 *
 *  \param TBase the basic type used to store geometric coordinates.
 */
/*----------------------------------------------------------------------------*/
	class EXPORT_GMDS GeomManager {

public:

	/*------------------------------------------------------------------------*/
	/** \brief  creation of a geometric volume
	 */
	virtual GeomVolume* newVolume() =0;

	/*------------------------------------------------------------------------*/
	/** \brief  creation of a geometric surface
	 */
	virtual GeomSurface* newSurface() =0;

	/*------------------------------------------------------------------------*/
	/** \brief  creation of a geometric curve
	 */
	virtual GeomCurve* newCurve() =0;

	/*------------------------------------------------------------------------*/
	/** \brief  creation of a geometric point
	 */
	virtual GeomPoint* newPoint() =0;

	/*------------------------------------------------------------------------*/
	/** \brief  Get the number of points of the model.
	 *
	 *	\return the number of points.
	 */
	virtual TInt getNbPoints() const =0;

	/*------------------------------------------------------------------------*/
	/** \brief  Get the number of curves of the model.
	 *
	 *	\return the number of curves.
	 */
	virtual TInt getNbCurves() const =0;

	/*------------------------------------------------------------------------*/
	/** \brief  Get the number of surfaces of the model.
	 *
	 *	\return the number of surfaces.
	 */
	virtual TInt getNbSurfaces() const =0;

	/*------------------------------------------------------------------------*/
	/** \brief  Get the number of volumes of the model.
	 *
	 *	\return the number of volumes.
	 */
	virtual TInt getNbVolumes() const =0;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the points of the model.
	 *
	 *  \param points the points of the model.
	 */
	virtual void getPoints(std::vector<GeomPoint*>& points) const =0;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the curves of the model.
	 *
	 *  \param curves the curves of the model.
	 */
	virtual void getCurves(std::vector<GeomCurve*>& curves) const =0;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the surface of the model.
	 *
	 *  \param surfaces the surfaces of the model.
	 */
	virtual void getSurfaces(std::vector<GeomSurface*>& surfaces)const=0;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the volumes of the model.
	 *
	 *  \param volumes the volumes of the model.
	 */
	virtual void getVolumes(std::vector<GeomVolume*>& volumes) const =0;
};
/*----------------------------------------------------------------------------*/
} // namespace geom
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_GEOM_GEOMMANAGER_H_ */

