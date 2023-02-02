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
 * FacetedGeomManager.h
 *
 *  Created on: 1 juil. 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_GEOM_FACETEDGEOMMANAGER_H_
#define GMDS_GEOM_FACETEDGEOMMANAGER_H_
/*----------------------------------------------------------------------------*/
#include <stack>
/*----------------------------------------------------------------------------*/
// GMDS File headers
/*----------------------------------------------------------------------------*/
#include "GMDS/Utils/CommonTypes.h"
#include "GMDS/Utils/Exception.h"

#include "GMDS/Math/VectorDyn.h"

#include "GMDS/CAD/GeomCurve.h"
#include "GMDS/CAD/GeomManager.h"
#include "GMDS/CAD/GeomServiceAbstractFactory.h"
#include "GMDS/CAD/FacetedVolume.h"
#include "GMDS/CAD/FacetedSurface.h"
#include "GMDS/CAD/FacetedCurve.h"
#include "GMDS/CAD/FacetedPoint.h"
#include "GMDS/CAD/FacetedTriangulationService.h"

#include "GMDS/IG/Edge.h"
#include "GMDS/IG/Node.h"
#include "GMDS/IG/IGMeshDoctor.h"
#include "GMDS/IG/IGMesh.h"

//#include "GMDSIO/VTKReader.h"
//#include "GMDSIO/VTKWriter.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace geom{
/*----------------------------------------------------------------------------*/
/** \class  FacetedGeomManager
 *  \brief  This class creates all the model entitites and services relative to
 *  		the faceted representation of the geometry
 *
 */
/*----------------------------------------------------------------------------*/
	class  EXPORT_GMDS FacetedGeomManager :
		public GeomManager, public GeomServiceAbstractFactory
{

public:
	/*------------------------------------------------------------------------*/
	/** \brief  Default constructor
	 */
	FacetedGeomManager();

	/*------------------------------------------------------------------------*/
	/** \brief  Destructor
	 */
	virtual ~FacetedGeomManager();

	/*------------------------------------------------------------------------*/
	/** \brief  creation of a geometric volume
	 */
	virtual GeomVolume* newVolume();

	/*------------------------------------------------------------------------*/
	/** \brief  creation of a geometric surface
	 */
	virtual GeomSurface* newSurface();

	/*------------------------------------------------------------------------*/
	/** \brief  creation of a geometric curve
	 */
	virtual GeomCurve* newCurve();

	/*------------------------------------------------------------------------*/
	/** \brief  creation of a geometric point
	 */
	virtual GeomPoint* newPoint();

	/*------------------------------------------------------------------------*/
	/** \brief  Get the number of points of the model.
	 *
	 *	\return the number of points.
	 */
	TInt getNbPoints() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Get the number of curves of the model.
	 *
	 *	\return the number of curves.
	 */
	TInt getNbCurves() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Get the number of surfaces of the model.
	 *
	 *	\return the number of surfaces.
	 */
	TInt getNbSurfaces() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Get the number of volumes of the model.
	 *
	 *	\return the number of volumes.
	 */
	TInt getNbVolumes() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the points of the model.
	 *
	 *  \param points the points of the model.
	 */
	void getPoints(std::vector<GeomPoint*>& points) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the curves of the model.
	 *
	 *  \param curves the curves of the model.
	 */
	void getCurves(std::vector<GeomCurve*>& curves) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the surface of the model.
	 *
	 *  \param surfaces the surfaces of the model.
	 */
	void getSurfaces(std::vector<GeomSurface*>& surfaces)const;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the volumes of the model.
	 *
	 *  \param volumes the volumes of the model.
	 */
	void getVolumes(std::vector<GeomVolume*>& volumes) const;

	/*------------------------------------------------------------------------*/
	/** \brief  creation of a geometry triangulation service
	 */
	virtual GeomTriangulationService* newGeomTriangulationService();

	/*------------------------------------------------------------------------*/
	/** \brief  reorient (outward) the triangles describing the model
	 */
	void reorient();

	/*------------------------------------------------------------------------*/
	/** \brief  Gives access to the inner mesh view of the model.
	 *
	 *  Do not modify the referenced mesh content
	 */
	IGMesh& getMeshView();


	/*------------------------------------------------------------------------*/
	/** \brief  reinitializes the complete model from data stored in mesh_.
	 */
	void updateFromMesh();

	/*------------------------------------------------------------------------*/
        /** \brief  reinitializes the complete model from the skin of the mesh.
	 *
	 *  \param AMesh the mesh from which a model will be built.
	 *  \param ASingleSurface build the model with only one surface
         */
        void buildFromMesh(const bool ASingleSurface=false);

	/*------------------------------------------------------------------------*/
        /** \brief  reorient (outward) the triangles describing the volume
         */
        void reorient(FacetedVolume* AVol);

private:


	void OutwardNormal( Face* f);
	bool IsInsideFace(const int&, const math::VectorDyn&,
			const math::VectorDyn&, Face fi);
	math::VectorDyn GetBarycentricCoefs(math::VectorDyn& V1,
			math::VectorDyn& V2,
			math::VectorDyn& V3,
			math::VectorDyn& I,
			math::VectorDyn& normal);


	IGMesh  mesh_;

	std::vector<FacetedVolume* > volumes_;
	std::vector<FacetedSurface*> surfaces_;
	std::vector<FacetedCurve*  > curves_;
	std::vector<FacetedPoint*  > points_;

	std::vector<GeomService*  > services_;
};
/*----------------------------------------------------------------------------*/
} // namespace geom
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_GEOM_FACETEDGEOMMANAGER_H_ */
/*----------------------------------------------------------------------------*/
