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
 * FacetedSurface.h
 *
 *  Created on: 27 juin 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_GEOM_FACETEDSURFACE_H_
#define GMDS_GEOM_FACETEDSURFACE_H_
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include "GMDS/CAD/GeomSurface.h"
#include "GMDS/IG/Face.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace geom{
/*----------------------------------------------------------------------------*/
class FacetedCurve;
class FacetedPoint;
class FacetedVolume;
/*----------------------------------------------------------------------------*/
/** \class FacetedSurface
 *  \brief This class implements the surface services that are required by the
 *  	   mesh to the geometrical model.
 *  \param TCoord the data type used to store geometrical data
 */
/*----------------------------------------------------------------------------*/
class EXPORT_GMDS FacetedSurface : public GeomSurface {

public:

	/*------------------------------------------------------------------------*/
	/** \brief  Default Constructor.
	 */
	FacetedSurface();

	/*------------------------------------------------------------------------*/
	/** \brief  Constructor.
	 *
	 *  \param AP 		the points adjacent to the surface
	 *  \param AP2 		the curves adjacent to the surface
	 *  \param APoints 	the triangles discretizing the surface
	 *  \param AName	the surface name
	 */
	FacetedSurface(std::vector<FacetedPoint* >& AP,
			std::vector<FacetedCurve* >& AC,
			std::vector<Face>& ADiscret,
			const std::string& AName="Unknown surface");

	/*------------------------------------------------------------------------*/
	/** \brief  Destructor.
	 */
	virtual ~FacetedSurface();


	/*------------------------------------------------------------------------*/
	/** \brief  computes the area of the entity.
	 */
	virtual TCoord computeArea() const;

	/*------------------------------------------------------------------------*/
	/** \brief  computes the bounding box
	 *
	 *	\param minXYZ The minimum coordinate of the bounding box.
	 *	\param maxXYZ The maximum coordinate of the bounding box.
	 */
	virtual void computeBoundingBox(TCoord minXYZ[3], TCoord maxXYZ[3]) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Add an adjacent volume
	 *
	 *  \param AVol the new adjacent volume to add
	 */
	virtual void add(FacetedVolume* AVol);
	// in order to avoid the hidden function warning
	virtual void add(GeomVolume*) {throw GMDSException("FacetedSurface::add(GeomVolume*) should not be called");}
	virtual void add(GeomCurve*) {throw GMDSException("FacetedSurface::add(GeomCurve*) should not be called");}

	/*------------------------------------------------------------------------*/
	/** \brief  replace adjacent points
	 *
	 *  \param APnt the adjacent points.
	 */
	virtual void replace(std::vector<FacetedPoint* > APoints);

	/*------------------------------------------------------------------------*/
	/** \brief  replace adjacent curves
	 *
	 *  \param ACur the adjacent Curves.
	 */
	virtual void replace(std::vector<FacetedCurve* > ACurves);

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the adjacent points.
	 *
	 *  \param APnt the adjacent points.
	 */
	virtual void get(std::vector<GeomPoint*>& APnt) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the adjacent curves.
	 *
	 *  \param ACur the adjacent curves.
	 */
	virtual void get(std::vector<GeomCurve*>& ACur) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the adjacent volumes
	 *
	 *  \param AVol the adjacent volumes.
	 */
	virtual void get(std::vector<GeomVolume*>& AVol) const;

	/*------------------------------------------------------------------------*/
    /** \brief Move a point AP near the surface to the closest point on the
     * 		   surface.
	 *  \param AP
     */
//	virtual void project(math::Point<3,TCoord>& AP) const;

	/*------------------------------------------------------------------------*/
	/** \brief  computes normal at the closest point to AP in 3D. In this case,
	 * 			it is the normal to one of the triangles that compose the
	 * 			surface.
	 *
	 *  \param AP the point
	 *
	 *  \param AV 	normal vector at the closest point of AP on this
	 */
	virtual void computeNormal(	const math::Point& AP, math::Vector& AV) const;

	/*------------------------------------------------------------------------*/
    /** \brief Get the closest point from AP on the surface
	 *  \param AP a 3D point
	 *
	 *  \return the closest point of APoint on the surface
     */

	virtual math::Point	closestPoint(const math::Point& AP) const;

	/*------------------------------------------------------------------------*/
    /** \brief Get center of the surface.
	 *
	 *  \return the center of the surface
     */
	virtual math::Point	getCenter() const;

	/*------------------------------------------------------------------------*/
	int getId() const {return id_;}

	/*------------------------------------------------------------------------*/
	/** \brief  Returns a copy of the internal mesh representation
	 *
	 *  \param AFaces a vector of mesh faces
	 */
	void getMeshFaces(std::vector<Face>& AFaces) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Returns a triangulation of the surface
	 *
	 *  \param ATri a triangulation
	 */
	virtual void getTriangulation(std::vector<math::Triangle >& ATri) const;

	void reorient(FacetedVolume* AVol);
	void propagateOrient(Face AFace, int AMarkTreatedFaces, int AMarkFacesInverted, IGMesh* AMesh);
	bool checkSameOrientFace(Face AFaceRef, Face AFaceCheck);
	bool isOutwardDirection(FacetedVolume* AVol);
	void invertAllFaces();
	void invertFace(Face AFace);

private:

	int id_;
	static int next_id_;
	std::vector<FacetedPoint* > points_;
	std::vector<FacetedCurve* > curves_;
	std::vector<FacetedVolume* > volumes_;

	std::vector<Face > mesh_representation_;
};
/*----------------------------------------------------------------------------*/
} // namespace geom
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_GEOM_FACETEDSURFACE_H_ */
/*----------------------------------------------------------------------------*/
