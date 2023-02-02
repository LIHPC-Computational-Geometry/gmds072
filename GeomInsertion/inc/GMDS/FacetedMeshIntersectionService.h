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
/*
 * FacetedMeshIntersectionService.h
 *
 *  Created on: 1 juil. 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_GEOM_FACETEDMESHINTERSECTIONSERVICE_H_
#define GMDS_GEOM_FACETEDMESHINTERSECTIONSERVICE_H_
/*----------------------------------------------------------------------------*/
#include <set>
#include <map>
#include <vector>
#include <gts.h>
/*----------------------------------------------------------------------------*/
#include "GMDS/CAD/FacetedVolume.h"
#include "GMDS/CAD/FacetedSurface.h"
#include "GMDS/CAD/FacetedCurve.h"

#include "GMDS/Math/Hexahedron.h"
/*----------------------------------------------------------------------------*/
#include "GeomMeshIntersectionService.h"
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
namespace geom {
/*----------------------------------------------------------------------------*/
// epsilon for numeric type in isInsideMC method
const int FacetedMeshIntersectionService_isInsideMC_EPS = 12;
// number of rays drawn in isInsideMC method
const int FacetedMeshIntersectionService_isInsideMC_NBDRAWS = 100;

const int FacetedMeshIntersectionService_project_EPS = 12;
const int FacetedMeshIntersectionService_intersect_curve = 13;
const int FacetedMeshIntersectionService_intersect_hex = 13;
/*----------------------------------------------------------------------------*/
/** \class  GeomModelAbstractFactory
 *  \brief  This interface gathers the factory methods required for every
 *  		geometric model.
 *
 *  \param TBase the basic type used to store geometric coordinates.
 */
/*----------------------------------------------------------------------------*/
class FacetedMeshIntersectionService: public GeomMeshIntersectionService {

public:

    /*------------------------------------------------------------------------*/
    /** \brief Constructor.
     *
     *  \param AMesh the mesh to build.
     */
	FacetedMeshIntersectionService();

    /*------------------------------------------------------------------------*/
    /** \brief  Destructor.	*/
	virtual ~FacetedMeshIntersectionService();

//	/*------------------------------------------------------------------------*/
//	/** \brief  computes the bounding box
//	 *
//	 *	\param minXYZ The minimum coordinate of the bounding box.
//	 *	\param maxXYZ The maximum coordinate of the bounding box.
//	 */
//	virtual void computeBoundingBox(std::vector<GeomSurface<TBase>* >& ASurfaces,
//						TBase minXYZ[3], TBase maxXYZ[3]) const;

	/*------------------------------------------------------------------------*/
	/** \brief Initialize the
	 *
	 *  \param ASurfaces the surfaces of the model
	 */
	virtual void initialization(std::vector<GeomSurface* >& ASurfaces);

	/*------------------------------------------------------------------------*/
	/** \brief check whether a hexahedra intersects the volume surface
	 *
	 *  \param ASurfaces the surfaces of the model
	 *  \param AHex a hexahedra
	 *
	 *  \return a boolean
	 */
	virtual bool intersects(std::vector<GeomSurface* >& ASurfaces,
						gmds::math::Hexahedron& AHex) const;

	/*------------------------------------------------------------------------*/
	/** \brief check whether a region intersects the volume surface
	 *
	 *  \param ASurfaces the surfaces of the model
	 *  \param ARegion a region
	 *
	 *  \return a boolean
	 */
	virtual bool intersects(std::vector<GeomSurface* >& ASurfaces,
						gmds::Region& ARegion) const;

	/*------------------------------------------------------------------------*/
	/** \brief  check whether a region intersects the volume surface and stores the
	 *		   	intersected triangles of the surfaces triangulations.
	 *
	 *  \param ASurfaces the surfaces of the model
	 *  \param ARegion a region
	 *  \param ATriangles a vector of the intersected triangles
	 *
	 *  \return a boolean
	 */
	virtual bool intersects(std::vector<GeomSurface* >& ASurfaces,
						gmds::Region& ARegion,
						std::vector<gmds::math::Triangle>& Atriangles) const;

	virtual bool intersects(
						gmds::Region& ARegion,
						std::vector<gmds::math::Triangle>& Atriangles) const;

	virtual bool intersects(
							GNode* ATree,
							std::vector<GeomSurface* >& ASurfaces,
							gmds::Region& ARegion,
							std::vector<gmds::math::Triangle>& Atriangles) const;

	virtual void intersects(
							GNode* ATree,
							std::map<gmds::TCellID, std::vector<gmds::math::Triangle*> >& AInOutTriangles) const;
	virtual void intersects(
							GNode* ATree,
							std::map<gmds::TCellID, std::map<gmds::geom::GeomSurface*,std::vector<gmds::math::Triangle*> > >& AInOutSurfTriangles) const;

	/*------------------------------------------------------------------------*/
	/** \brief check whether a region intersects a triangle
	 *
	 *  \param ATri a triangle
	 *  \param ARegion a region
	 *
	 *  \return a boolean
	 */
	virtual bool intersects(gmds::math::Triangle& ATri,
						gmds::Region& ARegion) const;

//	/*------------------------------------------------------------------------*/
//	/** \brief check whether a hexahedra intersects twice the volume surface,
//	 * 		   with two surfaces that are not neighbors.
//	 *
//	 *  \param ASurfaces the surfaces of the model
//	 *  \param AHex a hexahedra
//	 *
//	 *  \return a boolean
//	 */
	virtual bool intersectsTwiceTheModel(std::vector<GeomSurface*>& ASurfaces,
						std::map<gmds::geom::GeomSurface*, std::set<gmds::geom::GeomSurface*> > ASurfacesNeighbors,
						gmds::Region& ARegion) const;
//
//	/*------------------------------------------------------------------------*/
//	/** \brief check whether a hexahedra intersects the curve
//	 *
//	 *  \param ACurve a curve of the model
//	 *  \param AHex a hexahedra
//	 *
//	 *  \return a boolean
//	 */
//	virtual bool intersects(GeomCurve<TBase>& ACurve,
//						GEPETO::Hexahedron<TBase>& AHex) const;
//
//	/*------------------------------------------------------------------------*/
//	/** \brief check whether a hexahedron is inside the model. Based on the number of
//	 * 		intersections between lines passing through the hexahedron with
//	 *      with seven rays tested.
//	 *
//	 *  \param ASurfaces the surfaces of the model
//	 *  \param AHex a hexahedron
//	 *
//	 *  \return a boolean
//	 */
//	virtual bool isInside(std::vector<GeomSurface<TBase>* >& ASurfaces,
//						GEPETO::Hexahedron<TBase>& AHex) const;
//
//	/*------------------------------------------------------------------------*/
//	/** \brief check whether a hexahedra is inside the model.
//	 * 		   Based on a Monte-Carlo algorithm with multiple rays tested.
//	 * 		   Takes longer than the previous method but is more robust.
//	 *
//	 *  \param ASurfaces the surfaces of the model
//	 *  \param AHex a hexahedron
//	 *
//	 *  \return a boolean
//	 */
//	virtual bool isInsideMC(
//			std::vector<GeomSurface<TBase>* >& ASurfaces,
//			std::vector<GEPETO::Vector<3,TBase > >& ARays,
//			GEPETO::Point<3, TBase >  & APoint) const;
//
	/*------------------------------------------------------------------------*/
	/** \brief check whether a hexahedra is inside the model.
	 * 		   Based on a Monte-Carlo algorithm with multiple rays tested.
	 * 		   Takes longer than the previous method but is more robust.
	 *
	 *  \param ASurfaces the surfaces of the model
	 *  \param AHex a hexahedron
	 *
	 *  \return a boolean
	 */
	virtual bool isInsideMC(std::vector<GeomSurface* >& ASurfaces,
						gmds::math::Hexahedron& AHex) const;

	/*------------------------------------------------------------------------*/
	/** \brief check whether a region is inside the model.
	 * 		   Based on a Monte-Carlo algorithm with multiple rays tested.
	 * 		   Takes longer than the previous method but is more robust.
	 *
	 *  \param ASurfaces the surfaces of the model
	 *  \param ARegion a region
	 *
	 *  \return a boolean
	 */
	virtual bool isInsideMC(std::vector<GeomSurface* >& ASurfaces,
						gmds::Region& ARegion) const;
//
//	/*------------------------------------------------------------------------*/
//	/** \brief check whether a hexahedra is mainly inside the model.
//	 *
//	 *  \param ASurfaces the surfaces of the model
//	 *  \param AHex a hexahedron
//	 *
//	 *  \return a boolean
//	 */
//	virtual bool isMainlyInside(
//			std::vector<GEPETO::Triangle<3,TBase > >& ATriangles,
//			GEPETO::Hexahedron<TBase>& AHex) const;
//
	/*------------------------------------------------------------------------*/
	/** \brief check whether a hexahedra is mainly inside the model, using
	 * 		   a "left of nearest triangle" criteria.
	 *
	 *  \param ATriangles the triangles that intersect the region
	 *  \param AHex a hexahedron
	 *  \param ARatio the ratio of the region that is inside
	 *
	 *  \return a boolean
	 */
	virtual bool isMainlyInsideMC(
			std::vector<gmds::math::Triangle>& ATriangles,
			gmds::math::Hexahedron& AHex,
			double& ARatio) const;

	/*------------------------------------------------------------------------*/
	/** \brief check whether a region is mainly inside the model, using
	 * 		   a "left of nearest triangle" criteria.
	 *
	 *  \param ATriangles the triangles that intersect the region
	 *  \param ARegion a region
	 *  \param ARatio the ratio of the region that is inside
	 *
	 *  \return a boolean
	 */
	virtual bool isMainlyInsideMC(
			std::vector<gmds::math::Triangle>& ATriangles,
			gmds::Region& ARegion,
			double& ARatio) const;

	/*------------------------------------------------------------------------*/
	/** \brief check whether a point is inside; based on the normal orientation
	 * 		   of the triangles.
	 *
	 *  \param ATriangles the triangles that intersect the region
	 *  \param APoint a point
	 *
	 *  \return a boolean
	 */
	virtual bool isInside(
			std::vector<gmds::math::Triangle>& ATriangles,
			gmds::math::Point& APoint) const;

//	/*------------------------------------------------------------------------*/
//	/** \brief project a point on the surface of the model, orthogonal projection
//	 *
//	 *  \param ASurfaces the surfaces of the model
//	 *  \param APoint a point
//	 *
//	 *  \return a point
//	 */
//	virtual GEPETO::Point<3,TBase> project(const std::vector<GeomSurface<TBase>* >& ASurfaces,
//						 const GEPETO::Point<3,TBase>& APoint) const;
//
//	/*------------------------------------------------------------------------*/
//	/** \brief project a point on the surface of the model, orthogonal projection
//	 *
//	 *  \param ASurfaces the surfaces of the model
//	 *  \param APoint a point
//	 *  \param ASurfaceDestID the id of the surface on which the point is projected
//	 *
//	 *  \return a point
//	 */
//	virtual GEPETO::Point<3,TBase> project(const std::vector<GeomSurface<TBase>* >& ASurfaces,
//						 const GEPETO::Point<3,TBase>& APoint, int& ASurfaceDestID) const;
//
//	/*------------------------------------------------------------------------*/
//	/** \brief project a point on a curve, orthogonal projection
//	 *
//	 *  \param ACurve the curve of the model
//	 *  \param APoint a point
//	 *
//	 *  \return a point
//	 */
//	virtual GEPETO::Point<3,TBase> project(GeomCurve<TBase>& ACurve,
//						 const GEPETO::Point<3,TBase>& APoint) const;
//
//	/*------------------------------------------------------------------------*/
//	/** \brief project a point on the surface of the model, projection along a vector
//	 *
//	 *  \param ASurfaces the surfaces of the model
//	 *  \param APoint a point
//	 *  \param AVector
//	 *
//	 *  \return a point
//	 */
//	virtual GEPETO::Point<3,TBase> project(const std::vector<GeomSurface<TBase>* >& ASurfaces,
//						 const GEPETO::Point<3,TBase>& APoint,
//						 const GEPETO::Vector<3,TBase>& AVector) const;
//
//	/*------------------------------------------------------------------------*/
//	/** \brief project a point on the surface of the model, projection along a vector
//	 *
//	 *  \param ASurfaces the surfaces of the model
//	 *  \param APoint a point
//	 *  \param AVector
//	 *  \param ASurfaceDestID the id of the surface on which the point is projected
//	 *
//	 *  \return a point
//	 */
//	virtual GEPETO::Point<3,TBase> project(const std::vector<GeomSurface<TBase>* >& ASurfaces,
//						 const GEPETO::Point<3,TBase>& APoint,
//						 const GEPETO::Vector<3,TBase>& AVector,
//						 int& ASurfaceDestID) const;



	bool intersect(gmds::math::Hexahedron& AHex, gmds::math::Triangle& ATri) const;


	GNode* buildAABBSurfacesTriangulationTree(
			std::vector<GeomSurface*>& ASurfaces,
			std::map<GeomSurface*,GNode*>& ASurfacesTriangulationTrees);

	void buildAABBSurfacesTriangulationTree(
                        std::vector<GeomSurface*>& ASurfaces);

	void project(gmds::math::Point& APoint, GNode* ATree) const;

	virtual void project(
                        gmds::geom::GeomSurface* ASurf,
                        gmds::math::Point& APoint) const;

	static GtsPoint* FacetedMeshIntersectionService_triangle_project(
					GtsPoint* AP,
					gpointer bounded)
	{
		gmds::math::Point p(AP->x,AP->y,AP->z);
		gmds::math::Triangle* tri = (gmds::math::Triangle*) bounded;
		gmds::math::Point newP = tri->project(p);
		GtsPoint* gtsp = gts_point_new(gts_point_class (),
                                newP.X(),newP.Y(),newP.Z());
		
		return gtsp;
	}

	/*------------------------------------------------------------------------*/
        /** \brief Compute the normal at a point. GeomMeshIntersectionService must be a 
         *         closed set of surfaces (basically it must be the boundary of a volume).
         *
         *  \param APoint      the point from which we compute the normal 
         *  \param AGeomEntity the GeomEntity of which we compute the normal
         *  \param AOutward    a boolean that indicates whether we want the outward or the inward normal 
         *  \param ANormal     the computed normal
         *
        */
        virtual void computeNormal(
                        const gmds::math::Point APoint,
                        gmds::geom::GeomEntity* AGeomEntity,
                        const bool AOutward,
                        gmds::math::Vector& ANormal) const;	

protected :

	/*------------------------------------------------------------------------*/
	/** \brief initialize random generator
	 *
	 */
	virtual void random_init() const;

	/*------------------------------------------------------------------------*/
	/** \brief returns a random double ranging between 0 and 1.
	 *
	 *  \return a random double
	 */
	virtual double random_value() const;

	// Axis-Aligned Bounding Box tree for the surfaces triangles
        GNode* aabbSurfacesTrianglesTree_;

        // Axis-Aligned Bounding Box tree for the surfaces triangles;
        // one tree for each surface
        std::map<gmds::geom::GeomSurface*,GNode*> aabbSurfacesTrianglesTrees_;

	mutable std::map<gmds::geom::GeomSurface*,std::vector<gmds::math::Triangle> > m_surfacesTriangulation;
	mutable std::map<gmds::geom::GeomSurface*,std::vector<gmds::math::Point> > m_surfacesBoundingBox;

};
/*----------------------------------------------------------------------------*/
} // namespace geom
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_GEOM_FACETEDMESHINTERSECTIONSERVICE_H_ */
/*----------------------------------------------------------------------------*/
