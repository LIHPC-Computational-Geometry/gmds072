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
/** \file    FacetedCurve.h
 *  \author  F. LEDOUX
 *  \date    30/05/2011
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_GEOM_FACETEDCURVE_H_
#define GMDS_GEOM_FACETEDCURVE_H_
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include "GMDS/Utils/CommonTypes.h"
#include "GMDS/Utils/Exception.h"
#include "GMDS/CAD/GeomCurve.h"
#include "GMDS/CAD/FacetedPoint.h"
#include "GMDS/IG/Edge.h"
#include "GMDS/IG/Node.h"
#include "GMDS/Math/Point.h"
#include "GMDS/Math/Vector.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace geom{
/*----------------------------------------------------------------------------*/
class FacetedSurface;
/*----------------------------------------------------------------------------*/
/** \class FacetedCurve
 *  \brief This class implements the curve services that are required by the
 *  	   mesh to the geometrical model.
 */
/*----------------------------------------------------------------------------*/
class EXPORT_GMDS FacetedCurve: public GeomCurve {

public:
	/*------------------------------------------------------------------------*/
	/** \brief  Default Constructor
	 */
	FacetedCurve();

	/*------------------------------------------------------------------------*/
	/** \brief  Constructor. A geometric curve is built as an ordered
	 * 			collection of points.
	 * 			The points vector must contain
	 * 		 	number of segments + 1 nodes : when the curve is a loop,
	 * 			the last point is also the first one.
	 *  \param AP1 		the first end point
	 *  \param AP2 		the second end point
	 *  \param APoints 	the node used to discretize the curve
	 *  				(including end points).
	 */
	FacetedCurve( FacetedPoint* AP1,
				  FacetedPoint* AP2,
				  std::vector<Node>& APoints,
				  std::vector<Edge>& AEdges,
				  const std::string& AName = "Unknown curve");

	/*------------------------------------------------------------------------*/
	/** \brief  Destructor
	 */
	virtual ~FacetedCurve();

	/*------------------------------------------------------------------------*/
	/** \brief  Add an adjacent surface
	 *
	 *  \param ASurf the new adjacent surface to add
	 */
	virtual void add(FacetedSurface* ASurf);
	// in order to avoid the hidden function warning
	virtual void add(GeomPoint*){throw GMDSException("FacetedCurve::add(GeomPoint*) should not be called");}
	virtual void add(GeomSurface*){throw GMDSException("FacetedCurve::add(GeomSurface*) should not be called");}

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the adjacent points.
	 *
	 *  \param APnt the adjacent points.
	 */
	virtual void get(std::vector<GeomPoint*>& APnt) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the adjacent surfaces.
	 *
	 *  \param ASurf the adjacent surfaces
	 */
	virtual void get(std::vector<GeomSurface*>& ASurf)const ;


	/*------------------------------------------------------------------------*/
	/** \brief  Access to the adjacent volumes
	 *
	 *  \param AVol the adjacent volumes.
	 */
	virtual void get(std::vector<GeomVolume*>& AVol) const;

	/*------------------------------------------------------------------------*/
        /** \brief  Length of the curve
         */
	double length() const;
	
	/*------------------------------------------------------------------------*/
    /** \brief Move a point AP near the surface to the closest point on the
     * 		   surface.
	 *  \param AP
     */
	virtual void project(math::Point& AP) const;


	/*------------------------------------------------------------------------*/
    /** \brief Move a point AP near the surface to the closest point on the
     * 		   surface. Also fills the tangent vector.
	 *  \param AP
	 *  \param AV the tangent vector
     */
	virtual void project(math::Point& AP, math::Vector& AV) const;

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
	/** \brief  Access to the first end point
	 *
	 */
	GeomPoint*  getFirstPoint() const {return p1_;}

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the second end point
	 *
	 */
	GeomPoint*  getSecondPoint() const {return p2_;}

    /*------------------------------------------------------------------------*/
    /** \brief  Return whether the curve is a loop or not
     *
     *  \return a boolean
     */
    virtual bool isALoop() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Compute the dihedral angle (max for each edge of the curve)
	 *
	 *  \return the dihedral angle
	 */
	TCoord computeDihedralAngle() const;

	/*------------------------------------------------------------------------*/
    /** \brief Compute the tangent vector at one of the endpoints.
	 *  \param AP one of the endpoints
	 *  \param AV the tangent vector
     */
	virtual void computeVector(const GeomPoint& AP, math::Vector& AV) const;


	/*------------------------------------------------------------------------*/
	/** \brief  Access to the surface on the left;
	 * 			left is defined as triangles oriented
	 * 			in the same direction as the curve's segments
	 *
	 */
	virtual GeomSurface* getLeftSurface() const;
	virtual GeomSurface* getRightSurface() const;

	virtual unsigned int getNbSegments() const;
	virtual math::Segment  getSegment(const unsigned int& AISegment) const;

	virtual void getMultiplePoints(const int& ANbPoints, math::Point* APoints) const;

	virtual TCoord computeDistanceHaussdorf(const math::Segment& ASegment) const;
	virtual TCoord computeDistanceHaussdorfSubCurve(const int& ANbPoints,
			const int& AIPoint, const math::Segment& ASegment) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Returns a copy of the internal mesh representation nodes
	 *
	 *  \param ANodes a vector of mesh nodes
	 */
	void getMeshNodes(std::vector<Node>& ANodes) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Returns a copy of the internal mesh representation edges
	 *
	 *  \param AEdges a vector of mesh edges
	 */
	void getMeshEdges(std::vector<Edge>& AEdges) const;

	int getId() const {return id_;}
private:

	int id_;
	static int next_id_;

	FacetedPoint* p1_;
	FacetedPoint* p2_;

	std::vector<FacetedSurface*> surfaces_;

	std::vector<Node > mesh_representation_;
	std::vector<Edge > mesh_representation_edges_;

	std::vector<math::Segment> geom_representation_;

};
/*----------------------------------------------------------------------------*/
} // namespace geom
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_GEOM_FACETEDCURVE_H_ */
/*----------------------------------------------------------------------------*/
