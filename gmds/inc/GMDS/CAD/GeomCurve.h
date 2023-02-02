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
/** \file    GeomCurve.h
 *  \author  F. LEDOUX
 *  \date    09/21/2010
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_GEOM_GEOMCURVE_H_
#define GMDS_GEOM_GEOMCURVE_H_
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include "GMDS/Utils/Exception.h"
#include "GMDS/Utils/CommonTypes.h"

#include "GMDS/CAD/GeomEntity.h"
#include "GMDS/CAD/GeomPoint.h"
#include "GMDS/CAD/GeomSurface.h"
#include "GMDS/CAD/GeomVolume.h"
#include "GMDS/Math/Segment.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace geom{
/*----------------------------------------------------------------------------*/
/** \class GeomCurve
 *  \brief This class describe the services that are required by the
 *  	   mesh to the geometrical model. As a consequence, this interface only
 *  	   contains query methods.
 */
/*----------------------------------------------------------------------------*/
	class EXPORT_GMDS GeomCurve : public GeomEntity {
public:


	/*------------------------------------------------------------------------*/
	/** \brief  Constructor
	 */
	GeomCurve(const std::string& AName = "Unknown curve")
	:GeomEntity(AName){;}

	/*------------------------------------------------------------------------*/
	/** \brief  provides the dimension of the geometrical entity.
	 */
	int getDim() const {return 1;}

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the adjacent points.
	 *
	 *  \param APnt the adjacent points.
	 */
	virtual void get(std::vector<GeomPoint*>& APnt) const  =0;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the adjacent surfaces.
	 *
	 *  \param ASurf the adjacent surfaces
	 */
	virtual void get(std::vector<GeomSurface*>& ASurf) const  =0;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the adjacent volumes
	 *
	 *  \param AVol the adjacent volumes.
	 */
	virtual void get(std::vector<GeomVolume*>& AVol) const =0;

	/*------------------------------------------------------------------------*/
        /** \brief  Length of the curve
         */
        virtual double length() const =0;

	/*------------------------------------------------------------------------*/
	/** \brief  Compute the dihedral angle (max of each edge of the curve)
	 *
	 *  \return the dihedral angle
	 */
	virtual TCoord computeDihedralAngle() const =0;

	/*------------------------------------------------------------------------*/
    /** \brief Compute the tangent vector at one of the endpoints.
	 *  \param AP one of the endpoints
	 *  \param AV the tangent vector
     */
	virtual void computeVector(const GeomPoint& AP,math::Vector& AV) const =0;

	/*------------------------------------------------------------------------*/
    /** \brief Move a point AP near the surface to the closest point on the
     * 		   surface.
	 *  \param AP
     */
	virtual void project(math::Point& AP) const = 0;

	virtual void getMultiplePoints(const int& ANbPoints,math::Point* APoints) const =0;

	virtual TCoord computeDistanceHaussdorf(const math::Segment& ASegment) const =0;
	virtual TCoord computeDistanceHaussdorfSubCurve(const int& ANbPoints,
			const int& AIPoint, const math::Segment& ASegment) const =0;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the first end point
	 *
	 */
	virtual GeomPoint*  getFirstPoint() const =0;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the second end point
	 *
	 */
	virtual GeomPoint*  getSecondPoint() const =0;

    /*------------------------------------------------------------------------*/
    /** \brief  Return whether the curve is a loop or not
     *
     *  \return a boolean
     */
    virtual bool isALoop() const =0;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the surface on the left;
	 * 			left is defined as triangles oriented
	 * 			in the same direction as the curve's segments
	 *
	 */
	virtual GeomSurface* getLeftSurface() const =0;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the surface on the left;
	 * 			left is defined as triangles oriented
	 * 			in the same direction as the curve's segments
	 *
	 */
	virtual GeomSurface* getRightSurface() const =0;

    /*------------------------------------------------------------------------*/
    /** \brief  Add an adjacent point
     *
     *  \param APoint the new adjacent surface to add
     */
    virtual void add(GeomPoint* APoint){throw GMDSException("GeomCurve::add Not yet implemented!");}

    /*------------------------------------------------------------------------*/
    /** \brief  Add an adjacent surface
     *
     *  \param ASurf the new adjacent surface to add
     */
    virtual void add(GeomSurface* ASurf){throw GMDSException("GeomCurve::add Not yet implemented!");}

    /*------------------------------------------------------------------------*/
    /** \brief  Count the number of surfaces in common and store them
     *
     *  \param ACurve the curve this curve is compared to
     *  \param ASurfaces the resulting common surfaces
     *
     *  \return the number of surfaces in common
     */

    virtual int commonSurfaces(GeomCurve* ACurve, std::vector<GeomSurface*>& ASurfaces) const;

};
/*----------------------------------------------------------------------------*/
} // namespace geom
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_GEOM_GEOMCURVE_H_ */

