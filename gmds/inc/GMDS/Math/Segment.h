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
 * Segment.h
 *
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_SEGMENT_H_
#define GMDS_MATH_SEGMENT_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <vector>
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/CommonTypes.h>
#include <GMDS/Math/Point.h>
#include <GMDS/Math/Vector.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace math{
/*----------------------------------------------------------------------------*/
class Plane;
/*----------------------------------------------------------------------------*/
/** \class Segment
 *  \brief class implementing a geometrical segment
 */
/*----------------------------------------------------------------------------*/
class EXPORT_GMDS Segment
{
	private:

	/* points */
	Point pnts_[2];

	/*------------------------------------------------------------------------*/
	/* derived attributes */
	mutable bool isVunit_;
	mutable Vector vunit_;
	// TODO dans le cas du segment... vaudrait mieux pas !!??


	/*------------------------------------------------------------------------*/
	/** \brief  reset the "evaluators" stocked are false as center,
	 * sphere including, plane including
	 *
	 */
	void _reset();

	public:

		/*------------------------------------------------------------------------*/
	/** \brief  constructor.
	 *
	 */
	Segment();

	/*------------------------------------------------------------------------*/
	/** \brief  constructor.
	 *
	 * \param AP1 the first point of the segment
	 * \param AP2 the second point of the segment
	 */
	Segment(const Point&AP1, const Point&AP2);

	/*------------------------------------------------------------------------*/
	/** \brief  constructor.
	 *
	 * \param AP1 the first point of the segment
	 * \param AP2 the second point of the segment
	 */
	void set(const Point &AP1, const Point  &AP2);

	/*------------------------------------------------------------------------*/
	/** \brief  Overloaded operator=
	 */
	virtual Segment& operator=(const Segment&);


	/*------------------------------------------------------------------------*/
	/** \brief  Getter for the segment points
	 *
	 * \param AIndex an integer, either 0 ot 1
	 *
	 * \return a point
	 */
	Point getPoint(const int AIndex) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Getter for having the vector corresponding to the segment
	 *
	 * \return a vector
	 */
	Vector getDir() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Return the unit vector defining the segment direction
	 *
	 * \return  the vector
	 */
	Vector& getUnitVector() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Return the center of the segment.
	 *
	 * \return  the center
	 */
	Point computeCenter() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Return the length of the segment.
	 *
	 * \return  the length value
	 */
	TCoord computeLength2() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Return the length of the segment.
	 *
	 * \return  the length value
	 */
	TCoord computeLength() const;


	/*------------------------------------------------------------------------*/
	/** \brief  Overloaded operator==
	 */
	bool operator==(const Segment&) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Overloaded operator!=
	 */
	bool operator!=(const Segment&) const;

	/*------------------------------------------------------------------------*/
	/** \brief  predicate indicating if a point is on a segment
	 *
	 * \param AP a point
	 */
	bool isIn(const Point& AP, const bool AProper = false) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Compute the distance of a point AP to a segment
	 *
	 * \param AP a point
	 *
	 * \return the distance between AP and AS
	 */
	TCoord distance(const Point& AP) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Compute the distance of a segment ASegment to a segment
	 *
	 * \param ASegmùent a segment
	 *
	 * \return the distance between ASegment and this segment
	 */
	TCoord distanceInf(const Segment& ASegment) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Orthogonal projection of a point AP onto a segment
	 *
	 * \param AP a point
	 *
	 * \return the orthogonal projection of AP
	 */
	Point project(const Point& AP) const;


	/*------------------------------------------------------------------------*/
	/** \brief  predicate indicating if a segment and a plane intersect each
	 * 			other.
	 *			With AProper=true, if the segment lies in the plane we get
	 * 			GEOM_UNDEF
	 *
	 * 			Warning this predicate is only meaningful in 3D
	 *
	 * \param AP a plane
	 * \param AProper indicates if limit cases must be GEOM_YES or GEOM_UNDEF
	 *
	 * \return GEOM_YES if they intersect, GEOM_NO if they don't and GEOM_UNDEF
	 * 		   if we cannot conclude
	 */
	bool intersect(const Plane& AP, const bool AProper = false) const;
	bool intersect(const Plane& AP, Point &API, const bool AProper = false) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Orthogonal projection of a point AP onto a segment
         *
         * \param AP a point
         *
         * \return the orthogonal projection of AP
         */
        bool intersect2D(const Segment& ASeg, const bool AProper = false) const;

       
	/*------------------------------------------------------------------------*/
        /** \brief 2D Intersection with a segment
         *
         * \param[IN]  AS the segment we want to get the intersection with
	 * \param[OUT] AP the intersection point if it exists
         *
         * \return true if (*this) intersects AS
         */
        bool intersect2D(const Segment& AS, Point& AP) const;

    bool intersect3D(const Segment& AS, Point& AP, double& AParamSeg,
                     double& AParamThis) const;
    
	/*------------------------------------------------------------------------*/
        /** \brief  Overloaded operator<< for output
         */
        friend EXPORT_GMDS std::ostream& operator<<(std::ostream&, const Segment&);
};
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_SEGMENT_H_ */
/*----------------------------------------------------------------------------*/
