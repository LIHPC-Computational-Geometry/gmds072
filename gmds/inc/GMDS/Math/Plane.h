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
 * Plane.h
 *
 *  Created on: 5 octobre 2011
 *      Author: collegial
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_PLANE_H_
#define GMDS_MATH_PLANE_H_
/*----------------------------------------------------------------------------*/
// Gepeto File Headers
#include <GMDS/Utils/CommonTypes.h>
#include <GMDS/IG/Face.h>
#include <GMDS/IG/Node.h>
#include <GMDS/Math/Point.h>
#include <GMDS/Math/Vector.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace math{
/*----------------------------------------------------------------------------*/
class Segment;
class Triangle;
/*----------------------------------------------------------------------------*/
/** \class Plane
 *  \brief template class implementing a geometrical plane in space
 */
/*----------------------------------------------------------------------------*/
	class EXPORT_GMDS Plane
{
private:

	/* point */
	Point pnt_;

	/* vector normal to the plane */
	Vector normal_;

	/* unit vector normal to the plane */
	mutable bool isNormalUnit_;
	mutable Vector normalUnit_;

public:

	/*------------------------------------------------------------------------*/
	/** \brief  constructor.
	 *
	 */
	Plane();

	/*------------------------------------------------------------------------*/
	/** \brief  constructor.
	 *
	 * \param AP point of the plane
	 * \param AN normal vector
	 */
	Plane(const Point&AP, const Vector& AN);

	/*------------------------------------------------------------------------*/
	/** \brief  constructor.
	 *
	 * \param AP1 a point of the plane
	 * \param AP2 a point of the plane
	 * \param AP3 a point of the plane
	 */
	Plane(const Point&AP1, const Point &AP2, const Point &AP3);

	/*------------------------------------------------------------------------*/
	/** \brief  constructor.
	 *
	 * \param AT A 3D Triangle
	 */
	Plane(const Face &AF);

        /*------------------------------------------------------------------------*/
        /** \brief  constructor.
         *
         * \param AT A 3D Triangle
         */
        Plane(const Triangle& AT);

	/*------------------------------------------------------------------------*/
	/** \brief  Overloaded operator=
	 */
	virtual Plane& operator= (const Plane&);

	/*------------------------------------------------------------------------*/
	/** \brief  Overloaded operator==
	 */
	bool operator==(const Plane&) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Overloaded operator!=
	 */
	bool operator!=(const Plane&) const;

	/*------------------------------------------------------------------------*/
	/** \brief  constructor.
	 *
	 * \param AP point of the plane
	 * \param AN normal vector
	 */
	void set(const Point &AP, const Vector & AN) {
		pnt_ = AP;
		normal_ = AN;
		isNormalUnit_ = false;
	}

	/*------------------------------------------------------------------------*/
	/** \brief  constructor.
	 *
	 * \param AP1 a point of the plane
	 * \param AP2 a point of the plane
	 * \param AP3 a point of the plane
	 */
	void set(const Point &AP1, const Point &AP2, const Point &AP3) {
		pnt_ = AP1;
		Vector v1(AP1, AP2);
		Vector v2(AP1, AP3);
		normal_ = v1.cross(v2);
		isNormalUnit_ = false;
	}

	/*------------------------------------------------------------------------*/
	/** \brief  provides the coeff AA, AB, AC and AD for the parametric
	 * 			representation of the plane where AAx,+ABy+ACz=AD is the plane
	 * 			equation
	 *
	 * \param AA first component
	 * \param AB first component
	 * \param AC first component
	 * \param AD first component
	 */
	void getEquationCoeffs(TCoord& AA, TCoord& AB, TCoord& AC, TCoord& AD) const
	{
		AA = normal_.X();
		AB = normal_.Y();
		AC = normal_.Z();
		Vector v(pnt_.X(), pnt_.Y(), pnt_.Z());
		AD = normal_.dot(v);
	}

	/*------------------------------------------------------------------------*/
	/** \brief  predicate indicating if a point is on a plane
	 *
	 * \param AP  a point
	 */
	bool isIn(const Point& AP) const;

	/*------------------------------------------------------------------------*/
	/** \brief  predicate indicating if this point is on the left of the plane.
	 *
	 * \param AP a point
	 */
	bool isStrictlyOnLeft(const Point& AP) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Compute the distance of a point AP to a plane APlane
	 *			Warning, only meaningful in 3D.
	 *
	 * \param AP 	 a point
	 *
	 * \return the distance between AP and AR
	 */
	TCoord distance(const Point& AP) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Orthogonal projection of a point AP onto a plane.
	 *
	 * \param AP 	 a point
	 *
	 * \return the orthogonal projection of AP onto APlane
	 */
	Point project(const Point& AP) const;

	/*------------------------------------------------------------------------*/
	/** \brief  predicate indicating if a segment and a plane intersect each
	 * 			other.
	 * \param AS a segment
	 */
	bool intersect(const Segment& AS, const bool AProper = false) const;
	bool intersect(const Segment& AS, Point &PI, const bool AProper = false) const;



	/*------------------------------------------------------------------------*/
	/** \brief  predicate indicating if two planes intersect each other.
	 *
	 * \param AP a plane
	 * \param AProper indicates if limit cases must be true or false
	 *
	 */
	bool intersect(const Plane& AP, const bool AProper = false) const;



	/*------------------------------------------------------------------------*/
	/** \brief  Getter for the plane point
	 *
	 * \return a point
	 */
	const Point& getPoint() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Getter for the plane normal
	 *
	 * \return a vector
	 */
	const Vector& getNormal() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Getter for the plane normal unit vector
	 *
	 * \return a vector
	 */
	const Vector& getNormalUnit() const;

};
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_PLANE_H_ */
/*----------------------------------------------------------------------------*/
