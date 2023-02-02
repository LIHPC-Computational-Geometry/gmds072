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
 * Triangle.h
 *
 *  Created on: 3 juil. 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_TRIANGLE_H_
#define GMDS_MATH_TRIANGLE_H_
/*----------------------------------------------------------------------------*/
// gmds file headers
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/CommonTypes.h>
#include <GMDS/Math/Point.h>
#include <GMDS/Math/Vector.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace math{
/*----------------------------------------------------------------------------*/
class Segment;
class Ray;
class Plane;
/*----------------------------------------------------------------------------*/
/** \class Triangle
 *  \brief Defines a 3D point
 */
/*----------------------------------------------------------------------------*/
	class EXPORT_GMDS Triangle
{

public:

	/*------------------------------------------------------------------------*/
	/** \brief  constructor
	 *
	 */
	Triangle();

	/*------------------------------------------------------------------------*/
	/** \brief  constructor
	 *
	 * \param AP1 a point of the triangle
	 * \param AP2 a point of the triangle (after AP1)
	 * \param AP3 a point of the triangle (after AP2, before AP1)
	 *
	 */
	Triangle(const Point& AP1, const Point& AP2, const Point& AP3);

	/*------------------------------------------------------------------------*/
	/** \brief  constructor
	 *
	 * \param AT a triangle
	 *
	 */
	Triangle(const Triangle& AT);

	/*------------------------------------------------------------------------*/
	/** \brief  destructor
	 */
	virtual ~Triangle();

    /*------------------------------------------------------------------------*/
    /** \brief  Getter for the triangle point
     *
     * \param AIndex an integer
     *
     * \return a point
     */
    const Point& getPoint(const TInt& AIndex) const;
    /*------------------------------------------------------------------------*/
    /** \brief  Computes the area of the triangle
     */
    double area() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Return the number of points
         *
         * \return the number of points
         */
        int getNbPoints() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Computes if the triangle is not degenerated (i.e. flat)
	 *
	 * \return 	false if the triangle points are aligned
	 * 			else true
	 */
	bool isGood() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Computes the normal of the triangle
         *
         * \return the normal of the triangle
         */
        Vector getNormal() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Computes the center of the triangle
         *
         * \return the center of the triangle
         */
        Point getCenter() const;

	/*------------------------------------------------------------------------*/ 
	/** \brief  Compute the scaled jacobian of the quadrilateral
         *
         * \return the scaled jacobian
         */
        double computeScaledJacobian2D() const;

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the scaled jacobian of the quadrilateral,
	 *          normalized between [-1., 1.]
         *
         * \return the scaled jacobian
         */
        double computeNormalizedScaledJacobian2D() const;

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the mean edge length of the triangle.
         *
         * \return the mean edge length
         */
        double computeMeanEdgeLength() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Computes the bounding box of the triangle.
         *
         * \param AMinXYZ the lower front left coordinates
         * \param AMaxXYZ the upper back right coordinates
         */
        void computeBoundingBox(TCoord AMinXYZ[3], TCoord AMaxXYZ[3]) const;

        /*------------------------------------------------------------------------*/
        /** \brief  predicate indicating if a triangle and a triangle intersect each
         *          other.
         * \param ATri a triangle
         */
        bool intersect(const Triangle& ATri, const bool AProper = false) const;

        /*------------------------------------------------------------------------*/
        /** \brief  predicate indicating if a triangle and a triangle intersect each
         *          other, in 2D.
         * \param ATri a triangle
         */
        bool intersect2D(const Triangle& ATri, const bool AProper = false) const;

        /*------------------------------------------------------------------------*/
        /** \brief  predicate indicating if a triangle and a segment intersect each
         *          other.
         * \param ASeg a segment
         */
        bool intersect(const Segment& ASeg, const bool AProper = false) const;

        /*------------------------------------------------------------------------*/
        /** \brief  predicate indicating if a triangle and a segment intersect each
         *          other, in 2D.
         * \param ASeg a segment
         */
        bool intersect2D(const Segment& ASeg, const bool AProper = false) const;

        /*------------------------------------------------------------------------*/
        /** \brief  predicate indicating if a triangle and a ray intersect each
         *          other.
         * \param ARay a ray
         */
        bool intersect(const Ray& ARay, const bool AProper = false) const;

	
        /*------------------------------------------------------------------------*/
        /** \brief  predicate indicating if a point is in a triangle
         *  
         * \param AP a point
         */
	bool isIn(const Point& AP) const;

        /*------------------------------------------------------------------------*/
        /** \brief  predicate indicating if a point is strictly in a triangle
         *
         * \param AP a point
         */
	bool isStrictlyIn(const Point& AP) const;

	/*------------------------------------------------------------------------*/
        /** \brief  Return a plane built from the triangle
         *
         * \return a plane
         */
	Plane getPlaneIncluding() const;

        /*------------------------------------------------------------------------*/
        /** \brief  Return the square distance between a point and a triangle
         *
         * \return a scalar
         */
        Point project(const Point& APoint) const;


        /*------------------------------------------------------------------------*/
        /** \brief  Return the square distance between a point and a triangle
         *
         * \return a scalar
         */
        TCoord distance2(const Point& APoint) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Return the distance between a point and a triangle
         *
         * \return a scalar
         */
        TCoord distance(const Point& APoint) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Overloaded operator<< for output
         */
        friend std::ostream& operator<<(std::ostream&, const Triangle&);


protected:
	Point m_pnts[3];

};
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_TRIANGLE_H_ */
/*----------------------------------------------------------------------------*/
