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
 * Line.h
 *
 *  Created on: sept. 3, 2015
 *      Author: franck ledoux
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_LINE_H_
#define GMDS_MATH_LINE_H_
/*----------------------------------------------------------------------------*/
// Gepeto File Headers
#include <GMDS/Utils/CommonTypes.h>
#include <GMDS/Math/Point.h>
#include <GMDS/Math/Vector.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace math{
/*----------------------------------------------------------------------------*/
class Segment;
class Plane;
class Triangle;
        /*----------------------------------------------------------------------------*/
        /** \class Line
         *  \brief Geometrical class implementing a line defining by 2 points
         */
        /*----------------------------------------------------------------------------*/
        class EXPORT_GMDS Line
        {
        public:
            
            /*------------------------------------------------------------------------*/
            /** \brief  constructor.
             *
             * \param AP1 first point of the line
             * \param AP2 second point of the line
             */
            Line(const Point& AP1, const Point& AP2):m_p1(AP1),m_p2(AP2) {;}
           
	    /*------------------------------------------------------------------------*/
            /** \brief  constructor.
             *
             * \param AP a point on the line
             * \param AVec direction vector of the line
             */
            Line(const Point& AP, const Vector& AVec):m_p1(AP),m_p2(AP+AVec) {;}
 
            /*------------------------------------------------------------------------*/
            /** \brief Copy constructor.
             *
             * \param AL another line
             */
            Line(const Line& AL):m_p1(AL.m_p1),m_p2(AL.m_p2) {;}
            
            /*------------------------------------------------------------------------*/
            /** \brief  Overloaded operator=
             */
            virtual Line& operator= (const Line&);
            
            /*------------------------------------------------------------------------*/
            /** \brief  Getter for the first point
             *
             * \return a point
             */
            const Point& getFirstPoint() const;
            
            /*------------------------------------------------------------------------*/
            /** \brief  Getter for the first point
             *
             * \return a point
             */
            const Point& getSecondPoint() const;
            
            /*------------------------------------------------------------------------*/
            /** \brief 2D Intersection with another line
             *
             * \param[IN]  AL the Line we want to get the intersection with
             * \param[OUT] AP the intersection point if it exists
             * \param[OUT] AParam the a parameter such that AP = a AS[0]+ (1-a) AS[1]
             *
             * \return true if (*this) intersects AL
             */
            bool intersect2D(const Line& AL, Point& AP, double& AParam) const;
           
	    /*------------------------------------------------------------------------*/
            /** \brief 3D Intersection with a plane
             *
             * \param[IN]  APlane the line we want to get the intersection with
             * \param[OUT] AP the intersection point if it exists
             * \param[OUT] AParam the a parameter such that AP = a AS[0]+ (1-a) AS[1]
             *
             * \return true if (*this) intersects APlane
             */
            bool intersect3D(const Plane& APlane, Point& AP, double& AParam) const;

	    /*------------------------------------------------------------------------*/
            /** \brief 3D Intersection with a triangle
             *
             * \param[IN]  ATri the triangle we want to get the intersection with
             * \param[OUT] AP the intersection point if it exists
             * \param[OUT] AParam the a parameter such that AP = a AS[0]+ (1-a) AS[1]
             *
             * \return true if (*this) intersects ATri
             */
            bool intersect3D(const Triangle& ATri, Point& AP, double& AParam) const;
 
            /*------------------------------------------------------------------------*/
            /** \brief 2D distance to a segment. The distance is the minimal orthogonal
             *         distance to @AS
             *
             * \param[IN]  AS  the segment we want to get the distance to
             * \param[OUT] AP1 the point of *this where the min. distance is computed
             * \param[OUT] AP2 the point of AS where the minimal distance is computed
             *
             * \return true if (*this) intersects AS
             */
            TCoord distance2D(const Segment& AS, Point& AP1, Point& AP2) const;

            /*------------------------------------------------------------------------*/
            /** \brief Compute the orthogonal projection of AP on (*this)
             *
             * \param[IN]  AP the point we want to project onto (*this)
             *
             * \return the projected point
             */
            Point project(const Point& AP) const;
            
        private:
            
            /* first point */
            Point m_p1;
            /* second point */
            Point m_p2;
            
        };
        /*----------------------------------------------------------------------------*/
    } // namespace math
    /*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_RAY_H_ */
/*----------------------------------------------------------------------------*/
