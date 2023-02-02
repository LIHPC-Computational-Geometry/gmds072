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
 * Point.h
 *
 *  Created on: 6 févr. 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_POINT_H_
#define GMDS_MATH_POINT_H_
/*----------------------------------------------------------------------------*/
// gmds file headers
#include <GMDS/Utils/CommonTypes.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
    /*----------------------------------------------------------------------------*/
    namespace math{
        /*----------------------------------------------------------------------------*/
        /** \class Point
         *  \brief Defines a 3D point
         */
        /*----------------------------------------------------------------------------*/
        class EXPORT_GMDS Point
        {
            
        public:
            
            Point(const TCoord& AX=0.0, const TCoord& AY=0.0, const TCoord& AZ=0.0);
            
            virtual ~Point();
            
            inline TCoord X() const {return m_x;}
            inline TCoord Y() const {return m_y;}
            inline TCoord Z() const {return m_z;}
            
            inline TCoord& X() {return m_x;}
            inline TCoord& Y() {return m_y;}
            inline TCoord& Z() {return m_z;}
            
            inline void setX(const TCoord AVal){m_x = AVal;}
            inline void setY(const TCoord AVal){m_y = AVal;}
            inline void setZ(const TCoord AVal){m_z = AVal;}
            
            inline void setXYZ(const TCoord& AX, const TCoord& AY, const TCoord& AZ)
            {m_x=AX; m_y=AY; m_z=AZ; }
            
            /*------------------------------------------------------------------------*/
            /** \brief  Compute the distance of a point AP to a point
             *
             * \param AP a point
             *
             * \return the distance between AP and this
             */
            TCoord distance(const Point& AP) const;
            
            /*------------------------------------------------------------------------*/
            /** \brief  Compute the square distance of a point AP to a point
             *
             * \param AP a point
             *
             * \return the distance between AP and this
             */
            TCoord distance2(const Point& AP) const;
            
            /*------------------------------------------------------------------------*/
            /** \brief  predicate indicating if three points are colinear
             *
             * \param AP2 a second point
             * \param AP3 a third point
             */
            bool areColinear(const Point &AP2, const Point& AP3) const;
            
            /*------------------------------------------------------------------------*/
            /** \brief  predicate indicating if four points are coplanar
             *
             * 			Warning this predicate is only meaningful in 3D
             *
             *  \param AP2 a second point
             *  \param AP3 a third point
             *  \param AP4 a fourth point
             */
            bool areCoplanar(const Point& AP2, const Point& AP3, const Point& AP4) const;
            
            /*------------------------------------------------------------------------*/
            /** \brief  predicate indicating if a point is on the left of the line formed
             *          by the two other points
             *
             *  \param AP1 a point
             *  \param AP2 a second point
             */
            bool isStrictlyOnLeft2D(const Point& AP1, const Point& AP2) const;
            
            /*------------------------------------------------------------------------*/
            /** \brief  Overloaded operator==
             */
            bool operator==(const Point&) const;
            
            /*------------------------------------------------------------------------*/
            /** \brief  Overloaded operator!=
             */
            bool operator!=(const Point&) const;
            
            /*------------------------------------------------------------------------*/
            /** \brief  Overloaded operator+ to create a new point from 2 points
             */
            friend EXPORT_GMDS Point operator+(const Point&, const Point&);
            
            /*------------------------------------------------------------------------*/
            /** \brief  Overloaded operator- to create a new point from 2 points
             */
            friend EXPORT_GMDS Point operator-(const Point&, const Point&);
            
            /*------------------------------------------------------------------------*/
            /** \brief  Overloaded operator* to create a new point
             */
            friend EXPORT_GMDS Point operator*(const TCoord&, const Point&);
            friend EXPORT_GMDS Point operator*(const Point&, const TCoord&);
            
            /*------------------------------------------------------------------------*/
            /** \brief  Overloaded operator<< for output
             */
            friend EXPORT_GMDS std::ostream& operator<<(std::ostream&, const Point&);
            
            static void computeBarycentric(const math::Point& AT1,
                                           const math::Point& AT2,
                                           const math::Point& AT3,
                                           const math::Point& AT4,
                                           const math::Point& AP,
                                           TCoord& A1,
                                           TCoord& A2,
                                           TCoord& A3,
                                           TCoord& A4);
            
            static void computeBarycentric(const math::Point& AT1,
                                           const math::Point& AT2,
                                           const math::Point& AT3,
                                           const math::Point& AP,
                                           TCoord& AX, TCoord& AY, TCoord& AZ);
            
            
            static void computeBarycentric(const std::vector< math::Point>& AT,
                                           const math::Point& AP,
                                           std::vector<TCoord>& ACoeff);
            
            static Point massCenter(const std::vector<Point>& AT);
            
        protected:
            static void computeBarycentric2D(const math::Point& AT1, const math::Point& AT2,
                                             const math::Point& AT3,	const math::Point& AP,
                                             TCoord& AX, TCoord& AY, TCoord& AZ);
        protected:
            TCoord m_x;
            TCoord m_y;
            TCoord m_z;
        };
        /*--------------------------------------------------------------------*/
    } // namespace math
    /*------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
/** \brief Provides a weak order between two points \p AP1 and \p AP2 in
 *         order to apply STL algorithms on containers of gmds::math::Point.
 *         3D points are ordered along X then Y, then Z. It gives a strict
 *         weak ordering.
 *
 * \param[in] AP1 a first point
 * \param[in] AP2 a second point
 *
 * \return true if AP1 <= AP2, false otherwise
 */
EXPORT_GMDS bool operator<(const gmds::math::Point& AP1,
                           const gmds::math::Point& AP2);
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_POINT_H_ */
/*----------------------------------------------------------------------------*/
