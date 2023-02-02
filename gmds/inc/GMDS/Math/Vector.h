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
 * Vector.h
 *
 *  Created on: 3 juin 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_VECTOR_H_
#define GMDS_MATH_VECTOR_H_
/*----------------------------------------------------------------------------*/
// gmds file headers
/*----------------------------------------------------------------------------*/
#include <GMDS/Math/Constants.h>
#include <GMDS/Math/Point.h>
/*----------------------------------------------------------------------------*/
#include<cmath>
#include<string.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
    /*----------------------------------------------------------------------------*/
    namespace math{
        /*----------------------------------------------------------------------------*/
        /** \class Vector
         *  \brief Defines a 3D Vector
         */
        class EXPORT_GMDS Vector {
        public:
            
            /*------------------------------------------------------------------------*/
            /** \brief Constructor.
             */
            Vector(const TCoord& AX = 0, const TCoord& AY = 0, const TCoord& AZ = 0);
            
            /*------------------------------------------------------------------------*/
            /** \brief Vector constructor from 2 points
             */
            Vector(const Point& AP1, const Point& AP2);
            /*------------------------------------------------------------------------*/
            /** \brief Vector constructor from 1 point. It gives vector OP with O the
             * 		   origin (0.0,...,0.0) in the nD space.
             */
            Vector(const Point& AP);
            
            /*------------------------------------------------------------------------*/
            /** \brief Constructor from a C array
             */
            Vector(TCoord AParam[3]);
            
            /*------------------------------------------------------------------------*/
            /** \brief  Copy constructor.
             */
            Vector(const Vector& AV);
            
            /*------------------------------------------------------------------------*/
            /** \brief  Named accessors to X, Y and Z coordinates
             */
            inline TCoord X() const { return m_tab[0]; }
            inline TCoord Y() const { return m_tab[1]; }
            inline TCoord Z() const { return m_tab[2]; }
            
            /*------------------------------------------------------------------------*/
            /** \brief  Overloaded operator=
             */
            inline Vector& operator= (const Vector& AV){
                if (AV == *this)
                    return *this;
                
                memcpy(&m_tab[0], &AV.m_tab[0], 3 * sizeof(TCoord));
                return *this;
            }
            
            
            /*------------------------------------------------------------------------*/
            /** \brief  Overloaded operator==
             */
            inline bool operator== (const Vector& AV) const {
                return (m_tab[0] == AV.m_tab[0] &&
                        m_tab[1] == AV.m_tab[1] &&
                        m_tab[2] == AV.m_tab[2]);
            }
            
            /*------------------------------------------------------------------------*/
            /** \brief  Overloaded operator!=
             */
            bool operator!= (const Vector &AV) const{
                return !(operator ==(AV));
            }
            
            /*------------------------------------------------------------------------*/
            /** \brief  Destructor.
             */
            virtual ~Vector();
            
            /*------------------------------------------------------------------------*/
            /** \brief normalize
             */
            inline void normalize(){
                TCoord n = norm();
                if (n != 0.0){
                    m_tab[0] = m_tab[0] / n;
                    m_tab[1] = m_tab[1] / n;
                    m_tab[2] = m_tab[2] / n;
                }
            }
            
            Vector getNormalize() const{
                Vector v(*this);
                v.normalize();
                return v;
            }
            
            /*------------------------------------------------------------------------*/
            /** \brief update 1 component
             */
            inline void set(const int AIndex, const TCoord AVal)
            {
                if (AIndex >= 0 && AIndex < 3)
                    m_tab[AIndex] = AVal;
            }
            inline TCoord get(const int AIndex) const
            {
                
                return m_tab[AIndex];
            }
            
            /*------------------------------------------------------------------------*/
            /** \brief update X component
             */
            inline void setX(const TCoord AVal)
            {
                m_tab[0] = AVal;
            }
            /*------------------------------------------------------------------------*/
            /** \brief update Y component
             */
            inline void setY(const TCoord AVal)
            {
                m_tab[1] = AVal;
            }
            /*------------------------------------------------------------------------*/
            /** \brief update Z component
             */
            inline void setZ(const TCoord AVal)
            {
                m_tab[2] = AVal;
            }
            /*------------------------------------------------------------------------*/
            /** \brief compute the square norm of the vector
             */
            inline TCoord norm2() const
            {
                return m_tab[0] * m_tab[0] + m_tab[1] * m_tab[1] + m_tab[2] * m_tab[2];
            }
            
            /*------------------------------------------------------------------------*/
            /** \brief compute the L2 norm of the vector
             */
            inline TCoord norm()   const { return sqrt(norm2()); }
            inline TCoord normL2() const { return sqrt(norm2()); }
            
            /*------------------------------------------------------------------------*/
            /** \brief compute the L1 norm of the vector
             */
            inline TCoord normL1() const {
                return fabs(m_tab[0]) + fabs(m_tab[1]) + fabs(m_tab[2]);
            }
            
            /*------------------------------------------------------------------------*/
            /** \brief compute the Linf norm of the vector
             */
            inline TCoord normLinf() const{
                TCoord v0 = fabs(m_tab[0]);
                TCoord v1 = fabs(m_tab[1]);
                TCoord v2 = fabs(m_tab[2]);
                if (v0 > v1 && v0 > v2)
                    return v0;
                else if (v1 > v2 && v1 > v0)
                    return v1;
                else
                    return v2;
            }
            
            /*------------------------------------------------------------------------*/
            /** \brief compute the Lp norm of the vector
             *
             *  \param AP the value of p in Lp
             */
            inline TCoord normLp(const int& AP) const{
                return pow(pow(m_tab[0], AP) + pow(m_tab[1], AP) + pow(m_tab[2], AP), AP);
            }
            /*------------------------------------------------------------------------*/
            /** \brief compute the angle between two vectors from this to AV
             *
             *  \param AV a vector
             */
            inline TCoord angle(const Vector& AV) const{
                Vector v1(AV);
                v1.normalize();
                Vector v2(*this);
                v2.normalize();
                if(v1.isZero() || v2.isZero()) {
                    throw GMDSException("Vector::angle one of the vector is zero.");
                }
                double dotProduct = v2.dot(v1);
                
                if(dotProduct>1.) {
                    return 0.;
                } else if(dotProduct<-1.) {
                    return gmds::math::Constants::PI;
                }
                return acos(dotProduct);
            }
            
            /*------------------------------------------------------------------------*/
            /** \brief compute the oriented angle from this to AV. This operation
             *         requires to give a reference vector, that must be orthogonal to
             *         the plane containing (*this) and AV. By default, we put it
             *         equals to Oz(0,0,1)
             *
             *  \param AV a vector
             */
            inline TCoord orientedAngle(const Vector& AV,
                                        const Vector& AOrtho=Vector(0,0,1)) const
            {
                Vector ref = cross(AOrtho);
                double a = (AV.dot(ref)>0)?-angle(AV): angle(AV);
                return a;
            }
            
            /*------------------------------------------------------------------------*/
            /** \brief compute the oriented angle from this to AV. This operation
             *         requires to give a reference vector, that must be orthogonal to
             *         the plane containing (*this) and AV. By default, we put it
             *         equals to Oz(0,0,1)
             *
             *  \param AV a vector
             */
            inline TCoord angleIn02PI(const Vector& AV,
                                      const Vector& AOrtho=Vector(0,0,1)) const
            {
                Vector ref = cross(AOrtho);
                double oriented_angle = orientedAngle(AV,AOrtho);
                double a=0;
                if(oriented_angle >0.0)
                    a = oriented_angle;
                else
                    a = Constants::PI2+oriented_angle;
                
                return a;
            }
            
            
            /*------------------------------------------------------------------------*/
            /** \brief provides the max absolute value component index
             */
            inline	TInt getMaxAbsComponentIndex() const {
                TInt index = 0;
                for (int i = 1; i < 3; ++i) {
                    if (fabs(m_tab[i]) > fabs(m_tab[index])) {
                        index = i;
                    }
                }
                return index;
            }
            /*------------------------------------------------------------------------*/
            inline TCoord const& operator[](const TInt& i) const{
                return m_tab[i];
            }
            /*------------------------------------------------------------------------*/
            inline TCoord& operator[](const TInt& i) {
                return m_tab[i];
            }
            
            /*------------------------------------------------------------------------*/
            /** \brief Return if the vector is zero
             */
            bool isZero() const {
                if (m_tab[0] == 0. && m_tab[1] == 0. && m_tab[2] == 0.) {
                    return true;
                }
                else {
                    return false;
                }
            }
            
            /*------------------------------------------------------------------------*/
            /** \brief  Dot product
             *
             *  \param AV a second vector
             *
             *  \return this.AV
             */
            inline TCoord dot(const Vector& AV) const {
                return m_tab[0] * AV.m_tab[0] + m_tab[1] * AV.m_tab[1] + m_tab[2] * AV.m_tab[2];
            }
            
            /*------------------------------------------------------------------------*/
            /** \brief  Cross product
             *
             *  \param AV a second vector
             *
             *  \return this x AV
             */
            inline Vector cross(const Vector& AV) const {
                return Vector( m_tab[1] * AV.m_tab[2] - m_tab[2] * AV.m_tab[1],
                              -m_tab[0] * AV.m_tab[2] + m_tab[2] * AV.m_tab[0],
                              m_tab[0] * AV.m_tab[1] - m_tab[1] * AV.m_tab[0]);
            }
            /*------------------------------------------------------------------------*/
            /** \brief  opposite of a vector
             *
             *  \return (-x,-y,-z)
             */
            inline Vector opp() const {
                return Vector(-m_tab[0], -m_tab[1], -m_tab[2]);
            }
            
            /*------------------------------------------------------------------------*/
            /** \brief sum of all the components of the vector
             */
            TCoord sumComponents() const;
            
            /*------------------------------------------------------------------------*/
            /** \brief  indicates if two vector are colinear
             *
             *  \param AV a second vector
             *
             *  \return true if vector are colinear, no otherwise
             */
            bool isColinear(const Vector&) const;
            
            /*------------------------------------------------------------------------*/
            /** \brief  indicates if two vector are orthogonal
             *
             *  \param AV a second vector
             *
             *  \return GEOM_YES if vector are orthogonal, GEOM_NO if they are
             *  		not, GEOM_UNDEF otherwise
             */
            bool isOrthogonal(const Vector&) const;
            
            /*------------------------------------------------------------------------*/
            /** \brief operator power
             */
            inline Vector operator^(const int AN){
                return Vector(pow(m_tab[0], AN), pow(m_tab[1], AN), pow(m_tab[2], AN));
            }
            
            /*------------------------------------------------------------------------*/
            /** \brief  Overloaded operator- to get the difference between two vectors
             */
            friend EXPORT_GMDS  Vector operator-(const Vector&, const Vector&);
            
            /*------------------------------------------------------------------------*/
            /** \brief  Overloaded operator+ to get the sum of two vectors
             */
            friend EXPORT_GMDS Vector operator+(const Vector&, const Vector&);
            
            /*------------------------------------------------------------------------*/
            /** \brief  Overloaded operator+ to create a new point from a point and a
             * 			vector
             */
            friend EXPORT_GMDS Point operator+(const Point&, const Vector&);
            
            /*------------------------------------------------------------------------*/
            /** \brief  Overloaded operator- to create a new point from a point and a
             * 			vector
             */
            friend EXPORT_GMDS Point operator-(const Point&, const Vector&);
            /*------------------------------------------------------------------------*/
            /** \brief  Overloaded operator* the product of a vector by a scalar value
             */
            friend EXPORT_GMDS Vector operator*(const TCoord&, const Vector&);
            friend EXPORT_GMDS Vector operator*(const Vector&, const TCoord&);
            
            /*------------------------------------------------------------------------*/
            /** \brief  Overloaded operator/ of a vector by a scalar value
             */
            friend	EXPORT_GMDS Vector operator/(const Vector&, const TCoord&);
            
            /*------------------------------------------------------------------------*/
            /** \brief  Overloaded operator<< for output
             */
            friend EXPORT_GMDS  std::ostream& operator<<(std::ostream&, const Vector&);
            
        protected:
            TCoord m_tab[3];
        };
        /*----------------------------------------------------------------------------*/
    } // namespace math
    /*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_VECTOR_H_ */
/*----------------------------------------------------------------------------*/
