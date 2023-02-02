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
 * AxisAngleRotation
 *
 *  Created on: December 2, 2015
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_AXIS_ANGLE_ROTATION_H_
#define GMDS_MATH_AXIS_ANGLE_ROTATION_H_
/*----------------------------------------------------------------------------*/
// gmds file headers
/*----------------------------------------------------------------------------*/
#include <GMDS/Math/Matrix.h>
#include <GMDS/Math/Quaternion.h>
#include <GMDS/Math/Chart.h>
/*----------------------------------------------------------------------------*/
#include<cmath>
#include<string.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
    /*------------------------------------------------------------------------*/
    namespace math{
        /*--------------------------------------------------------------------*/
        /** \class Vector
         *  \brief Defines a 3D Vector
         */
        class EXPORT_GMDS AxisAngleRotation {
        public:
            
            /*---------------------------------------------------------------*/
            /** \brief Default Constructor with null vector and null angle
             */
            AxisAngleRotation();
            
            
            /*---------------------------------------------------------------*/
            /** \brief Constructor.
             *
             * \param[in] AV the rotation vector
             * \param[in] AA the rotation angle
             */
            AxisAngleRotation(const Vector3d& AV, const double AA);
            /*---------------------------------------------------------------*/
            /** \brief Constructor.
             *
             * \param[in] AV a vector defining the rotation axis and the angle (via
             *           its norm
             */
            AxisAngleRotation(const Vector3d& AV);
            
            /*---------------------------------------------------------------*/
            /** \brief Constructor from a quaternion
             *
             * \param AQ a unit quaternion corresponding to a rotation
             */
            AxisAngleRotation(const Quaternion& AQ);
            
            
            
            /*---------------------------------------------------------------*/
            /** \brief Constructor from a chart \p AC. It gives the rotation
             *         bringing (OX,OY,OZ) onto the (X,Y,Z) axis of \p AC.
             *
             * \param[in] AC a chart
             */
            AxisAngleRotation(const Chart& AC);
            
            
            /*---------------------------------------------------------------*/
            /**  \brief Constructs an axis-angle rotation which transforms 
             *          \p AFrom axis into \p ATo axis
             *
             * \param[in] AFrom a first vector
             * \param[in] ATo   a second vector
             */
            AxisAngleRotation(const Vector3d& AFrom, const Vector3d& ATo);
            
            /*---------------------------------------------------------------*/
            /** \brief Provides a quaternion representation of this rotation.
             */
            Quaternion quaternion() const;
            
            /*---------------------------------------------------------------*/
            /** \brief Returns the rotation vector
             */
            const Vector3d& axis() const {return m_axis;}
            Vector3d& axis() {return m_axis;}
            
            
            /*---------------------------------------------------------------*/
            /** \brief Returns the rotation angle
             */
            const double angle() const {return m_axis.norm();}
            
            /*---------------------------------------------------------------*/
            /** \brief Build the inverse rotation, that is the an axis-angle
             *         rotation with opposite rotation angle */
            AxisAngleRotation inverse() const
            { return AxisAngleRotation(m_axis,m_axis.norm()); }
            
            /*---------------------------------------------------------------*/
            /** \brief Build the identity axis-angle rotation */
            AxisAngleRotation identity() const
            { return AxisAngleRotation(Vector3d(1,0,0),0); }
            
            /*---------------------------------------------------------------*/
            /** \brief Build a corresponding 3x3 rotation matrix */
            Matrix<3,3,double> toRotationMatrix() const;
            Vector3d toRotationAxis(const int AIndex) const;
            
            /*---------------------------------------------------------------*/
            /** \brief Build the chart corresponding to applying *this to
             *         the reference Chart(OX, OY,OZ)*/
            Chart toChart() const;
            /*---------------------------------------------------------------*/
            /**  \brief Gets an axis-angle rotation which transforms Z axis into
             *          \p AV axis
             *
             * \param[in]  AV the constrained axis
             * \return an axis-angle rotation
             */
            static AxisAngleRotation alignZ(const Vector3d& AV);
            /*---------------------------------------------------------------*/
            /**  \brief Gets an axis-angle rotation which transforms Y and Z 
             *          axis into \p AV1 and \p AVE respectively
             *
             * \param[in]  AY the axis Y must be aligned with
             * \param[in]  AZ the axis Z must be aligned with
             * \return an axis-angle rotation
             */

            static AxisAngleRotation alignYZ(const Vector3d& AY,
                                             const Vector3d& AZ);

        protected:
            /** rotation axis */
            Vector3d m_axis;
        };
        
        Vector3d operator*(const AxisAngleRotation& AR,
                           const Vector3d &AV);
        AxisAngleRotation operator*(const AxisAngleRotation& AR0,
                                    const AxisAngleRotation& AR1);
        /*----------------------------------------------------------------------------*/
    } // namespace math
    /*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_AXIS_ANGLE_ROTATION_H_ */
/*----------------------------------------------------------------------------*/
