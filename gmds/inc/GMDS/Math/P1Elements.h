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
 * P1ShapeFunctions.h
 *
 *  Created on: April 19, 2016
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_P1_ELEMENTS_H_
#define GMDS_MATH_P1_ELEMENTS_H_
/*----------------------------------------------------------------------------*/
// gmds file headers
/*----------------------------------------------------------------------------*/
#include <GMDS/Math/Matrix.h>
#include <GMDS/Math/VectorND.h>
#include <GMDS/Math/Point.h>
/*----------------------------------------------------------------------------*/
#include<cmath>
#include<string.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
    /*------------------------------------------------------------------------*/
    namespace math{
        
        class ShapeFunction{
            
            /*---------------------------------------------------------------*/
            /** \brief gives the value of the shape function for the point
             *         (\p AX, \p AY, \p AZ) living in the physical space.
             *
             * \param[in] AX X coordinate of the physical point
             * \param[in] AY Y coordinate of the physical point
             * \param[in] AZ Z coordinate of the physical point, default value
                          id 0 for the 2D case
             */
            virtual double value(const double AX,
                                 const double AY,
                                 const double AZ=0) const =0;
            
            virtual Vector3d grad(const double AX,
                                  const double AY,
                                  const double AZ=0) const =0;
            
public:
            virtual ~ShapeFunction(){;};
            
        };
        
        class ShapeFunctionP1Tetrahedron :public ShapeFunction {
          public:
           ~ShapeFunctionP1Tetrahedron() {};
 
        };
        class ShapeFunctionP1Triangle2D :public ShapeFunction {
          public:
           ~ShapeFunctionP1Triangle2D() {};    
        };
        
        /*--------------------------------------------------------------------*/
        /** \class P1Element
         *  \brief Describe the services provided by a P1 Element built from
         *         a mesh cell using usual reference elements to perform
         *         computations. In a P1 element, points live in the physical
         *         elements while nodes live in the reference element.
         */
        class EXPORT_GMDS P1Element {
        protected:
            
            /*---------------------------------------------------------------*/
            /** \brief Gives access to the \p AI-th shape function of this 
             *         element
             *
             * \param[in] APnts the quadrature points defining the element
             */
            ShapeFunction* shapeFunction(const int AI);
            /*---------------------------------------------------------------*/
            /** \brief Constructor.
             *
             * \param[in] APnts the quadrature points defining the element
             */
            P1Element(const std::vector<Point>& APnts);
            
        public:
            /*---------------------------------------------------------------*/
            /** \brief Destructor.
             */
            virtual ~P1Element();

            Matrix<3,3,double> stiffness();
            Matrix<3,3,double> mass();

            /*---------------------------------------------------------------*/
            /** \brief Set the tolerance used for some computation. The default
             *         value is 1.e-6
             *
             * \param[in] ATol the new tolerance to use
             */
            static void setTolerance (const double tol);
            /*---------------------------------------------------------------*/
            /** \brief Get the tolerance used for computations on P1 elements
             */
            static double getTolerance ();
            /*---------------------------------------------------------------*/
            /** \brief Acces to point \p AIndex of the element.
             *
             * \param[in] AIndex the index of the point we want to access to
             */
            Point getPoint(const int AIndex) const;
            /*---------------------------------------------------------------*/
            /** \brief Acces to dimension of the element.
             */
            virtual int getDimension() = 0;
            /*---------------------------------------------------------------*/
            /** \brief Acces to the number of nodes/points of the element
             */
            virtual int getNbPoints() = 0;
            
            /*---------------------------------------------------------------*/
            /** \brief Return the (u,v,w) values of the node in the reference
             *         element corresponding to the AIndex point of the element.
             *
             * \param[in] AIndex the index of the point we want to access to
             *
             * \return the vector of (u,v,w) cooordinates
             */
            virtual Vector3d getNode(const int AIndex) const = 0;
            
            virtual void getShapeFunction(int num, double u, double v, double w, double &s) = 0;
            virtual void getGradShapeFunction(int num, double u, double v, double w, double s[3]) = 0;

        protected:
            // points are stored by components
            double *m_x, *m_y, *m_z;
            static double m_tolerance;
        public:

            
             };

        /*----------------------------------------------------------------------------*/
    } // namespace math
    /*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_P1_ELEMENTS_H_ */
/*----------------------------------------------------------------------------*/
