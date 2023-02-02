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
/*-----------------------------------------------------------------*/
/*
 * Chart.h
 *
 *  Created on: 03/19/2015
 *      Author: ledouxf
 */
/*-----------------------------------------------------------------*/
#ifndef GMDS_MATH_CHART_H_
#define GMDS_MATH_CHART_H_
/*-----------------------------------------------------------------*/
#include <vector>
/*-----------------------------------------------------------------*/
#include "GMDS/Math/Vector.h"
#include "GMDS/Math/VectorND.h"
#include "GMDS/Math/Matrix.h"
/*-----------------------------------------------------------------*/
#include <iostream>
/*-----------------------------------------------------------------*/
namespace gmds{
    /*-----------------------------------------------------------------*/
    namespace math{
        class Quaternion;
        /*------------------------------------------------------------*/
        /* \class Chart
         * \brief A char is a right-handed basis in R3
         */
        class EXPORT_GMDS Chart
        {
        public:
            /*------------------------------------------------------------*/
            /* \brief Defaut constructor
             *  Defines a chart with (1,0,0), (0,1,0), (0,0,1)
             */
            Chart();
            /*------------------------------------------------------------*/
            /* \brief Constructor. Warning there is no check about the parameter
             * values
             */
            Chart(const Vector& AX, const Vector& AY,const Vector& AZ);
            
            
            /*------------------------------------------------------------*/
            /* \brief Constructor. Warning there is no check about the parameter
             * values
             */
            Chart(const Vector3d& AX, const Vector3d& AY,const Vector3d& AZ);
            /*------------------------------------------------------------*/
            /* \brief Constructor. Warning there is no check about the parameter
             * values. Last vector is computed on-the-fly as AX x AY
             */
            Chart(const Vector3d& AX, const Vector3d& AY);
            Chart(const Vector&   AX, const Vector&   AY);
            
            /*------------------------------------------------------------*/
            /*  \brief A constructor of chart from quaternion.
             * This constructor is used for conversion purposes.
             */
            explicit Chart(const Quaternion& AQ);
            
            Vector X() const { return Vector(m_v[0].X(),m_v[0].Y(),m_v[0].Z());}
            Vector Y() const { return Vector(m_v[1].X(),m_v[1].Y(),m_v[1].Z());}
            Vector Z() const { return Vector(m_v[2].X(),m_v[2].Y(),m_v[2].Z());}
            Vector get(const int AIndex) const;
            
            Vector3d VX() const { return m_v[0];}
            Vector3d VY() const { return m_v[1];}
            Vector3d VZ() const { return m_v[2];}
            
            Vector3d operator[](const int AIndex) const;
            /*------------------------------------------------------------*/
            /*  \brief Return the six vector corresponding to the chirral
             *  rotation group of this
             */
            std::vector<Vector> image()const;
            
            /*------------------------------------------------------------*/
            /*  \brief Return the 3x3 matrix whose each column is one chart
             *         vector
             *
             * \return a matrix view of *this
             */
            Matrix<3,3,double> toMatrix() const;
            
            /*------------------------------------------------------------*/
            /*  \brief Return the 3x3 rotation matrix that transforms 
             *         *this to chart \p AC
             *
             * \param[in] AC a chart we want to get the rotation matrix to
             *
             * \return a rotation matrix bringing $this to \p AC
             */
            Matrix<3,3,double> computeRotationTo(const Chart& AC) const;
            

            /*------------------------------------------------------------*/
            /*  \brief align the charyt with AV following the minimal rotation
             */
            void align(const Vector& AV);
            
            /*------------------------------------------------------------*/
            /*  \brief Give the matching of the vectors between this and
             *         AChart. AMatching[i]=j means the i^th vector of this
             *         matchs the j^th vector of AChart. By matching, we
             *         mean the most aligned with (max absolute dot product)
             *
             * \param[in] AChart a chart to compare with
             * \param[in] AMatching a 3-size tabular storing the computed matching
             */
            void matchVectors(const Chart&AChart, int (&AMatching)[3]) const;
            
            /*------------------------------------------------------------*/
            /* \brief Nested class storing the mapping from one chart to
             *        another one in a compact manner.
             * 
             * \details Considering two charts C1 and C2, a mapping object
             *          indicates which vectors of C2 corresponds to each
             *          vector of C1
             */
            class Mapping {
            public:
                /*-------------------------------------------------------*/
                /** \brief Default constructor providing the identity 
                 *         mapping
                 */
                Mapping();
                        
                /*-------------------------------------------------------*/
                /** \brief Constructs a new mapping from \p AC1 to \p AC2
                 *
                 * \param[in] AC1 an origin chart
                 * \param[in] AC2 a destination chart
                 */
                Mapping(const Chart& AC1, const Chart& AC2);
                
                /*--------------------------------------------------------*/
                /** \brief Gets the direction vector
                 *
                 * \return a integer vector
                 */
                
                const Vector3i& getDirections() const ;
                /*--------------------------------------------------------*/
                /** \brief Gets the index permutation vector
                 *
                 * \return a integer vector
                 */
                const Vector3i& getPermutations() const;
                
                /*-------------------------------------------------------*/
                /** \brief Provide the inverse mapping
                 *
                 * \return the inverse mapping of *this
                 */
                Mapping inverse() const;
                /*-------------------------------------------------------*/
                /** \brief Check if we have the identity mapping
                 */
                bool isIdentity() const;
                
                /*-------------------------------------------------------*/
                /** \brief Computes the composition of two mappings, the
                 *         resulting mapping consists in applyng \p AM 
                 *         then (*this)
                 *
                 * \param[in] AM another mapping
                 * \return the mapping (*this)*\p AM
                 */
                Mapping operator*(const Mapping& AM) const;
                /*-------------------------------------------------------*/
                /** \brief Transport a solution vector by *this. This 
                 *         operation must be understand in the context of
                 *         a chart changement. If a vector \p AV is 
                 *         expressed in a chart c1 and (*this) is a mapping
                 *         from c1 to c2 then (*this)*\p AV is the 
                 *         expression of \p AV in the basis defined by c2.
                 *
                 * \param[in] AV a vector
                 *
                 * \return the vector \p AV expressed in *this
                 */
                Vector3d operator*(const Vector3d& AV) const;
                
                /*------------------------------------------------------------*/
                /*  \brief Return the 3x3 rotation matrix corresponding to
                 *         this mapping
                 *
                 * \return a rotation matrix made of 0 and (-)1
                 */
                Matrix<3,3,double> toMatrix() const;
            public: //WARNING TMP
                Vector3i m_map;
                Vector3i m_dir;
            } ;
            
        private:
            /** three orthogonal vectors */
            Vector3d m_v[3];
            
        };
        /*-----------------------------------------------------------------*/
        EXPORT_GMDS  std::ostream & operator<<(std::ostream & op_g,
                                               const Chart & op_d);
        /*-----------------------------------------------------------------*/
    }
}
/*-----------------------------------------------------------------*/
#endif /* GMDS_MATH_CHART_H_ */
/*-----------------------------------------------------------------*/
