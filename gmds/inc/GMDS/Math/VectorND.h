/*------------------------------------------------------------------------*/
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
/*------------------------------------------------------------------------*/
//  VectorND.h
//
//  GMDSSuite
//
//  Created by F. Ledoux on 02/11/2015. Inspirated by the generic
//  vectors of Eigen and and Geogram (INRIA)
//
//
/*------------------------------------------------------------------------*/
#ifndef GMDS_MATH_VECTOR_ND_H_
#define GMDS_MATH_VECTOR_ND_H_
/*------------------------------------------------------------------------*/
// STL Headers
#include <cassert>
/*------------------------------------------------------------------------*/
// GMDS Headers
#include <GMDS/Math/Vector.h>
/*------------------------------------------------------------------------*/
namespace gmds{
    /*--------------------------------------------------------------------*/
    namespace math{
        /*--------------------------------------------------------------------*/
        /** \struct CrossNDPolicy
         *  \brief  Template structure providing a tailored process for
         *          computing the cross product between two N-dim vectors
         *
         *  \tparam TDim	vector dimension
         *  \tparam TType	value type
         */
        /*-------------------------------------------------------------------*/
        template<int TDim, typename TType> struct CrossNDPolicy {
            static void perform(const TType t1[TDim],
                                const TType t2[TDim],
                                TType(&i)[TDim]) {
                throw GMDSException("Not yet implemented");
            }

        };
        /*-------------------------------------------------------------------*/
        template<typename TType> struct CrossNDPolicy<3,TType>  {
            static void perform(const TType t1[3],
                                const TType t2[3],
                                TType(&i)[3])
            {
                i[0] = t1[1] * t2[2] - t1[2] * t2[1];
                i[1] =-t1[0] * t2[2] + t1[2] * t2[0];
                i[2] = t1[0] * t2[1] - t1[1] * t2[0];
            }
        };

        /*----------------------------------------------------------------*/
        /**
         * \brief Mathematical vector of dimension \p TDim and where each
         *        component is of type \p TType. Usual mathematical
         *        operations are provided to handle them in a simple manner.
         *
         * \tparam TDim  vector dimension
         * \tparam TType component type
         */
        template <int TDim, class TType> class VectorND {
            
        public:
            
            /** alias on this vector type */
            typedef VectorND<TDim, TType> vector_type;
            /** alias on the type of components */
            typedef TType value_type;
            /** alias on the vector dimension */
            static const int dimension=TDim;
            
            /*-----------------------------------------------------------*/
            /** \brief Default constructor (all component are zero)
             */
            VectorND() {
                for(int i=0; i<dimension; i++) {
                    m_data[i] = value_type(0);
                }
            }
            
            /*-----------------------------------------------------------*/
            /** \brief Constructor from a C-style array \p ATab. The first
             *         dimension component of \p ATab will be copied into
             *         *this. Size of \p ATab must so be greater than
             *         dimension (no check). Component of ATab ar of \p T 
             *         type, whic must can be cast into value_type.
             *
             * \tparam T the type of components in \p v
             *
             * \param[in] ATab an array of \p T - type values
             */
            template <class T> explicit VectorND(const T* ATab) {
                for(int i = 0; i < dimension; i++) {
                    m_data[i] = value_type(ATab[i]);
                }
            }
            /*-----------------------------------------------------------*/
            /** \brief overloaded operator=.
             *
             * \param[in]  AV a same type generic vector
             * \return     a reference onto *this
             */
            vector_type& operator= (const vector_type& AV) {
                memcpy(m_data, AV.data(), dimension * sizeof(value_type));
                return *this;
            }
            
            /*-----------------------------------------------------------*/
            /** \brief Gets read-write access to the vector data
             *
             * \return a pointer onto the first item of the vector
             */
            value_type* data() {
                return m_data;
            }
            
            
            /*-----------------------------------------------------------*/
            /** \brief Gets read-only access to the vector data
             *
             * \return a pointer onto the first item of the vector
             */
            const value_type* data() const {
                return m_data;
            }
            /*-----------------------------------------------------------*/
            /** \brief Gets read-write access to a vector component
             * 
             * \param[in]  AI component index
             * \return     a reference onto the \p i th component
             */
             value_type& operator[] (const int AI) {
                assert(AI>=0 && AI < dimension);
                return m_data[AI];
            }
            
            /*-----------------------------------------------------------*/
            /** \brief Gets read-only access to a vector component
             *
             * \param[in]  AI component index
             * \return     a reference onto the \p i th component
             */
            const value_type& operator[] (const int AI) const {
                assert(AI>=0 && AI < dimension);
                return m_data[AI];
            }
            
            /*-----------------------------------------------------------*/
            /** \brief Gets the vector dimension
             * \return the value of \p TDIM
             */
            int dim() const {
                return dimension;
            }
            /*-----------------------------------------------------------*/
            /** \brief Gets the squared L2 norm of *this
             */
            value_type norm2() const {
                value_type r = value_type(0.);
                for(int i = 0; i < dimension; i++) {
                    r += m_data[i] * m_data[i];
                }
                return r;
            }
            
            /*-----------------------------------------------------------*/
            /** \brief Gets the  L2 norm of *this
             */
            value_type norm() const {
                return sqrt(norm2());
            }
            
            /*-----------------------------------------------------------*/
            /** \brief Normalizes this vector if its norm is not zero
             */
            void normalize() {
                value_type n = norm();
                if(n>1e-24) {
                    for(int i = 0; i < dimension; i++) {
                        m_data[i] /= n;
                    }
                }
            }

            /*-----------------------------------------------------------*/
            /** \brief Adds vector \p AV to *this
             *
             * \param[in] AV a same type vector
             * \return    a reference to *this
             */
            vector_type& operator+= (const vector_type& AV) {
                for(int i = 0; i < dimension; i++) {
                    m_data[i] += AV.m_data[i];
                }
                return *this;
            }
            /*-----------------------------------------------------------*/
            /** \brief Substracts vector \p AV to *this
             *
             * \param[in] AV a same type vector
             * \return    a reference to *this
             */
            vector_type& operator-= (const vector_type& AV) {
                for(int i = 0; i < dimension; i++) {
                    m_data[i] -= AV.m_data[i];
                }
                return *this;
            }
            
            /*-----------------------------------------------------------*/
            /** \brief Multiplies all component of *this by the scalar
             *         \p AS
             *
             * \tparam    T  scalar type compatible with AType
             * \param[in] AS a \p T-type calar
             * \return    a reference to *this
             */
            template <class T>
            vector_type& operator*= (const T AS) {
                for(int i = 0; i < dimension; i++) {
                    m_data[i] *= value_type(AS);
                }
                return *this;
            }
            
            /*-----------------------------------------------------------*/
            /** \brief Divides all component of *this by the scalar \p AS
             *
             * \tparam    T  scalar type compatible with AType
             * \param[in] AS a \p T-type calar
             * \return    a reference to *this
             */
            template <class T>
            vector_type& operator/= (const T AS) {
                for(int i = 0; i < dimension; i++) {
                    m_data[i] /= value_type(AS);
                }
                return *this;
            }
            
            /*-----------------------------------------------------------*/
            /** \brief Addition of 2 vectors
             *
             * \param[in]  AV another vector
             * \return     Vector \p *this + \p AV
             */
            vector_type operator+ (const vector_type& AV) const {
                vector_type r(*this);
                for(int i = 0; i < dimension; i++) {
                    r.m_data[i] += AV.m_data[i];
                }
                return r;
            }
            /*-----------------------------------------------------------*/
            /** \brief Difference of 2 vectors
             *
             * \param[in]  AV another vector
             * \return     Vector \p *this - \p AV
             */
            vector_type operator- (const vector_type& AV) const {
                vector_type r(*this);
                for(int i = 0; i < dimension; i++) {
                    r.m_data[i] -= AV.m_data[i];
                }
                return r;
            }
            
            /*-----------------------------------------------------------*/
            /** \brief Multiplies a vector \p AV by a scalar \p AS (of a
             *         compatible type \p T)
             *
             * \tparam    T scalar type
             * \param[in] AS a scalar value
             * \return    Vector \p *this * \p AS
             */
            template <class T>
            vector_type operator* (const T AS) const {
                vector_type r(*this);
                for(int i = 0; i < dimension; i++) {
                    r.m_data[i] *= value_type(AS);
                }
                return r;
            }
            /*-----------------------------------------------------------*/
            /** \brief Divides a vector \p AV by a scalar \p AS (of a
             *         compatible type \p T)
             *
             * \tparam    T scalar type
             * \param[in] AS a scalar value
             * \return    Vector \p *this / \p AS
             */
            template <class T>
            vector_type operator/ (const T AS) const {
                vector_type r(*this);
                for(int i = 0; i < dimension; i++) {
                    r.m_data[i] /= value_type(AS);
                }
                return r;
            }
            
            
            /*-----------------------------------------------------------*/
            /** \brief Provides the opposite vector
             *
             * \return Vector -\p *this
             */
            vector_type operator- () const {
                vector_type r;
                for(int i = 0; i < dimension; i++) {
                    r.m_data[i] = -m_data[i];
                }
                return r;
            }
            /*-----------------------------------------------------------*/
            /** \brief Gets the dot product of this vector with \p AV
             *
             * \param[in] AV a vector
             * \return the dot product (\p *this . \p AV)
             */
            value_type  dot( const vector_type& AV) const{
                value_type r = 0;
                for(int i = 0; i < dimension; i++) {
                    r += m_data[i] * AV.m_data[i];
                }
                return r;
            }
            /*-----------------------------------------------------------*/
            /** \brief Gets the cross product of this vector with \p AV
             *
             * \param[in] AV a vector
             * \return the cross product (\p *this x \p AV)
             */
            vector_type  cross( const vector_type& AV) const{
                
                value_type c_tab[dimension];
                CrossNDPolicy<dimension, value_type>::perform(m_data,
                                                              AV.m_data,
                                                              c_tab);
                
                return VectorND<dimension,value_type>(c_tab);
            }
            /*-----------------------------------------------------------*/
            /** \brief Compute the angle between 0 and 360 degrees
             *
             * \param[in] AV a vector
             * \return the angle between \p *this and \p AV
             */
            value_type angle( const vector_type& AV) const{
                
                vector_type v1 = *this;
                vector_type v2 = AV;
                v1.normalize();
                v2.normalize();
                value_type d = v1.dot(v2);
                value_type c = v1.cross(v2).norm();
                return std::atan2(c, d)*math::Constants::INVPIDIV180;
            }
            /*-----------------------------------------------------------*/
            /** \brief Indicate if all the component of *this are 0.0
             * \return true if zero, false otherwise
             */
            bool isNull() const {
                for(int i=0; i<TDim; i++){
                    if (m_data[i]!=0.0)
                        return false;
                }
                return true;
            }

        protected:
            TType m_data[TDim];

            
        }; //class VectorND
        
        
        /*--------------------------------------------------------------*/
        /** \brief Multiplies a scalar by a vector like done in the class
         *         VectorND
         *
         * \tparam TScalar scalar type
         * \tparam TDim    vector dimension
         * \tparam TVec    vector component type
         *
         * \param[in] AS a scalar value of type \p TScalar
         * \param[in] AV a vector of type dimension \p TDim and component
         *                 type \p TVec
         *
         * \return Vector \p AS * \p AV
         */
        template <typename TScalar, int TDim, typename TVec>
        VectorND<TDim, TVec> operator*(const TScalar AS,
                                            const VectorND<TDim, TVec>& AV)
        {
            VectorND<TDim, TVec> r;
            for(int i = 0; i < TDim; i++) {
                r[i] = TVec(AS) * AV[i];
            }
            return r;
        }
        
        /*--------------------------------------------------------------*/
        /** \brief Specific implementations of VectorND
         */
        class Vector2d: public VectorND<2, double>{
        public:
            /*-----------------------------------------------------------*/
            /** \brief Specific constructor from 2 components
             *
             * \param AX first  component
             * \param AY second component
             */
            Vector2d(const double& AX, const double& AY)
            :VectorND<2, double>()
            {
                m_data[0] = AX;
                m_data[1] = AY;
            }
            
            /*-----------------------------------------------------------*/
            /** \brief Gets read-only access to the X component
             *
             * \return X Component
             */
            const double& X() const { return m_data[0];}
            
            /*-----------------------------------------------------------*/
            /** \brief Gets read-only access to the Y component
             *
             * \return Y Component
             */
            const double& Y() const { return m_data[1];}
        };

        class Vector3d: public VectorND<3, double>{
        public:
            
            /*-----------------------------------------------------------*/
            /** \brief Default Constructor
             *
             * \param AV a vector object
             */
            Vector3d():VectorND<3, double>(){;}

            /*-----------------------------------------------------------*/
            /** \brief Specific constructor from 3 components
             *
             * \param AX first  component
             * \param AY second component
             * \param AZ third  component
             */
            Vector3d(const double& AX, const double& AY, const double&AZ)
            :VectorND<3, double>()
            {
                m_data[0] = AX;
                m_data[1] = AY;
                m_data[2] = AZ;
            }
            /*-----------------------------------------------------------*/
            /** \brief Constructor from deprecated Vector class
             *
             * \param AV a vector object
             */
            Vector3d(const Vector& AV):VectorND<3, double>(){
                m_data[0] = AV.X();
                m_data[1] = AV.Y();
                m_data[2] = AV.Z();
            }
            
            /*-----------------------------------------------------------*/
            /** \brief Constructor the vector going from \p AP1 to \p AP2.
             *
             * \param AP1 first point
             * \param AP2 second point
             */
            Vector3d(const Point& AP1, const Point& AP2):VectorND<3, double>(){
                m_data[0] = AP2.X() -AP1.X();
                m_data[1] = AP2.Y() -AP1.Y();
                m_data[2] = AP2.Z() -AP1.Z();
            }
            
            /*-----------------------------------------------------------*/
            /** \brief Constructor from mother class
             *
             * \param AV a vector object
             */
            Vector3d(const VectorND<3, double>& AV)
            :VectorND<3, double>(AV){;}
            
            /*-----------------------------------------------------------*/
            /** \brief Gets read-only access to the X component
             *
             * \return X Component
             */
            const double& X() const { return m_data[0];}
            
            /*-----------------------------------------------------------*/
            /** \brief Gets read-only access to the Y component
             *
             * \return Y Component
             */
            const double& Y() const { return m_data[1];}
            
            /*-----------------------------------------------------------*/
            /** \brief Gets read-only access to the Z component
             *
             * \return Z Component
             */
            const double& Z() const { return m_data[2];}
            
        };
        
        class Vector4d: public VectorND<4, double>{
        public:
            /*-----------------------------------------------------------*/
            /** \brief Specific constructor from 4 components
             *
             * \param A1 Component 1
             * \param A2 Component 2
             * \param A3 Component 3
             * \param A4 Component 4
             */
            Vector4d(const double& A1, const double& A2,
                     const double& A3, const double& A4)
            :VectorND<4, double>()
            {
                m_data[0] = A1;
                m_data[1] = A2;
                m_data[2] = A3;
                m_data[3] = A4;
            }
        };
        
        typedef VectorND<9, double> Vector9d;
        typedef VectorND<3, double> Vector3i;
        
        
    } //namespace math
    /*-----------------------------------------------------------------*/
}//namespace gmds
/*---------------------------------------------------------------------*/
/** \brief Writes a vector to a stream
 *
 * \param[in] AStream the output stream
 * \param[in] AV      the vector to write
 * \return a reference to the output stream \p AStream
 */
template <int TDim, class TType> std::ostream&
operator<< (std::ostream& AStream,
            const gmds::math::VectorND<TDim,TType>& AV)
{
    AStream<<"(";
    for(int i = 0; i < TDim-1; i++) {
        AStream<< AV[i]<<", ";
    }
    AStream<<AV[TDim-1]<<")";
    return AStream;
}
/*---------------------------------------------------------------------*/
#endif /*GMDS_MATH_VECTOR_ND_H_*/
/*---------------------------------------------------------------------*/
