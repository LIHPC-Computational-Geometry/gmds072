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
 * Point.cpp
 *
 *  Created on: 6 févr. 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/Math/Point.h>
#include <GMDS/Math/Plane.h>
#include <GMDS/Math/Matrix.h>
#include <GMDS/Math/Vector.h>
#include <GMDS/Math/VectorND.h>
/*-----------------------------------------------------------------------------*/
#include <cmath>
/*-----------------------------------------------------------------------------*/
namespace gmds{
    /*-------------------------------------------------------------------------*/
    namespace math{
        /*---------------------------------------------------------------------*/
        Point::Point(const TCoord& AX, const TCoord& AY, const TCoord& AZ)
        : m_x(AX),m_y(AY),m_z(AZ)
        {}
        /*---------------------------------------------------------------------*/
        Point::~Point()
        {}
        /*---------------------------------------------------------------------*/
        TCoord Point::distance(const Point& AP) const{
            return sqrt(distance2(AP));
        }
        /*---------------------------------------------------------------------*/
        TCoord Point::distance2(const Point& AP) const{
            return Vector(*this,AP).norm2();
        }
        /*---------------------------------------------------------------------*/
        Point operator+(const Point& AP1, const Point& AP2){
            return Point(
                         AP1.m_x + AP2.m_x,
                         AP1.m_y + AP2.m_y,
                         AP1.m_z + AP2.m_z);
        }
        /*---------------------------------------------------------------------*/
        Point operator-(const Point& AP1, const Point& AP2){
            return Point(
                         AP1.m_x - AP2.m_x,
                         AP1.m_y - AP2.m_y,
                         AP1.m_z - AP2.m_z);
        }
        /*---------------------------------------------------------------------*/
        bool Point::operator==(const Point& AP) const
        {
            if (&AP == this)
                return true;
            
            return AP.m_x==m_x && AP.m_y==m_y && AP.m_z==m_z;
        }
        /*---------------------------------------------------------------------*/
        bool Point::operator!= (const Point& AP) const
        {
            if (&AP == this)
                return false;
            
            return AP.m_x!=m_x || AP.m_y!=m_y || AP.m_z!=m_z;
        }
        /*---------------------------------------------------------------------*/
        bool Point::areColinear(const Point& AP2, const Point& AP3) const
        {
            Point AP1 = *this;
            Vector v1(AP1, AP2);
            Vector v2(AP1, AP3);
            Vector v = v1.cross(v2);
            
            return (isZero(v[0]) && isZero(v[1]) && isZero(v[2]) );
        }
        /*---------------------------------------------------------------------*/
        bool Point::
        areCoplanar(const Point& AP2, const Point& AP3, const Point& AP4) const{
            
            Point AP1 = *this;
            Vector v12(AP1,AP2);
            Vector v13(AP1,AP3);
            Vector v14(AP1,AP4);
      //      std::cout<<"Coplanar? "<<v12.dot(v13.cross(v14))<<std::endl;
            return isZero(v12.dot(v13.cross(v14)));
            
//            // is this more robust than computing the volume?
//            if(isZero(v41.norm2()) || isZero(v42.norm2()) || isZero(v43.norm2()) )
//                return true;
//            
//            Plane pl(AP1,AP2,AP3);
//            v41.normalize();
//            Vector normal = pl.getNormal();
//            normal.normalize();
//            TCoord dotProduct = v41.dot(normal);
//            return isZero(dotProduct);
        }
        /*---------------------------------------------------------------------*/
        bool
        Point::isStrictlyOnLeft2D(const Point& AP1, const Point& AP2) const
        {
            
            double r = (AP2.X() - AP1.X()) * (this->Y() - AP1.Y()) - (this->X() - AP1.X()) * (AP2.Y() - AP1.Y());
            
            return (r > 0.0) ? true : false;
        }
        /*---------------------------------------------------------------------*/
        void Point::computeBarycentric2D(
                                         const math::Point& AT1,
                                         const math::Point& AT2,
                                         const math::Point& AT3,
                                         const math::Point& AP,
                                         TCoord& AX, TCoord& AY, TCoord& AZ)
        {
            TCoord det;
            
            if (AT1.areColinear(AT2,AT3))
            {
                throw GMDSException("flat triangle in barycentric computation");
            }
            TCoord x1 = AT1.X();
            TCoord y1 = AT1.Y();
            TCoord x2 = AT2.X();
            TCoord y2 = AT2.Y();
            TCoord x3 = AT3.X();
            TCoord y3 = AT3.Y();
            TCoord x = AP.X();
            TCoord y = AP.Y();
            
            Matrix<2, 2,double> T;
            T.set(0, 0, x1 - x3);
            T.set(0, 1, x2 - x3);
            T.set(1, 0, y1 - y3);
            T.set(1, 1, y2 - y3);
            
            det = T.det();
            
            if (isZero(det)) {
                throw GMDSException("error in barycentric computation");
            }
            
            TCoord a1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3));
            TCoord a2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3));
            //TCoord a3 = det - a1 - a2;
            
            AX = a1 / det;
            AY = a2 / det;
            AZ = 1.0 - AX - AY;
        }
        
        /*---------------------------------------------------------------------*/
        void Point::computeBarycentric(const math::Point& AT0,
                                       const math::Point& AT1,
                                       const math::Point& AT2,
                                       const math::Point& AT3,
                                       const math::Point& AP,
                                       TCoord& A0, TCoord& A1,
                                       TCoord& A2, TCoord& A3) {
            double val_A[4][4] = {
                {AT0.X(), AT1.X(), AT2.X(), AT3.X()},
                {AT0.Y(), AT1.Y(), AT2.Y(), AT3.Y()},
                {AT0.Z(), AT1.Z(), AT2.Z(), AT3.Z()},
                {1.0    , 1.0    , 1.0    , 1.0    }};
            
            double val_b[4] = {AP.X(),AP.Y(),AP.Z(),1.0};
            VectorND<4, double>   b(val_b);
            Matrix<4, 4, double>  A(val_A);
            
            //We solve AX=b
            VectorND<4, double> x = A.solve(b);
            
            A0=x[0]; A1=x[1]; A2=x[2]; A3=x[3];
        }
        
        /*---------------------------------------------------------------------*/
        void Point::computeBarycentric(const std::vector<math::Point>& AT,
                                       const math::Point& AP,
                                       std::vector<TCoord>& ACoeff) {
            if(AT.size()==4){
                double val_A[4][4] = {
                    {AT[0].X(), AT[1].X(), AT[2].X(), AT[3].X()},
                    {AT[0].Y(), AT[1].Y(), AT[2].Y(), AT[3].Y()},
                    {AT[0].Z(), AT[1].Z(), AT[2].Z(), AT[3].Z()},
                    {1.0      , 1.0      , 1.0      , 1.0      }};
                
                double val_b[4] = {AP.X(),AP.Y(),AP.Z(),1.0};
                VectorND<4, double>   b(val_b);
                Matrix<4, 4, double>  A(val_A);
                //We solve AX=b
                VectorND<4, double> x = A.solve(b);
                ACoeff.resize(4);
                ACoeff[0]=x[0];
                ACoeff[1]=x[1];
                ACoeff[2]=x[2];
                ACoeff[3]=x[3];
            }
            else if (AT.size()==3){
                double val_A[3][3] = {
                    {AT[0].X(), AT[1].X(), AT[2].X()},
                    {AT[0].Y(), AT[1].Y(), AT[2].Y()},
                    {AT[0].Z(), AT[1].Z(), AT[2].Z()}};
                
                double val_b[3] = {AP.X(),AP.Y(),AP.Z()};
                VectorND<3, double>   b(val_b);
                Matrix<3, 3, double>  A(val_A);
                //We solve AX=b
                VectorND<3, double> x = A.solve(b);
                ACoeff.resize(3);
                ACoeff[0]=x[0];
                ACoeff[1]=x[1];
                ACoeff[2]=x[2];
                
                
            }
            else
                throw GMDSException("Point::computeBarycentric only for 3 or 4 points");
        }
        /*---------------------------------------------------------------------*/
        Point Point::massCenter(const std::vector<math::Point>& AP) {
            double x=0, y=0, z=0;
            for(int p=0; p<AP.size(); p++){
                x+=AP[p].X();
                y+=AP[p].Y();
                z+=AP[p].Z();
            }
            return Point(x/AP.size(),y/AP.size(),z/AP.size());
        }
        /*---------------------------------------------------------------------*/
        void Point::computeBarycentric(
                                       const math::Point& AT1,
                                       const math::Point& AT2,
                                       const math::Point& AT3,
                                       const math::Point& AP,
                                       TCoord& AX, TCoord& AY, TCoord& AZ)
        {
            math::Point p0 = AT1;
            math::Point p1 = AT2;
            math::Point p2 = AT3;
            
            if(p0.areColinear(p1,p2))
            {
                throw GMDSException("flat triangle in barycentric computation");
            }
            if(!p0.areCoplanar(p1,p2,AP))
            {
                //throw GMDSException("Coplanarity is mandatory to compute barycentric coordinates of a point into a 3D triangle");
            }
            
            math::Vector v1(p1.X()-p0.X(), p1.Y()-p0.Y(), p1.Z()-p0.Z());
            math::Vector v3(p2.X()-p0.X(), p2.Y()-p0.Y(), p2.Z()-p0.Z());
            
            Vector normal = (v1.cross(v3));
            normal.normalize();
            
            int maxIndex = normal.getMaxAbsComponentIndex();
            
            if(maxIndex==2) {
                // we can project on plane Oxy
                Point p(AP.X() ,AP.Y());
                Point a(AT1.X(),AT1.Y());
                Point b(AT2.X(),AT2.Y());
                Point c(AT3.X(),AT3.Y());
                computeBarycentric2D(a,b,c,p,AX,AY,AZ);
            } else if(maxIndex==1) {
                // we can project on plane Oxz
                Point p(AP.X() ,AP.Z());
                Point a(AT1.X(),AT1.Z());
                Point b(AT2.X(),AT2.Z());
                Point c(AT3.X(),AT3.Z());
                computeBarycentric2D(a,b,c,p,AX,AY,AZ);
            } else {
                // we can project on plane Oyz
                Point p(AP.Y() ,AP.Z());
                Point a(AT1.Y(),AT1.Z());
                Point b(AT2.Y(),AT2.Z());
                Point c(AT3.Y(),AT3.Z());
                computeBarycentric2D(a,b,c,p,AX,AY,AZ);
            }
        }
        
        /*---------------------------------------------------------------------*/
        Point operator*(const TCoord& AK, const Point& AP){
            return Point(AK*AP.m_x, AK*AP.m_y, AK*AP.m_z);
        }
        /*---------------------------------------------------------------------*/
        Point operator*(const Point& AP, const TCoord& AK){
            return Point(AK*AP.m_x, AK*AP.m_y, AK*AP.m_z);
        }
        /*---------------------------------------------------------------------*/
        std::ostream& operator<<(std::ostream& AStr, const Point& AP){
            AStr<<"("<<AP.m_x<<", "<<AP.m_y<<", "<<AP.m_z<<")";
            return AStr;
        }
        
        /*---------------------------------------------------------------------*/
    } // namespace math
    /*-------------------------------------------------------------------------*/
} // namespace gmds
bool operator<(const gmds::math::Point& AP1, const gmds::math::Point& AP2){
    return ((AP1.X()<AP2.X()) ||
            (AP1.X()==AP2.X() && AP1.Y()<AP2.Y()) ||
            (AP1.X()==AP2.X() && AP1.Y()==AP2.Y() && AP1.Z()<AP2.Z()) );
}
/*----------------------------------------------------------------------------*/
