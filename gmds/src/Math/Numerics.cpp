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
 * Numerics.cpp
 *
 *  Created on: March 01, 2015
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
// GMDS header files
#include <GMDS/Math/Constants.h>
#include <GMDS/Math/Vector.h>
#include <GMDS/Math/Segment.h>
#include <GMDS/Math/Numerics.h>
/*----------------------------------------------------------------------------*/
// Eigen header files
#include <Eigen/Eigenvalues>
/*----------------------------------------------------------------------------*/
namespace gmds{
    /*----------------------------------------------------------------------------*/
    namespace math{
        /*----------------------------------------------------------------------------*/
        TCoord modulo(const TCoord AVal, const TCoord AMod){
            double v = AVal;
            while(v >= AMod || v < 0.) {
                if(v > 0)
                    v -= AMod;
                else
                    v += AMod;
            }
            return v;
        }
        /*----------------------------------------------------------------------------*/
        TCoord modulo2PI(const TCoord AVal){
            return modulo(AVal,Constants::PI2);
        }
        /*----------------------------------------------------------------------------*/
        TCoord min2(const TCoord AV1, const TCoord AV2){
            double val = AV1;
            if(AV2<val)
                val=AV2;
            return val;
        }
        /*----------------------------------------------------------------------------*/
        TCoord max2(const TCoord AV1, const TCoord AV2){
            double val = AV1;
            if(AV2>val)
                val=AV2;
            return val;
        }
        /*----------------------------------------------------------------------------*/
        TCoord min3(const TCoord AV1, const TCoord AV2, const TCoord AV3){
            double val = AV1;
            if(AV2<val)
                val=AV2;
            if(AV3<val)
                val=AV3;
            return val;
        }
        /*----------------------------------------------------------------------------*/
        TCoord max3(const TCoord AV1, const TCoord AV2, const TCoord AV3){
            double val = AV1;
            if(AV2>val)
                val=AV2;
            if(AV3>val)
                val=AV3;
            return val;
        }
        /*----------------------------------------------------------------------------*/
        bool near(const TCoord AV1, TCoord AV2, TCoord AEpsilon){
            return (fabs(AV1-AV2)<=AEpsilon);
        }
        /*----------------------------------------------------------------------------*/
        Matrix<3,3,double> stiffnessMatrix2D(const Point& AP1, const Point& AP2,
                                             const Point& AP3)
        {
            math::Matrix<3,3,double> s;
            double x1 = AP1.X(), y1 = AP1.Y();
            double x2 = AP2.X(), y2 = AP2.Y();
            double x3 = AP3.X(), y3 = AP3.Y();
            
            double x23 = x2-x3;
            double x31 = x3-x1;
            double x12 = x1-x2;
            double y23 = y2-y3;
            double y31 = y3-y1;
            double y12 = y1-y2;
            
            Vector v12(AP1,AP2);
            Vector v13(AP1,AP3);
            double area2 = x1*y2 - x1*y3 - x2*y1 + x2*y3 + x3*y1 -x3*y2;
            //      double area2 = v12.cross(v13).norm(); //get Area Triangle(AP1,AP2,AP3)
            double area4 = 2*area2;
            double inv_area4 = 1.0/area4;
            // First row
            s.set(0, 0, (y23*y23 + x23*x23)*inv_area4);
            s.set(0, 1, (y23*y31 + x23*x31)*inv_area4);
            s.set(0, 2, (y23*y12 + x23*x12)*inv_area4);
            // Second row
            s.set(1, 0, (y23*y31 + x23*x31)*inv_area4);
            s.set(1, 1, (y31*y31 + x31*x31)*inv_area4);
            s.set(1, 2, (y31*y12 + x31*x12)*inv_area4);
            // Second row
            s.set(2, 0, (y23*y12 + x23*x12)*inv_area4);
            s.set(2, 1, (y31*y12 + x31*x12)*inv_area4);
            s.set(2, 2, (y12*y12 + x12*x12)*inv_area4);
            
            return s;
        }
 
        /*----------------------------------------------------------------------------*/
        int solve2ndDegreePolynomial(const double& AA,
                                        const double& AB,
                                        const double& AC,
                                        std::vector<double>& AX)
        {
            AX.clear();
            int nb_solutions = 0;
            
            double a = AA;
            double b = AB;
            double c = AC;
            
            double det = b*b - 4*a*c;
            
            if (det < 0.0) {
                std::cout << "No solution for: "<<a<<" X2 + "<<b<<" X + "<<c <<" with Det= "<<det<< std::endl;
                nb_solutions = 0;
            }
            else if (det==0) {
                nb_solutions = 1;
                double x = -b/(2*a);
                AX.push_back(x);
            }
            else {
                nb_solutions = 2;
                double x1 = (-b+sqrt(det))/(2*a);
                double x2 = (-b-sqrt(det))/(2*a);
                AX.push_back(x1);
                AX.push_back(x2);
                
            }
            return nb_solutions;
        }
        /*----------------------------------------------------------------------------*/
        double dihedralAngle(const Point& A,
                             const Point& B,
                             const Point& C,
                             const Point& D)
        {
            math::Vector3d u(A,B);
            math::Vector3d v(A,C);
            math::Vector3d w(A,D);
            math::Vector3d uv = u.cross(v);
            math::Vector3d uw = u.cross(w);
            uv.normalize();
            uw.normalize();
            auto cos_angle =uv.dot(uw);
            if(cos_angle>1)
                cos_angle=1;
            else if(cos_angle<0){
                cos_angle=0;
            }
            return std::acos(cos_angle);
            
        }
        /*----------------------------------------------------------------------------*/
        double cotangentWeight(const Point& A,
                               const Point& B,
                               const Point& C,
                               const Point& D)
        {
            math::Vector3d cd(C,D);
            double opp_angle = dihedralAngle(C,D,A,B);
            return (cd.norm()*std::atan(opp_angle))/6.0;
        }
        /*----------------------------------------------------------------------------*/
        void computeLeastSquarePlane(const std::vector<Point>& AP,
                                     Point& APlanePnt,
                                     Vector3d& APlaneNormal)
        {
            Point mc = Point::massCenter(AP);
            
            
            Eigen::Matrix3d A;
            
            for(auto i=0; i<3; i++)
                for(auto j=0; j<3; j++)
                    A(i,j)=0;
            
            for(int p=0; p<AP.size(); p++){
                double x = AP[p].X()-mc.X();
                double y = AP[p].Y()-mc.Y();
                double z = AP[p].Z()-mc.Z();
                A(0,0) += x*x; A(0,1) += x*y; A(0,2) += x*z;
                A(1,0) += x*y; A(1,1) += y*y; A(1,2) += y*z;
                A(2,0) += x*z; A(2,1) += y*z; A(2,2) += z*z;
            }
            Eigen::EigenSolver<Eigen::MatrixXd> solver(A);
            
            int sev=0;
            double val = solver.eigenvalues()(0).real();
            for(auto i=1;i<3;i++){
                double val_i =solver.eigenvalues()(i).real();
                if(val_i<val){
                    sev = i;
                    val = val_i;
                }
            }
            
            double a = solver.eigenvectors().col(sev)(0).real();
            double b = solver.eigenvectors().col(sev)(1).real();
            double c = solver.eigenvectors().col(sev)(2).real();
            
            
            APlanePnt = mc;
            APlaneNormal = math::Vector3d(a,b,c);
        }

        /*----------------------------------------------------------------------------*/
    } // namespace math
    /*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
