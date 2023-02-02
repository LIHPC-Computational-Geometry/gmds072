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
 * Chart.cpp
 *  Created on: 03/19/2015
 *      Author: ledouxf
 */
/*-----------------------------------------------------------------*/
#include "GMDS/Math/Chart.h"
#include "GMDS/Math/Quaternion.h"
#include "GMDS/Math/AxisAngleRotation.h"
/*-----------------------------------------------------------------*/
#include <cmath>
/*-----------------------------------------------------------------*/
namespace gmds {
    /*-----------------------------------------------------------------*/
    namespace math {
        /*-----------------------------------------------------------------*/
        Chart::Chart()
        {
            m_v[0] = Vector3d(1,0,0);
            m_v[1] = Vector3d(0,1,0);
            m_v[2] = Vector3d(0,0,1);
        }
        /*-----------------------------------------------------------------*/
        Chart::Chart(const Vector& AX, const Vector& AY,const Vector& AZ)
        {
            
            m_v[0] = Vector3d(AX.X(),AX.Y(),AX.Z());
            m_v[1] = Vector3d(AY.X(),AY.Y(),AY.Z());
            m_v[2] = Vector3d(AZ.X(),AZ.Y(),AZ.Z());
        }
        
        /*-----------------------------------------------------------------*/
        Chart::Chart(const Vector3d& AX, const Vector3d& AY,
                     const Vector3d& AZ)
        {
            m_v[0] = AX;
            m_v[1] = AY;
            m_v[2] = AZ;
        }
        /*-----------------------------------------------------------------*/
        Chart::Chart(const Vector3d& AX, const Vector3d& AY)
        {
            m_v[0] = AX;
            m_v[1] = AY;
            m_v[2] = AX.cross(AY);
        }
        /*-----------------------------------------------------------------*/
        Chart::Chart(const Vector& AX, const Vector& AY)
        {
            m_v[0] = AX;
            m_v[1] = AY;
            m_v[2] = AX.cross(AY);
        }
        /*-----------------------------------------------------------------*/
        Chart::Chart(const Quaternion & AQ)
        {
            m_v[0][0] = AQ.X()*AQ.X() + AQ.I() * AQ.I() - AQ.J()*AQ.J() - AQ.K()*AQ.K();
            m_v[0][1] = 2.0 * AQ.X() * AQ.K() + 2.0 * AQ.I() * AQ.J();
            m_v[0][2] = 2.0 * AQ.K() * AQ.I() - 2.0 * AQ.X() * AQ.J();
            
            m_v[1][0] = 2.0 * AQ.I() * AQ.J() - 2.0 * AQ.X() * AQ.K();
            m_v[1][1] = AQ.X() * AQ.X() - AQ.I() * AQ.I() + AQ.J() * AQ.J() - AQ.K() * AQ.K();
            m_v[1][2] = 2.0 * AQ.X() * AQ.I() + 2.0 * AQ.K() * AQ.J();
            
            m_v[2][0] = 2.0 * AQ.X() * AQ.J() + 2.0 * AQ.I() * AQ.K();
            m_v[2][1] = 2.0 * AQ.K() * AQ.J() - 2.0 * AQ.X() * AQ.I();
            m_v[2][2] = AQ.X() * AQ.X() - AQ.I() * AQ.I() - AQ.J() * AQ.J() + AQ.K() * AQ.K();
            
            m_v[0].normalize();
            m_v[1].normalize();
            m_v[2].normalize();
            
        }
        /*-----------------------------------------------------------------*/
        Vector Chart::get(const int AIndex) const
        {
            if(AIndex<0 || AIndex>2)
                throw GMDSException("Chart - Wrong range index");
            
            return Vector(m_v[AIndex].X(),m_v[AIndex].Y(),m_v[AIndex].Z());
        }
        
        /*-----------------------------------------------------------------*/
        Vector3d Chart::operator[](const int AIndex) const
        {
            return m_v[AIndex];
        }
        
        /*-----------------------------------------------------------------*/
        std::vector<Vector> Chart::image() const
        {
            std::vector<Vector> res;
            res.reserve(6);
            res.push_back(math::Vector(m_v[0].X(),m_v[0].Y(),m_v[0].Z()));
            res.push_back(math::Vector(m_v[1].X(),m_v[1].Y(),m_v[1].Z()));
            res.push_back(math::Vector(m_v[2].X(),m_v[2].Y(),m_v[2].Z()));
            res.push_back(math::Vector(m_v[0].X(),m_v[0].Y(),m_v[0].Z()).opp());
            res.push_back(math::Vector(m_v[1].X(),m_v[1].Y(),m_v[1].Z()).opp());
            res.push_back(math::Vector(m_v[2].X(),m_v[2].Y(),m_v[2].Z()).opp());
            return res;
        }
        
        /*-----------------------------------------------------------------*/
        Matrix<3,3,double> Chart::toMatrix() const
        {
            double v[3][3]= {
                { m_v[0].X(), m_v[1].X(), m_v[2].X()},
                { m_v[0].Y(), m_v[1].Y(), m_v[2].Y()},
                { m_v[0].Z(), m_v[1].Z(), m_v[2].Z()} };
            return Matrix<3,3,double>(v);
        }
        
        /*-----------------------------------------------------------------*/
        Matrix<3,3,double> Chart::computeRotationTo(const Chart& AC) const
        {
            Matrix<3,3,double> m;
            Quaternion q_from(*this);
            Quaternion q_to(AC);
            Quaternion closest_q_to = q_to.closestImg(q_from);
            AxisAngleRotation r(closest_q_to);
            return r.toRotationMatrix();
        }
        /*-----------------------------------------------------------------*/
        void Chart::align(const Vector & AV)
        {
            const vector<Vector> all_vectors = this->image();
            TCoord tmp;
            TCoord max = all_vectors.front().dot(AV);
            unsigned int i_max = 0;
            
            for (auto i = 1; i < all_vectors.size(); i++){
                tmp = all_vectors[i].dot(AV);
                if (tmp > max) {
                    max = tmp;
                    i_max = i;
                }
            }
            
            if (max != 1.0){
                //A rotation is mandatory to be align with AV
                Vector  invariant = AV.cross(all_vectors[i_max]);
                invariant.normalize();
                
                const TCoord cos_thet = max;
                TCoord sin_thet = 0.0;
                const TCoord acos_of_cos_thet = acos(cos_thet);
                
                if (acos_of_cos_thet > (-10.0) && (acos_of_cos_thet < 10.0))
                    sin_thet = sin(acos_of_cos_thet);
                
                Vector line1, line2, line3;
                
                line1.setX(invariant[0] * invariant[0] + (1.0 - invariant[0] * invariant[0]) * cos_thet);
                line1.setY(invariant[0] * invariant[1] * (1.0 - cos_thet) + invariant[2] * sin_thet);
                line1.setZ(invariant[0] * invariant[2] * (1.0 - cos_thet) - invariant[1] * sin_thet);
                
                line2.setX(invariant[0] * invariant[1] * (1.0 - cos_thet) - invariant[2] * sin_thet);
                line2.setY(invariant[1] * invariant[1] + (1.0 - invariant[1] * invariant[1]) * cos_thet);
                line2.setZ(invariant[1] * invariant[2] * (1.0 - cos_thet) + invariant[0] * sin_thet);
                
                line3.setX(invariant[0] * invariant[2] * (1.0 - cos_thet) + invariant[1] * sin_thet);
                line3.setY(invariant[1] * invariant[2] * (1.0 - cos_thet) - invariant[0] * sin_thet);
                line3.setZ(invariant[2] * invariant[2] + (1.0 - invariant[2] * invariant[2]) * cos_thet);
                
                
                m_v[0][0] = all_vectors[0].dot(line1);
                m_v[0][1] = all_vectors[0].dot(line2);
                m_v[0][2] = all_vectors[0].dot(line3);
                
                m_v[1][0] = all_vectors[1].dot(line1);
                m_v[1][1] = all_vectors[1].dot(line2);
                m_v[1][2] = all_vectors[1].dot(line3);
                
                m_v[2][0] = all_vectors[2].dot(line1);
                m_v[2][1] = all_vectors[2].dot(line2);
                m_v[2][2] = all_vectors[2].dot(line3);
            }
        }
        /*-----------------------------------------------------------------*/
        void Chart::matchVectors(const Chart&AChart,
                                 int (&AMatching)[3]) const
        {
            math::Chart c1 = *this;
            math::Chart c2 = AChart;
            
            math::Vector c1_vec[3] = {c1.X(), c1.Y(), c1.Z()};
            math::Vector c2_vec[3] = {c2.X(), c2.Y(), c2.Z()};
            
            /* We keep int mind which vectors of c1 and c2 allows to get the
             * maximal dot product (in abs value). It means then that those
             * vector are the most aligned ones.
             */
            int iMax = 0, jMax = 0;
            // the dotMax value store the maximal dot product
            double dotMax = 0.0;
            double dotProducts[3][3];
            
            
            for (unsigned int i = 0; i < 3; i++){
                for (unsigned int j = 0; j < 3; j++){
                    double dot_ij = c1_vec[i].dot(c2_vec[j]);
                    if (dot_ij < 0.0){
                        dot_ij = -dot_ij;
                    }
                    dotProducts[i][j] = dot_ij;
                    if (dot_ij > dotMax){
                        dotMax = dot_ij;
                        iMax = i;
                        jMax = j;
                    }
                } //for (unsigned int j = 0; j < 3; j++){
            } //for (unsigned int i = 0; i < 3; i++){
            
            //We have the first matching
            AMatching[iMax] = jMax;
            //We look for the second matching
            int iSecond = 0, jSecond = 0;
            double dotSecond = 0.0;
            for (unsigned int i = 0; i < 3; i++){
                for (unsigned int j = 0; j < 3; j++){
                    if (i != iMax && j != jMax){
                        if (dotProducts[i][j] > dotSecond){
                            dotSecond = dotProducts[i][j];
                            iSecond = i;
                            jSecond = j;
                        }
                    }
                }
            }
            
            AMatching[iSecond] = jSecond;
            //We look for the third matching
            int iThird = 0, jThird = 0;
            for (unsigned int i = 0; i < 3; i++){
                if(i!=iMax && i!=iSecond)
                    iThird = i;
            }
            for (unsigned int j = 0; j < 3; j++){
                if(j!=jMax && j!=jSecond)
                    jThird = j;
            }
            AMatching[iThird] = jThird;
            
        }
        
        /*-----------------------------------------------------------------*/
        Chart::Mapping::Mapping(){
            int m[3]={0,1,2}, d[3]={1,1,1};
            m_map= Vector3i(m);
            m_dir= Vector3i(d);
        }
        /*-----------------------------------------------------------------*/
        Chart::Mapping::Mapping(const Chart& AC1, const Chart& AC2){
            int matching[3];
            AC1.matchVectors(AC2,matching);
            for(int i=0; i<3; i++){
                m_map[i] = matching[i];
                m_dir[i] = (AC1[i].dot(AC2[matching[i]])>0)?1:-1;
            }
        }
        /*-----------------------------------------------------------------*/
        const Vector3i& Chart::Mapping::getPermutations() const {
            return m_map;
        }
        /*-----------------------------------------------------------------*/
        const Vector3i& Chart::Mapping::getDirections() const {
            return m_dir;
        }
        /*-----------------------------------------------------------------*/
        Chart::Mapping Chart::Mapping::inverse() const {
            Mapping inv;
            for (int i = 0; i<3; i++){
                inv.m_map[m_map[i]] = i;
                inv.m_dir[m_map[i]] = m_dir[i];
            }
            return inv;
        }
        /*-----------------------------------------------------------------*/
        bool Chart::Mapping::isIdentity() const {
            for (int i = 0; i<3; i++) {
                if(m_map[i]!=i || m_dir[i]!=1)
                    return false;
            }
            return true;
        }
        /*-----------------------------------------------------------------*/
        Matrix<3,3,double> Chart::Mapping::toMatrix() const {
            Matrix<3,3,double> M;
            std::cout<<"map "<<m_map[0]<<" "
            <<m_map[1]<<" "
            <<m_map[2]<<std::endl;
            std::cout<<"dir "<<m_dir[0]<<" "
            <<m_dir[1]<<" "
            <<m_dir[2]<<std::endl;
            for (int i = 0; i<3; i++) {
                for (int j = 0; j<3; j++) {
                    if(m_map[i]==j)
                        M(j,i) = m_dir[i];
                    else
                        M(j,i) = 0.0;
                }
            }
            return M;
        }
        /*-----------------------------------------------------------------*/
        Chart::Mapping Chart::Mapping::operator*(const Mapping& AM) const {
            Mapping comp;
            for (int i = 0; i<3; i++) {
                comp.m_map[i] = m_map[AM.m_map[i]];
                comp.m_dir[i] = m_dir[AM.m_map[i]] * AM.m_dir[i];
            }
            return comp;
        }
        /*-----------------------------------------------------------------*/
        Vector3d Chart::Mapping::operator*(const Vector3d& AV) const {
            Vector3d v;
            for (int i = 0; i<3; i++) {
                v[i] = m_dir[i]*AV[m_map[i]];
            }
            return v;
        }

        /*-----------------------------------------------------------------*/
        ostream & operator << (ostream & AStream, const Chart & AC)
        {
            AStream<<"Chart"<<endl;
            AStream<<" x : " << AC.X() << endl;
            AStream<<" y : " << AC.Y() << endl;
            AStream<<" z : " << AC.Z();
            return AStream;
        }
        /*-----------------------------------------------------------------*/
    }//namespace math
    /*-----------------------------------------------------------------*/
} //namespace gmds

/*-----------------------------------------------------------------*/
