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
 * Quadrilateral.cpp
 *
 *  Created on: 26 oct. 2015
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/Math/Quadrilateral.h>
/*----------------------------------------------------------------------------*/
#include <GMDS/Math/Vector.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace math{
/*----------------------------------------------------------------------------*/
Quadrilateral::Quadrilateral()
{
	m_pnts[0] = Point(0,0,0);
	m_pnts[1] = Point(1,0,0);
	m_pnts[2] = Point(1,1,0);
	m_pnts[3] = Point(0,1,0);
}
/*----------------------------------------------------------------------------*/
Quadrilateral::Quadrilateral(const Point& AP1, const Point& AP2, const Point& AP3, const Point& AP4)
{
	m_pnts[0] = AP1;
	m_pnts[1] = AP2;
	m_pnts[2] = AP3;
	m_pnts[3] = AP4;
}
/*----------------------------------------------------------------------------*/
Quadrilateral::Quadrilateral(const Quadrilateral& AQ)
{
	m_pnts[0] = AQ.m_pnts[0];
	m_pnts[1] = AQ.m_pnts[1];
	m_pnts[2] = AQ.m_pnts[2];
	m_pnts[3] = AQ.m_pnts[3];
}
/*----------------------------------------------------------------------------*/
Quadrilateral::~Quadrilateral(){;}
/*----------------------------------------------------------------------------*/
const Point& Quadrilateral::getPoint(const TInt& AIndex) const
{
	return m_pnts[AIndex];
}
/*----------------------------------------------------------------------------*/
int Quadrilateral::getNbPoints() const
{
        return 4;
}
/*----------------------------------------------------------------------------*/
Point
Quadrilateral::getCenter() const {
	Point pt((m_pnts[0]+m_pnts[1]+m_pnts[2]+m_pnts[3])*(1./4.));
	return pt;
}
/*----------------------------------------------------------------------------*/
double
Quadrilateral::computeScaledJacobian2D() const
{        
	int neighbors[4][2] =
        {       
                {1,3},
                {2,0},
                {3,1},
                {0,2}
        };
        
        double scaledJ[4]; 
        for (int iVertex = 0; iVertex < 4; iVertex++) {
                
                Matrix<2,2,double> A = this->jacobian2D(iVertex);
                
                int i0    = neighbors[iVertex][0];
                int i1    = neighbors[iVertex][1];
                double l0 = math::Vector(m_pnts[iVertex]- m_pnts[i0]).norm2();
                double l1 = math::Vector(m_pnts[iVertex]- m_pnts[i1]).norm2();
                double l01 = sqrt(l0*l1);
                double det = A.det();
                if (l01 == 0 || det == 0) 
                        scaledJ[iVertex] = 0.;
                else    
                        scaledJ[iVertex] = det/l01;
        }
        double scaledJmin = HUGE_VALF; 
        for(int iVertex = 0; iVertex < 4; iVertex++) {
                if(scaledJ[iVertex] < scaledJmin) {
                        scaledJmin = scaledJ[iVertex];
                }
        }
        
        return scaledJmin;
}
/*----------------------------------------------------------------------------*/
double
Quadrilateral::computeNormalizedScaledJacobian2D() const
{
        double scaledJ = this->computeScaledJacobian2D();

	if(scaledJ >  1.) scaledJ =  1.;
	if(scaledJ < -1.) scaledJ = -1.;

	return scaledJ;
}
/*----------------------------------------------------------------------------*/
double
Quadrilateral::computeMeanEdgeLength() const
{
  double sumLength = 0;
  sumLength += m_pnts[0].distance(m_pnts[1]);
  sumLength += m_pnts[1].distance(m_pnts[2]);
  sumLength += m_pnts[2].distance(m_pnts[3]);
  sumLength += m_pnts[3].distance(m_pnts[0]);

  sumLength /= 4.;  
  return sumLength;
}
/*---------------------------------------------------------------------------*/
math::Matrix<2,2,double>
Quadrilateral::jacobian2D(const int iVertex) const
{
        TCoord matValues[2][2];

        switch(iVertex) {
        case 0:
                matValues[0][0] = m_pnts[1].X()-m_pnts[0].X();
                matValues[1][0] = m_pnts[1].Y()-m_pnts[0].Y();
                matValues[0][1] = m_pnts[3].X()-m_pnts[0].X();
                matValues[1][1] = m_pnts[3].Y()-m_pnts[0].Y();
                break;
        case 1:
                matValues[0][0] = m_pnts[1].X()-m_pnts[0].X();
                matValues[1][0] = m_pnts[1].Y()-m_pnts[0].Y();
                matValues[0][1] = m_pnts[2].X()-m_pnts[1].X();
                matValues[1][1] = m_pnts[2].Y()-m_pnts[1].Y();
                break;
        case 2:
                matValues[0][0] = m_pnts[2].X()-m_pnts[3].X();
                matValues[1][0] = m_pnts[2].Y()-m_pnts[3].Y();
                matValues[0][1] = m_pnts[2].X()-m_pnts[1].X();
                matValues[1][1] = m_pnts[2].Y()-m_pnts[1].Y();
                break;
        case 3:
		matValues[0][0] = m_pnts[2].X()-m_pnts[3].X();
                matValues[1][0] = m_pnts[2].Y()-m_pnts[3].Y();
                matValues[0][1] = m_pnts[3].X()-m_pnts[0].X();
                matValues[1][1] = m_pnts[3].Y()-m_pnts[0].Y();
                break;
	default:
                throw GMDSException("Quadrilateral::jacobian invalid vertex number.");
                break;
        }

        return math::Matrix<2,2,double> (matValues);
}
/*----------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& AStr, const Quadrilateral& AQuad){
        AStr<<"Quadrilateral"<<std::endl;
        for(int iPoint=0; iPoint<4; iPoint++) {
                AStr<<std::endl<<AQuad.getPoint(iPoint);
        }
        return AStr;
}
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/

