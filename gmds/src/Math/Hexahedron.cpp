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
 * Hexahedron.cpp
 *
 *  Created on: 16 oct 2014
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/Math/Hexahedron.h>
/*----------------------------------------------------------------------------*/
#include <GMDS/Math/Vector.h>
#include <GMDS/Math/Triangle.h>
#include <GMDS/Math/Pyramid.h>
#include <GMDS/Math/Matrix.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace math{
/*----------------------------------------------------------------------------*/
Hexahedron::Hexahedron()
{
	m_pnts[0] = Point(0,0,0);
	m_pnts[1] = Point(1,0,0);
	m_pnts[2] = Point(1,1,0);
	m_pnts[3] = Point(0,1,0);
	m_pnts[4] = Point(0,0,1);
	m_pnts[5] = Point(1,0,1);
	m_pnts[6] = Point(1,1,1);
	m_pnts[7] = Point(0,1,1);
}
/*----------------------------------------------------------------------------*/
Hexahedron::Hexahedron(const Point& AP1, const Point& AP2, const Point& AP3, const Point& AP4,
		const Point& AP5, const Point& AP6, const Point& AP7, const Point& AP8)
{
	m_pnts[0] = AP1;
	m_pnts[1] = AP2;
	m_pnts[2] = AP3;
	m_pnts[3] = AP4;
	m_pnts[4] = AP5;
	m_pnts[5] = AP6;
	m_pnts[6] = AP7;
	m_pnts[7] = AP8;
}
/*----------------------------------------------------------------------------*/
Hexahedron::Hexahedron(Point APoints[8])
{
        m_pnts[0] = APoints[0];
        m_pnts[1] = APoints[1];
        m_pnts[2] = APoints[2];
        m_pnts[3] = APoints[3];
        m_pnts[4] = APoints[4];
        m_pnts[5] = APoints[5];
        m_pnts[6] = APoints[6];
        m_pnts[7] = APoints[7];
}
/*----------------------------------------------------------------------------*/
Hexahedron::Hexahedron(const std::vector<Point>& APoints)
{
	if(APoints.size() != 8) {
                throw GMDSException("Hexahedron::Hexahedron APoints not of size 6.");
        }

        m_pnts[0] = APoints[0];
        m_pnts[1] = APoints[1];
        m_pnts[2] = APoints[2];
        m_pnts[3] = APoints[3];
        m_pnts[4] = APoints[4];
        m_pnts[5] = APoints[5];
        m_pnts[6] = APoints[6];
        m_pnts[7] = APoints[7];
}
/*----------------------------------------------------------------------------*/
Hexahedron::Hexahedron(const Hexahedron& AHex)
{
	m_pnts[0] = AHex.m_pnts[0];
	m_pnts[1] = AHex.m_pnts[1];
	m_pnts[2] = AHex.m_pnts[2];
	m_pnts[3] = AHex.m_pnts[3];
	m_pnts[4] = AHex.m_pnts[4];
	m_pnts[5] = AHex.m_pnts[5];
	m_pnts[6] = AHex.m_pnts[6];
	m_pnts[7] = AHex.m_pnts[7];
}
/*----------------------------------------------------------------------------*/
Hexahedron::~Hexahedron(){;}
/*----------------------------------------------------------------------------*/
const Point& Hexahedron::getPoint(const TInt& AIndex) const
{
	return m_pnts[AIndex];
}
/*----------------------------------------------------------------------------*/
const Point Hexahedron::getCenter() const
{
	TCoord coordX = 0.;
	TCoord coordY = 0.;
	TCoord coordZ = 0.;
	
	for(int iPoint=0; iPoint<8; iPoint++) {
		coordX += m_pnts[iPoint].X();
		coordY += m_pnts[iPoint].Y();
		coordZ += m_pnts[iPoint].Z();
	}
	coordX /= 8.;
	coordY /= 8.;
	coordZ /= 8.;	

        return Point(coordX,coordY,coordZ);
}
/*----------------------------------------------------------------------------*/
const double
Hexahedron::getVolume() const
{
        
        Point center(getCenter());
	Pyramid p0(m_pnts[0],m_pnts[1],m_pnts[2],m_pnts[3],center);
	Pyramid p1(m_pnts[7],m_pnts[6],m_pnts[5],m_pnts[4],center);
	Pyramid p2(m_pnts[2],m_pnts[1],m_pnts[5],m_pnts[6],center);
	Pyramid p3(m_pnts[0],m_pnts[3],m_pnts[7],m_pnts[4],center);
	Pyramid p4(m_pnts[1],m_pnts[0],m_pnts[4],m_pnts[5],center);
	Pyramid p5(m_pnts[3],m_pnts[2],m_pnts[6],m_pnts[7],center);
	
	double vol = p0.getVolume() + p1.getVolume() + p2.getVolume() + p3.getVolume() + p4.getVolume() + p5.getVolume();
	return vol;
}
/*----------------------------------------------------------------------------*/
double
Hexahedron::computeScaledJacobian() const
{
	const int neighbors[8][3] = 
	{
		{1,3,4},
		{2,0,5},
		{3,1,6},
		{0,2,7},
		{7,5,0},
		{4,6,1},
		{5,7,2},
		{6,4,3}
	};
	
	double scaledJ[8];
	for (int iVertex = 0; iVertex < 8; iVertex++) {
		
		Matrix<3,3,double> A = this->jacobian(iVertex);
		
		int i0    = neighbors[iVertex][0];
		int i1    = neighbors[iVertex][1];
		int i2    = neighbors[iVertex][2];
		double l0 = math::Vector(m_pnts[iVertex]- m_pnts[i0]).norm2();
		double l1 = math::Vector(m_pnts[iVertex]- m_pnts[i1]).norm2();
		double l2 = math::Vector(m_pnts[iVertex]- m_pnts[i2]).norm2();
		double l012 = sqrt(l0*l1*l2);
		double det = A.det();
		if (l012 == 0 || det == 0)
			scaledJ[iVertex] = 0.;
		else
			scaledJ[iVertex] = det/l012;
    	}
	double scaledJmin = HUGE_VALF;
	for(int iVertex = 0; iVertex < 8; iVertex++) {
		if(scaledJ[iVertex] < scaledJmin) {
			scaledJmin = scaledJ[iVertex];
		}
	}

	if(scaledJmin >  1.) scaledJmin =  1.;
	if(scaledJmin < -1.) scaledJmin = -1.;

	return scaledJmin;	
}
/*----------------------------------------------------------------------------*/
double
Hexahedron::computeNormalizedScaledJacobian() const
{
  return this->computeScaledJacobian();
}
/*----------------------------------------------------------------------------*/
double
Hexahedron::computeMeanRatio() const
{
	throw GMDSException("Hexahedron::computeMeanRatio not available.");

	const int neighbors[8][3] =
        {       
                {1,3,4},
                {2,0,5},
                {3,1,6},
                {0,2,7},
                {7,5,0},
                {4,6,1},
                {5,7,2},
                {6,4,3}
        };
       
	double meanRatio = 0.;
 
        for (int iVertex = 0; iVertex < 8; iVertex++) {
                
                Matrix<3,3,double> A = this->jacobian(iVertex);
                
                int i0    = neighbors[iVertex][0];
                int i1    = neighbors[iVertex][1];
                int i2    = neighbors[iVertex][2];
                double l0 = math::Vector(m_pnts[iVertex]- m_pnts[i0]).norm2();
                double l1 = math::Vector(m_pnts[iVertex]- m_pnts[i1]).norm2();
                double l2 = math::Vector(m_pnts[iVertex]- m_pnts[i2]).norm2();
                double frobA = l0+l1+l2;
                double det = A.det();

		if(frobA == 0.) {
			throw GMDSException("Hexahedron::computeMeanRatio frobA is zero.");
		}

        	meanRatio += (3. * std::cbrt(std::pow(det,2))) / frobA;
        }

	return (1./8.) * meanRatio;	
}
/*----------------------------------------------------------------------------*/
double
Hexahedron::computeMeanEdgeLength() const
{
  double sumLength = 0;
  sumLength += m_pnts[0].distance(m_pnts[1]);
  sumLength += m_pnts[1].distance(m_pnts[2]);
  sumLength += m_pnts[2].distance(m_pnts[3]);
  sumLength += m_pnts[3].distance(m_pnts[4]);
  sumLength += m_pnts[4].distance(m_pnts[5]);
  sumLength += m_pnts[5].distance(m_pnts[6]);
  sumLength += m_pnts[6].distance(m_pnts[7]);
  sumLength += m_pnts[7].distance(m_pnts[4]);
  sumLength += m_pnts[0].distance(m_pnts[4]);
  sumLength += m_pnts[1].distance(m_pnts[5]);
  sumLength += m_pnts[2].distance(m_pnts[6]);
  sumLength += m_pnts[3].distance(m_pnts[7]);
  
  sumLength /= 12.;
  return sumLength;
}
/*----------------------------------------------------------------------------*/
bool
Hexahedron::intersect(const Triangle& ATri, const bool AProper) const
{
	Triangle T1(m_pnts[0], m_pnts[1], m_pnts[2]);
	bool inter = ATri.intersect(T1, AProper);
	if (inter) {
		return true;
	}

	Triangle T2(m_pnts[0], m_pnts[2], m_pnts[3]);
	inter = ATri.intersect(T2, AProper);
	if (inter) {
		return true;
	}

	Triangle T3(m_pnts[4], m_pnts[5], m_pnts[6]);
	inter = ATri.intersect(T3, AProper);
	if (inter) {
		return true;
	}

	Triangle T4(m_pnts[4], m_pnts[6], m_pnts[7]);
	inter = ATri.intersect(T4, AProper);
	if (inter) {
		return true;
	}

	Triangle T5(m_pnts[3], m_pnts[0], m_pnts[4]);
	inter = ATri.intersect(T5, AProper);
	if (inter) {
		return true;
	}

	Triangle T6(m_pnts[3], m_pnts[4], m_pnts[7]);
	inter = ATri.intersect(T6, AProper);
	if (inter) {
		return true;
	}

	Triangle T7(m_pnts[1], m_pnts[2], m_pnts[6]);
	inter = ATri.intersect(T7, AProper);
	if (inter) {
		return true;
	}

	Triangle T8(m_pnts[1], m_pnts[6], m_pnts[5]);
	inter = ATri.intersect(T8, AProper);
	if (inter) {
		return true;
	}

	Triangle T9(m_pnts[0], m_pnts[1], m_pnts[5]);
	inter = ATri.intersect(T9, AProper);
	if (inter) {
		return true;	
	}

	Triangle T10(m_pnts[0], m_pnts[5], m_pnts[4]);
	inter = ATri.intersect(T10, AProper);
	if (inter) {
		return true;
	}

	Triangle T11(m_pnts[2], m_pnts[3], m_pnts[7]);
	inter = ATri.intersect(T11, AProper);
	if (inter) {
		return true;
	}

	Triangle T12(m_pnts[2], m_pnts[7], m_pnts[6]);
	inter = ATri.intersect(T12, AProper);
	if (inter) {
		return true;
	}

	return false;	
}
/*---------------------------------------------------------------------------*/
math::Matrix<3,3,double>
Hexahedron::jacobian(const int iVertex) const
{
	TCoord matValues[3][3];

	switch(iVertex) {
	case 0:
		matValues[0][0] = m_pnts[1].X()-m_pnts[0].X();
		matValues[1][0] = m_pnts[1].Y()-m_pnts[0].Y();
		matValues[2][0] = m_pnts[1].Z()-m_pnts[0].Z();
                matValues[0][1] = m_pnts[3].X()-m_pnts[0].X();
                matValues[1][1] = m_pnts[3].Y()-m_pnts[0].Y();
                matValues[2][1] = m_pnts[3].Z()-m_pnts[0].Z();
                matValues[0][2] = m_pnts[4].X()-m_pnts[0].X();
                matValues[1][2] = m_pnts[4].Y()-m_pnts[0].Y();
                matValues[2][2] = m_pnts[4].Z()-m_pnts[0].Z();
		break;
        case 1:
                matValues[0][0] = m_pnts[1].X()-m_pnts[0].X();
                matValues[1][0] = m_pnts[1].Y()-m_pnts[0].Y();
                matValues[2][0] = m_pnts[1].Z()-m_pnts[0].Z();
                matValues[0][1] = m_pnts[2].X()-m_pnts[1].X();
                matValues[1][1] = m_pnts[2].Y()-m_pnts[1].Y();
                matValues[2][1] = m_pnts[2].Z()-m_pnts[1].Z();
                matValues[0][2] = m_pnts[5].X()-m_pnts[1].X();
                matValues[1][2] = m_pnts[5].Y()-m_pnts[1].Y();
                matValues[2][2] = m_pnts[5].Z()-m_pnts[1].Z();
                break;
        case 2:
                matValues[0][0] = m_pnts[2].X()-m_pnts[3].X();
                matValues[1][0] = m_pnts[2].Y()-m_pnts[3].Y();
                matValues[2][0] = m_pnts[2].Z()-m_pnts[3].Z();
                matValues[0][1] = m_pnts[2].X()-m_pnts[1].X();
                matValues[1][1] = m_pnts[2].Y()-m_pnts[1].Y();
                matValues[2][1] = m_pnts[2].Z()-m_pnts[1].Z();
                matValues[0][2] = m_pnts[6].X()-m_pnts[2].X();
                matValues[1][2] = m_pnts[6].Y()-m_pnts[2].Y();
                matValues[2][2] = m_pnts[6].Z()-m_pnts[2].Z();
                break;
        case 3:
                matValues[0][0] = m_pnts[2].X()-m_pnts[3].X();
                matValues[1][0] = m_pnts[2].Y()-m_pnts[3].Y();
                matValues[2][0] = m_pnts[2].Z()-m_pnts[3].Z();
                matValues[0][1] = m_pnts[3].X()-m_pnts[0].X();
                matValues[1][1] = m_pnts[3].Y()-m_pnts[0].Y();
                matValues[2][1] = m_pnts[3].Z()-m_pnts[0].Z();
                matValues[0][2] = m_pnts[7].X()-m_pnts[3].X();
                matValues[1][2] = m_pnts[7].Y()-m_pnts[3].Y();
                matValues[2][2] = m_pnts[7].Z()-m_pnts[3].Z();
                break;
        case 4:
                matValues[0][0] = m_pnts[5].X()-m_pnts[4].X();
                matValues[1][0] = m_pnts[5].Y()-m_pnts[4].Y();
                matValues[2][0] = m_pnts[5].Z()-m_pnts[4].Z();
                matValues[0][1] = m_pnts[7].X()-m_pnts[4].X();
                matValues[1][1] = m_pnts[7].Y()-m_pnts[4].Y();
                matValues[2][1] = m_pnts[7].Z()-m_pnts[4].Z();
                matValues[0][2] = m_pnts[4].X()-m_pnts[0].X();
                matValues[1][2] = m_pnts[4].Y()-m_pnts[0].Y();
                matValues[2][2] = m_pnts[4].Z()-m_pnts[0].Z();
                break;
        case 5:
                matValues[0][0] = m_pnts[5].X()-m_pnts[4].X();
                matValues[1][0] = m_pnts[5].Y()-m_pnts[4].Y();
                matValues[2][0] = m_pnts[5].Z()-m_pnts[4].Z();
                matValues[0][1] = m_pnts[6].X()-m_pnts[5].X();
                matValues[1][1] = m_pnts[6].Y()-m_pnts[5].Y();
                matValues[2][1] = m_pnts[6].Z()-m_pnts[5].Z();
                matValues[0][2] = m_pnts[5].X()-m_pnts[1].X();
                matValues[1][2] = m_pnts[5].Y()-m_pnts[1].Y();
                matValues[2][2] = m_pnts[5].Z()-m_pnts[1].Z();
                break;
        case 6:
                matValues[0][0] = m_pnts[6].X()-m_pnts[7].X();
                matValues[1][0] = m_pnts[6].Y()-m_pnts[7].Y();
                matValues[2][0] = m_pnts[6].Z()-m_pnts[7].Z();
                matValues[0][1] = m_pnts[6].X()-m_pnts[5].X();
                matValues[1][1] = m_pnts[6].Y()-m_pnts[5].Y();
                matValues[2][1] = m_pnts[6].Z()-m_pnts[5].Z();
                matValues[0][2] = m_pnts[6].X()-m_pnts[2].X();
                matValues[1][2] = m_pnts[6].Y()-m_pnts[2].Y();
                matValues[2][2] = m_pnts[6].Z()-m_pnts[2].Z();
                break;
        case 7:
                matValues[0][0] = m_pnts[6].X()-m_pnts[7].X();
                matValues[1][0] = m_pnts[6].Y()-m_pnts[7].Y();
                matValues[2][0] = m_pnts[6].Z()-m_pnts[7].Z();
                matValues[0][1] = m_pnts[7].X()-m_pnts[4].X();
                matValues[1][1] = m_pnts[7].Y()-m_pnts[4].Y();
                matValues[2][1] = m_pnts[7].Z()-m_pnts[4].Z();
                matValues[0][2] = m_pnts[7].X()-m_pnts[3].X();
                matValues[1][2] = m_pnts[7].Y()-m_pnts[3].Y();
                matValues[2][2] = m_pnts[7].Z()-m_pnts[3].Z();
                break;
	default:
		throw GMDSException("Hexahedron::jacobian invalid vertex number.");
		break;
	}	
	
	return math::Matrix<3,3,double> (matValues);
}
/*---------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& AStr, const Hexahedron& AHex){
	AStr<<"Hexahedron"<<std::endl;
	for(int iPoint=0; iPoint<8; iPoint++) {
		AStr<<std::endl<<AHex.getPoint(iPoint);
	}
        return AStr;
}
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/




