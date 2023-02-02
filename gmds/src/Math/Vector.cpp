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
 * Vector.cpp
 *
 *  Created on: 3 juin 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
// gmds file headers
#include <GMDS/Math/Vector.h>
#include <GMDS/Math/Numerics.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace math{
/*----------------------------------------------------------------------------*/
Vector::Vector(const TCoord& AX, const TCoord& AY, const TCoord& AZ)
{
	m_tab[0]=AX;
	m_tab[1]=AY;
	m_tab[2]=AZ;
}
/*----------------------------------------------------------------------------*/
Vector::Vector(const Point& AP1, const Point& AP2){
	m_tab[0]=AP2.X()-AP1.X();
	m_tab[1]=AP2.Y()-AP1.Y();
	m_tab[2]=AP2.Z()-AP1.Z();
}
/*----------------------------------------------------------------------------*/
Vector::Vector(const Point& AP)
{
	m_tab[0]=AP.X();
	m_tab[1]=AP.Y();
	m_tab[2]=AP.Z();
}
/*----------------------------------------------------------------------------*/
Vector::Vector(TCoord AParam[3])
{
	m_tab[0]=AParam[0];
	m_tab[1]=AParam[1];
	m_tab[2]=AParam[2];
}
/*----------------------------------------------------------------------------*/
Vector::Vector(const Vector& AV)
{
	m_tab[0] = AV.m_tab[0];
	m_tab[1] = AV.m_tab[1];
	m_tab[2] = AV.m_tab[2];
}
/*----------------------------------------------------------------------------*/
Vector::~Vector() {;}
/*----------------------------------------------------------------------------*/
TCoord Vector::sumComponents() const
{
	return m_tab[0] + m_tab[1] + m_tab[2];
}
/*----------------------------------------------------------------------------*/
bool Vector::isColinear(const Vector& AV) const
{
	Vector  v(*this), v2=AV;
	v.normalize();
	v2.normalize();
	TCoord result =fabs(v.dot(v2));
	return (near(result,1e+0));
}
/*----------------------------------------------------------------------------*/
bool Vector::isOrthogonal(const Vector& AV) const
{
	return (dot(AV)==0.0);
}
/*----------------------------------------------------------------------------*/
Vector operator-(const Vector& a, const Vector& b)
{
	return Vector(
			a.m_tab[0] - b.m_tab[0],
			a.m_tab[1] - b.m_tab[1],
			a.m_tab[2] - b.m_tab[2]);

}
/*----------------------------------------------------------------------------*/
Vector operator+(const Vector& a, const Vector& b)
{
	return Vector(
			a.m_tab[0] + b.m_tab[0],
			a.m_tab[1] + b.m_tab[1],
			a.m_tab[2] + b.m_tab[2]);
}
/*----------------------------------------------------------------------------*/
Point operator+(const Point& a, const Vector& b)
{
	return Point(
			a.X() + b.m_tab[0],
			a.Y() + b.m_tab[1],
			a.Z() + b.m_tab[2]);
}
/*----------------------------------------------------------------------------*/
Point operator-(const Point& a, const Vector& b)
{
	return Point(
			a.X() - b.m_tab[0],
			a.Y() - b.m_tab[1],
			a.Z() - b.m_tab[2]);
}
/*----------------------------------------------------------------------------*/
Vector operator*(const Vector& a, const TCoord& k)
{
	return Vector(k*a.m_tab[0],k*a.m_tab[1],k*a.m_tab[2]);
}
/*----------------------------------------------------------------------------*/
Vector operator*(const TCoord& k, const Vector& a)
{
	return Vector(k*a.m_tab[0],k*a.m_tab[1],k*a.m_tab[2]);
}
/*----------------------------------------------------------------------------*/
Vector operator/(const Vector& a, const TCoord& k)
{
	if(k==0)
		throw GMDSException("Division by 0 imposssible");
	return Vector(a.m_tab[0]/k,a.m_tab[1]/k,a.m_tab[2]/k);

}
/*----------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& stream, const Vector& vec)
{
	stream<<"["<<vec.m_tab[0]<<", "<<vec.m_tab[1]<<", "<<vec.m_tab[2]<<"]";
	return stream;
}
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/

