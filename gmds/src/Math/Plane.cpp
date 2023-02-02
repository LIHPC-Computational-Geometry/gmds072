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
 * Plane.cpp
 *
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/Math/Plane.h>
/*----------------------------------------------------------------------------*/
#include <GMDS/Math/Segment.h>
#include <GMDS/Math/Triangle.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace math{
/*----------------------------------------------------------------------------*/
Plane::Plane()
: isNormalUnit_(false)
{
	
}
/*----------------------------------------------------------------------------*/
Plane::Plane(const Point&AP, const Vector& AN) 
: pnt_(AP),normal_(AN),isNormalUnit_(false)
{

}
/*----------------------------------------------------------------------------*/
Plane::Plane(const Point&AP1, const Point &AP2, const Point &AP3)
: pnt_(AP1)
{
	Vector v1(AP1, AP2);
        Vector  v2(AP1, AP3);
        normal_ = v1.cross(v2);
        isNormalUnit_ = false;
}
/*----------------------------------------------------------------------------*/
Plane::Plane(const Face &AF)
{
	if(AF.getType()!=GMDS_TRIANGLE)
                throw GMDSException("A plane can be computed from triangular mesh faces only");
        std::vector<Node> nodes = AF.get<Node>();

        pnt_ = nodes[0].getPoint();
        Point p1 = nodes[1].getPoint();
        Point p2 = nodes[2].getPoint();
        Vector v1(pnt_, p1);
        Vector v2(pnt_, p2);
        normal_ = v1.cross(v2);
        isNormalUnit_ = false;
}
/*----------------------------------------------------------------------------*/
Plane::Plane(const Triangle& AT)
{

                pnt_ = AT.getPoint(0);
                Point p1 = AT.getPoint(1);
                Point p2 = AT.getPoint(2);
                Vector v1(pnt_, p1);
                Vector v2(pnt_, p2);
                normal_ = v1.cross(v2);
                isNormalUnit_ = false;
}
/*----------------------------------------------------------------------------*/
const Point& Plane::getPoint() const
{
	return pnt_;
}
/*----------------------------------------------------------------------------*/
const Vector& Plane::getNormal() const {
	return normal_;
}
/*----------------------------------------------------------------------------*/
const Vector& Plane::getNormalUnit() const {
	if(isNormalUnit_)
		return normalUnit_;

	isNormalUnit_ = true;
	normalUnit_ = normal_.getNormalize();
	return normalUnit_;
}

/*----------------------------------------------------------------------------*/
Plane& Plane::operator=(const Plane& APlane) {
	if (APlane == *this) {
		return *this;
	}
	set(APlane.getPoint(), APlane.getNormal());
	return *this;
}
/*----------------------------------------------------------------------------*/
bool Plane::operator==(const Plane& AP) const {
	if (&AP == this) {
		return true;
	}
	return AP.getPoint() == getPoint() &&
			(AP.getNormal() == getNormal() ||
					AP.getNormal().getNormalize() == getNormal().getNormalize());
}
/*----------------------------------------------------------------------------*/
bool Plane::operator!=(const Plane& AP) const {
	if (&AP == this) {
		return false;
	}
	return AP.getPoint() != getPoint() ||
			(AP.getNormal().getNormalize() != getNormal().getNormalize());
}
/*----------------------------------------------------------------------------*/
bool  Plane::isIn(const Point& AP) const
{
	return isZero(normal_.dot(Vector(pnt_, AP)));
}
/*----------------------------------------------------------------------------*/
bool Plane::isStrictlyOnLeft(const Point& AP) const
{
	return (normal_.dot(Vector(pnt_, AP)) < 0.0);
}
/*----------------------------------------------------------------------------*/
TCoord Plane::distance(const Point& AP) const
{
	Vector v(pnt_, AP);
	return fabs(getNormalUnit().dot(v));
}

/*----------------------------------------------------------------------------*/
Point Plane::project(const Point& AP) const
{
	Vector  v(pnt_, AP);
	TCoord num = v.dot(normal_);
	TCoord det = normal_.norm2();

	TCoord alpha = -num / det;

	return AP + alpha * normal_;
}

/*----------------------------------------------------------------------------*/
bool Plane::intersect(const Segment& AS, const bool AProper) const
{
	if (isIn(AS.getPoint(0)) && isIn(AS.getPoint(1)))
		return (AProper) ? false : true;

	if((this->isStrictlyOnLeft(AS.getPoint(0)) == true &&
			this->isStrictlyOnLeft(AS.getPoint(1)) == false) ||
			(this->isStrictlyOnLeft(AS.getPoint(0)) == false &&
					this->isStrictlyOnLeft(AS.getPoint(1)) == true))
		return true;

	return false;
}
/*----------------------------------------------------------------------------*/
bool Plane::intersect(const Segment& AS, Point &PI, const bool AProper) const
{
	Vector OP1 (AS.getPoint(0).X(), AS.getPoint(0).Y(), AS.getPoint(0).Z());
	Vector P1P2(AS.getPoint(0), AS.getPoint(1));

	TCoord a, b, c, d;
	getEquationCoeffs(a, b, c, d);

	TCoord t = (d - OP1.dot(normal_)) / (P1P2.dot(normal_));

	if (t < 0.0 || t > 1.0) {
		std::cout<<"plane "<<this->getPoint()<<std::endl;
		std::cout<<"plane "<<this->getNormal()<<std::endl;
		std::cout<<"segment "<<AS.getPoint(0)<<std::endl;
		std::cout<<"segment "<<AS.getPoint(1)<<std::endl;

	} else if (isZero(t)) {
		PI = AS.getPoint(0);
	} else if (isZero(t-1.0)) {
		PI = AS.getPoint(1);
	} else {
		PI = AS.getPoint(0) + t * P1P2;
	}
	return true;
}
/*----------------------------------------------------------------------------*/
bool Plane::intersect(const Plane& AP, const bool AProper) const
{
	//check parallelism
	if (isZero(normal_.cross(AP.normal_).norm2()))
	{
		//we have parallel planes
		if (isZero(normal_.dot(Vector(pnt_, AP.pnt_))))
			return (AProper) ? false : true; // the same plane !
		else
			return false; //different planes
	}
	else // now planes intersect each others
		return true;
}
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
