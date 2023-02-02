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
 * Ray.h
 *
 *  Created on: 22 oct. 2014
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_RAY_H_
#define GMDS_MATH_RAY_H_
/*----------------------------------------------------------------------------*/
// Gepeto File Headers
#include <GMDS/Utils/CommonTypes.h>
#include <GMDS/Math/Point.h>
#include <GMDS/Math/Vector.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace math{
/*----------------------------------------------------------------------------*/
class Segment;
class Plane;
class Triangle;
/*----------------------------------------------------------------------------*/
/** \class Ray
 *  \brief template class implementing a geometrical ray in space
 */
/*----------------------------------------------------------------------------*/
	class EXPORT_GMDS Ray
{
public:

	/*------------------------------------------------------------------------*/
	/** \brief  constructor.
	 *
	 * \param AP origin point of the ray
	 * \param ADir vector
	 */
	Ray(const Point& AP, const Vector& ADir)
	:m_pnt(AP),m_dir(ADir) {
	  m_isDirUnit = true;
	  m_dirUnit = m_dir.getNormalize();
	}

	/*------------------------------------------------------------------------*/
	/** \brief  constructor.
	 *
	 * \param AP1 a point of the plane
	 * \param AP2 a point of the plane
	 */
	Ray(const Point& AP1, const Point& AP2)
	:m_pnt(AP1),m_dir(Vector(AP1,AP2)) {
	  m_isDirUnit = true;
	  m_dirUnit = m_dir.getNormalize();
	}

        /*------------------------------------------------------------------------*/
        /** \brief  constructor.
         *
         * \param ARay a ray
         */
        Ray(const Ray& ARay)
        :m_pnt(ARay.m_pnt),m_dir(ARay.m_dir),
	  m_isDirUnit(ARay.m_isDirUnit),m_dirUnit(ARay.m_dirUnit) {
        }
	

	/*------------------------------------------------------------------------*/
	/** \brief  Overloaded operator=
	 */
	virtual Ray& operator= (const Ray&);

	/*------------------------------------------------------------------------*/
        /** \brief  Overloaded operator<< for output
         */
        friend EXPORT_GMDS std::ostream& operator<<(std::ostream&, const Ray&);

	/*------------------------------------------------------------------------*/
	/** \brief  Getter for the plane point
	 *
	 * \return a point
	 */
	const Point& getPoint() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Getter for the directional vector
	 *
	 * \return a vector
	 */
	const Vector& getDir() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Getter for the directional unit vector
	 *
	 * \return a vector
	 */
	const Vector& getDirUnit() const; 
    
    /*------------------------------------------------------------------------*/
    /** \brief 2D Intersection with a segment
     *
     * \param[in]  AS the segment we want to get the intersection with
     * \param[out] AP the intersection point if it exists
     * \param[out] AParam the a parameter such that AP = a AS[0]+ (1-a) AS[1]
     *
     * \return true if (*this) intersects AS
     */
    bool intersect2D(const Segment& AS, Point& AP, double& AParam) const;
    
    /*------------------------------------------------------------------------*/
    /** \brief 3D Intersection with a segment
     *
     * \param[in]  AS the segment we want to get the intersection with
     * \param[out] AP the intersection point if it exists
     * \param[out] AParamSeg the a parameter such that AP = a AS[0]+ (1-a) AS[1]
     * \param[out] AParamRay the a parameter along the ray
     *
     * \return true if (*this) intersects AS
     */
    bool intersect3D(const Segment& AS, Point& AP, double& AParamSeg,
                     double& AParamRay) const;

	/*------------------------------------------------------------------------*/
        /** \brief 2D Intersection with another ray
         *
         * \param[IN]  AR the ray we want to get the intersection with
	 * \param[OUT] AP the intersection point if it exists
         *
         * \return true if (*this) intersects AS
         */
        bool intersect2D(const Ray& AR, Point& AP) const;

	/*------------------------------------------------------------------------*/
        /** \brief 3D Intersection with a plane
         *
         * \param[IN]  APlane the plane we want to get the intersection with
         * \param[OUT] AP the intersection point if it exists
         *
         * \return true if (*this) intersects APlane
         */
        bool intersect3D(const Plane& APlane, Point& AP) const;

	/*------------------------------------------------------------------------*/
        /** \brief 3D Intersection with a triangle
         *
         * \param[IN]  ATri the triangle we want to get the intersection with
         * \param[OUT] AP the intersection point if it exists
         *
         * \return true if (*this) intersects ATri
         */
        bool intersect3D(const Triangle& ATri, Point& AP) const;

	/*------------------------------------------------------------------------*/
        /** \brief Compute the orthogonal projection of AP on (*this)
         *
         * \param[IN]  AP the point we want to project onto (*this)
         *
         * \return the projected point
         */
        Point project(const Point& AP) const;

private:

	/* point */
	Point m_pnt;

	/* directional vector  */
	Vector m_dir;

	/* unit directional vector */
	mutable bool m_isDirUnit;
	mutable Vector m_dirUnit;

};
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_RAY_H_ */
/*----------------------------------------------------------------------------*/
