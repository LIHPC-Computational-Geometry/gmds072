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
 * Hexahedron.h
 *
 *  Created on: 16 oct 2014
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATH_HEXAHEDRON_H_
#define GMDS_MATH_HEXAHEDRON_H_
/*----------------------------------------------------------------------------*/
#include <iostream>
/*----------------------------------------------------------------------------*/
// gmds file headers
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/CommonTypes.h>
#include <GMDS/Math/Point.h>
#include <GMDS/Math/Matrix.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace math{
/*----------------------------------------------------------------------------*/
class Triangle;
/*----------------------------------------------------------------------------*/
/** \class Hexahedron
 *  \brief Defines a hexahedron
 */
/*----------------------------------------------------------------------------*/
	class EXPORT_GMDS Hexahedron
{

public:

	/*------------------------------------------------------------------------*/
	/** \brief  constructor
	 *
	 */
	Hexahedron();

	/*------------------------------------------------------------------------*/
	/** \brief  constructor
  	 *
  	 * \param AP1 a point of the hex
  	 * \param AP2 a point of the hex
  	 * \param AP3 a point of the hex
  	 * \param AP4 a point of the hex
  	 * \param AP5 a point of the hex
  	 * \param AP6 a point of the hex
  	 * \param AP7 a point of the hex
  	 * \param AP8 a point of the hex
  	 *			    7--------------6
  	 * 			   /|		  /|
  	 * 			  / |	  	 / |
  	 * 			 4--------------5  |
  	 * 			 |  |		|  |
  	 * 			 |  3-----------|--2   
  	 * 			 | /  		| /   
  	 * 			 |/		|/
  	 * 			 0--------------1
  	 */
	Hexahedron(const Point& AP1, const Point& AP2, const Point& AP3, const Point& AP4,
		const Point& AP5, const Point& AP6, const Point& AP7, const Point& AP8);

        /*------------------------------------------------------------------------*/
        /** \brief  constructor
         *
         * \param APoints an array of points
         *
         */
	Hexahedron(Point APoints[8]);

	/*------------------------------------------------------------------------*/
        /** \brief  constructor
         *
         * \param APoints a vector of points
         *
         */
        Hexahedron(const std::vector<Point>& APoints);

	/*------------------------------------------------------------------------*/
	/** \brief  copy constructor
	 *
	 * \param AHex a hexahedron
	 *
	 */
	Hexahedron(const Hexahedron& AHex);

	/*------------------------------------------------------------------------*/
	/** \brief  destructor
	 */
	virtual ~Hexahedron();

        /*------------------------------------------------------------------------*/
        /** \brief  operator=
         */
        void operator=(const Hexahedron& AHex);


	/*------------------------------------------------------------------------*/
	/** \brief  Getter for the hexahedron point
	 *
	 * \param AIndex an integer
	 *
	 * \return a point
	 */
	const Point& getPoint(const TInt& AIndex) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Compute the center of the hexahedron
         *
         * \return a point
         */
        const Point getCenter() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Compute the signed volume of the hexahedron. It is computed as 
	 *          the sum of 6 pyramids.
	 *
	 * \return the signed volume of the hexahedron
	 */
	const TCoord getVolume() const;

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the scaled jacobian of the hexahedron.
         *
         * \return the scaled jacobian
         */
        double computeScaledJacobian() const;

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the scaled jacobian of the hexahedron
	 *          normalized between [-1., 1.]
         *
         * \return the scaled jacobian
         */
        double computeNormalizedScaledJacobian() const;

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the mean ratio of the hexahedron.
         *
         * \return the mean ratio
         */
        double computeMeanRatio() const;

	/*------------------------------------------------------------------------*/
        /** \brief  Compute the mean edge length of the hexahedron.
         *
         * \return the mean edge length
         */
        double computeMeanEdgeLength() const;

        /*------------------------------------------------------------------------*/
        /** \brief  predicate indicating if a triangle and a hex intersect each
 	 *          other.
         * \param AT a triangle
         */
        bool intersect(const Triangle& AT, const bool AProper = false) const;

        /*------------------------------------------------------------------------*/
        /** \brief  Overloaded operator<< for output
         */
        friend std::ostream& operator<<(std::ostream&, const Hexahedron&);

protected:
        /*------------------------------------------------------------------------*/
        /** \brief  Return the jacobian matrix at a vertex.
         *
         * \return the jacobian matrix
         */
    math::Matrix<3,3,double> jacobian(const int iVertex) const;
    
    
protected:
	Point m_pnts[8];

};
/*----------------------------------------------------------------------------*/
} // namespace math
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATH_HEXAHEDRON_H_ */
/*----------------------------------------------------------------------------*/
