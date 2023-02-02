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
 * IGMeshQualityEvaluation.h
 *
 *  Created on: 10 juin 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_IGMESHQUALITYEVALUATION_H_
#define GMDS_IGMESHQUALITYEVALUATION_H_
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/Node.h>
#include <GMDS/IG/Face.h>
#include <GMDS/IG/Region.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
/** \class IGMeshQualityEvaluation
 *
 *  \brief this class provides some quality criteria to evaluate mesh cells.
 *
 */
	class EXPORT_GMDS IGMeshQualityEvaluation
{
public:

	/*------------------------------------------------------------------------*/
	/** \brief Constructor
	 * \param AModel a mesh model defining available cells and  connectivities
	 */
	IGMeshQualityEvaluation();

	/*------------------------------------------------------------------------*/
	/** \brief Destructor
	 */
	virtual ~IGMeshQualityEvaluation();

	/*------------------------------------------------------------------------*/
	/** \brief compute the minimum angle criterion for a 3D face
	 */
	double minAngle(const Face& AF) const;


	/*------------------------------------------------------------------------*/
	/** \brief compute the maximum angle criterion for a 3D face
	 */
	double maxAngle(const Face& AF) const;

	/*------------------------------------------------------------------------*/
	/** \brief compute aspect ratio of a 3D face. The aspect ratio is equal to
	 * 		   the length of the longest edge divided by the length of the
	 * 		   shortest one.
	 */
	double aspectRatio(const Face& AF) const;

	/*------------------------------------------------------------------------*/
	/** \brief compute the Jacobian criterion for a region
	 */
	double jacobian(const Region& AR) const;

	/*------------------------------------------------------------------------*/
	/** \brief compute the scaled Jacobian criterion for a region
	 */
	double scaledJacobian(const Region& AR) const;

private:
	/*------------------------------------------------------------------------*/
	/** \brief compute the Jacobian criterion for a tetrahedral element
	 */
	double jacobianTet(const Region& AR) const;
	/*------------------------------------------------------------------------*/
	/** \brief compute the Jacobian criterion for a hexahedral element
	 */
	double jacobianHex(const Region& AR) const;

	/*------------------------------------------------------------------------*/
	/** \brief compute the Jacobian criterion for a tetrahedral element
	 */
	double scaledJacobianTet(const Region& AR) const;
	/*------------------------------------------------------------------------*/
	/** \brief compute the Jacobian criterion for a hexahedral element
	 */
	double scaledJacobianHex(const Region& AR) const;

	/*------------------------------------------------------------------------*/
	/** \brief compute the Jacobian criterion at node AN0 respectively to nodes
	 * 		   AN1, AN2 and AN3
	 */
	double jacobian(const Node& AN0,const Node& AN1,const Node& AN2,const Node& AN3) const;
};
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_IGMESHQUALITYEVALUATION_H_ */
/*----------------------------------------------------------------------------*/
