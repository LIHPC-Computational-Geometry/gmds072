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
/** \file    GeomEntity.h
 *  \author  F. LEDOUX
 *  \date    09/21/2010
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_GEOM_GEOMENTITY_H_
#define GMDS_GEOM_GEOMENTITY_H_
/*----------------------------------------------------------------------------*/
#include <string>
/*----------------------------------------------------------------------------*/
#include <GMDS/Math/Point.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace geom{
/*----------------------------------------------------------------------------*/
/** \class Entity
 *  \brief This class provides services thar are common to all the geometric
 *  	   entities (volume, surface, curve, point).
 */
/*----------------------------------------------------------------------------*/
	class EXPORT_GMDS GeomEntity {
public:

	/*------------------------------------------------------------------------*/
	/** \brief  Default constructor
	 */
	GeomEntity(const std::string& AName = "Unknown entity"):
		name_(AName),is_meshed_(false),isIncludedIn_(false),
		includedGeomEntity_(NULL){;}

	/*------------------------------------------------------------------------*/
	/** \brief  provides the name f the geometrical entity.
	 */
	std::string getName() const {return name_;}

	/*------------------------------------------------------------------------*/
	/** \brief  provides the dimension of the geometrical entity.
	 */
	virtual int getDim() const=0;

	/*------------------------------------------------------------------------*/
	/** \brief  computes the area of the entity.
	 */
	virtual TCoord computeArea() const= 0;

	/*------------------------------------------------------------------------*/
	/** \brief  computes the bounding box
	 *
	 *	\param minXYZ The minimum coordinate of the bounding box.
	 *	\param maxXYZ The maximum coordinate of the bounding box.
	 */
	virtual void computeBoundingBox(TCoord minXYZ[3], TCoord maxXYZ[3]) const=0;

	/*------------------------------------------------------------------------*/
        /** \brief Project the point AP unto the geometric entity.
 	 *
 	 *  \param AP the point to project
         */
    	virtual void project(gmds::math::Point& AP) const =0;

	/*------------------------------------------------------------------------*/
	/** \brief  indicates if this entity is already meshed.
	 *
	 *	\return \a true if it is meshed, \a false otherwise.
	 */
	virtual bool isMeshed() const {return is_meshed_;}

	/*------------------------------------------------------------------------*/
	/** \brief  set the mesh flag
	 *
	 *	\param AB the new value of the mesh flag
	 */
	void setMesh(const bool AB){is_meshed_=AB;}

	/*------------------------------------------------------------------------*/
	/** \brief  indicates if this entity is included in an other entity..
	 *
	 *	\return true if it is, false otherwise.
	 */
	virtual bool isIncludedIn() const {return isIncludedIn_;}

	/*------------------------------------------------------------------------*/
	/** \brief  returns the entity this entity is included in..
	 *
	 *	\return pointer to the .
	 */
	virtual GeomEntity* getIncludedGeomEntity() const {return includedGeomEntity_;}

protected:

	/* name of the entity if it has one */
	std::string name_;
	/* boolean determining whether this entity is meshed */
	bool is_meshed_;

	/* boolean determining whether this entity is included inside an entity
	 * of greater or equal dimension */
	bool isIncludedIn_;

	/* pointer to the entity in which this entity is included. Must of
	 * greater or equal dimension */
	GeomEntity* includedGeomEntity_;
};
/*----------------------------------------------------------------------------*/
} // namespace geom
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_GEOM_GEOMENTITY_H_ */

