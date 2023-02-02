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
/** \file    GeomVolume.h
 *  \author  F. LEDOUX
 *  \date    09/21/2010
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_GEOM_GEOMVOLUME_H_
#define GMDS_GEOM_GEOMVOLUME_H_
/*----------------------------------------------------------------------------*/
// GMDS File Headers
/*----------------------------------------------------------------------------*/
#include "GMDS/Utils/Exception.h"
#include "GMDS/Utils/CommonTypes.h"
#include "GMDS/CAD/GeomEntity.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace geom{
/*----------------------------------------------------------------------------*/
class GeomPoint;
class GeomCurve;
class GeomSurface;
/*----------------------------------------------------------------------------*/
/** \class GeomVolume
 *  \brief This class describe the services that are required by the
 *  	   mesh to the geometrical model. As a consequence, this interface only
 *  	   contains query methods.
 */
/*----------------------------------------------------------------------------*/
class EXPORT_GMDS GeomVolume : public GeomEntity {

public:


	/*------------------------------------------------------------------------*/
	/** \brief  Constructor
	 */
	GeomVolume(const std::string& AName = "Unknown volume")
	:GeomEntity(AName){;}

	/*------------------------------------------------------------------------*/
	/** \brief  provides the dimension of the geometrical entity.
	 */
	int getDim() const {return 3;}

	/*------------------------------------------------------------------------*/
        /** \brief Project the point AP unto the geometric entity.
         *
         *  \param AP the point to project
         */
        virtual void project(gmds::math::Point& AP) const{
                throw GMDSException("GeomVolume::project not implemented");
        };

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the adjacent points.
	 *
	 *  \param APnt the adjacent points.
	 */
	virtual void get(std::vector<GeomPoint* >& APnt) const=0;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the adjacent curves.
	 *
	 *  \param ACur the adjacent curves.
	 */
	virtual void get(std::vector<GeomCurve* >& ACur) const=0;

	/*------------------------------------------------------------------------*/
	/** \brief  Access to the adjacent surfaces.
	 *
	 *  \param ASurf the adjacent surfaces
	 */
	virtual void get(std::vector<GeomSurface* >& ASurf) const=0;

    /*------------------------------------------------------------------------*/
    /** \brief  Add an adjacent surface
     *
     *  \param ASurf the new adjacent surface to add
     */
    virtual void add(gmds::geom::GeomSurface* ASurf){throw GMDSException("GeomVolume::add Not yet implemented!");}

};
/*----------------------------------------------------------------------------*/
} // namespace geom
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_GEOM_GEOMVOLUME_H_ */

