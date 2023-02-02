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
 * FacetedTriangulationService.h
 *
 *  Created on: 1 juil. 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_GEOM_FACETEDTRIANGULATIONSERVICE_H_
#define GMDS_GEOM_FACETEDTRIANGULATIONSERVICE_H_
/*----------------------------------------------------------------------------*/
#include "GMDS/CAD/GeomTriangulationService.h"
#include "GMDS/CAD/FacetedVolume.h"
#include "GMDS/CAD/FacetedSurface.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace geom{
/*----------------------------------------------------------------------------*/
/** \class  GeomModelAbstractFactory
 *  \brief  This interface gathers the factory methods required for every
 *  		geometric model.
 *
 *  \param TBase the basic type used to store geometric coordinates.
 */
/*----------------------------------------------------------------------------*/

	class EXPORT_GMDS FacetedTriangulationService : public GeomTriangulationService {

public:

	/*------------------------------------------------------------------------*/
	/** \brief  perform the triangulation process for a volume
	 *
	 *  \param AVol the volume we want to get a triangulation
	 *  \param ATri the set of triangles to get as a result
	 */
	virtual void perform(GeomVolume& AVol,
						 std::vector<gmds::Face >& ATri);
	/*------------------------------------------------------------------------*/
	/** \brief  perform the triangulation process for a surface
	 *
	 *  \param ASurf the surface we want to get a triangulation
	 *  \param ATri  the set of triangles to get as a result
	 */
	virtual void perform(GeomSurface& ASurf, std::vector<gmds::Face>& ATri);

};
/*----------------------------------------------------------------------------*/
} // namespace geom
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_GEOM_FACETEDTRIANGULATIONSERVICE_H_ */
/*----------------------------------------------------------------------------*/
