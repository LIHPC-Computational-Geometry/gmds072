/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * This software is a computer program whose purpose is to provide a set of
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
 * data to be ensured and, more generally, to use and operate it in the 
 * same conditions as regards security. 
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \file    MeshInsertDetail.h
 *  \author  legoff
 *  \date    23/01/2015
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MESHINSERTDETAIL_H_
#define GMDS_MESHINSERTDETAIL_H_
/*----------------------------------------------------------------------------*/
#include "GMDS/IG/IGMesh.h"
#include "GMDS/CAD/GeomManager.h"
#include "GeomMeshIntersectionService.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
/**
 *  \class MeshInsertDetail
 *
 *  \brief Interface
 */
/*----------------------------------------------------------------------------*/
class MeshInsertDetail {
public:

    /*------------------------------------------------------------------------*/
    /** \brief Constructor.
     *
     *  \param AMesh the mesh to build.
     */
	MeshInsertDetail(
			IGMesh& AMesh,
			gmds::geom::GeomManager& AManager,
			gmds::geom::GeomMeshIntersectionService& AService);

    /*------------------------------------------------------------------------*/
    /** \brief  Destructor.	*/
	virtual ~MeshInsertDetail();

//    /*------------------------------------------------------------------------*/
//	/** \brief  check the validity of the mesh model and the geometry model
//	 * 			for the desired algorithm.
//	 *
//	 * 			Currently empty; always returns true.
//	 *
//	 *  \return a boolean
//	 */
//	virtual bool checkValidity() =0;
//
//	/*------------------------------------------------------------------------*/
//	/** \brief  insert the geometric detail
//	 *
//	 */
//	virtual void autoInsert(gmds::geom::GeomVolume<TBase>* AVol) =0;
//
	/*------------------------------------------------------------------------*/
	/** \brief  export the
	 *
	 *  \param AFile the output file
	 */
	virtual void exportMeshVTK(const std::string& AFile);

	/*------------------------------------------------------------------------*/
	/** \brief  export the input model in the VTK format
	 *
	 *  \param AFile the output file
	 */
	virtual void exportModelVTK(const std::string& AFile);

	/*------------------------------------------------------------------------*/
	/** \brief  export the output mesh in the VTK format
	 *
	 *  \param AFile the output file
	 */
	virtual void exportMeshVTK(const std::string& AFile, const gmds::MeshModel& AMeshModel);

	virtual void exportCurvesEdgesVTK(const std::string& AFile) =0;

        /*------------------------------------------------------------------------*/
        /** \brief  execute a pillowing operation of the faces  
 	 * 	    classified on the geometric detail surfaces
	 *
	 *  \param AVol the volume around which the pillowing must be executed
	 */
	virtual void pillow(gmds::geom::GeomVolume* AVol);
	virtual void pillowExt(gmds::geom::GeomVolume* AVol);	

	virtual void pillowIntWithoutVerticesCurves(gmds::geom::GeomVolume* AVol);
	virtual void pillowExtWithoutVerticesCurves(gmds::geom::GeomVolume* AVol);

	virtual void pillowIntWithoutVerticesCurves_bis(gmds::geom::GeomVolume* AVol);
        virtual void pillowExtWithoutVerticesCurves_bis(gmds::geom::GeomVolume* AVol);

	virtual void pillowIntSolitaryRegions(gmds::geom::GeomVolume* AVol);
	virtual void pillowExtSolitaryRegions(gmds::geom::GeomVolume* AVol);

	/*------------------------------------------------------------------------*/
        /** \brief  remmoves the faces that are not associated to a surface
         */
        virtual void removeUnassociatedFaces();

	/*------------------------------------------------------------------------*/
        /** \brief  Display mesh quality mesurements.
         */
        virtual void displayMeshQuality();

protected:

	/* a mesh */
	IGMesh& mesh_;

	/* a geometric model */
	gmds::geom::GeomManager& manager_;

	/* service associated to the geometric model */
	gmds::geom::GeomMeshIntersectionService& service_;

//	/* association between GeomNode and mesh node*/
//	std::vector<id> geom2MeshNode_;
//	std::vector<std::vector<id> > geom2MeshEdge_;
//	std::vector<std::vector<id> > geom2MeshFace_;

};
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MESHINSERTDETAIL_H_ */
/*----------------------------------------------------------------------------*/
