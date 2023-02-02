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
/** \file    MeshSplitter.h
 *  \author  N. LE GOFF
 *  \date    18/10/2012
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MESHSPLITTER_H_
#define GMDS_MESHSPLITTER_H_
/*----------------------------------------------------------------------------*/
#include "GMDS/IG/IGMesh.h"
#include "GMDS/CAD/GeomManager.h"
#include "GeomMeshIntersectionService.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
/**
 *  \class MeshSplitter
 *
 *  \brief Interface
 */
/*----------------------------------------------------------------------------*/
class MeshSplitter {
public:

    	/*------------------------------------------------------------------------*/
    	/** \brief Constructor.
     	*
     	*  \param AMesh the mesh to split.
     	*/
	MeshSplitter(
			gmds::IGMesh& AMesh,
			gmds::geom::GeomManager& AManager,
			gmds::geom::GeomMeshIntersectionService& AService);

    	/*------------------------------------------------------------------------*/
    	/** \brief  Destructor.
     	*
     	*/
	virtual ~MeshSplitter();

    	/*------------------------------------------------------------------------*/
	/** \brief split the mesh.
	 *
	 */
	void splitMesh();

	/*------------------------------------------------------------------------*/
        /** \brief Refine in order to obtain a mesh with more edges adjacent to a node.
         *
         * \param ANode the node around which more edges are wanted.
         * \param AFace the face on which an edge will be added.
         */
        void addEdges2Node(gmds::Node ANode, gmds::Face AFace);

	/*------------------------------------------------------------------------*/
        /** \brief  Refine quads that have two edges associated to a same curve.
         */
        virtual void refineQuads2EdgesOnCurve();
	virtual void refineQuads2EdgesOnCurveImproved();

protected:

	/*------------------------------------------------------------------------*/
	/** \brief mark the regions that will be split.
	 *
	 */
	virtual void markRegionsToSplit(
			const int& AMarkRegionToSplit,
			const int& AMarkFaceToSplit,
			const int& AMarkEdgeToSplit,
			const int& AMarkNodeToSplit) =0;

	/*------------------------------------------------------------------------*/
        /** \brief mark the entities that will be split in order to add edges around
         *      a node.
         *
         */
        virtual void markEntitiesToSplitForAddEdges(
                        gmds::Node ANode,
                        gmds::Face AFace,
                        const int& AMarkRegionToSplit,
                        const int& AMarkFaceToSplit,
                        const int& AMarkEdgeToSplit,
                        const int& AMarkNodeToSplit) =0;

	/*------------------------------------------------------------------------*/
	/** \brief mark the regions that will be split in order to avoid bad
	 * 		   configurations depending on the refinement scheme used.
	 *
	 */
	virtual void avoidBadConfigurations(
			const int& AMarkRegionToSplit,
			const int& AMarkFaceToSplit,
			const int& AMarkEdgeToSplit,
			const int& AMarkNodeToSplit) =0;

	/*------------------------------------------------------------------------*/
	/** \brief split the mesh.
	 *
	 */
	virtual void split(
			const int& AMarkRegionToSplit,
			const int& AMarkFaceToSplit,
			const int& AMarkEdgeToSplit,
			const int& AMarkNodeToSplit) =0;

	/*------------------------------------------------------------------------*/
	/** \brief clean the mesh, ie remove all adjacency and entities, except N|R|R2N
	 *
	 */
	virtual void cleanMesh(
	const int& AMarkRegionToSplit,
	const int& AMarkFaceToSplit,
	const int& AMarkEdgeToSplit,
	const int& AMarkNodeToSplit);


	/*------------------------------------------------------------------------*/
	/* a mesh */
	gmds::IGMesh& mesh_;

	/* a geometric model */
	gmds::geom::GeomManager& manager_;

	/* service associated to the geometric model */
	gmds::geom::GeomMeshIntersectionService& service_;

};
/*----------------------------------------------------------------------------*/
} // end namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MESHSPLITTER_H_ */
/*----------------------------------------------------------------------------*/
