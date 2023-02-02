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
/** \file    Pillowing.h
 *  \author  legoff
 *  \date    22/07/2016
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_PILLOWING_H_
#define GMDS_PILLOWING_H_
/*----------------------------------------------------------------------------*/
#include "GMDS/IG/IG.h"
/*----------------------------------------------------------------------------*/
#include <vector>
#include <map>
#include <set>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
/** \class  Pillowing
 *  \brief  Class providing pillowing operation.
 *  
 *//*----------------------------------------------------------------------------*/
class Pillowing{
public:

    /*------------------------------------------------------------------------*/
    /** \brief Constructor.
     *
     *  \param AMesh the mesh where pillowing operations are performed.
     */
	Pillowing(IGMesh& AMesh);

    /*------------------------------------------------------------------------*/
    /** \brief  Destructor.	*/
	virtual ~Pillowing();

    /*------------------------------------------------------------------------*/
    /** \brief  Pillow the set of regions ARegion along the surfaces defined by
     * 			the set of faces in AFaces.
     *
     * 			This operation has no meaning in 2D.
     *
     * 	\param ARegions The regions to be pillowed
     * 	\param ANewRegions The regions that has been added (forming the pillow layer)
     *  \param AFaces   The faces to be inflated
     *  \param AADisplacement The created nodes will be created with a small displacement
     *			from the original node
     */
	void pillow(
			const std::vector<gmds::Face>& AFaces,
			const std::vector<Region>& ARegions,
			std::vector<Region>& ANewRegions,
			bool ADisplacement=false);

    /*------------------------------------------------------------------------*/
    /** \brief  Pillow the set of regions 
     *
     *                  This operation has no meaning in 2D.
     *
     *  \param ARegions The regions to be pillowed
     *  \param ANewRegions The regions that has been added (forming the pillow layer)
     *  \param AADisplacement The created nodes will be created with a small displacement
     *                  from the original node
     */
        void pillow(
                        std::vector<Region>& ARegions,
                        std::vector<Region>& ANewRegions,
                        bool ADisplacement=false);

    /*------------------------------------------------------------------------*/
    /** \brief  Get the set of newly created internal nodes after pillowing
     *
     */
        std::set<Node> getNewInternalNodes();

protected:

	void pillow (
			std::vector<gmds::FakeFace>& AFakeFaces,
                        std::vector<Region>& ARegions,
                        std::vector<Region>& ANewRegions,
                        bool ADisplacement);

	/*------------------------------------------------------------------------*/
        /** \brief  Check whether the mesh model meets the requirements for the pillowing
         *
         */
	bool checkModel();

	/*------------------------------------------------------------------------*/
	/** \brief  Identify the different shrink sets of regions
 	 *
	 */
	void identifyShrinkSets(const std::vector<gmds::Region>& ARegions);

	void propagateRecurse(gmds::Region ARegion, uint ASetIndex, std::set<gmds::Region>& AUnTreatedRegions, std::set<gmds::Region>& ATreatedRegions);

	/*------------------------------------------------------------------------*/
        /** \brief  Build the fake faces and identify those on the boundary of the shrink sets
         *
         */
	void retrieveBoundary(const std::vector<gmds::Region>& ARegions);

	/*------------------------------------------------------------------------*/
        /** \brief  Create the new nodes
         *
         */
        void createNewNodes();

	/*------------------------------------------------------------------------*/
        /** \brief  Create the new cells
         *
         */
        void createNewCells(std::vector<gmds::Region>& ANewRegions);

	/*------------------------------------------------------------------------*/
        /** \brief  Set the new nodes for the affected regions
         *
         */
        void setNewNodes();
	

	/* tells if the new nodes should be created with a small distance from their antecedent */
	bool m_displacement;	

	/* a mesh */
	IGMesh& m_mesh;

	/* shrink sets */
	std::map<gmds::Region, uint> m_regionsSets;

	/* actual faces and their corresponding fake faces */
	std::map<gmds::Face, gmds::FakeFace> m_faces2FakeFaces;

	/* fake faces and their region adjacencies */
	std::map<gmds::FakeFace, std::vector<gmds::Region> > m_fakeFaces2Regions;

	/* regions and their fake faces */
        std::map<gmds::Region, std::vector<gmds::FakeFace> > m_regions2FakeFaces;

	/* fake faces along which we pillow, and their corresponding regions */
	std::map<gmds::FakeFace, gmds::Region> m_pillowedFakeFaces;

	/* stores whether the fake faces along which we pillow is outward oriented relative to its region */
	std::map<gmds::FakeFace, bool> m_pillowedFakeFacesIsOutward;

	/* nodes of the boundary and their adjacent fake faces */
	std::map<gmds::TCellID, std::vector<gmds::FakeFace> > m_node2FakeFaces;

	/* old nodes to new nodes mapping */
	std::map<uint, std::map<gmds::Node, gmds::Node> > m_nodesOld2New;

	std::set<gmds::Node> m_newNodes;
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_PILLOWING_H_ */
/*----------------------------------------------------------------------------*/
