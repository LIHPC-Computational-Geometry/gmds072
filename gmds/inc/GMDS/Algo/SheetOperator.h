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
/** \file    SheetOperator.h
 *  \author  F. LEDOUX
 *  \date    11/06/2008
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_SHEETOPERATOR_H_
#define GMDS_SHEETOPERATOR_H_
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IG.h>
/*----------------------------------------------------------------------------*/
#include <vector>
#include <map>
#include <set>
/*----------------------------------------------------------------------------*/
namespace gmds{

/*----------------------------------------------------------------------------*/
/** \class  SheetOperator
 *  \brief  Class gathering sheet operations performed on 2D quad and 3D hex
 *  		meshes
 *//*----------------------------------------------------------------------------*/
class SheetOperator{
public:

    /*------------------------------------------------------------------------*/
    /** \brief Constructor.
     *
     *  \param AMesh the mesh where sheet operations are performed.
     */
	SheetOperator(IGMesh& AMesh);

    /*------------------------------------------------------------------------*/
    /** \brief  Destructor.	*/
	virtual ~SheetOperator();

    /*------------------------------------------------------------------------*/
    /** \brief  Pillow the set of regions AVol along the surfaces defined by
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
	std::vector<Face > pillow(
			std::vector<Region>& ARegions,
			std::vector<Face  >& AFaces,
			std::vector<Region>& ANewRegions,
			bool ADisplacement=false);

    /*------------------------------------------------------------------------*/
    /** \brief  Get the set of newly created internal nodes after pillowing
     *
     */
        std::set<Node> getNewInternalNodes();

protected:
	void updateGeomClassification(Face& AFace);

	/* a mesh */
	IGMesh& m_mesh;

	/* newly created internal nodes after pillowing */
	std::set<gmds::Node> m_newNodesL;
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_SHEETOPERATOR_H_ */
/*----------------------------------------------------------------------------*/
