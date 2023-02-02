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
/** \file    IReader_def.h
 *  \author  F. LEDOUX
 *  \date    03/17/2009
 */
/*----------------------------------------------------------------------------*/
template<typename TMesh>
class IReader{
public:

    /*------------------------------------------------------------------------*/
    /** \brief  Destructor.	*/
	virtual ~IReader();

protected:


    /*------------------------------------------------------------------------*/
    /** \brief  Constructor.
     *
     *  \param AMesh the mesh we want to write into a file.
     */
	IReader(TMesh& AMesh);

    /*------------------------------------------------------------------------*/
    /** \brief  Prepare the mesh node container to receive node ids such that
     * 			the maximum node id is AID.
     *
     *  \param AID maximum id we add in the mesh
     */
	void specifyMaxNodeID(const TInt& AID);

    /*------------------------------------------------------------------------*/
    /** \brief  Prepare the mesh edge container to receive edge ids such that
     * 			the maximum edge id is AID.
     *
     *  \param AID maximum id we add in the mesh
     */
	void specifyMaxEdgeID(const TInt& AID);

    /*------------------------------------------------------------------------*/
    /** \brief  Prepare the mesh face container to receive face ids such that
     * 			the maximum face id is AID.
     *
     *  \param AID maximum id we add in the mesh
     */
	void specifyMaxFaceID(const TInt& AID);

    /*------------------------------------------------------------------------*/
    /** \brief  Prepare the mesh region container to receive region ids such
     * 			that the maximum region id is AID.
     *
     *  \param AID maximum id we add in the mesh
     */
	void specifyMaxRegionID(const TInt& AID);

    /*------------------------------------------------------------------------*/
    /** \brief  Update the mesh id containers. This operations is NECESSARY to
     * 			keep valid meshes.
     */
	void updateMeshIDContainers();

	/*------------------------------------------------------------------------*/
    /** \brief  Add a Node into the mesh. In 2D, only AX and AY are used.
     *
     *  \param AX X coordinate
     *  \param AY Y coordinate
     *  \param AZ Z coordinate
     *  \param AGID global id we want to assign to give to the new cell
     */

	void newNode(const TCoord& AX, const TCoord& AY, const TCoord& AZ,
				 const TCellID& AGID);

	void newNode(const TCoord& AX, const TCoord& AY, const TCellID& AGID);

	/*------------------------------------------------------------------------*/
	/** \brief  Add an edge defined by two vertice ids.
	 *
	 *  \param AV1 first  Node id
	 *  \param AV2 second Node id
	 *  \param AGID global id we want to assign to give to the new cell
	 */
	void newEdge(const TCellID& AN1, const TCellID& AN2,
				 const TCellID& AGID);

    /*------------------------------------------------------------------------*/
    /** \brief  Add a triangle defined by three vertex ids.
     *
     *  \param AN1 Node 1 id
     *  \param AN2 Node 2 id
     *  \param AN3 Node 3 id
     *  \param AGID global id we want to assign to give to the new cell
	 */
	void newTriangle(const TCellID& AN1, const TCellID& AN2,
					 const TCellID& AN3, const TCellID& AGID);
    /*------------------------------------------------------------------------*/
    /** \brief  Add a quad defined by four ordered vertices ids.
     *
     *  \param AN1 Node 1 id
     *  \param AN2 Node 2 id
     *  \param AN3 Node 3 id
     *  \param AN4 Node 4 id
     *  \param AGID global id we want to assign to give to the new cell
	 */
	void newQuad(const TCellID& AN1, const TCellID& AN2,
				 const TCellID& AN3, const TCellID& AN4, const TCellID& AGID);

	/*------------------------------------------------------------------------*/
    /** \brief  Add a polygon defined by an ordered collection of vertex ids.
     *
     *  \param ANodes a collection of vertex ids
     *  \param AGID global id we want to assign to give to the new cell
	 */
	void newPolygon(std::vector<TCellID>& ANodes, const TCellID& AGID);

	/*------------------------------------------------------------------------*/
    /** \brief  Add a face defined by an ordered collection of vertex ids.
     *
     *  \param ANodes a collection of vertex ids
     *  \param AGID global id we want to assign to give to the new cell
	 */
	void newFace(std::vector<TCellID>& AIDs, const TCellID& AGID);

	/*------------------------------------------------------------------------*/
    /** \brief  Add a tetrahedron defined by four vertices ids
     *
     *  \param AN1 Node 1 id
     *  \param AN2 Node 2 id
     *  \param AN3 Node 3 id
     *  \param AN4 Node 4 id
     *  \param AGID global id we want to assign to give to the new cell
	 */
	void newTet(const TCellID& AN1, const TCellID& AN2,
				const TCellID& AN3, const TCellID& AN4, const TCellID& AGID);

	/*------------------------------------------------------------------------*/
    /** \brief  Add a pyramid defined by 5 vertices id whose the fourth first
     * 			vertices define the square face.
     *
     *			              5
     *
     * 				  	2 ----------- 3
     * 				   /		     /
     *			      /             /
     *			     1 ----------- 4
     *
     *  \param AN1 Node 1 id
     *  \param AN2 Node 2 id
     *  \param AN3 Node 3 id
     *  \param AN4 Node 4 id
     *  \param AN5 Node 5 id
     *  \param AGID global id we want to assign to give to the new cell
	 */
	void newPyramid(const TCellID& AN1, const TCellID& AN2,
					const TCellID& AN3, const TCellID& AN4,
					const TCellID& AN5, const TCellID& AGID);

	/*------------------------------------------------------------------------*/
    /** \brief  Add a hexahedron defined by eight vertices ids whose order is
     *
     * 					2 ----------- 3
     * 				   /|            /|
     *			      / |           / |
     *			     1 ----------- 4  |
     *			     |  |          |  |
     * 				 |	6 ---------|- 7
     * 				 | /		   | /
     *			     |/            |/
     *			     5 ----------- 8
     *
     *  \param AN1 Node 1 id
     *  \param AN2 Node 2 id
     *  \param AN3 Node 3 id
     *  \param AN4 Node 4 id
     *  \param AN5 Node 5 id
     *  \param AN6 Node 6 id
     *  \param AN7 Node 7 id
     *  \param AN8 Node 8 id
     *  \param AGID global id we want to assign to give to the new cell
     *
	 */
	void newHex(const TCellID& AN1, const TCellID& AN2,
				const TCellID& AN3, const TCellID& AN4,
				const TCellID& AN5, const TCellID& AN6,
				const TCellID& AN7, const TCellID& AN8,
				const TCellID& AGID);

protected:

	/* a mesh */
	TMesh& mesh_;

};
/*----------------------------------------------------------------------------*/
