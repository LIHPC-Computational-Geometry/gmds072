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
 * EdgeContainer.h
 *
 *  Created on: 19 mai 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_EDGECONTAINER_H_
#define GMDS_EDGECONTAINER_H_
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/Node.h>
#include <GMDS/IG/Edge.h>
#include <GMDS/IG/Face.h>
#include <GMDS/Utils/SmartBitVector.h>
#include <GMDS/Utils/SmartVector.h>
#ifdef GMDS_PARALLEL
#include <GMDS/Parallel/DistributedCellData.h>
#endif //GMDS_PARALLEL

/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
class IGMesh;
/*----------------------------------------------------------------------------*/
class EXPORT_GMDS EdgeContainer {

	friend class Node;
	friend class Edge;
	friend class Face;
	friend class Region;
	friend class IGMesh;
public:

	EdgeContainer(IGMesh* AMesh);
	virtual ~EdgeContainer();

	/*------------------------------------------------------------------------*/
	/** \brief Add an edge in the container
	 *  \param AN1 first node id
	 *  \param AN2 second node id
	 */
	Edge add(const TCellID& AN1,const TCellID& AN2);

	/*------------------------------------------------------------------------*/
	/** \brief Indicate if this container contains a cell of id AID
	 *
	 *  \param AID a cell id
	 */
	bool has(const TCellID& AID);
	/*------------------------------------------------------------------------*/
	/** \brief Update of the mesh model
	 *
	 *  \param AModel the new model
	 */
	inline void setModel(const MeshModel& AModel) {
		m_model = AModel;
	}

	/*------------------------------------------------------------------------*/
	/** \brief add the containers for accessing to ADim-dimensional cells
	 *
	 *  \param ADim the dimension of the cells we want to store the adjacency
	 *  			to
	 */
	void addConnectivityContainers(const TInt ADim);

	/*------------------------------------------------------------------------*/
	/** \brief remove the containers for accessing to ADim-dimensional cells
	 *
	 *  \param ADim the dimension of the cells we want to suppress the
	 *  			adjacency to
	 */
	void removeConnectivityContainers(const TInt ADim);

	TInt getNbElements() const {return m_edge_ids.size();}
	TCellID getMaxID()   const {return m_edge_ids.top()-1;}
	
	class EXPORT_GMDS iterator{
	public:
		iterator();
		iterator(EdgeContainer* nc);
		Edge value() const;
		void next();
		bool isDone() const;
	private:
		EdgeContainer* m_container;
		SmartBitVector::Iterator m_iterator;
	};

	iterator getIterator();

	/*------------------------------------------------------------------------*/
	/** \brief Get the node infos for the edge of id AID
	 *
	 * \param AID the id of the edge we want to get infos
	 * \param ANbNodes the number of nodes of the edge
	 */
	void getNodesData(const TCellID& AID, int& ANbNodes);

	/*------------------------------------------------------------------------*/
	/** \brief Get the edge infos for the edge of id AID
	 * \param AID the id of the edge we want to get infos
	 * \param ANbEdges the number of adjacent edges
	 */
	void getEdgesData(const TCellID& AID, int& ANbEdges);

	/*------------------------------------------------------------------------*/
	/** \brief Get the face infos for the edge of id AID
	 * \param AID the id of the edge we want to get infos
	 * \param ANbFaces the number of adjacent faces
	 */
	void getFacesData(const TCellID& AID, int& ANbFaces);


	/*------------------------------------------------------------------------*/
	/** \brief Get the region infos for the edge of id AID
	 *
	 * \param AID the id of the edge we want to get infos
	 * \param ANbRegions the number of regions adjacent to the edge
	 */
	void getRegionsData(const TCellID& AID, int& ANbRegions);

	inline void remove(TInt index)
	{
		m_edge_ids.unselect(index);
		if(m_model.has(E2N))
			m_E2N->remove(index);
		if(m_model.has(E2E))
			m_E2E->remove(index);
		if(m_model.has(E2F))
			m_E2F->remove(index);
		if(m_model.has(E2R))
			m_E2R->remove(index);
	}
	void clear();
	void resize(const TInt);

	/*------------------------------------------------------------------------*/
	/** \brief  This method is necessary when you want to regularize the
	 * 			container after several assignements. Indeed assignement method
	 * 			has the advantage to check nothing before insertion. The problem
	 * 			is then that the container is no more coherent. This method fix
	 * 			the container.
	 */
	void update();

	/*------------------------------------------------------------------------*/
	/** \brief serialize (*this) in AStr
	 *
	 * \param AStr an output streammap
	 */
	void serialize(std::ostream& AStr);

	/*------------------------------------------------------------------------*/
	/** \brief unserialize (*this) from AStr
	 *
	 * \param AStr an input stream
	 */
	void unserialize(std::istream& AStr);

	/*------------------------------------------------------------------------*/
        /** \brief return the storage capacity
         *
         * \return the storage capacity
         */
        TInt capacity() {return m_edge_ids.capacity();}

protected:
	Edge buildEdge(const TInt);
	/*------------------------------------------------------------------------*/
	/** \brief Add an edge in the container
	 *  \param AN1 first node id
	 *  \param AN2 second node id
	 *  \param AID the id of the edge
	 */
	Edge add(const TCellID& AN1,const TCellID& AN2, const TCellID& AID);
protected:

	IGMesh* m_mesh;
	/** supported mesh model */
	MeshModel m_model;

	SmartBitVector m_edge_ids;
#ifdef GMDS_PARALLEL
	IndexedVector<DistributedCellData> m_distributed_data;
#endif //GMDS_PARALLEL


	SmartVector<TabCellID<2> >* m_E2N;
	SmartVector<TabCellID<size_undef> >* m_E2E;
	SmartVector<TabCellID<size_undef> >* m_E2F;
	SmartVector<TabCellID<size_undef> >* m_E2R;

};
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_EDGECONTAINER_H_ */
/*----------------------------------------------------------------------------*/




