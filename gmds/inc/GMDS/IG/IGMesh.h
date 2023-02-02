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
/* IGMesh_def.h
 *
 *  Created on: 20 mai 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_IGMESH_H_
#define GMDS_IGMESH_H_
/*----------------------------------------------------------------------------*/
#include <list>
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/CommonTypes.h>
#include <GMDS/Utils/SmartBitVector.h>
#include <GMDS/Utils/IndexedVector.h>
#include <GMDS/Utils/VariableManager.h>
#include <GMDS/Utils/Marks32.h>
#include <GMDS/IG/NodeContainer.h>
#include <GMDS/IG/EdgeContainer.h>
#include <GMDS/IG/FaceContainer.h>
#include <GMDS/IG/RegionContainer.h>
#include <GMDS/IG/CellGroup_def.h>

#include <GMDS/IO/IReader.h>
#ifdef GMDS_PARALLEL
#include <GMDS/Parallel/DistributedCellData.h>
#endif //GMDS_PARALLEL
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
template<class T> class VTKWriter;
class IGMeshDoctor;
namespace geom{
class GeomEntity;
}
/*----------------------------------------------------------------------------*/
/** \class IGMesh
 *
 *  \brief this class represents meshes as general incidence graphs. It is possible
 *  to select the cells and connectivities.
 *
 */
class EXPORT_GMDS IGMesh
{
	friend class IReader<IGMesh>;
public:


	/*------------------------------------------------------------------------*/
	//friend class to access to private data
	friend class VTKWriter<IGMesh>;
	friend class IGMeshDoctor;
	/*------------------------------------------------------------------------*/
	/** \brief friend class to access to protected/private data
	 */
	friend class Region;
	friend class Face;
	friend class Edge;
	friend class Node;

	/*------------------------------------------------------------------------*/
	/** \brief Definition of local iterators
	 */
	typedef  NodeContainer::iterator node_iterator;
	typedef  EdgeContainer::iterator edge_iterator;
	typedef  FaceContainer::iterator face_iterator;
	typedef  RegionContainer::iterator region_iterator;

	typedef Node  	node;
	typedef Edge  	edge;
	typedef Face  	face;
	typedef Region  region;

	typedef CellGroup<Node  > cloud;
	typedef CellGroup<Edge  > line;
	typedef CellGroup<Face  > surface;
	typedef CellGroup<Region> volume;

	typedef std::list<cloud>::iterator 	clouds_iterator;
	typedef std::list<line>::iterator 	 	lines_iterator;
	typedef std::list<surface>::iterator 	surfaces_iterator;
	typedef std::list<volume>::iterator  	volumes_iterator;

	/*------------------------------------------------------------------------*/
	/** \brief Constructor
	 * \param AModel a mesh model defining available cells and  connectivities
	 */
	IGMesh(MeshModel model);

	/*------------------------------------------------------------------------*/
	/** \brief Destructor
	 */
	virtual ~IGMesh();

	/*------------------------------------------------------------------------*/
	/** \brief Accessor to the mesh model
	 *
	 * \return the mesh model
	 */
	MeshModel getModel() const;


	/*------------------------------------------------------------------------*/
	/** \brief Remove all the cells, groups and variables stored in this mesh
	 */
	void clear();

	/*------------------------------------------------------------------------*/
	/** \brief Change the mesh model
	 *
	 * \param the mesh model
	 * \param a boolean that chooses whether to create new entites/adjacencies
	 */
	void changeModel(const MeshModel& AModel, const bool& ACallDoctor=true);

	/*------------------------------------------------------------------------*/
	/** \brief Accessor to the mesh dimension
	 *
	 * \return the mesh dimension
	 */
	inline TInt getDim(){return m_model.getDim();}

	/*------------------------------------------------------------------------*/
	/** \brief return the number of cells in the mesh (per dimension)
	 */
	inline TInt getNbNodes  () const {return m_nodes_container->getNbElements();}
	inline TInt getNbEdges  () const {return m_edges_container->getNbElements();}
	inline TInt getNbFaces  () const {return m_faces_container->getNbElements();}
	inline TInt getNbRegions() const {return m_regions_container->getNbElements();}


	inline TInt getNbTriangles() const {return m_faces_container->m_T2N->size();}
	inline TInt getNbQuadrilaterals() const {return m_faces_container->m_Q2N->size();}

	/*------------------------------------------------------------------------*/
	/** \brief return the max id for dimension ADim
	 *
	 * \param ADim is the dimension of cells we want to retrieve the max id
	 */
	TCellID getMaxLocalID(const TInt& ADim) const;

	/*------------------------------------------------------------------------*/
	/** \brief retrieve a cell from its id
	 *
	 * \param AID the id of the cell we look for
	 */
	template<typename T> EXPORT_GMDS T get(const TCellID& AID);

	/*------------------------------------------------------------------------*/
        /** \brief retrieve the point from its node id.
         *
         * \param AID the id of the node of the point we look for
         */
        EXPORT_GMDS math::Point getPoint(const TCellID& AID);

	/*------------------------------------------------------------------------*/
	/** \brief provdes a STL vector containing all the T-type cells of the mesh
	 *
	 * \param AVec the vector of cells
	 */
	template<typename T> EXPORT_GMDS void getAll(std::vector<T>& AVec);

	/*------------------------------------------------------------------------*/
	/** \brief Tell if there is a T-type cell having id AID
	 *
	 * \param AID the id of the cell we look for
	 */
	template<typename T>  EXPORT_GMDS bool has(const TCellID& AID);

	/*------------------------------------------------------------------------*/
	/** \brief Factory method to create a node
	 *
	 * \param AX X coordinate
	 * \param AY Y coordinate
	 * \param AZ Z coordinate
	 *
	 * \return a node object that encapsulates access to the mesh node
	 */
	Node newNode(const TCoord& AX=0.0, const TCoord& AY=0.0, const TCoord AZ=0.0);
	Node newNode(const math::Point& APnt);

	/*------------------------------------------------------------------------*/
	/** \brief Factory method to create an edge
	 *
	 * \param AN1 Node 1
	 * \param AN2 Node 2
	 * \return an edgece object that encapsulates access to the mesh edge
	 */
	Edge newEdge(const Node& AN1, const Node& AN2);
	Edge newEdge(const TCellID& AN1, const TCellID& AN2);

	/*------------------------------------------------------------------------*/
	/** \brief Factory method to create a triangle
	 *
	 * \param AN1 Node 1
	 * \param AN2 Node 2
	 * \param AN3 Node 3
	 * \return a face object that encapsulates access to the mesh face
	 */
	Face newTriangle(const Node& AN1, const Node& AN2, const Node& AN3);
	Face newTriangle(const TCellID& AN1=NullID, const TCellID& AN2=NullID, const TCellID& AN3=NullID);

	/*------------------------------------------------------------------------*/
	/** \brief Factory method to create a quad
	 *
	 * \param AN1 Node 1
	 * \param AN2 Node 2
	 * \param AN3 Node 3
	 * \param AN4 Node 4
	 * \return a face object that encapsulates access to the mesh face
	 */
	Face newQuad(const Node& AN1, const Node& AN2, const Node& AN3, const Node& AN4);
	Face newQuad(const TCellID& AN1=NullID, const TCellID& AN2=NullID,
			const TCellID& AN3=NullID, const TCellID& AN4=NullID);
	/*------------------------------------------------------------------------*/
		/** \brief Factory method to create a polygon
		 *
		 * \param ANodes ordered vector of nodes
		 * \return a face object that encapsulates access to the mesh face
		 */
	Face newPolygon(const std::vector<Node>& ANodes);
	Face newPolygon(const std::vector<TCellID>& ANodes);
	/*------------------------------------------------------------------------*/
	/** \brief Factory method to create a face
	 *
	 * \param ANodes ordered vector of nodes
	 * \return a face object that encapsulates access to the mesh face
	 */
	Face newFace(const std::vector<Node>& ANodes);
	Face newFace(const std::vector<TCellID>& ANodes);
	/*------------------------------------------------------------------------*/
	/** \brief Factory method to create a tetrahedral element
	 *
	 * \param AN1 Node 1
	 * \param AN2 Node 2
	 * \param AN3 Node 3
	 * \param AN4 Node 4
	 * \return a region object that encapsulates access to the mesh region
	 */
	Region newTet(const Node& AN1, const Node& AN2, const Node& AN3,
			const Node& AN4);
	Region newTet(const TCellID& AN1=NullID, const TCellID& AN2=NullID,
			const TCellID& AN3=NullID, const TCellID& AN4=NullID);

	/*------------------------------------------------------------------------*/
	/** \brief Factory method to create a hexahedral element
	 *
	 * \param AN1 Node 1
	 * \param AN2 Node 2
	 * \param AN3 Node 3
	 * \param AN4 Node 4
	 * \param AN5 Node 5
	 * \param AN6 Node 6
	 * \param AN7 Node 7
	 * \param AN8 Node 8
	 * \return a region object that encapsulates access to the mesh region
	 */
	Region newHex(const Node& AN1, const Node& AN2, const Node& AN3,
				const Node& AN4, const Node& AN5, const Node& AN6,
				const Node& AN7, const Node& AN8);

	Region newHex(const TCellID& AN1=NullID, const TCellID& AN2=NullID,
			const TCellID& AN3=NullID, const TCellID& AN4=NullID,
			const TCellID& AN5=NullID, const TCellID& AN6=NullID,
			const TCellID& AN7=NullID, const TCellID& AN8=NullID);
	/*------------------------------------------------------------------------*/
	/** \brief Factory method to create a pyramid element. The last node is
	 * 		the top of the pyramid
	 *
	 * \param AN1 Node 1
	 * \param AN2 Node 2
	 * \param AN3 Node 3
	 * \param AN4 Node 4
	 * \param AN5 Node 5
	 * \return a region object that encapsulates access to the mesh region
	 */
	Region newPyramid(const Node& AN1, const Node& AN2, const Node& AN3,
			const Node& AN4, const Node& AN5);
	Region newPyramid(const TCellID& AN1=NullID, const TCellID& AN2=NullID,
			const TCellID& AN3=NullID, const TCellID& AN4=NullID,
			const TCellID& AN5=NullID);

	/*------------------------------------------------------------------------*/
	/** \brief Factory method to create a prism3 element.
	 * \param AN1 Node 1
	 * \param AN2 Node 2
	 * \param AN3 Node 3
	 * \param AN4 Node 4
	 * \param AN5 Node 5
	 * \param AN6 Node 6
	 * \return a region object that encapsulates access to the mesh region
	 */
	Region newPrism3(const Node& AN1, const Node& AN2, const Node& AN3,
			const Node& AN4, const Node& AN5, const Node& AN6);
	Region newPrism3(const TCellID& AN1=NullID, const TCellID& AN2=NullID,
			const TCellID& AN3=NullID, const TCellID& AN4=NullID,
			const TCellID& AN5=NullID, const TCellID& AN6=NullID);

	/*------------------------------------------------------------------------*/
	/** \brief Deletion of the node AN from the mesh
	 *  \param AN ANode
	 */
	inline void deleteNode(Node AN){
		m_nodes_container->remove(AN.getID());
		m_node_variable_manager.removeEntry(AN.getID());
	}
	/*------------------------------------------------------------------------*/
	/** \brief Deletion of a node from the mesh
	 *  \param AID A node id
	 */
	inline void deleteNode(TCellID n){
		m_nodes_container->remove(n);
		m_node_variable_manager.removeEntry(n);
	}

	/*------------------------------------------------------------------------*/
	/** \brief Deletion of the edge AE from the mesh
	 * \param AE An edge
	 */
	inline void deleteEdge(Edge e){
		m_edges_container->remove(e.getID());
		m_edge_variable_manager.removeEntry(e.getID());
	}
	/*------------------------------------------------------------------------*/
	/** \brief Deletion of an edge
	 *  \param AID An edge id
	 */
	inline void deleteEdge(TCellID e){
		m_edges_container->remove(e);
		m_edge_variable_manager.removeEntry(e);
	}

	/*------------------------------------------------------------------------*/
	/** \brief Deletion of a face
	 * \param AF AFace
	 */
	inline void deleteFace(Face f){
		m_faces_container->remove(f.getID());
		m_face_variable_manager.removeEntry(f.getID());
	}
	/*------------------------------------------------------------------------*/
	/** \brief Deletion of the face AF from the mesh
	 * \param AF A face id
	 */
	inline void deleteFace(TCellID f){
		m_faces_container->remove(f);
		m_face_variable_manager.removeEntry(f);
	}
	/*------------------------------------------------------------------------*/
	/** \brief Deletion of a region
	 * \param AR A region
	 */
	inline void deleteRegion(Region AR){
		m_regions_container->remove(AR.getID());
		m_region_variable_manager.removeEntry(AR.getID());
	}
	/*------------------------------------------------------------------------*/
	/** \brief Deletion of a region
	 *  \param AR A region id
	 */
	inline void deleteRegion(TCellID AR){
		m_regions_container->remove(AR);
		m_region_variable_manager.removeEntry(AR);
	}

	/*------------------------------------------------------------------------*/
	/** \brief Factory method to get an iterator on the first node
	 *
	 * \return a node iterator
	 */
	node_iterator nodes_begin(){return m_nodes_container->getIterator();}

	/*------------------------------------------------------------------------*/
	/** \brief Factory method to get an iterator on the first node
	 *
	 * \return a node iterator
	 */
	edge_iterator edges_begin(){return m_edges_container->getIterator();}

	/*------------------------------------------------------------------------*/
	/** \brief Factory method to get an iterator on the first face
	 *
	 * \return a face iterator
	 */
	face_iterator faces_begin(){return m_faces_container->getIterator();}

	/*------------------------------------------------------------------------*/
	/** \brief Factory method to get an iterator on the first region
	 *
	 * \return a region iterator
	 */
	region_iterator regions_begin(){return m_regions_container->getIterator();}

	/*------------------------------------------------------------------------*/
    /** \brief  Reserve a mark of the mesh for the cell type T
     *
     *  \return A mark number
     */
	template<typename T> EXPORT_GMDS TInt getNewMark();

	/*------------------------------------------------------------------------*/
    /** \brief Get the marks of a cell 
     *
     *  \return A mark 
     */
	template<typename T> EXPORT_GMDS Marks32 getMarks(const T& ACell);

	/*------------------------------------------------------------------------*/
    /** \brief  Free mark AMarkNumber which was previously reserved with
     * 			getNewMark()
     *
     *  \param AMarkNumber  A mark number
     */
	template<typename T> EXPORT_GMDS void freeMark(const TInt AMarkNumber);

	/*------------------------------------------------------------------------*/
    /** \brief  Invert the mark value for all the cell in the mesh in O(1). This
     * 			is useful for unmark all the cells of a mesh at the end of an
     * 			algorithm
     *
     *  \param AMarkNumber A mark number
     */
	template<typename T> EXPORT_GMDS void negateMaskMark(const TInt AMarkNumber);

	/*------------------------------------------------------------------------*/
    /** \brief  Invert the mark value for all the cell in the mesh in O(n) where
     * 			n is the number of cells in the mesh. It is the only way to
     * 			ensure to have a uniform state for this mark.
     *
     *  \param AMarkNumber A mark number
     */
	template<typename T> EXPORT_GMDS void unmarkAll(const TInt AMarkNumber);

	/*------------------------------------------------------------------------*/
    /** \brief  Test if ACell is marked with mark AMarkNumber.
     *
     *  \param ACell		a cell
     *  \param AMarkNumber 	A mark number
     *
     *  \return A boolean providing the value of mark AMarkNumber for cell ACell
     */
    bool isMarked(const Node& ACell, int AMarkNumber) const;
    bool isMarked(const Edge& ACell, int AMarkNumber) const;
    bool isMarked(const Face& ACell, int AMarkNumber) const;
    bool isMarked(const Region& ACell, int AMarkNumber) const;
	template<typename T> EXPORT_GMDS bool isMarked(const TCellID& ACellID, int AMarkNumber) const;

	/*------------------------------------------------------------------------*/
    /** \brief  Update value of mark AMarkNumber for cell ACell.
     *
     *  \param ACell		A cell
     *  \param AMarkNumber 	A mark number
     *  \param AState		The new value of mark AMarkNumber for ACell
     */
    void markTo(const Node& ACell, int AMarkNumber, bool AState);
    void markTo(const Edge& ACell, int AMarkNumber, bool AState);
    void markTo(const Face& ACell, int AMarkNumber, bool AState);
    void markTo(const Region& ACell, int AMarkNumber, bool AState);
	template<typename T>  EXPORT_GMDS void markTo(const TCellID& ACellID, int AMarkNumber, bool AState);

	/*------------------------------------------------------------------------*/
    /** \brief  Mark cell ACell with mark AMarkNumber.
     * 			Equivalent to markTo(ADart, AMarkNumber, true).
     *
     *  \param ACell		A cell
     *  \param AMarkNumber 	A mark number
     */
    void mark(const Node& ACell, int AMarkNumber);
    void mark(const Edge& ACell, int AMarkNumber);
    void mark(const Face& ACell, int AMarkNumber);
    void mark(const Region& ACell, int AMarkNumber);
	template<typename T> EXPORT_GMDS void mark(const TCellID& ACellID, int AMarkNumber);
	/*------------------------------------------------------------------------*/
    /** \brief  Unmark cell ACell with mark AMarkNumber.
     * 			Equivalent to markTo(ADart, AMarkNumber, false).
     *
     *  \param ACell 		A cell
     *  \param AMarkNumber 	A mark number
     */
    void unmark(const Node& ACell, int AMarkNumber);
    void unmark(const Edge& ACell, int AMarkNumber);
    void unmark(const Face& ACell, int AMarkNumber);
    void unmark(const Region& ACell, int AMarkNumber);
    template<typename T> EXPORT_GMDS void unmark(const TCellID& ACellID, int AMarkNumber);

    /*------------------------------------------------------------------------*/
    /** \brief  Create a new variable attached to a generic cell type
     * 			(GMDS_NODE, GMDS_EDGE, GMDS_FACE or GNDS_REGION)
     *
     *  \param AType the cell type
     *  \param AName the name of the variable. If this name is already used, the
     *  	   variable is not created and an exception is thrown
     *
     *  \return A pointer on the variable
     */
    template<typename T> EXPORT_GMDS Variable<T>* newVariable(ECellType AType,
							      const std::string& AName);

    /*------------------------------------------------------------------------*/
    /** \brief  Returns whether the variable attached to a generic cell type
     * 			(GMDS_NODE, GMDS_EDGE, GMDS_FACE or GNDS_REGION) exists
     *
     *  \param AType the cell type
     *  \param AName the name of the queried variable.
     *
     *  \return A boolean
     */
    bool doesVariableExist(ECellType AType, const std::string& AName);

    /*------------------------------------------------------------------------*/
    /** \brief  Access to variable attached to a generic cell type
     * 			(GMDS_NODE, GMDS_EDGE, GMDS_FACE or GNDS_REGION)
     *
     *  \param AType the cell type
     *  \param AName the name of the variable. If this name does not exist, the
     *  	   an exception is thrown
     *
     *  \return A pointer on the variable. This pointer can be null if the
     *  		specified type is wrong
     */
    template<typename T> EXPORT_GMDS Variable<T>* getVariable(ECellType AType,
							      const std::string& AName);

    /*------------------------------------------------------------------------*/
    /** \brief  Access to variable attached to a generic cell type
     * 			(GMDS_NODE, GMDS_EDGE, GMDS_FACE or GNDS_REGION)
     *
     *  \param AType the cell type
     *  \param AName the name of the variable. If this name does not exist, the
     *  	   an exception is thrown
     *
     *  \return A pointer on the variable. This pointer can be null if the
     *  		specified type is wrong
     */
    std::vector<VariableItf*> getAllVariables(ECellType AType);


    /*------------------------------------------------------------------------*/
    /** \brief  Delete a variable attached to a cell type (GMDS_NODE, GMDS_EDGE,
     * 			GMDS_FACE or GNDS_REGION)
     *
     *  \param AType the cell type
     *  \param AName the name of the variable to be deleted
     */
    void deleteVariable(ECellType AType, const std::string& AName);

    /** \brief  Delete a variable attached to a cell type (GMDS_NODE, GMDS_EDGE,
     * 			GMDS_FACE or GNDS_REGION)
     *
     *  \param AType the cell type
     *  \param AVar a pointer on the variable to be deleted
     */
    void deleteVariable(ECellType AType, VariableItf* AVar);


    /*------------------------------------------------------------------------*/
    /** \brief  Create a new cloud inside the mesh. This cloud allows
     * 			users to gather nodes having common properties. Be careful when
     * 			you remove a node from the mesh. If this cell belongs to a
     * 			surface, this cloud will keep an invalid pointer onto this
     * 			node.
     *
     * 	\param AName the name of the new cloud
     *
     *  \return a cloud
     */
    cloud& newCloud(const std::string& AName);

    /*------------------------------------------------------------------------*/
    /** \brief  Delete a cloud already available in the mesh.
     *
     * 	\param ACloud the cloud to delete
     */
    void deleteCloud(cloud& ACloud);

    /*------------------------------------------------------------------------*/
    /** \brief  return the cloud named AName if it exists, an exception
     * 			otherwise.
     *
     * 	\param AName the name of the cloud
     *
     *  \return a cloud
     */
    cloud& getCloud(const std::string& AName);

    /*------------------------------------------------------------------------*/
    /** \brief  return the AIndex-indexed cloud if it exists, an exception
     * 			otherwise.
     *
     * 	\param AIndex the index of the cloud
     *
     *  \return a cloud
     */
    cloud& getCloud(const unsigned int AIndex);

    /*------------------------------------------------------------------------*/
    /** \brief  return the number of clouds stored in the mesh.
     *
     *  \return the number of clouds
     */
    unsigned int getNbClouds() const;

    /*------------------------------------------------------------------------*/
    /** \brief  return an iterator on the first cloud of the mesh
     *
     *  \return an iterator located on the first cloud of the mesh
     */
    clouds_iterator clouds_begin();

    /*------------------------------------------------------------------------*/
    /** \brief  return an iterator on the last cloud of the mesh
     *
     *  \return an iterator located on the last cloud of the mesh
     */
    clouds_iterator clouds_end();

    /*------------------------------------------------------------------------*/
    /** \brief  Create a new line inside the mesh. This line allows
     * 			users to gather edges having common properties. Be careful when
     * 			you remove a line from the mesh. If this cell belongs to a
     * 			line, this line will keep an invalid pointer onto this
     * 			edge.
     *
     * 	\param AName the name of the new line
     *
     *  \return a line
     */
    line& newLine(const std::string& AName);

    /*------------------------------------------------------------------------*/
    /** \brief  Delete a line already available in the mesh.
     *
     * 	\param ALine the line to delete
     */
    void deleteLine(line& ALine);

    /*------------------------------------------------------------------------*/
    /** \brief  return the line named AName if it exists, an exception
     * 			otherwise.
     *
     * 	\param AName the name of the line
     *
     *  \return a line
     */
    line& getLine(const std::string& AName);

    /*------------------------------------------------------------------------*/
    /** \brief  return the AIndex-indexed line if it exists, an exception
     * 			otherwise.
     *
     * 	\param AIndex the index of the line
     *
     *  \return a line
     */
    line& getLine(const unsigned int AIndex);

    /*------------------------------------------------------------------------*/
    /** \brief  return the number of lines stored in the mesh.
     *
     *  \return the number of lines
     */
    unsigned int getNbLines() const;

    /*------------------------------------------------------------------------*/
    /** \brief  return an iterator on the first line of the mesh
     *
     *  \return an iterator located on the first line of the mesh
     */
    lines_iterator lines_begin();

    /*------------------------------------------------------------------------*/
    /** \brief  return an iterator on the last line of the mesh
     *
     *  \return an iterator located on the last line of the mesh
     */
    lines_iterator lines_end();

    /*------------------------------------------------------------------------*/
    /** \brief  Create a new surface inside the mesh. This surface allows
     * 			users to gather faces having common properties. Be careful when
     * 			you remove a face from the mesh. If this cell belongs to a
     * 			surface, this surface will keep an invalid pointer onto this
     * 			face.
     *
     * 	\param AName the name of the new surface
     *
     *  \return a surface
     */
    surface& newSurface(const std::string& AName);

    /*------------------------------------------------------------------------*/
    /** \brief  Delete a surface already available in the mesh.
     *
     * 	\param ASurf the surface to delete
     */
    void deleteSurface(surface& ASurf);

    /*------------------------------------------------------------------------*/
    /** \brief  return the surface named AName if it exists, an exception
     * 			otherwise.
     *
     * 	\param AName the name of the surface
     *
     *  \return a surface
     */
    surface& getSurface(const std::string& AName);

    /*------------------------------------------------------------------------*/
    /** \brief  return the AIndex-indexed surface if it exists, an exception
     * 			otherwise.
     *
     * 	\param AIndex the index of the surface
     *
     *  \return a surface
     */
    surface& getSurface(const unsigned int AIndex);
    /*------------------------------------------------------------------------*/
    /** \brief  return an iterator on the first surface of the mesh
     *
     *  \return an iterator located on the first surface of the mesh
     */
    surfaces_iterator surfaces_begin();

    /*------------------------------------------------------------------------*/
    /** \brief  return an iterator on the last surface of the mesh
     *
     *  \return an iterator located on the last surface of the mesh
     */
    surfaces_iterator surfaces_end();


    /*------------------------------------------------------------------------*/
    /** \brief  return the number of surfaces stored in the mesh.
     *
     *  \return the number of surfaces
     */
    unsigned int getNbSurfaces() const;

    /*------------------------------------------------------------------------*/
    /** \brief  Create a new volume inside the mesh. This volume allows
     * 			users to gather regions having common properties. Be careful when
     * 			you remove a region from the mesh. If this cell belongs to a
     * 			volume, this volume will keep an invalid pointer onto the
     * 			region.
     *
     * 	\param AName the name of the new volume
     *
     *  \return a volume
     */
    volume& newVolume(const std::string& AName);


    /*------------------------------------------------------------------------*/
    /** \brief  Delete a volume already available in the mesh.
     *
     * 	\param AVol the volume to delete
     */
    void deleteVolume(volume& AVol);

    /*------------------------------------------------------------------------*/
    /** \brief  return the volume named AName if it exists, an exception
     * 			otherwise.
     *
     * 	\param AName the name of the volume
     *
     *  \return a volume
     */
    volume& getVolume(const std::string& AName);


    /*------------------------------------------------------------------------*/
    /** \brief  return the AIndex-indexed volume if it exists, an exception
     * 			otherwise.
     *
     * 	\param AIndex the index of the volume
     *
     *  \return a volume
     */
    volume& getVolume(const unsigned int AIndex);

    /*------------------------------------------------------------------------*/
    /** \brief  return an iterator on the first volume of the mesh
     *
     *  \return an iterator located on the first volume of the mesh
     */
    volumes_iterator volumes_begin();

    /*------------------------------------------------------------------------*/
    /** \brief  return an iterator on the last volume of the mesh
     *
     *  \return an iterator located on the last volume of the mesh
     */
    volumes_iterator volumes_end();

    /*------------------------------------------------------------------------*/
    /** \brief  Returs the nodes common to two faces AF1 and AF2
     *
     * 	\param AF1  a first face
     * 	\param AF2  a second face
     *  \param AVec a vector of nodes incident to both AF1 and AF2
     */
    std::vector<Node> getCommonNodes(const Face& AF1, const Face& AF2);
    void getCommonNodes(const Face& AF1, const Face& AF2, std::vector<Node>& Avec);

    /*------------------------------------------------------------------------*/
    /** \brief  return the number of volumes stored in the mesh.
     *
     *  \return the number of volumes
     */
    unsigned int getNbVolumes() const;

    /*------------------------------------------------------------------------*/
    /** \brief  Initialize the geometry classification. In other words, this
     * 			method attaches a geometry variable for each cell dimension.
     * 			After that, each cell is responsible of the value of the
     * 			geometric classification.
     */
    void initializeGeometryClassification();

    /*------------------------------------------------------------------------*/
    /** \brief  Returns whether the classification variable for dimension ADim
     * 			exists.
     *
     *  \param ADim the dimension of cells we want to know whether the classification
     *  			exists.
     */
    bool doesGeometricClassificationExist(const int ADim);

    /*------------------------------------------------------------------------*/
    /** \brief  Provide the classification variable for dimension ADim
     *
     *  \param ADim the dimension of cells we want the classification
     */
    Variable<geom::GeomEntity*> * getGeometricClassification(const int ADim);

    /*------------------------------------------------------------------------*/
    /** \brief  Provide the geometric entity associated to ACell
     *
     *  \param ACell the mesh cell we want to get the classification
     */
    geom::GeomEntity* getGeometricClassification(const Node ACell);
    geom::GeomEntity* getGeometricClassification(const Edge ACell);
    geom::GeomEntity* getGeometricClassification(const Face ACell);
    geom::GeomEntity* getGeometricClassification(const Region ACell);
    /*------------------------------------------------------------------------*/
    /** \brief  Update the geometric entity associated to ACell
     *
     *  \param ACell the mesh cell we want to modify the classification
     */
    void setGeometricClassification(const Node ACell, geom::GeomEntity* e);
    void setGeometricClassification(const Edge ACell, geom::GeomEntity* e);
    void setGeometricClassification(const Face ACell, geom::GeomEntity* e);
    void setGeometricClassification(const Region ACell, geom::GeomEntity* e);



    /*------------------------------------------------------------------------*/
    /** \brief  Serializes the mesh into stream AStr
     *
     *  \param AStr an output stream where the mesh data is written
     */
    void serialize(std::ostream& AStr) ;
    /*------------------------------------------------------------------------*/
    /** \brief  Unserializes a mesh from stream AStr
     *
     *  \param AStr an input stream where the mesh data is read
     */
    void unserialize(std::istream& AStr);


#ifdef GMDS_PARALLEL

    TInt getPartID() const;

    /*------------------------------------------------------------------------*/
    /** \brief  Set a cell as being (or not) a master cell. If this cell is not
     * 			shared, this operation makes it shared too.
     *
     * \param ACell 	the cell to transform as being a master cell
     * \para  AMaster	transforms the cell as a master (true) or not (false)
     */
    template<typename T> void setMaster(const TCellID& AID, const bool& AMaster = true);
    void setMaster(const Node  & ACell, const bool& AMaster = true);
    void setMaster(const Edge  & ACell, const bool& AMaster = true);
    void setMaster(const Face  & ACell, const bool& AMaster = true);
    void setMaster(const Region& ACell, const bool& AMaster = true);

    /*------------------------------------------------------------------------*/
    /** \brief  Check if a cell is shared (master or slave)
     *
     *  \param ACell a mesh cell
     *
     *  \return true if the cell is shared
     */
    template<typename T> bool isShared(const TCellID& AID) const;
    bool isShared(const Node  & ACell) const;
    bool isShared(const Edge  & ACell) const;
    bool isShared(const Face  & ACell) const;
    bool isShared(const Region& ACell) const;

    /*------------------------------------------------------------------------*/
    /** \brief  Check if a cell is a master cell
     *
     *  \param ACell a mesh cell
     *
     *  \return true if the cell is a master cell
     */
    template<typename T> bool isMaster(const TCellID& AID) const;
    bool isMaster(const Node  & ACell) const;
    bool isMaster(const Edge  & ACell) const;
    bool isMaster(const Face  & ACell) const;
    bool isMaster(const Region& ACell) const;


    /*------------------------------------------------------------------------*/
    /** \brief  Check if a cell is a slave cell
     *
     *  \param ACell a mesh cell
     *
     *  \return true if the cell is a slave cell
     */
    template<typename T> bool isSlave(const TCellID& AID) const;
    bool isSlave(const Node  & ACell) const;
    bool isSlave(const Edge  & ACell) const;
    bool isSlave(const face  & ACell) const;
    bool isSlave(const Region& ACell) const;


    /*------------------------------------------------------------------------*/
    /** \brief  Check if a cell is a slave cell of a master on partition APID.
     *
     *  \param ACell a local mesh cell
     *  \param APID  id of a partition
     *
     *  \return true if the cell master of ACell is on partition APID
     */
    template<typename T>
      bool isSlaveOf(const TCellID& AID, const TInt& APID) const;
    bool isSlaveOf(Node  & ACell, const TInt& APID) const;
    bool isSlaveOf(Edge  & ACell, const TInt& APID) const;
    bool isSlaveOf(Face  & ACell, const TInt& APID) const;
    bool isSlaveOf(Region& ACell, const TInt& APID) const;

    /*------------------------------------------------------------------------*/
    /** \brief make a cell a slaver of the master cell defined by (AMasterID,
     * 		   AMasterPartition)
     *
     *  \param ACell 		the cell that must be shared
     *  \param AMasterID 	the local of the master in AMasterPart
     *  \param AMasterPart 	the part id where the master is located
     */
    template<typename T>
      void setSlave(const TCellID& ACell, const TCellID& AMasterID, const TInt& AMasterPart);
    void setSlave(Node&    ACell, const TCellID& AMasterID, const TInt& AMasterPart);
    void setSlave(Edge&    ACell, const TCellID& AMasterID, const TInt& AMasterPart);
    void setSlave(Face&    ACell, const TCellID& AMasterID, const TInt& AMasterPart);
    void setSlave(Region&  ACell, const TCellID& AMasterID, const TInt& AMasterPart);

    /*------------------------------------------------------------------------*/
    /** \brief make a cell not a slaver knowing the part of its master cell
     *
     *  \param ACell 		the cell that must be shared
     *  \param AMasterPart 	the part id where the master is located
     */
    template<typename T>
      void unsetSlave(const TCellID& ACell, const TInt& AMasterPart);
    void unsetSlave(Node&   ACell, const TInt& AMasterPart);
    void unsetSlave(Edge&   ACell, const TInt& AMasterPart);
    void unsetSlave(Face&   ACell, const TInt& AMasterPart);
    void unsetSlave(Region& ACell, const TInt& AMasterPart);

    /*------------------------------------------------------------------------*/
    /** \brief make a cell not a slaver without knowing the part of its
     * 		   master cell
     *
     *  \param ACell 		the cell that must be shared
     */
    template<typename T> void unsetSlave(const TCellID& ACell);
    void unsetSlave(Node&   ACell);
    void unsetSlave(Edge&   ACell);
    void unsetSlave(Face&   ACell);
    void unsetSlave(Region& ACell);

    /*------------------------------------------------------------------------*/
    /** \brief Get master data. If ACell is not a slave cell, this operation
     * 		   returns false
     *
     *  \param ACell  the cell that we want to collect the data
     *  \param AMPart the part the  master belongs to
     * 	\param AMID   the id of the master in AMPart
     *
     *  \return true if ACell is a slave cell, false otherwise
     */
    template<typename T>
      bool getMasterData(const TCellID&  AID, TInt& AMPart, TCellID& AMID) ;
    bool getMasterData(const Node&   ACell, TInt& AMPart, TCellID& AMID) ;
    bool getMasterData(const Edge&   ACell, TInt& AMPart, TCellID& AMID) ;
    bool getMasterData(const Face&   ACell, TInt& AMPart, TCellID& AMID) ;
    bool getMasterData(const Region& ACell, TInt& AMPart, TCellID& AMID) ;

    /*------------------------------------------------------------------------*/
    /** \brief Get interface with partition APartID.
     *
     *  \param  ACells  the pair of ids we want to get with
     *  							 (local slave id, distant master id)
     *  \param  APartID the partition we want the shared cells with
     *  \param  ADim	the dimension of the cells we work on
     *
     *  \return the collection of shared cells with this partition
     */
    void getInterface(SmartVector<std::pair<TCellID,TCellID> >& ACells,
		      const TInt& APartID, const int& ADim) const;

    /*------------------------------------------------------------------------*/
    /** \brief Clear slave interfaces.
     * */
    void clearInterface(const TInt& APartitionID, const TInt& ADimCell);
    void clearInterface(const TInt& ADimCell);


    void printShareInfo(std::ostream& str);

    /*------------------------------------------------------------------------*/
    /** \brief Add the part APart as being a part containing slave ADim-cell of
     * 		   whose master is on this mesh.
     * \param ADim  the dimension of slave cell
     * \param APart the id of the part containing the slave cell
     * */
    void addSlavePart(const int ADim, const int APart);

    /*------------------------------------------------------------------------*/
    /** \brief Get the parts containing slaves
     * \param ADim  the dimension of slave cell
     */
    std::set<TInt> getSlaveParts(const int ADim) const;
    std::set<TInt> getSlaveParts() const;

#endif //GMDS_PARALLEL

 protected:
    /*------------------------------------------------------------------------*/
    /** \brief  Resize node container when the user exactly know the ids of
     * 		    cells he wants to create. Warning, in this case, the user must
     * 			take care of creating a coherent data structure. Otherwise,
     * 			segmentation faults might occur (not meet yet but it could
     * 			happen).
     *
     * 			Just before resizing all the ids are removed and all the cells
     * 			are erased in their respective memory allocators.
     *
     *  \param  AMaxID the max id stored for ADim-cells
     */
    void clearAndResizeNodeIDContainer(const TInt AMaxID);

    /*------------------------------------------------------------------------*/
    /** \brief  Resize edge container when the user exactly know the ids of
     * 		    cells he wants to create. Warning, in this case, the user must
     * 			take care of creating a coherent data structure. Otherwise,
     * 			segmentation faults might occur (not meet yet but it could
     * 			happen).
     *
     * 			Just before resizing all the ids are removed and all the cells
     * 			are erased in their respective memory allocators.
     *
     *  \param  AMaxID the max id stored for ADim-cells
     */
    void clearAndResizeEdgeIDContainer(const TInt AMaxID);

    /*------------------------------------------------------------------------*/
    /** \brief  Resize face container when the user exactly know the ids of
     * 		    cells he wants to create. Warning, in this case, the user must
     * 			take care of creating a coherent data structure. Otherwise,
     * 			segmentation faults might occur (not meet yet but it could
     * 			happen).
     *
     * 			Just before resizing all the ids are removed and all the cells
     * 			are erased in their respective memory allocators.
     *
     *  \param  AMaxID the max id stored for ADim-cells
     */
    void clearAndResizeFaceIDContainer(const TInt AMaxID);

    /*------------------------------------------------------------------------*/
    /** \brief  Resize region container when the user exactly know the ids of
     * 		    cells he wants to create. Warning, in this case, the user must
     * 			take care of creating a coherent data structure. Otherwise,
     * 			segmentation faults might occur (not meet yet but it could
     * 			happen).
     *
     * 			Just before resizing all the ids are removed and all the cells
     * 			are erased in their respective memory allocators.
     *
     *  \param  AMaxID the max id stored for ADim-cells
     */
    void clearAndResizeRegionIDContainer(const TInt AMaxID);

    /*------------------------------------------------------------------------*/
    /** \brief  Update the mesh id containers. This operations is NECESSARY to
     * 			keep valid meshes when we let the user to insert cells with
     * 			specified ids. Each time a new cell with its specified id is
     * 			added, the corresponding id container is not updated to keep
     * 			a valid index of its free spaces. As a consequence, this
     * 			operation is necessary as being a final step to get a valid
     * 			container.
     */
    void updateIDContainers();

    /*------------------------------------------------------------------------*/
    /** \brief  Add a Node into the mesh. In 2D, only AX and AY are used.
     *
     * 			This operation must be called by reader class only. It is
     * 			particularly true in a distributed memory context. Otherwise,
     * 		    you encounter some troubles with the global ids.
     *
     * 			The size of the id container must be specified before adding
     * 		    nodes with ids. At the end, an update of the id container must
     * 			be done to have a coherent container.
     *
     *  \param AX X coordinate
     *  \param AY Y coordinate
     *  \param AZ Z coordinate
     *  \param AGID global id we want to assign to give to the new cell
     *
     *  \return a pointer on the new Node
     */
    Node newNodeWithID(const TCoord& AX, const TCoord& AY, const TCoord& AZ,
		       const TCellID& AGID);

    /*------------------------------------------------------------------------*/
    /** \brief  Add an edge defined by two vertice ids
     * 			This operation must be called by reader class only. It is
     * 			particularly true in a distributed memory context. Otherwise,
     * 		    you encounter some troubles with the global ids.
     *
     * 			The size of the id container must be specified before adding
     * 		    nodes with ids. At the end, an update of the id container must
     * 			be done to have a coherent container.
     *
     *  \param AV1 first  Node id
     *  \param AV2 second Node id
     *  \param AGID global id we want to assign to give to the new cell
     *
     *  \return a pointer on the new edge
     */
    Edge newEdgeWithID(const TCellID& AV1, const TCellID& AV2,
		       const TCellID& AGID);


    /*------------------------------------------------------------------------*/
    /** \brief  Add a triangle defined by three vertex ids.
     * 			This operation must be called by reader class only. It is
     * 			particularly true in a distributed memory context. Otherwise,
     * 		    you encounter some troubles with the global ids.
     *
     * 			The size of the id container must be specified before adding
     * 		    nodes with ids. At the end, an update of the id container must
     * 			be done to have a coherent container.
     *
     *  \param AN1 Node 1 id
     *  \param AN2 Node 2 id
     *  \param AN3 Node 3 id
     *  \param AGID global id we want to assign to give to the new cell
     *
     *  \return a pointer on the new face
     */
    Face newTriangleWithID(const TCellID& AN1,const TCellID& AN2,
			   const TCellID& AN3, const TCellID& AGID);


    /*------------------------------------------------------------------------*/
    /** \brief  Add a quad defined by four ordered vertices ids.
     * 			This operation must be called by reader class only. It is
     * 			particularly true in a distributed memory context. Otherwise,
     * 		    you encounter some troubles with the global ids.
     *
     * 			The size of the id container must be specified before adding
     * 		    nodes with ids. At the end, an update of the id container must
     * 			be done to have a coherent container.
     *
     *  \param AN1 Node 1 id
     *  \param AN2 Node 2 id
     *  \param AN3 Node 3 id
     *  \param AN4 Node 4 id
     *  \param AGID global id we want to assign to give to the new cell
     *
     *  \return a pointer on the new face
     */
    Face newQuadWithID(const TCellID& AN1, const TCellID& AN2,
		       const TCellID& AN3, const TCellID& AN4, const TCellID& AGID);

    /*------------------------------------------------------------------------*/
    /** \brief  Add a polygon defined by an ordered collection of vertex ids.
     * 			This operation must be called by reader class only. It is
     * 			particularly true in a distributed memory context. Otherwise,
     * 		    you encounter some troubles with the global ids.
     *
     * 			The size of the id container must be specified before adding
     * 		    nodes with ids. At the end, an update of the id container must
     * 			be done to have a coherent container.
     *
     *  \param ANodes a collection of vertex ids
     *  \param AGID global id we want to assign to give to the new cell
     *
     *  \return a pointer on the new face
     */
    Face newPolygonWithID(const std::vector<TCellID>& ANodes, const TCellID& AGID);

    /*------------------------------------------------------------------------*/
    /** \brief  Add a face defined by an ordered collection of vertex ids.
     * 			Contrary to the polygon builder, the most appropriate face is
     * 			built in this case.
     *
     * 			This operation must be called by reader class only. It is
     * 			particularly true in a distributed memory context. Otherwise,
     * 		    you encounter some troubles with the global ids.
     *
     * 			The size of the id container must be specified before adding
     * 		    nodes with ids. At the end, an update of the id container must
     * 			be done to have a coherent container.
     *
     *  \param ANodes a collection of vertex ids
     *  \param AGID global id we want to assign to give to the new cell
     *
     *  \return a pointer on the new face
     */
    Face newFaceWithID(const std::vector<TCellID>& AIDs, const TCellID& AGID);



    /*------------------------------------------------------------------------*/
    /** \brief  Add a tetrahedron defined by four vertices ids
     * 			This operation must be called by reader class only. It is
     * 			particularly true in a distributed memory context. Otherwise,
     * 		    you encounter some troubles with the global ids.
     *
     * 			The size of the id container must be specified before adding
     * 		    nodes with ids. At the end, an update of the id container must
     * 			be done to have a coherent container.
     *
     *  \param AN1 Node 1 id
     *  \param AN2 Node 2 id
     *  \param AN3 Node 3 id
     *  \param AN4 Node 4 id
     *  \param AGID global id we want to assign to give to the new cell
     *
     *  \return a pointer on the new region
     */
    Region newTetWithID(const TCellID& AN1, const TCellID& AN2,
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
     * 			This operation must be called by reader class only. It is
     * 			particularly true in a distributed memory context. Otherwise,
     * 		    you encounter some troubles with the global ids.
     *
     * 			The size of the id container must be specified before adding
     * 		    nodes with ids. At the end, an update of the id container must
     * 			be done to have a coherent container.
     *
     *  \param AN1 Node 1 id
     *  \param AN2 Node 2 id
     *  \param AN3 Node 3 id
     *  \param AN4 Node 4 id
     *  \param AN5 Node 5 id
     *  \param AGID global id we want to assign to give to the new cell
     *
     *  \return a pointer on the new region
     */
    Region newPyramidWithID(const TCellID& AN1, const TCellID& AN2,
			    const TCellID& AN3, const TCellID& AN4,
			    const TCellID& AN5, const TCellID& AGID);



    /*------------------------------------------------------------------------*/
    /** \brief  Add a prism defined by 6 vertices.
     *
     * 			This operation must be called by reader class only. It is
     * 			particularly true in a distributed memory context. Otherwise,
     * 		    you encounter some troubles with the global ids.
     *
     * 			The size of the id container must be specified before adding
     * 		    nodes with ids. At the end, an update of the id container must
     * 			be done to have a coherent container.
     *
     *  \param AN1 Node 1 id
     *  \param AN2 Node 2 id
     *  \param AN3 Node 3 id
     *  \param AN4 Node 4 id
     *  \param AN5 Node 5 id
     *  \param AN6 Node 6 id
     *  \param AGID global id we want to assign to give to the new cell
     *
     *  \return a pointer on the new region
     */
    Region newPrism3WithID(const TCellID& AN1, const TCellID& AN2,
			   const TCellID& AN3, const TCellID& AN4,
			   const TCellID& AN5, const TCellID& AN6,
			   const TCellID& AGID);



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
     * 			This operation must be called by reader class only. It is
     * 			particularly true in a distributed memory context. Otherwise,
     * 		    you encounter some troubles with the global ids.
     *
     * 			The size of the id container must be specified before adding
     * 		    nodes with ids. At the end, an update of the id container must
     * 			be done to have a coherent container.
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
     *  \return a pointer on the new region
     */
    Region newHexWithID(const TCellID& AN1, const TCellID& AN2,
			const TCellID& AN3, const TCellID& AN4,
			const TCellID& AN5, const TCellID& AN6,
			const TCellID& AN7, const TCellID& AN8,
			const TCellID& AGID);

    /*------------------------------------------------------------------------*/
    /** \brief Change the mesh model and build the new entities/adjacnecies
     *
     * \param the mesh model
     */
    void changeModelWithDoctor(const MeshModel& AModel);

    /*------------------------------------------------------------------------*/
    /** \brief Change the mesh model but DO NOT build the new entities/adjacnecies
     *
     * \param the mesh model
     */
    void changeModelWithoutDoctor(const MeshModel& AModel);


 protected:


    /** implemented mesh model */
    MeshModel m_model;

    /** Cells container */
    NodeContainer*   m_nodes_container  ;
    EdgeContainer*   m_edges_container  ;
    FaceContainer*   m_faces_container  ;
    RegionContainer* m_regions_container;

    VariableManager m_node_variable_manager;
    VariableManager m_edge_variable_manager;
    VariableManager m_face_variable_manager;
    VariableManager m_region_variable_manager;

    std::list<cloud>   m_clouds;
    std::list<line>    m_lines;
    std::list<surface> m_surfaces;
    std::list<volume>  m_volumes;

    /** geometrical classification. Useful in some cases,
     *  WARNING it must be activated in the mesh interface.
     */
    Variable<geom::GeomEntity*>* classification[4];

    /** Boolean marks management
     */
    Variable<Marks32>* m_marks[4];

    /* marks currently used*/
    Marks32 m_usedMarks_nodes;
    Marks32 m_usedMarks_edges;
    Marks32 m_usedMarks_faces;
    Marks32 m_usedMarks_regions;
    /* mask of the marks (altered by negateMaskMark()) */
    Marks32 m_maskMarks_nodes;
    Marks32 m_maskMarks_edges;
    Marks32 m_maskMarks_faces;
    Marks32 m_maskMarks_regions;
    /* indicates free marks*/
    TInt m_marks_nodes	[32];
    TInt m_marks_edges	[32];
    TInt m_marks_faces	[32];
    TInt m_marks_regions[32];
    /* number of used marks*/
    TInt m_nbUsedMarks_nodes;
    TInt m_nbUsedMarks_edges;
    TInt m_nbUsedMarks_faces;
    TInt m_nbUsedMarks_regions;

#ifdef GMDS_PARALLEL

    /** the part id of this mesh in distributed parallel context */
    TInt m_part_id;

    /** As this meh can contain master cells, it must known where slave cells
     *  are. This is the purpose of the two following containers.*/
    std::set<TInt> m_slave_cells_parts_ids[4]; //by cell dimension
    std::set<TInt> m_slave_all_part_ids;		//for all cells

    /** For each cell dimension, we associate to a local slave cell (the key), a pair p
     *  where p.first is the part number of the mesh containing the master cell and
     *  p.second is the master id locally to its part.
     */
    std::map<TCellID, std::pair<TInt,TCellID> > m_slave_mapping[4];

    /** In m_interfaces[i][k] are stored the i-cells that are shared (slave or master) with the
     *  part k.
     */
    std::map<TInt, std::vector<TCellID> > m_interfaces[4];
#endif //GMDS_PARALLEL

#ifdef _DEBUG
    TInt m_maxNbUsedMarks_nodes;
    TInt m_maxNbUsedMarks_edges;
    TInt m_maxNbUsedMarks_faces;
    TInt m_maxNbUsedMarks_regions;
#endif // _DEBUG

};
 template<typename T>
   Variable<T>* IGMesh::newVariable(ECellType AType, const std::string& AName)
   {
     Variable<T> *v;
     std::vector<int> ids;
     switch(AType){
     case GMDS_NODE:{
       node_iterator it = nodes_begin();
       for(;!it.isDone();it.next())
	 ids.push_back(it.value().getID());
       v = m_node_variable_manager.newVariable<T>(AName, m_nodes_container->capacity(),&ids);
     }
       break;
     case GMDS_EDGE:{
       edge_iterator it = edges_begin();
       for(;!it.isDone();it.next())
	 ids.push_back(it.value().getID());
       v = m_edge_variable_manager.newVariable<T>(AName, m_edges_container->capacity(),&ids);
     }
       break;
     case GMDS_FACE:{
       face_iterator it = faces_begin();
       for(;!it.isDone();it.next())
	 ids.push_back(it.value().getID());
       v = m_face_variable_manager.newVariable<T>(AName, m_faces_container->capacity(),&ids);
     }
       break;
     case GMDS_REGION:{
       region_iterator it = regions_begin();
       for(;!it.isDone();it.next())
	 ids.push_back(it.value().getID());
       v = m_region_variable_manager.newVariable<T>(AName, m_regions_container->capacity(),&ids);
     }
       break;
     default:
       throw GMDSException("Unmanaged type of cell -> impossible to create a variable");
     }

     return v;
   }
 /*----------------------------------------------------------------------------*/
 template<typename T>
   Variable<T>* IGMesh::getVariable(ECellType AType, const std::string& AName)
   {
     Variable<T> *v;
     std::vector<int> ids;
     switch(AType){
     case GMDS_NODE:{
       v = m_node_variable_manager.getVariable<T>(AName);
     }
       break;
     case GMDS_EDGE:{
       v = m_edge_variable_manager.getVariable<T>(AName);
     }
       break;
     case GMDS_FACE:{
       v = m_face_variable_manager.getVariable<T>(AName);
     }
       break;
     case GMDS_REGION:{
       v = m_region_variable_manager.getVariable<T>(AName);
     }
       break;
     default:
       throw GMDSException("Unmanaged type of cell -> impossible to access to a variable");
     }
     return v;
   }
 /*----------------------------------------------------------------------------*/
 template<> EXPORT_GMDS void IGMesh::getAll<Node>(std::vector<Node>& AVec);
 template<> EXPORT_GMDS void IGMesh::getAll<Edge>(std::vector<Edge>& AVec);
 template<> EXPORT_GMDS void IGMesh::getAll<Face>(std::vector<Face>& AVec);
 template<> EXPORT_GMDS void IGMesh::getAll<Region>(std::vector<Region>& AVec);

 /*----------------------------------------------------------------------------*/
 template<> EXPORT_GMDS bool IGMesh::has<Node>(const TCellID& AID);
 template<> EXPORT_GMDS bool IGMesh::has<Edge>(const TCellID& AID);
 template<> EXPORT_GMDS bool IGMesh::has<Face>(const TCellID& AID);
 template<> EXPORT_GMDS bool IGMesh::has<Region>(const TCellID& AID);
 /*----------------------------------------------------------------------------*/
 template<> EXPORT_GMDS Node IGMesh::get<Node>(const TCellID& AID);
 template<> EXPORT_GMDS Edge IGMesh::get<Edge>(const TCellID& AID);
 template<> EXPORT_GMDS Face IGMesh::get<Face>(const TCellID& AID);
 template<> EXPORT_GMDS Region IGMesh::get<Region>(const TCellID& AID);

 /*----------------------------------------------------------------------------*/
 template<> EXPORT_GMDS void IGMesh::freeMark<Node>(int AMarkNumber);
 template<> EXPORT_GMDS void IGMesh::freeMark<Edge>(int AMarkNumber);
 template<> EXPORT_GMDS void IGMesh::freeMark<Face>(int AMarkNumber);
 template<> EXPORT_GMDS void IGMesh::freeMark<Region>(int AMarkNumber);

 /*----------------------------------------------------------------------------*/
 template<> EXPORT_GMDS void IGMesh::unmarkAll<Node>(const TInt AMarkNumber);
 template<> EXPORT_GMDS void IGMesh::unmarkAll<Edge>(const TInt AMarkNumber);
 template<> EXPORT_GMDS void IGMesh::unmarkAll<Face>(const TInt AMarkNumber);
 template<> EXPORT_GMDS void IGMesh::unmarkAll<Region>(const TInt AMarkNumber);

 /*----------------------------------------------------------------------------*/
 template<> EXPORT_GMDS void IGMesh::unmark<Node>(const TCellID& ACellID, int AMarkNumber);
 template<> EXPORT_GMDS void IGMesh::unmark<Edge>(const TCellID& ACellID, int AMarkNumber);
 template<> EXPORT_GMDS void IGMesh::unmark<Face>(const TCellID& ACellID, int AMarkNumber);
 template<> EXPORT_GMDS void IGMesh::unmark<Region>(const TCellID& ACellID, int AMarkNumber);

 /*----------------------------------------------------------------------------*/
 template<> EXPORT_GMDS void IGMesh::mark<Node>(const TCellID& ACellID, int AMarkNumber);
 template<> EXPORT_GMDS void IGMesh::mark<Edge>(const TCellID& ACellID, int AMarkNumber);
 template<> EXPORT_GMDS void IGMesh::mark<Face>(const TCellID& ACellID, int AMarkNumber);
 template<> EXPORT_GMDS void IGMesh::mark<Region>(const TCellID& ACellID, int AMarkNumber);
 /*----------------------------------------------------------------------------*/
 template<> EXPORT_GMDS void IGMesh::
   markTo<Node>(const TCellID& ACellID, int AMarkNumber, bool AState);
 template<> EXPORT_GMDS void IGMesh::
   markTo<Edge>(const TCellID& ACellID, int AMarkNumber, bool AState);
 template<> EXPORT_GMDS void IGMesh::
   markTo<Face>(const TCellID& ACellID, int AMarkNumber, bool AState);
 template<> EXPORT_GMDS void IGMesh::
   markTo<Region>(const TCellID& ACellID, int AMarkNumber, bool AState);

 /*----------------------------------------------------------------------------*/
 template<> EXPORT_GMDS
   bool IGMesh::isMarked<Node>(const TCellID& ACellID, int AMarkNumber) const;
 template<> EXPORT_GMDS
   bool IGMesh::isMarked<Edge>(const TCellID& ACellID, int AMarkNumber) const;
 template<> EXPORT_GMDS
   bool IGMesh::isMarked<Face>(const TCellID& ACellID, int AMarkNumber) const;
 template<> EXPORT_GMDS
   bool IGMesh::isMarked<Region>(const TCellID& ACellID, int AMarkNumber) const;
 /*----------------------------------------------------------------------------*/
 template<> EXPORT_GMDS void IGMesh::negateMaskMark<Node>(const TInt AMarkNumber);
 template<> EXPORT_GMDS void IGMesh::negateMaskMark<Edge>(const TInt AMarkNumber);
 template<> EXPORT_GMDS void IGMesh::negateMaskMark<Face>(const TInt AMarkNumber);
 template<> EXPORT_GMDS void IGMesh::negateMaskMark<Region>(const TInt AMarkNumber);
 /*----------------------------------------------------------------------------*/
 template<> EXPORT_GMDS TInt IGMesh::getNewMark<Node>();
 template<> EXPORT_GMDS TInt IGMesh::getNewMark<Edge>();
 template<> EXPORT_GMDS TInt IGMesh::getNewMark<Face>();
 template<> EXPORT_GMDS TInt IGMesh::getNewMark<Region>();
 /*----------------------------------------------------------------------------*/
 template<> EXPORT_GMDS Marks32 IGMesh::getMarks<Node>(const Node& ACell);
 template<> EXPORT_GMDS Marks32 IGMesh::getMarks<Edge>(const Edge& ACell);
 template<> EXPORT_GMDS Marks32 IGMesh::getMarks<Face>(const Face& ACell);
 template<> EXPORT_GMDS Marks32 IGMesh::getMarks<Region>(const Region& ACell);
 /*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_IGMESH_H_ */
/*----------------------------------------------------------------------------*/
