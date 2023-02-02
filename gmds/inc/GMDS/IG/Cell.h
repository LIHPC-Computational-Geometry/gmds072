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
/** \file    Cell.h
 *  \author  F. LEDOUX
 *  \date    January 6, 2014
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_CELL_H_
#define GMDS_CELL_H_
/*----------------------------------------------------------------------------*/
// gmds file headers
#include <GMDS/Utils/CommonTypes.h>
#include <GMDS/Utils/Exception.h>
#include "GMDS/Math/Point.h"
/*----------------------------------------------------------------------------*/
// STL file headers
#include <vector>
#include <iostream>
/*----------------------------------------------------------------------------*/
namespace gmds{

  class IGMesh;
  class Node;
  class Edge;
  class Face;
  class Region;
  /*----------------------------------------------------------------------------*/
  /** \class Cell
   *  \brief Defines the functions that are common to any type of cells (nodes,
   *  	   edges, faces, regions) in an Incidence Graph Representation.
   */
  /*----------------------------------------------------------------------------*/
  class EXPORT_GMDS Cell
  {
  public:
      /*------------------------------------------------------------------------*/
      /** \struct Data
       *  \brief nested structure to represent a cell by its dim and it id
       */
      /*------------------------------------------------------------------------*/
      struct Data{
          int dim;
          TCellID id;
          
          Data(const int ADim=-1, const TCellID AID=0):dim(ADim),id(AID){;}
      };
      
    /*------------------------------------------------------------------------*/
    /** \brief  Accessor to the global id. There exists a unique id for every
     * 			cell of a specific dimension but a face and an edge can have the
     * 			same id.
     *
     *  \return an id
     */
    TCellID getID() const;
    /*------------------------------------------------------------------------*/
    /** \brief  Accesor to the cell dim.
     */
    virtual int getDim() const = 0;
    /*------------------------------------------------------------------------*/
    /** \brief  Accesor to the cell type.
     */
    ECellType getType() const;

    /*------------------------------------------------------------------------*/
    /** \brief Accessor th the number of incident nodes, edges, faces and
     * 		   adjacent regions
     */
    virtual TInt getNbNodes()   const = 0;
    virtual TInt getNbEdges()   const = 0;
    virtual TInt getNbFaces()   const = 0;
    virtual TInt getNbRegions() const = 0;

    /*------------------------------------------------------------------------*/
    /** \brief  Compute the center of the cell
     *
     * \return the center of the region
     */
    virtual math::Point center() const = 0;

    /*------------------------------------------------------------------------*/
    /** \brief  Compute the bounding box
     *
     * \param minXYZ the minimum corner of the bounding box
     * \param maxXYZ the maximum corner of the bounding box
     */
    void computeBoundingBox(TCoord minXYZ[3], TCoord maxXYZ[3]) const;

    /*------------------------------------------------------------------------*/
    /** \brief  Accessor to the incident cells. Only the non-null cells are
     * 			provided.
     *
     * 			T can be Node, Edge, Face or Region
     */
    template<class T> EXPORT_GMDS std::vector<T> get() const;
    template<class T> EXPORT_GMDS void get(std::vector<T>&) const;

    /*------------------------------------------------------------------------*/
    /** \brief  Accessor to the ids of the incident cells.
     * 			Only the non-null cells are provided.
     *
     * 			T can be Node, Edge, Face or Region
     */
    template<class T> EXPORT_GMDS std::vector<TCellID> getIDs() const;
    template<class T> EXPORT_GMDS void getIDs(std::vector<TCellID>&) const;

    /*------------------------------------------------------------------------*/
    /** \brief  Accessor to the incident nodes including null cells. In this case
     * 			cells are provided in the storage order
     */
    template<class T> EXPORT_GMDS std::vector<T> getAll() const;
    template<class T> EXPORT_GMDS void getAll(std::vector<T>&) const;

    /*------------------------------------------------------------------------*/
    /** \brief  Accessor to the ids of incident nodes including null cells. In this case
     * 			cells are provided in the storage order
     */
    template<class T> EXPORT_GMDS std::vector<TCellID> getAllIDs() const;
    template<class T> EXPORT_GMDS void getAllIDs(std::vector<TCellID>&) const;


    /*------------------------------------------------------------------------*/
    /** \brief  update adjacency of type T to be equal to ACells.
     *  \param ACells cells to be added
     */
    template<class T> EXPORT_GMDS bool has(TCellID AId);
    template<class T> EXPORT_GMDS bool has(T& AElt);

  /*------------------------------------------------------------------------*/
    /** \brief  update adjacency of type T to be equal to ACells.
     *  \param ACells cells to be added
     */
    template<class T> EXPORT_GMDS void set(const std::vector<T>& ACells);


    template<class T> EXPORT_GMDS void set(const std::vector<TCellID>& ACells);
    /*------------------------------------------------------------------------*/
    /** \brief  update adjacency of type T to insert the new element AElt.
     *  \param AElt the element to be added
     */
    template<class T> EXPORT_GMDS void add(T& AElt);

    /*------------------------------------------------------------------------*/
    /** \brief  update adjacency of type T to insert the new element AElt.
     *  \param AElt the element to be added
     */
    template<class T> EXPORT_GMDS void add(TCellID AElt);

    /*------------------------------------------------------------------------*/
    /** \brief  update adjacency of type T to remove the element AElt.
     *  \param AElt the element to be removed
     */
    template<class T> EXPORT_GMDS  void remove(T& AElt);

    /*------------------------------------------------------------------------*/
    /** \brief  update adjacency of type T to remove the element AElt.
     *  \param AElt the element to be removed
     */
    template<class T> EXPORT_GMDS void remove(TCellID AElt);

    /*------------------------------------------------------------------------*/
    /** \brief  update adjacency of type T to remove all the elements of type T
     */
    template<class T> EXPORT_GMDS  void removeAll();

    /*------------------------------------------------------------------------*/
    /** \brief  replace an incident T-typed cell AC1 by cell AC2 in the
     * 			incident elements
     *
     *  \param AC1 the cell to be replaced
     *  \param AC2 the new cell
     */
    template<class T>  EXPORT_GMDS void replace(T& AC1, T& AC2);

    /*------------------------------------------------------------------------*/
    /** \brief  replace the incident T-typed cell of ID AID11 by AID2 in the
     * 			incident elements
     *
     *  \param AID1 the id of the cell to be replaced
     *  \param AID2 the id of the new cell
     */
    template<class T>  EXPORT_GMDS void replace(TCellID AC1, TCellID AC2);

#ifdef GMDS_PARALLEL
    /*------------------------------------------------------------------------*/
    /** \brief  Gives the partition index of the mesh owner
     *
     *  \return Partition index
     */
    TInt getPartID() const;
#endif //GMDS_PARALLEL
  protected:

    /*------------------------------------------------------------------------*/
    /** \brief  Constructor.
     */
    Cell(IGMesh* AMesh, const ECellType& AType, const TCellID& AID);

    /*------------------------------------------------------------------------*/
    /** \brief  Destructor.
     */
    virtual ~Cell(){ ; }
    /*------------------------------------------------------------------------*/
    /** \brief  Serializes the cell data into stream AStr
     *
     *  \param AStr an output stream where the cell data is written
     */
    void serializeCellData(std::ostream& AStr) const;
    /*------------------------------------------------------------------------*/
    /** \brief  Unserializes the cell data from stream AStr
     *
     *  \param AStr an input stream where the cell data is read from
     */
    void unserializeCellData(std::istream& AStr);

    /*------------------------------------------------------------------------*/
    /** \brief  Accessor to the incident cells. Only the non-null cells are
     * 			provided. T can be Node, Edge, Face or Region.
     */
    virtual void delegateGet(std::vector<Node>&   ACells) const = 0;
    virtual void delegateGet(std::vector<Edge>&   ACells) const = 0;
    virtual void delegateGet(std::vector<Face>&   ACells) const = 0;
    virtual void delegateGet(std::vector<Region>& ACells) const = 0;

    virtual void delegateGetNodeIDs(std::vector<TCellID>& ACells) const = 0;
    virtual void delegateGetEdgeIDs(std::vector<TCellID>& ACells) const = 0;
    virtual void delegateGetFaceIDs(std::vector<TCellID>& ACells) const = 0;
    virtual void delegateGetRegionIDs(std::vector<TCellID>& ACells) const = 0;

    virtual void delegateGetAll(std::vector<Node>&   ACells) const = 0;
    virtual void delegateGetAll(std::vector<Edge>&   ACells) const = 0;
    virtual void delegateGetAll(std::vector<Face>&   ACells) const = 0;
    virtual void delegateGetAll(std::vector<Region>& ACells) const = 0;

    virtual void delegateGetAllNodeIDs(std::vector<TCellID>& ACells) const = 0;
    virtual void delegateGetAllEdgeIDs(std::vector<TCellID>& ACells) const = 0;
    virtual void delegateGetAllFaceIDs(std::vector<TCellID>& ACells) const = 0;
    virtual void delegateGetAllRegionIDs(std::vector<TCellID>& ACells) const = 0;

    virtual void delegateSetNodeIDs(const std::vector<TCellID>& ACells) = 0;
    virtual void delegateSetEdgeIDs(const std::vector<TCellID>& ACells) = 0;
    virtual void delegateSetFaceIDs(const std::vector<TCellID>& ACells) = 0;
    virtual void delegateSetRegionIDs(const std::vector<TCellID>& ACells) = 0;


    virtual void delegateNodeAdd(TCellID AElt) = 0;
    virtual void delegateEdgeAdd(TCellID AElt) = 0;
    virtual void delegateFaceAdd(TCellID AElt) = 0;
    virtual void delegateRegionAdd(TCellID AElt) = 0;

    virtual void delegateNodeRemove(TCellID AElt) = 0;
    virtual void delegateEdgeRemove(TCellID AElt) = 0;
    virtual void delegateFaceRemove(TCellID AElt) = 0;
    virtual void delegateRegionRemove(TCellID AElt) = 0;

    virtual void delegateNodeReplace(TCellID AID1, TCellID AID2) = 0;
    virtual void delegateEdgeReplace(TCellID AID1, TCellID AID2) = 0;
    virtual void delegateFaceReplace(TCellID AID1, TCellID AID2) = 0;
    virtual void delegateRegionReplace(TCellID AID1, TCellID AID2) = 0;

    template<class T> EXPORT_GMDS std::vector<TCellID> convertCellToID(const std::vector<T>& ACells);
  protected:

    /** mesh containing *this*/
    IGMesh* m_owner;

    /** cell type (node, edge, triangle, quad, polygon, tetrahedron, etc.)*/
    ECellType m_type;

    /** cell id locally to a part*/
    TCellID m_id;

    //
  };
  /*----------------------------------------------------------------------------*/
  // IMPLEMENTATION
  /*----------------------------------------------------------------------------*/
  template<class T> struct Type2Dim{};
  template<> struct Type2Dim<Node>  { static const int val = 0; };
  template<> struct Type2Dim<Edge>  { static const int val = 1; };
  template<> struct Type2Dim<Face>  { static const int val = 2; };
  template<> struct Type2Dim<Region>{ static const int val = 3; };
  /*----------------------------------------------------------------------------*/
  template<class T> std::vector<T> Cell::get() const{
    std::vector<T> vec;
    get(vec);
    return vec;
  }
  /*----------------------------------------------------------------------------*/
  template<class T> void Cell::get(std::vector<T>& ACells) const {
    delegateGet(ACells);
  }
  /*----------------------------------------------------------------------------*/
  template<class T> std::vector<TCellID> Cell::getIDs() const{
    std::vector<TCellID> vec;
    getIDs<T>(vec);
    return vec;
  }
  /*----------------------------------------------------------------------------*/
  template<class T> std::vector<T> Cell::getAll() const{
    std::vector<T> vec;
    getAll(vec);
    return vec;
  }
  /*----------------------------------------------------------------------------*/
  template<class T> void Cell::getAll(std::vector<T>& ACells) const {
    delegateGetAll(ACells);
  }
  /*----------------------------------------------------------------------------*/
  template<class T> std::vector<TCellID> Cell::getAllIDs() const{
    std::vector<TCellID> vec;
    getAllIDs<T>(vec);
    return vec;
  }
  /*----------------------------------------------------------------------------*/
  template<class T> void Cell::add(T& AElt) {
    add<T>(AElt.getID());
  }
  /*----------------------------------------------------------------------------*/
  template<class T> bool Cell::has(T& AElt) {
    return has<T>(AElt.getID());
  }
  /*----------------------------------------------------------------------------*/
  template<class T> bool Cell::has(TCellID AElt) {
    std::vector<TCellID> cellsIDs;
    getIDs<T>(cellsIDs);
    for(int iCell=0; iCell<cellsIDs.size(); iCell++) {
      if(cellsIDs[iCell] == AElt) {
        return true;
      }
    }
    return false;
  }
  /*----------------------------------------------------------------------------*/
  template<class T> void Cell::remove(T& AElt) {
    remove<T>(AElt.getID());
  }
  /*----------------------------------------------------------------------------*/
  template<class T> void Cell::removeAll() {
    std::vector<TCellID> cellsIDs;
    getIDs<T>(cellsIDs);
    for(int iCell=0; iCell<cellsIDs.size(); iCell++) {
      remove<T>(cellsIDs[iCell]);
    }
  }
  /*----------------------------------------------------------------------------*/
  template<class T> void Cell::replace(T& AElt1, T& AElt2) {
    replace<T>(AElt1.getID(), AElt2.getID());
  }
  /*----------------------------------------------------------------------------*/
  template<class T> std::vector<TCellID>
    Cell::convertCellToID(const std::vector<T>& ACells)
    {
      std::vector<TCellID> cellIDs;
      cellIDs.resize(ACells.size());
      for (unsigned int i = 0; i < ACells.size(); i++){
	cellIDs[i] = ACells[i].getID();
      }
      return cellIDs;
    }

  template<> EXPORT_GMDS void Cell::set<Node>(const std::vector<Node>& ACells);
  template<> EXPORT_GMDS void Cell::set<Edge>(const std::vector<Edge>& ACells);
  template<> EXPORT_GMDS void Cell::set<Face>(const std::vector<Face>& ACells);
  template<> EXPORT_GMDS void Cell::set<Region>(const std::vector<Region>& ACells);

  template<> EXPORT_GMDS void Cell::getIDs<Node>(std::vector<TCellID>& ACells) const;
  template<> EXPORT_GMDS void Cell::getIDs<Edge>(std::vector<TCellID>& ACells) const;
  template<> EXPORT_GMDS void Cell::getIDs<Face>(std::vector<TCellID>& ACells) const;
  template<> EXPORT_GMDS void Cell::getIDs<Region>(std::vector<TCellID>& ACells) const;

  template<> EXPORT_GMDS void Cell::getAllIDs<Node>(std::vector<TCellID>& ACells) const;
  template<> EXPORT_GMDS void Cell::getAllIDs<Edge>(std::vector<TCellID>& ACells) const;
  template<> EXPORT_GMDS void Cell::getAllIDs<Face>(std::vector<TCellID>& ACells) const;
  template<> EXPORT_GMDS void Cell::getAllIDs<Region>(std::vector<TCellID>& ACells) const;

  template<> EXPORT_GMDS void Cell::add<Node>(TCellID AElt);
  template<> EXPORT_GMDS void Cell::add<Edge>(TCellID AElt);
  template<> EXPORT_GMDS void Cell::add<Face>(TCellID AElt);
  template<> EXPORT_GMDS void Cell::add<Region>(TCellID AElt);

  template<> EXPORT_GMDS void Cell::remove<Node>(TCellID AElt);
  template<> EXPORT_GMDS void Cell::remove<Edge>(TCellID AElt);
  template<> EXPORT_GMDS void Cell::remove<Face>(TCellID AElt);
  template<> EXPORT_GMDS void Cell::remove<Region>(TCellID AElt);

  template<> EXPORT_GMDS void Cell::replace<Node>(TCellID AID1, TCellID AID2);
  template<> EXPORT_GMDS void Cell::replace<Edge>(TCellID AID1, TCellID AID2);
  template<> EXPORT_GMDS void Cell::replace<Face>(TCellID AID1, TCellID AID2);
  template<> EXPORT_GMDS void Cell::replace<Region>(TCellID AID1, TCellID AID2);

  /*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /*GMDS_CELL_H_*/
/*----------------------------------------------------------------------------*/

