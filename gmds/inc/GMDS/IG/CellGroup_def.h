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
/** \file    CellGroup_def.h
 *  \author  F. LEDOUX
 *  \date    11/17/2008
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_CELLGROUP_DEF_H_
#define GMDS_CELLGROUP_DEF_H_
/*----------------------------------------------------------------------------*/
/** \class CellGroup
 *  \brief Define a group of cells having the same dimension (nodes, edges,
 *  	   faces, regions). These groups are used as containers for node clouds,
 *  	   lines, surfaces by the mesh class.
 *
 *  	   WARNING: if a cell contained in the cell group is removed, the group
 *  			    keep pointing on it (can induce an invalid pointer).
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/CommonTypes.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
class IGMesh;
/*----------------------------------------------------------------------------*/
  template<typename TCellType>
  class EXPORT_GMDS CellGroup{
public:

	/*------------------------------------------------------------------------*/
    /** \brief  Constructor.
     */
	CellGroup(IGMesh* AMesh, const std::string& AName="")
	: mesh_(AMesh), name_(AName){;}

	/*------------------------------------------------------------------------*/
    /** \brief  Copy constructor.
     */
	CellGroup(const CellGroup<TCellType>& ACGroup)
	: mesh_(ACGroup.mesh_), name_(ACGroup.name_), cells_(ACGroup.cells_){;}


	bool operator==(const CellGroup<TCellType>& ACGroup)
	{
		return  (mesh_=ACGroup.mesh_) &&
				(name_==ACGroup.name_) &&
				(cells_==ACGroup.cells_);
	}

	/*------------------------------------------------------------------------*/
    /** \brief  Destructor.
     */
	virtual ~CellGroup(){;}

	/*------------------------------------------------------------------------*/
    /** \brief  Add a cell into the group. Warning, A cell can be twice inside
     * 			a group.
     *
     *  \param ACell the cell to add
     */
	void add(TCellType ACell);

	/*------------------------------------------------------------------------*/
    /** \brief  Del a cell from the group. Warning, A cell can be twice inside
     * 			a group.
     *
     *  \param ACell the cell to add
     */
	void del(TCellType ACell);


	/*------------------------------------------------------------------------*/
    /** \brief  Indicates if ACell belongs to this.
     *
     *  \param ACell the cell to find
     */
	bool has(TCellType ACell) const
	{
		const unsigned int size = cells_.size();
		for(unsigned int i=0; i<size;i++)
			if(cells_[i]==ACell.getID())
				return true;

		return false;
	}


	/*------------------------------------------------------------------------*/
    /** \brief Indicates if the group is empty.
     *
     *  \return a boolean whose value is true if the the group is empty
     */
	bool empty() const {return cells_.empty();}

    /*------------------------------------------------------------------------*/
    /** \brief  Get cells.
     *
     *  \return the cells contained in this group
     */
	std::vector<TCellType> cells();
    /*------------------------------------------------------------------------*/
    /** \brief  Get the ids of the cells.
     *
     *  \return the ids of the cells contained in this group
     */
	std::vector<TCellID>& cellIDs() {return cells_;}

    /*------------------------------------------------------------------------*/
    /** \brief  Get group name.
     *
     *  \return a STL string which is the group name
     */
	std::string name() const {return name_;}

    /*------------------------------------------------------------------------*/
    /** \brief  Set the group name.
     *
     *  \param AName the new name
     */
	void setName(const std::string& AName) {name_=AName;}

    /*------------------------------------------------------------------------*/
    /** \brief  Gives the number of cells composing this.
     *
     *  \return the number of cells
     */
	int size() const {return cells_.size();}

    /*------------------------------------------------------------------------*/
    /** \brief  Provides Overload the bracket operator
     *
     *  \param i the index of the elt you want to get
     *
     *  \return parameter in location i
     */
	TCellType operator[] (const int& i) ;

    /*------------------------------------------------------------------------*/
    /** \brief  Computes the bounding box of the CellGroup
     *
     *  \param AMinCoords array that will be filled with the lower
     *   coordinates of the bounding box. Must be allocated prior to
     *   the call to this function.
     *  \param AMaxCoords array that will be filled with the lower
     *   coordinates of the bounding box. Must be allocated prior to
     *   the call to this function.
     *
     */
	void computeBoundingBox (TCoord* AMinCoords, TCoord* AMaxCoords);


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

private:

	/** mesh containing the cells of (*this) */
	IGMesh* mesh_;

	/** name of the group */
	std::string name_;

	/** ids of the cells contained in this group*/
	std::vector<TCellID> cells_;
};
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_CELLGROUP_DEF_H_ */
/*----------------------------------------------------------------------------*/
