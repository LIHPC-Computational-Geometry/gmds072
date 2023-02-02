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
/** \file    CellGroup.h
 *  \author  N. Le Goff
 *  \date    16/09/2013
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_CELLGROUP_H_
#define GMDS_CELLGROUP_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <vector>
#include <string>
#include <math.h>
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/CommonTypes.h>
#include <GMDS/IG/IGMesh.h>
#include <GMDS/IG/Cell.h>
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/CellGroup_def.h>
/*----------------------------------------------------------------------------*/
namespace gmds {
/*----------------------------------------------------------------------------*/
template<typename TCellType>
void CellGroup<TCellType>::add(TCellType ACell)
{
	cells_.push_back(ACell.getID());
}
/*----------------------------------------------------------------------------*/
template<typename TCellType>
void CellGroup<TCellType>::del(TCellType ACell)
{
	const unsigned int size = cells_.size();
	bool notfound = true;
	unsigned int k=0;
	TCellID cell_id = ACell.getID();
	for(unsigned int i=0; i<size && notfound;i++)
	{
		if(cells_[i]==cell_id)
		{
			notfound = false;
			k=i;
		}
	}
	if(!notfound){
		cells_[k]=cells_[size-1];
		cells_.pop_back();
	}

}
/*----------------------------------------------------------------------------*/
template<typename TCellType>
void CellGroup<TCellType>::computeBoundingBox (TCoord* AMinCoords, TCoord* AMaxCoords)
{
	if(mesh_->getModel().has(N)) {
		unsigned int geomDimension = mesh_->getDim();

		AMinCoords[0] = HUGE_VALF;
		if(geomDimension>1) {
			AMinCoords[1] = HUGE_VALF;
		}
		if(geomDimension>2) {
			AMinCoords[2] = HUGE_VALF;
		}
		AMaxCoords[0] = -HUGE_VALF;
		if(geomDimension>1) {
			AMaxCoords[1] = -HUGE_VALF;
		}
		if(geomDimension>2) {
			AMaxCoords[2] = -HUGE_VALF;
		}

		if(this->empty()) {
			return;
		}

		int cellDim = Type2Dim<TCellType>::val;

		// check whether the cells know the nodes
		if(0 == cellDim) {

		} else if(1 == cellDim) {
			if(!mesh_->getModel().has(E2N)) {
				throw GMDSException("CellGroup::computeBoundingBox cannot be called "
									"when cells (GMDS_EDGE) do not know nodes.");
			}
		} else if(2 == cellDim) {
			if(!mesh_->getModel().has(F2N)) {
				throw GMDSException("CellGroup::computeBoundingBox cannot be called "
									"when cells (GMDS_FACE) do not know nodes.");
			}
		} else if(3 == cellDim) {
			if(!mesh_->getModel().has(R2N)) {
				throw GMDSException("CellGroup::computeBoundingBox cannot be called "
									"when cells (GMDS_REGION) do not know nodes.");
			}
		}

		if(0 == cellDim) {
			for(unsigned int iCell=0; iCell<this->size(); iCell++) {

				Node current_node = mesh_->get<Node>(cells_[iCell]);

				TCoord xCoord = current_node.X();
				TCoord yCoord = current_node.Y();
				TCoord zCoord = current_node.Z();

				if(AMinCoords[0] > xCoord) {
					AMinCoords[0] = xCoord;
				}
				if(AMaxCoords[0] < xCoord) {
					AMaxCoords[0] = xCoord;
				}
				if(geomDimension>1) {
					if(AMinCoords[1] > yCoord) {
						AMinCoords[1] = yCoord;
					}
					if(AMaxCoords[1] < yCoord) {
						AMaxCoords[1] = yCoord;
					}
				}
				if(geomDimension>2) {
					if(AMinCoords[2] > zCoord) {
						AMinCoords[2] = zCoord;
					}
					if(AMaxCoords[2] < zCoord) {
						AMaxCoords[2] = zCoord;
					}
				}
			}
		} else { //if(0 == cellDim)

		for(unsigned int iCell=0; iCell<this->size(); iCell++) {
			std::vector<Node> nodes;
			if(cellDim==1){
				nodes = mesh_->get<Edge>(cells_[iCell]).template get<Node>();
			}
			else if(cellDim==2)
				nodes = mesh_->get<Face>(cells_[iCell]).template get<Node>();
			else
				nodes = mesh_->get<Region>(cells_[iCell]).template get<Node>();

			for(unsigned int iNode=0; iNode<nodes.size(); iNode++) {

				Node current_node = nodes[iNode];

				TCoord xCoord = current_node.X();
				TCoord yCoord = current_node.Y();
				TCoord zCoord = current_node.Z();

				if(AMinCoords[0] > xCoord) {
					AMinCoords[0] = xCoord;
				}
				if(AMaxCoords[0] < xCoord) {
					AMaxCoords[0] = xCoord;
				}
				if(geomDimension>1) {
					if(AMinCoords[1] > yCoord) {
						AMinCoords[1] = yCoord;
					}
					if(AMaxCoords[1] < yCoord) {
						AMaxCoords[1] = yCoord;
					}
				}
				if(geomDimension>2) {
					if(AMinCoords[2] > zCoord) {
						AMinCoords[2] = zCoord;
					}
					if(AMaxCoords[2] < zCoord) {
						AMaxCoords[2] = zCoord;
					}
				}
			}
		}

		}

	} else {
		throw GMDSException("CellGroup::computeBoundingBox cannot be called "
				"when there are no nodes in the mesh.");
	}
}
/*----------------------------------------------------------------------------*/
template<typename TCellType>
std::vector<TCellType> CellGroup<TCellType>::cells()
{
	std::vector<TCellType> return_cells;
	return_cells.resize(cells_.size());
	for(unsigned int i=0;i<cells_.size();i++)
	{
		TCellID val = cells_[i];
		return_cells[i] = mesh_->get<TCellType>(val);
	}
	return return_cells;
}
/*----------------------------------------------------------------------------*/
template<typename TCellType>
TCellType CellGroup<TCellType>::operator[] (const int& i)
{
	return mesh_->get<TCellType>(cells_[i]);
}
/*----------------------------------------------------------------------------*/
template<typename TCellType>
void CellGroup<TCellType>::serialize(std::ostream& AStr)
{
	const int name_size= name_.size();
	const int nb_values = cells_.size();
	const int total_size = sizeof(int) 				/* 1) total size */
						+ sizeof(int)  				/* 2) name size*/
						+name_size*sizeof(char) 	/* 3) name */
						+sizeof(int)				/* 4) nb values*/
						+ nb_values*sizeof(TCellID);/* 5) values */


	AStr.write((char*)&total_size,sizeof(int)); 			/* fill (1) */
	AStr.write((char*)&name_size,sizeof(int));				/* fill (2) */
	AStr.write(name_.c_str(),sizeof(char)*name_.size());	/* fill (3) */
	AStr.write((char*)&nb_values,sizeof(int));				/* fill (4) */
	AStr.write((char*)&cells_[0],nb_values*sizeof(TCellID));	/* fill (5) */

}
/*----------------------------------------------------------------------------*/
template<typename TCellType>
void CellGroup<TCellType>::unserialize(std::istream& AStr)
{
	int total_size = 0;
	int name_size  = 0;
	int nb_values = 0;
	AStr.read((char*)&total_size,sizeof(int)); /* read (1) */
	AStr.read((char*)&name_size,sizeof(int));  /* read (2) */
	char *n = new char[name_size];
	AStr.read(n,sizeof(char)*name_size);   	 /* read (3) */
	name_.clear();
	name_.assign(n,n+name_size);
	delete[] n;


	AStr.read((char*)&nb_values,sizeof(int));	/* read (4) */
	cells_.clear();
	cells_.resize(nb_values);

	for(int i=0;i<nb_values;i++){
		TCellID elt;
		AStr.read((char*)&elt  ,sizeof(TCellID));
		cells_[i]=elt;

	}
}
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_CELLGROUP_H_ */
/*----------------------------------------------------------------------------*/
