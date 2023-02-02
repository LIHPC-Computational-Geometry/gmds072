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
 * SmartVector.cpp
 *
 *  Created on: 19 juin 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/SmartVector.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
template<> SmartVector<TCellID>::
SmartVector(const int capacity):
top_(0),size_(0)
{
	vec_.resize(capacity,NullID);
	mark_.resize(capacity,false);
}
/*----------------------------------------------------------------------------*/
template<> void SmartVector<TCellID>::resize(const TInt& ASize)
{
	vec_.resize(ASize,NullID);
	mark_.resize(ASize,false);

	update();
}
/*----------------------------------------------------------------------------*/
template<> void SmartVector<TabCellID<size_undef> >::
serialize(std::ostream& stream)
{
	int container_size = top_;
	stream.write((char*)&container_size,sizeof(int));
	for(int i=0;i<top_;i++)
	{
		TabCellID<size_undef> elt = vec_[i];
		int m = mark_[i];
		elt.serialize(stream);
		stream.write((char*)&m,sizeof(int));
	}
}

/*----------------------------------------------------------------------------*/
template<> void SmartVector<TabCellID<size_undef> >::
unserialize(std::istream& stream)
{
	int nb_items=0;

	stream.read((char*)&nb_items,sizeof(int));
	clear();
	resize(nb_items);
	top_=nb_items;
	for(int i=0;i<nb_items;i++){
		TabCellID<size_undef> elt;
		int m;
		elt.unserialize(stream);
		stream.read((char*)&m  ,sizeof(int));
		vec_[i]=elt;
		mark_[i]=m;
		if(m)
			size_++;
	}
	update();
}
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
