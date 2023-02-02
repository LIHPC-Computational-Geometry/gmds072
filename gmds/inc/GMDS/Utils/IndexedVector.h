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
/** \file    IndexedVector.h
 *  \author  F. LEDOUX
 *  \date    02/06/2014
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_INDEXEDVECTOR_H_
#define GMDS_INDEXEDVECTOR_H_
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/CommonTypes.h>
#include <vector>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/IndexedVector_def.h>
/*----------------------------------------------------------------------------*/
template<typename T> IndexedVector<T>::
IndexedVector(const int capacity)
{
	m_vec.resize(capacity);
}
/*----------------------------------------------------------------------------*/
template<typename T>
IndexedVector<T>::IndexedVector(const IndexedVector<T>& v)
: m_vec(v.m_vec)
{}
/*----------------------------------------------------------------------------*/
template<typename T> IndexedVector<T>::
~IndexedVector()
{
	m_vec.clear();
}
/*----------------------------------------------------------------------------*/
template<typename T>
TInt IndexedVector<T>::capacity() const
{return m_vec.size();}
/*----------------------------------------------------------------------------*/
template<typename T>
void IndexedVector<T>::clear()
{
	m_vec.clear();
	m_vec.resize(2);
}
/*----------------------------------------------------------------------------*/
template<typename T>
void IndexedVector<T>::resize(const TInt& ASize)
{
	m_vec.resize(ASize);
}
/*----------------------------------------------------------------------------*/
template<typename T>
T const& IndexedVector<T>::operator[]  (const TInt& AIndex) const
{
#ifdef __DEBUG__
	if(AIndex>=capacity())
		throw GMDSException("Bad index in IndexedVector<T>::operator[]");
#endif //__DEBUG__

	return m_vec[AIndex];
}
/*----------------------------------------------------------------------------*/
template<typename T>
T& IndexedVector<T>::operator[]  (const TInt& AIndex)
{
#ifdef __DEBUG__
	if(AIndex>=capacity())
		throw GMDSException("Bad index in IndexedVector<T>::operator[]");
#endif //__DEBUG__
	return m_vec[AIndex];
}
/*----------------------------------------------------------------------------*/
template<typename T>
void IndexedVector<T>::assign(T& AElt, const TInt& AIndex)
{
	if(AIndex>=m_vec.size()) {
		TInt prev_size = m_vec.size();
		m_vec.resize(prev_size*2);
	}
	m_vec[AIndex]=AElt;
}
/*----------------------------------------------------------------------------*/
template<typename T>
void IndexedVector<T>::serialize(std::ostream& stream)
{
	int container_size = m_vec.size();
	stream.write((char*)&container_size,sizeof(int));

	for(int i=0;i<container_size;i++)
	{
		T elt = m_vec[i];
		stream.write((char*)&elt,sizeof(T));
	}
}
/*----------------------------------------------------------------------------*/
template<typename T>
void IndexedVector<T>::unserialize(std::istream& stream)
{
	int nb_items=0;

	stream.read((char*)&nb_items,sizeof(int));
	clear();
	resize(nb_items);
	for(int i=0;i<nb_items;i++){
		T elt;
		stream.read((char*)&elt  ,sizeof(T));
		m_vec[i]=elt;
	}
}
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_INDEXEDVECTOR_H_ */
/*----------------------------------------------------------------------------*/
