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
 * SmartBitVector.cpp
 *
 *  Created on: 6 févr. 2014
 *      Author: ledouxf
 */

/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/SmartBitVector.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
SmartBitVector::SmartBitVector(const int capacity):
m_top(0),m_size(0)
{
	m_bits.resize(capacity,false);
	m_capacity = capacity;
}
/*----------------------------------------------------------------------------*/
SmartBitVector::SmartBitVector(const SmartBitVector& v)
: m_top(v.m_top),m_size(v.m_size), m_capacity(v.m_capacity),
  m_free_stack(v.m_free_stack), m_bits(v.m_bits)
{}
/*----------------------------------------------------------------------------*/
SmartBitVector::
~SmartBitVector()
{}
/*----------------------------------------------------------------------------*/
TInt SmartBitVector::size() const
{return m_size;}
/*----------------------------------------------------------------------------*/

TInt SmartBitVector::top() const
{return m_top;}
/*----------------------------------------------------------------------------*/
TInt SmartBitVector::capacity() const
{return m_bits.size();}
/*----------------------------------------------------------------------------*/
void SmartBitVector::clear()
{
	m_free_stack.clear();

	m_bits.clear();
	m_bits.resize(GChunkSize,false);

	m_top=0;
	m_size=0;
	m_capacity=GChunkSize;
}
/*----------------------------------------------------------------------------*/
void SmartBitVector::resize(const TInt& ASize)
{
	m_bits.resize(ASize,false);
	m_capacity=ASize;
	update();
}/*----------------------------------------------------------------------------*/
void SmartBitVector::update()
{
	/* the free stack is made empty */
	m_free_stack.clear();

	TInt free=0;
	/* we traverse the m_bits tabular to know which items are free and we link
	 * them together in the m_bits tabular.
	 */
	m_size=0;

	TInt i =m_bits.size()-1;
	while(i>=0 && m_bits[i]==0){
		i--;
	}

	/* in this case, all the items are free */
	if(i==-1){
		m_top=0;
		return;
	}
	/* we are on the greatest used item */
	m_top=i+1;
	free = i;

	for(;i>=0;i--){
		if(m_bits[i]==0)
			m_free_stack.push_back(i);
		else
			m_size++;
	}

}
/*----------------------------------------------------------------------------*/
bool  SmartBitVector::operator[]  (const TInt& AIndex) const
{
#ifdef __DEBUG__
	if(AIndex>=m_capacity || m_bits[AIndex]==0)
		throw GMDSException("Bad index in SmartBitVector::operator[]");
#endif //__DEBUG__

	bool val = m_bits[AIndex];
	return val;
}
/*----------------------------------------------------------------------------*/
bool SmartBitVector::empty() const
{
	return (m_size==0);
}
/*----------------------------------------------------------------------------*/
TInt SmartBitVector::selectNewBit()
{
	TInt index = getFreeIndex();
	assign(index);
	return index;
}
/*----------------------------------------------------------------------------*/
void SmartBitVector::unselect(const TInt& AIndex)
{
#ifdef __DEBUG__
	if(AIndex>=m_capacity || AIndex<0)
		throw GMDSException("Bad index in SmartBitVector::unselect()");
#endif //__DEBUG__
	/* the item AIndex is already free */
	if(m_bits[AIndex]==false)
		return;

	m_bits[AIndex]=false;
	m_size--;

	if(AIndex==m_top)
		m_top--;
	else
		m_free_stack.push_back(AIndex);

}
/*----------------------------------------------------------------------------*/
bool SmartBitVector::isAvailable(const TInt& AIndex) const
{

	if(AIndex>=m_capacity)
		return false;
	else
		return !(m_bits[AIndex]);
}
/*----------------------------------------------------------------------------*/
bool SmartBitVector::isOutOfContainer(const TInt& AIndex) const
{
	return (AIndex>=m_capacity);
}
/*----------------------------------------------------------------------------*/
void SmartBitVector::fillAll()
{
	TInt cap = m_capacity;
	for(int i=0;i<cap;i++)
		m_bits[i]=1;
	update();

}
/*----------------------------------------------------------------------------*/
void SmartBitVector::display()
{
	for(unsigned int i=0;i<m_bits.size();i++)
		std::cout<<m_bits[i]<<" ";
	std::cout<<std::endl;
	std::cout<<"top: "<<m_top
			 <<" - size: "<<m_size
			 <<" - capacity:"<<m_capacity<<std::endl;
}
/*----------------------------------------------------------------------------*/
void SmartBitVector::compact(std::vector<int>& AMove)
{
//	TInt first_not_free=0;
//	TInt last_used = m_top;
//
//	AMove.clear();
//	/* +1 is used in case of size_=0*/
//	AMove.reserve(m_top+1);
//	AMove.resize(m_top+1);
//
//	while(last_used>first_not_free){
//		/* increase first free*/
//		while(m_bits[first_not_free]==1){
//			AMove[first_not_free]= first_not_free;
//			first_not_free++;
//		}
//		/* decrease last used*/
//		while(m_bits[last_used]==0){
//			AMove[last_used]=last_used;
//			last_used--;
//		}
//
//		if(last_used>first_not_free){
//			vec_[first_not_free] = vec_[last_used];
//			mark_[first_not_free]=1;
//			mark_[last_used]=0;
//			AMove[last_used]= first_not_free;
//		}
//		//else, we have a compact collection
//
//	}
//
//	/* new dimension of the container */
//	m_top=first_not_free;
//	std::vector<T> new_vec;
//	new_vec.reserve(top_);
//	new_vec.resize(top_);
//	std::vector<bool> new_mark;
//	new_mark.resize(top_,false);
//
//	new_vec.assign(&vec_[0],&vec_[top_]);
//
//	for(int i=0;i<top_;i++)
//		new_mark[i] = mark_[i];
//
//	vec_.swap(new_vec);
//	mark_.swap(new_mark);
//
//	/* the free stack is made empty */
//	free_stack_.reserve(1);
//	free_stack_.resize(0);
}
/*----------------------------------------------------------------------------*/
void SmartBitVector::serialize(std::ostream& stream)
{
	int container_size = m_top;
	stream.write((char*)&container_size,sizeof(int));

	for(int i=0;i<m_top;i++)
	{
		int m = m_bits[i];
		stream.write((char*)&m,sizeof(int));
	}
}
/*----------------------------------------------------------------------------*/
void SmartBitVector::unserialize(std::istream& stream)
{
	int nb_items=0;

	stream.read((char*)&nb_items,sizeof(int));
	clear();
	resize(nb_items);
	m_top=nb_items;
	for(int i=0;i<nb_items;i++){
		int m;
		stream.read((char*)&m  ,sizeof(int));
		m_bits[i]=m;
		if(m)
			m_size++;
	}
	update();
}
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
