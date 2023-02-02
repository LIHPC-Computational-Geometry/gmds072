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
/** \file    SmartVector.t.h
 *  \author  F. LEDOUX
 *  \date    04/06/2009
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_SMARTVECTOR_H_
#define GMDS_SMARTVECTOR_H_
/*----------------------------------------------------------------------------*/
#include <vector>
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/CommonTypes.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/SmartVector_def.h>
/*----------------------------------------------------------------------------*/
template<typename T> SmartVector<T>::
SmartVector(const int capacity):
top_(0),size_(0)
{
	vec_.resize(capacity);
	mark_.resize(capacity,false);
}
/*----------------------------------------------------------------------------*/
template<typename T>
SmartVector<T>::SmartVector(
	std::vector<bool> AMarks):
mark_(AMarks)
{
	vec_.resize(AMarks.size());

	update();
}
/*----------------------------------------------------------------------------*/
template<typename T>
SmartVector<T>::SmartVector(const SmartVector<T>& v)
:top_(v.top_),size_(v.size_),
 free_stack_(v.free_stack_), vec_(v.vec_), mark_(v.mark_)
{}
/*----------------------------------------------------------------------------*/
template<typename T> SmartVector<T>::
~SmartVector()
{}
/*----------------------------------------------------------------------------*/
template<typename T>
TInt SmartVector<T>::size() const
{return size_;}
/*----------------------------------------------------------------------------*/
template<typename T>
TInt SmartVector<T>::top() const
{return top_;}
/*----------------------------------------------------------------------------*/
template<typename T>
TInt SmartVector<T>::capacity() const
{return vec_.size();}
/*----------------------------------------------------------------------------*/
template<typename T>
void SmartVector<T>::clear()
{
	free_stack_.clear();

	vec_.clear();
	vec_.resize(GChunkSize);

	mark_.clear();
	mark_.resize(GChunkSize,false);

	top_=0;
	size_=0;
}
/*----------------------------------------------------------------------------*/
template<typename T>
void SmartVector<T>::resize(const TInt& ASize)
{
	vec_.resize(ASize);
	mark_.resize(ASize,false);

	update();
}
/*----------------------------------------------------------------------------*/
template<typename T>
void SmartVector<T>::mark(const std::vector<int>& ref)
{
	for(unsigned int i=0; i<ref.size();i++)
		mark_[ref[i]]=true;
	update();
}
/*----------------------------------------------------------------------------*/
template<typename T>
void SmartVector<T>::update()
{
	/* the free stack is made empty */
	free_stack_.clear();

	TInt free=0;
	/* we traverse the marks_ tabular to know which items are free and we link
	 * them together in the vec_ tabular.
	 */
	size_=0;

	TInt i =mark_.size()-1;
	while(i>=0 && mark_[i]==false){
		i--;
	}

	/* in this case, all the items are free */
	if(i==-1){
		top_=0;
		return;
	}
	/* we are on the greatest used item */
	top_=i+1;
	free = i;

	for(;i>=0;i--){
		if(mark_[i]==false)
			free_stack_.push_back(i);
		else
			size_++;
	}

}
/*----------------------------------------------------------------------------*/
template<typename T>
bool SmartVector<T>::empty() const
{
	return (size_==0);
}
/*----------------------------------------------------------------------------*/
template<typename T>
TInt SmartVector<T>::selectNewIndex()
{
	TInt free;
	TInt cap=capacity();

	if(free_stack_.empty()){
		free= top_;
		top_++;
		if(top_>=cap){
			vec_.resize((cap+1)*2);
			mark_.resize((cap+1)*2,false);
		}
	}
	else{
		free = free_stack_.back();
		free_stack_.pop_back();
	}
	return free;
}
/*----------------------------------------------------------------------------*/
template<typename T>
void SmartVector<T>::assign(T AElt, const TInt& AIndex)
{
	vec_[AIndex]=AElt;
	if(mark_[AIndex]==false)
	{
		mark_[AIndex]=true;
		size_++;
	}

}
/*----------------------------------------------------------------------------*/
template<typename T>
void SmartVector<T>::add(T AElt)
{
	assign(AElt,selectNewIndex());
}
/*----------------------------------------------------------------------------*/
template<typename T>
bool SmartVector<T>::find(T AElt) //const
{
	//const_handle_iterator it = this->begin();

	Iterator it = this->begin();

	while(!it.isDone()){
		if(it.currentItem()==AElt)
			return true;
		it.next();
	}
	return false;
}
/*----------------------------------------------------------------------------*/
template<typename T>
void SmartVector<T>::remove(const TInt& AIndex)
{
#ifdef __DEBUG__
	if(AIndex>=capacity() || AIndex<0)
		throw GMDSException("Bad index in SmartVector<T>::remove()");
#endif //__DEBUG__
	/* the item AIndex is already free */
	if(mark_[AIndex]==false)
		return;
	T t = T();
	vec_[AIndex]=t;
	mark_[AIndex]=false;
	size_--;

	if(AIndex==top_)
		top_--;
	else
		free_stack_.push_back(AIndex);

}
/*----------------------------------------------------------------------------*/
template<typename T>
void SmartVector<T>::removeElement(const T& AElt)
{
	/* the item AIndex is already free */
	Iterator it = this->begin();
	Iterator ite = this->end();
	while(it!=ite){
		if(*it==AElt){
			remove(it);
			return;
		}
	}
}
/*----------------------------------------------------------------------------*/
template<typename T>
void SmartVector<T>::remove(Iterator& AIt)
{
	//necessary downcast. Design issue ?
	Iterator* local_it = dynamic_cast<Iterator*>(AIt.operator->());
	remove(local_it->current_);
}
/*----------------------------------------------------------------------------*/
template<typename T>
bool SmartVector<T>::isAvailable(const TInt& AIndex) const
{

	if(AIndex>=capacity())
		return false;
	else if (mark_[AIndex]==true)
		return false;

	return true;
}
/*----------------------------------------------------------------------------*/
template<typename T>
bool SmartVector<T>::isOutOfContainer(const TInt& AIndex) const
{
	return (AIndex>=capacity());
}
/*----------------------------------------------------------------------------*/
template<typename T>
void SmartVector<T>::fillAll()
{
	TInt cap = capacity();
	for(int i=0;i<cap;i++)
		mark_[i]=true;
	update();

}
/*----------------------------------------------------------------------------*/
template<typename T>
void SmartVector<T>::display()
{
	for(unsigned int i=0;i<vec_.size();i++)
		std::cout<<vec_[i]<<" ";
	std::cout<<std::endl;
	for(unsigned int i=0;i<mark_.size();i++)
		std::cout<<mark_[i]<<" ";
	std::cout<<std::endl;
	std::cout<<"top: "<<top_
			 <<" - size: "<<size_
			 <<" - capacity:"<<capacity()<<std::endl;
}
/*----------------------------------------------------------------------------*/
template<typename T>
void SmartVector<T>::compact()
{
	TInt first_not_free=0;
	TInt last_used = top_;

	while(last_used>first_not_free){
		/* increase first free*/
		while(mark_[first_not_free]==true)
			first_not_free++;

		/* decrease last used*/
		while(mark_[last_used]==false)
			last_used--;

		if(last_used>first_not_free){
			vec_[first_not_free] = vec_[last_used];
			mark_[first_not_free]=true;
			mark_[last_used]=false;
		}
		//else, we have a compact collection
	}

	/* new dimension of the container */
	top_=first_not_free;
	std::vector<T> new_vec;
	new_vec.resize(top_);
	std::vector<bool> new_mark;
	new_mark.resize(top_,false);

	new_vec.assign(vec_.begin(),vec_.end());

	for(int i=0;i<top_;i++)
		new_mark[i] = mark_[i];

	vec_.swap(new_vec);
	mark_.swap(new_mark);

	/* the free stack is made empty */
	free_stack_.clear();
}
/*----------------------------------------------------------------------------*/
template<typename T>
void SmartVector<T>::compactWithMemory(std::vector<int>& AMove)
{
	TInt first_not_free=0;
	TInt last_used = top_;

	AMove.clear();
	/* +1 is used in case of size_=0*/
	AMove.resize(top_+1);

	while(last_used>first_not_free){
		/* increase first free*/
		while(mark_[first_not_free]==true){
			AMove[first_not_free]= first_not_free;
			first_not_free++;
		}
		/* decrease last used*/
		while(mark_[last_used]==false){
			AMove[last_used]=last_used;
			last_used--;
		}

		if(last_used>first_not_free){
			vec_[first_not_free] = vec_[last_used];
			mark_[first_not_free]=true;
			mark_[last_used]=false;
			AMove[last_used]= first_not_free;
		}
		//else, we have a compact collection

	}

	/* new dimension of the container */
	top_=first_not_free;
	std::vector<T> new_vec;
	new_vec.resize(top_);
	std::vector<bool> new_mark;
	new_mark.resize(top_,false);

	new_vec.assign(&vec_[0],&vec_[top_]);

	for(int i=0;i<top_;i++)
		new_mark[i] = mark_[i];

	vec_.swap(new_vec);
	mark_.swap(new_mark);

	/* the free stack is made empty */
	free_stack_.clear();
}
/*----------------------------------------------------------------------------*/
template<typename T>
void SmartVector<T>::serialize(std::ostream& stream)
{
	int container_size = top_;
	stream.write((char*)&container_size,sizeof(int));
	for(int i=0;i<top_;i++)
	{
		T elt = vec_[i];
		int m = mark_[i];
		stream.write((char*)&elt,sizeof(T));
		stream.write((char*)&m,sizeof(int));
	}
}
/*----------------------------------------------------------------------------*/
template<typename T>
void SmartVector<T>::unserialize(std::istream& stream)
{
	int nb_items=0;

	stream.read((char*)&nb_items,sizeof(int));
	clear();
	resize(nb_items);
	top_=nb_items;
	for(int i=0;i<nb_items;i++){
		T elt;
		int m;
		stream.read((char*)&elt  ,sizeof(T  ));
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
#endif /* GMDS_SMARTVECTOR_H_ */
/*----------------------------------------------------------------------------*/
