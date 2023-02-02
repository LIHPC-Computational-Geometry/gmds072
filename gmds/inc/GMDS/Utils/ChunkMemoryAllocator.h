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
/** \file    ChunkMemoryAllocator.h
 *  \author  F. LEDOUX
 *  \date    January 6, 2014
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_CHUNKMEMORYALLOCATOR_H_
#define GMDS_CHUNKMEMORYALLOCATOR_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <cstddef>
#include <map>
#include <string>
#include <malloc.h>
#include <stdexcept>
/*----------------------------------------------------------------------------*/
// POSIX File Header
#include <pthread.h>
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/Exception.h>
#include <GMDS/Utils/CommonTypes.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
/** \class IAllocator
 *  \brief Interface defining the two functions needed to be defined for an
 * 		   allocator.
 */
/*----------------------------------------------------------------------------*/
class IAllocator
{
public:

    /*------------------------------------------------------------------------*/
    /** \brief  Default destructor.
*/
	virtual ~IAllocator() {;}


    /*------------------------------------------------------------------------*/
    /** \brief  Provide the address of a ASize space.
     *
     *  \return a pointer on the first free space
*/
	virtual void* allocate(const int& ASize) =0;

    /*------------------------------------------------------------------------*/
    /** \brief  Make the space taken by the element APtr free.
     *
     *  \param APtr       a pointer on a memory address managed by this
*/
	virtual void deallocate(void* APtr) =0;
};
class IDAllocator
{
public:

    /*------------------------------------------------------------------------*/
    /** \brief  Default destructor.
*/
	virtual ~IDAllocator() {;}

    /*------------------------------------------------------------------------*/
    /** \brief  Provide the address of a ASize space.
     *
     *  \return a pointer on the first free space
*/
	virtual TCellID* allocate(const int& ASize) =0;

    /*------------------------------------------------------------------------*/
    /** \brief  Make the space taken by the element APtr free.
     *
     *  \param APtr       a pointer on a memory address managed by this
*/
	virtual void deallocate(TCellID* APtr) =0;


};
/*----------------------------------------------------------------------------*/
/** \class FreeElementForLIDVectorAllocator
 */
/*----------------------------------------------------------------------------*/
class FreeElementForLIDVectorAllocator
{
public:
	FreeElementForLIDVectorAllocator(TCellID* ABuf, const TInt& ASize,
									 FreeElementForLIDVectorAllocator* ANext)
	: buf(ABuf), size(ASize), next(ANext) {;}


	TCellID* buf;
	TInt size;
	FreeElementForLIDVectorAllocator* next;
};
/*----------------------------------------------------------------------------*/
/** \class LIDVectorAllocator
 *  \brief Manages the memory allocation/desallocation for a collection of
 *  	   local ids.
 *
 * 		   This allocator is a chunk-based allocator which gathers vectors of
 * 		   local ids. A vector is stored with in the first pointed space, a
 * 		   pointer on the memory space after the last item of the vector.
 *
 * 		   The chunks are allocated when necessary. They can also be
 * 		   de-allocated if they are not used. Free elements are known and
 * 		   managed with a map of linked list. When a vector is erased, its
 * 		   space in memory is added at the list head of the adequate length.
 * 		   When a vector is added, it uses the memory space pointed by the
 * 		   head list (with the right size or more). When the vector size is
 * 		   increased, a new space is found and its containt is moved.
 *
 *  \param TNbElements	number of ids
 */
/*----------------------------------------------------------------------------*/
template<int TNbElements>
class LIDVectorAllocator : public IDAllocator
{
public:


    /*------------------------------------------------------------------------*/
    /** \brief  Default Constructor.
     */
	LIDVectorAllocator();

    /*------------------------------------------------------------------------*/
    /** \brief  Default destructor.
     */
	~LIDVectorAllocator();

    /*------------------------------------------------------------------------*/
    /** \brief  Provide the address of the first free address allowing to store
     * 			a series of ASize ids.
     *
     *  \param ASize the size of the series
     *
     *  \return a pointer on the first free space
     */
	TCellID* allocate(const TInt& ASize);

    /*------------------------------------------------------------------------*/
    /** \brief  Make the space taken by the vector APtr free.
     *
     *  \param AIC a pointer on an IC instance
     */
	void deallocate(TCellID* APtr);

	/*------------------------------------------------------------------------*/
    /** \brief Make the collection empty.
     */
	void clear() {removeAllChunks();}

	/*------------------------------------------------------------------------*/
    /** \brief  Display a description of the memory allocator for the current
     * 		    thread.
     */
	void trace_memory(const std::string& AMessage);

	/*------------------------------------------------------------------------*/
    /** \brief  Compact the allocator.
     */
	void compact(std::map<TCellID*,TCellID*>& AMove);

	void print();
private:

    /*------------------------------------------------------------------------*/
    /** \struct MemoryChunk
     * 	\brief  A Memory chunk containing TNbElements of size TSize and a
     * 			pointer to the next memory chunk.
     */
	struct MemoryChunk
  	{
		/** Next Memory Chunk */
	    MemoryChunk* next;
	    /** memory space*/
	    TCellID buf[TNbElements];
  	};

    /*------------------------------------------------------------------------*/
    /** \brief Add a memory chunk if the list of free element is empty.
     */
  	void addChunk();

  	/*------------------------------------------------------------------------*/
    /** \brief Remove all the memory chunks when they are empty (ref count=0).
     */
  	void removeAllChunks();

  	/*------------------------------------------------------------------------*/
    /** \brief  Remove all the empty chunks
     */
  	void removeEmptyChunks();

  	/*------------------------------------------------------------------------*/
    /** \brief  Check if AChunk is empty
     *
     *  \param  AChunk the chunk to check
     *  \return true if the chunk is empty
     */
  	bool isAnEmptyChunk(MemoryChunk* AChunk) const;

  	/*------------------------------------------------------------------------*/
    /** \brief  Initialize the free elements chained list
     */
  	void initFreeElts();

  	/*------------------------------------------------------------------------*/
    /** \brief  clear the free elements list
     */
  	void clearFreeElts();

  	/*------------------------------------------------------------------------*/
    /** \brief  add a free elt in the list of free elements.
     *
     *  \param ALoc  the location of the new free element
     *  \param ASize the size of the new free element
     */
  	void newFreeElt(TCellID* ALoc, const TInt& ASize);

  	/*------------------------------------------------------------------------*/
    /** \brief Provide the address of a sufficient space for storing a ASize
     * 		   vector. If there is not enough space, 0 is returned.
     *
     *  \param ASize the size of the vector
     */
  	TCellID* getSpace(const TInt& ASize);

  	/** List of memory chunks in the global heap*/
  	MemoryChunk*  first_chunk_;
  	/** Number of stored elements in the global heap*/
  	TInt        nb_stored_elements_;

  	/** Free elements in the global heap. They are classified according to the
  	 *  available free space they provide. */
  	std::map<TInt,FreeElementForLIDVectorAllocator*> free_elements_;
};
/*----------------------------------------------------------------------------*/
/** \class ChunkCollection
 *  \brief Provides a collection of elements of type T and manages the
 * 		   memory allocation/desallocation of these elements by gathering them
 * 		   into chunks of TNbElements elements.
 *
 * 		   It gives iterators to traverse used elements. A free space contains a
 * 		   pointer to the next free space + an extra value indicating the space
 * 		   is free.
 *
 * 		   The chunks are allocated when necessary. They can also be
 * 		   de-allocated if they are not used. Free elements are known and
 * 		   managed with a linked list. When an element is erased, its space in
 * 		   memory is added at the list head. When an element is added, it uses
 * 		   the memory space pointed bt the head list.
 *
 * 		   allocate and deallocate operators are necessary to use this class as
 * 		   an allocator for overloading operator new and operator delete of
 * 		   classes we want to store inside chunk collections.
 *
 *		   WARNING: size(T) must be greater than 2 bytes in 32 bits and 4 bytes
 * 		   in 64 bytes to store the free list information. Moreover, the first
 * 		   attribute of type T must be a pointer. Invalid pointer address are
 *         used to detect free items in the collection.
 *
 *  \param T			type of element
 *  \param TNbElements	number of elements
 */
/*----------------------------------------------------------------------------*/
template<typename T, int TNbElements>
class ChunkCollection: public IAllocator
{
public:


	/*------------------------------------------------------------------------*/
    /** \brief  Default Constructor.
     */
	ChunkCollection();

	/*------------------------------------------------------------------------*/
    /** \brief  Default Destructor.
     */
	 virtual ~ChunkCollection();

    /*------------------------------------------------------------------------*/
    /** \brief  Provide the address of the first free element. Warnig ASize is
     * 			useless here since we store elements of size sizeof(T). It is
     * 			just here to overload traditional CC+ operator new.
     *
     *  \return a pointer on the first free space
     */
	void* allocate(const int& ASize);

    /*------------------------------------------------------------------------*/
    /** \brief Make the space taken by the element APtr free.
     *
     *  \param APtr a pointer on a memory address managed by this.
     */
	void deallocate(void* APtr);

    /*------------------------------------------------------------------------*/
    /** \brief Make the collection empty.
     */
	void clear() {removeAllChunks();}

	void print();

	int nbElements() {return nb_stored_elements_;}

	/*------------------------------------------------------------------------*/
    /** \brief  Serializes the cell data into stream AStr
     *
     *  \param AStr an output stream where the cell data is written
     */
	void serialize(std::ostream& AStr);

	/*------------------------------------------------------------------------*/
    /** \brief  Unserializes the cell data from stream AStr
     *
     *  \param AStr an input stream where the cell data is read from
     */
	void unserialize(std::istream& AStr);

	/*------------------------------------------------------------------------*/
    /** \brief  Makes the collection compact
     */
	void compact();

	/*------------------------------------------------------------------------*/
    /** \brief  Provides the number of chunks used to store the current
     * 			collection of T elements
     */
	int getNbChunks() const;
private:

    /*------------------------------------------------------------------------*/
    /** \struct MemoryChunk
     * 	\brief  A Memory chunk containing TNbElements of size sizeof(T) and a
     * 			pointer to the next memory chunk.
*/
	struct MemoryChunk
  	{
		/** Next Memory Chunk */
	    MemoryChunk* next;

	    /** Previous Memory Chunk */
	    MemoryChunk* prev;

	    /** memory space*/
	    char buf[sizeof(T)*TNbElements];
  	};

public:

	/*------------------------------------------------------------------------*/
    /** \class iterator.
     *
     * 	\brief nested class providing an iterator which traverses an allocator.
*/
	class iterator
	{
	public:

		friend class ChunkCollection<T,TNbElements>;

		/*--------------------------------------------------------------------*/
	    /** \brief Default constructor.
	     *
	     *  \param ACollection the collection traversed by this.
	     */
		/*--------------------------------------------------------------------*/
		iterator(ChunkCollection<T,TNbElements>* AAlloc=0)
		: allocator_(AAlloc),chunk_((AAlloc)?AAlloc->first_chunk_:0), index_in_chunk_(0){;}

		/*--------------------------------------------------------------------*/
	    /** \brief Copy constructor.
	     */
		/*--------------------------------------------------------------------*/
		iterator(const iterator& AIt)
		: allocator_(AIt.allocator_),chunk_(AIt.chunk_),
		  index_in_chunk_(AIt.index_in_chunk_){;}

		/*--------------------------------------------------------------------*/
	    /** \brief Equality for iterators.
	     *
	     * 	\param AIt an other iterator
	     */
		/*--------------------------------------------------------------------*/
		bool operator== (const iterator& AIt) const {
			//iterators are equal if they point to the same element
			return ( allocator_ == AIt.allocator_ &&  chunk_ == AIt.chunk_ &&
					 index_in_chunk_ == AIt.index_in_chunk_);
		}

		/*--------------------------------------------------------------------*/
	    /** \brief Inequality for iterators.
	     *
	     * 	\param AIt an other iterator
	     */
		/*--------------------------------------------------------------------*/
		bool operator!= (const iterator& AIt) const {
			//iterators are inequal if they point to different elements
			return ( allocator_ != AIt.allocator_ || chunk_ != AIt.chunk_ ||
					 index_in_chunk_ != AIt.index_in_chunk_);
		}

		/*--------------------------------------------------------------------*/
	    /** \brief Pointer dereference operator.
	     */
		/*--------------------------------------------------------------------*/
		T* operator* () {
#ifdef _DEBUG
			/* WARNING: we just verify we point into the tab addresses. We accept
			 * the location after the last element for end() iterator */
			if( index_in_chunk_<0|| index_in_chunk_>TNbElements )
				throw GMDSException();
#endif //DEBUG
			return reinterpret_cast<T*>(&(chunk_->buf[index_in_chunk_*sizeof(T)]));
		}

		/*--------------------------------------------------------------------*/
	    /** \brief Prefix increment. Move forward to get to the next "real" element.
	     *         Some items in the allocator are free. We don't have to provide
	     * 		   them with this iterator.
	     */
		/*--------------------------------------------------------------------*/
		iterator& operator++() {
			// move to the first real successor of index_

			//traversal of the chunks
			while(chunk_){
				//traversal of current chunk items
				char* c;
				do{
					index_in_chunk_++;
					c = *reinterpret_cast<char**>(chunk_->buf+index_in_chunk_*sizeof(T));
				}
				while ((c == (char*)((unsigned long)c|1)) && index_in_chunk_!=TNbElements);

				// we reached the end of the chunk
				if(index_in_chunk_ == TNbElements)
				{

					if(chunk_->next!=0)
					{
						chunk_ = chunk_->next;
						index_in_chunk_=-1;
					}
					else // end()
						return *this;
				}
				else
					return *this;
			};

			return *this;
		}

		/*--------------------------------------------------------------------*/
	    /** \brief Postfix increment. Move forward to get to the next "real"
	     * 		   element. Some items in the allocator are free. We don't have
	     * 		   to provide them with this iterator.
	     */
		/*--------------------------------------------------------------------*/
		iterator& operator++(int) {
				iterator& it_temp= *this;
				this->operator ++();
				return it_temp;
		}
		/*--------------------------------------------------------------------*/
	    /** \brief Gives the beginning of the container.
	     */
		/*--------------------------------------------------------------------*/
		iterator& begin() {
			// move to the first real successor of index_
			if(chunk_==0)
			{
				index_in_chunk_ = TNbElements;
				return  *this;
			}

			char* c = *reinterpret_cast<char**>(chunk_->buf);
			char* c1 = (char*)((unsigned long)c|1);
			// it is a real free element
			if (c == c1)
				return operator ++();
			else
				return *this;
		}
		/*--------------------------------------------------------------------*/
	    /** \brief Gives the end of the container.
	     */
		/*--------------------------------------------------------------------*/
		iterator& end() {
			chunk_ = allocator_->last_chunk_;
			index_in_chunk_ = TNbElements;
			return *this;
		}

	private:

		ChunkCollection<T,TNbElements>* allocator_;

		MemoryChunk* chunk_;

		int index_in_chunk_;
	};


	/*------------------------------------------------------------------------*/
	// iterators
	/*------------------------------------------------------------------------*/
    /** \brief Returns a read/write iterator that points to the first element
     * 		   in the collection.  Iteration is done in ordinary element order.
*/
	iterator begin()
	{
		return iterator(this).begin();
	}
	/*------------------------------------------------------------------------*/
    /** \brief  Returns a read/write iterator that points on past the last
     * 			element in the collection. Iteration is done in ordinary
     * 			element order.
*/
   iterator end() {return iterator(this).end();}


private:

    /*------------------------------------------------------------------------*/
    /** \brief  Add a memory chunk if the list of free element is empty.
     */
  	void addChunk();

  	/*------------------------------------------------------------------------*/
    /** \brief  Remove all the memory chunks when they are empty (ref count=0)
     */
  	void removeAllChunks();

  	/*------------------------------------------------------------------------*/
    /** \brief  Remove all the empty chunks
     */
  	void removeEmptyChunks();

  	/*------------------------------------------------------------------------*/
    /** \brief  Check if AChunk is empty
     *
     *  \param  AChunk the chunk to check
     *  \return true if the chunk is empty
     */
  	bool isAnEmptyChunk(MemoryChunk* AChunk) const;

  	/*------------------------------------------------------------------------*/
    /** \brief  Initialize the free elements chained list
     *
     */
  	void initFreeElts();

private:
  	/** pointers to the free elements in the global heap. The first free
  	 *  element is pointed.*/
  	char* free_elements_;
  	/** First chunk of the list of memory chunks in the global heap*/
  	MemoryChunk*  first_chunk_;

  	/** Last chunk of the list of memory chunks in the global heap*/
  	MemoryChunk* last_chunk_;

  	/** Number of stored elements in the global heap*/
  	int        nb_stored_elements_;
};
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
template<int TNbElements>
LIDVectorAllocator<TNbElements>::LIDVectorAllocator()
: first_chunk_(0), nb_stored_elements_(0)
{}
/*----------------------------------------------------------------------------*/
template<int TNbElements>
LIDVectorAllocator<TNbElements>::~LIDVectorAllocator()
{
	removeAllChunks();
}
/*----------------------------------------------------------------------------*/
template<int TNbElements>
TCellID* LIDVectorAllocator<TNbElements>::allocate(const TInt& ASize)
{

	/* if there is not enough free space, a new chunk is added */
	TInt vector_size =ASize+1;//(ASize%2)?ASize+1:ASize+2;
	TCellID* result = getSpace(vector_size);
	if(!result)
	{
		addChunk();
		result = getSpace(vector_size);
	}

	/* we keep in mind that a new atomic memory block is used now */
	nb_stored_elements_ +=vector_size;

#ifdef _DEBUG
	//trace_memory("ALLOCATE");
#endif //DEBUG

	result[0]= ASize;

	return result;
}
/*----------------------------------------------------------------------------*/
template<int TNbElements>
void LIDVectorAllocator<TNbElements>::
deallocate(TCellID* APtr)
{
	/* APtr is supposed to be a pointer on the first items storing a series
	 * of local ids. Thus we can compute the size of the series from this item*/
	TInt size = APtr[0];
	TInt vector_size = size+1;//(size%2)?size+1:size+2;

	for(unsigned int i=0;i<vector_size;i++)
		APtr[i]=NullID;

	FreeElementForLIDVectorAllocator* new_list;

	/* the element which must be made free is simply added in the list of free
	 * elements of size size*/
	newFreeElt(APtr,vector_size);

	/* if all the elements are free, the memory chunks are free too */
	nb_stored_elements_ -= vector_size;
    if(!(nb_stored_elements_))
    	removeAllChunks();
#ifdef DEBUG
//	trace_memory("DEALLOCATE");
#endif //DEBUG
}
/*----------------------------------------------------------------------------*/
template<int TNbElements>
TCellID* LIDVectorAllocator<TNbElements>::
getSpace(const TInt& ASize)
{
	typename std::map<TInt,FreeElementForLIDVectorAllocator* >::iterator it;
	for(it=free_elements_.begin();it!=free_elements_.end();it++)
	{

		TInt size = it->first;

		if(size==ASize)
		{
			/* we catch the first free element whose size is ASize and we remove
			 * it from the list of free elements*/
			FreeElementForLIDVectorAllocator* elt = it->second;
			it->second = elt->next;

			/* if there is no more left space of this size, we remove the list
			 * from the free space lists*/
			if(it->second==0)
				free_elements_.erase(size);

			TCellID* byte = elt->buf;
			delete elt;
			return byte;
		}
		else if (size>ASize)
		{
			/* we catch the first free element and we remove
			 * it from the list of free elements*/
			FreeElementForLIDVectorAllocator* elt = it->second;
			it->second = elt->next;

			/* we have more space than needed. Thus we cut it it two spaces, one
			 * that we use and a smaller free space */
			int remainder = it->first-ASize;
			newFreeElt(elt->buf+ASize,remainder);

			/* if there is no more left space of this size, we remove the list
			 * from the free space lists*/
			if(it->second==0)
				free_elements_.erase(size);

			TCellID* byte = elt->buf;
			delete elt;
			return byte;
		}
	}

	/* we found nothing*/
	return 0;
}
/*----------------------------------------------------------------------------*/
template<int TNbElements>
void LIDVectorAllocator<TNbElements>::
addChunk()
{
  /* Allocation of a new memory chunk */
	MemoryChunk* new_chunk = new MemoryChunk();
		//(MemoryChunk*) malloc(sizeof(MemoryChunk));
	if(!new_chunk)
		throw std::bad_alloc();

	for(unsigned int i=0;i<TNbElements;i++)
		new_chunk->buf[i]=NullID;

  	/* the new chunk is added to the chunk list */
	new_chunk->next = first_chunk_;
	first_chunk_ = new_chunk;

	/* this chunk is seen as a huge free element*/
	newFreeElt(new_chunk->buf,TNbElements);
}
/*----------------------------------------------------------------------------*/
template<int TNbElements>
void LIDVectorAllocator<TNbElements>::
removeAllChunks()
{
	clearFreeElts();
	if(!first_chunk_)
		return;
	/* Traversal of the chunk list for making every chunk free */
	do{
		MemoryChunk* current_chunk = first_chunk_->next;
		delete(first_chunk_);
		//free(first_chunk_);
		first_chunk_ = current_chunk;
	} while(first_chunk_);
	nb_stored_elements_=0;
}
/*----------------------------------------------------------------------------*/
template<int TNbElements>
void LIDVectorAllocator<TNbElements>::
clearFreeElts()
{
	typename std::map<TInt,FreeElementForLIDVectorAllocator* >::iterator it;
	for(it=free_elements_.begin();it!=free_elements_.end();it++)
	{
		FreeElementForLIDVectorAllocator* elt = it->second;
		if(elt)
			delete elt;
	}

	free_elements_.clear();
}
/*----------------------------------------------------------------------------*/
template<int TNbElements>
void LIDVectorAllocator<TNbElements>::
print()
{
	/* Traversal of the chunk list for making every chunk free */
	std::cout<<"**********************************************"<<std::endl;
	std::cout<<std::endl;
	std::cout<<"**********************************************"<<std::endl;
	MemoryChunk* current_chunk = first_chunk_;
	if(current_chunk==0){
		std::cout<<"- empty -"<<std::endl;
		return;
	}
	int chunk_index=0;
	do{
		TCellID* buffer = current_chunk->buf;
		std::cout<<"Chunk "<<++chunk_index<<": ";
		std::cout<<"------------------------------------------------------"<<std::endl;
		for(int i =0; i<TNbElements;i++){
			std::cout<<" ";
			if(buffer[i]==NullID)
				std::cout<<"*";
			else
				std::cout<<buffer[i];
			if(i!=TNbElements-1)
				std::cout<<" - ";
		}
		std::cout<<std::endl;
		current_chunk = current_chunk->next;
	} while(current_chunk);
}

/*----------------------------------------------------------------------------*/
template<int TNbElements>
void LIDVectorAllocator<TNbElements>::
trace_memory(const std::string& AMessage)
{

	MemoryChunk* temp_chunk, *init_chunk = first_chunk_;
	TInt nb_stored_elements = nb_stored_elements_;

	temp_chunk = init_chunk;

	int nb_chunks=0;
	std::cout<<"============== "<<AMessage<<" =============="<<std::endl;
	std::cout<<"LIDVectorAllocator<T,"<<TNbElements<<">"<<std::endl;
	/* Traversal of the chunk list*/
	if(temp_chunk==0)
	{
		std::cout<<pthread_self()<<"NB Chunks= 0 - NB Stored Elts= 0 - NB Free Elts= 0"<<std::endl;
		return;
	}

	do{
		nb_chunks++;
		temp_chunk= temp_chunk->next;
	} while(temp_chunk);

	/* Number of free elements */
	typename std::map<TInt,FreeElementForLIDVectorAllocator* >::iterator it;

	std::map<TInt,TInt> nb_free_elts;

	for(it=free_elements_.begin();it!=free_elements_.end();it++)
	{
		FreeElementForLIDVectorAllocator* elt = it->second;

		if (elt!=0)
		{
			FreeElementForLIDVectorAllocator* temp_free_elts= elt;

			int nb_free=0;

			do{
				nb_free++;
				temp_free_elts= temp_free_elts->next;
			} while(temp_free_elts);

			nb_free_elts.insert(std::pair<TInt,TInt>(it->first,nb_free));
		}
	}
	std::cout<<"Thread: "<<pthread_self()<<" - NB Chunks= "<<nb_chunks;
	std::cout<<" - NB Stored Elts= "<<nb_stored_elements<<std::endl;

	std::map<TInt,TInt>::iterator it2;

	for(it2=nb_free_elts.begin();it2!=nb_free_elts.end();it2++)
		std::cout<<" - NB Free Elts of size "<<it2->first<<"= "<<it2->second<<std::endl;

}
/*----------------------------------------------------------------------------*/
template<int TNbElements>
void LIDVectorAllocator<TNbElements>::
compact(std::map<TCellID*,TCellID*>& AMove)
{
	AMove.clear();
	if(first_chunk_==0)
		return;
	/* We begin by traversing all the free elements to marks the corresponding
	 * space with NullID */
	typename std::map<TInt,FreeElementForLIDVectorAllocator* >::iterator it;
	for(it=free_elements_.begin();it!=free_elements_.end();it++)
	{
		FreeElementForLIDVectorAllocator* elt = it->second;

		if (elt!=0)	{
			FreeElementForLIDVectorAllocator* temp_free_elt= elt;

			do{
				TInt s = temp_free_elt->size;
				TCellID* b = temp_free_elt->buf;
				for(int i=0;i<s;i++)
					b[i] = NullID;
				temp_free_elt= temp_free_elt->next;
			} while(temp_free_elt);
		}
	}

	/* now all the free items in the memory chunks are equals to NullID. The second
	 * step consists in finding the first free item
	 * */
	MemoryChunk *current_chunk = first_chunk_, *free_chunk = first_chunk_;
	int index_in_chunk=-1, index_in_free_chunk=-1;

	bool not_found_first_free=true;
	while(not_found_first_free){
		//traversal of current chunk items
		TCellID d;
		do{
			index_in_chunk++;
			if (current_chunk->buf[index_in_chunk]==NullID){
				// it is a free item
				free_chunk = current_chunk;
				index_in_free_chunk = index_in_chunk;
				not_found_first_free=false;
			}
		}
		while (index_in_chunk!=TNbElements && not_found_first_free);

		// we reached the end of the chunk.
		if(index_in_chunk== TNbElements && current_chunk->next!=0){
			current_chunk= current_chunk->next;
			index_in_chunk=-1;
			index_in_free_chunk=-1;
		}
	};

	//traversal of the chunks to find used items and compact them
	while(current_chunk){
		//traversal of current chunk items

		while (index_in_chunk!=TNbElements-1){
			index_in_chunk++;

			if (current_chunk->buf[index_in_chunk]!=NullID){

				/* now we have to move an elt in order to compact the collection
				 * If both the free and current items are in the same chunk there
				 * is no problem. Otherwise, we have to ensure that the free chunk
				 * provides enough space to store the current elt */

				if(current_chunk!=free_chunk &&
				   current_chunk->buf[index_in_chunk]+1>TNbElements-index_in_free_chunk)
				{
					/* we do not have enough free space here, thus we have to
					 * available free space in the next chunk
					 */
					free_chunk = free_chunk->next;
					if(free_chunk==0)
						throw GMDSException("Unexpected behavior during LIDVectorAllocator compacting");
					index_in_free_chunk=-1;

					not_found_first_free=true;
					while(not_found_first_free){
						//traversal of current chunk items
						TCellID d;
						do{
							index_in_free_chunk++;
							if (free_chunk->buf[index_in_free_chunk]==NullID){
								not_found_first_free=false;
							}
						}
						while (index_in_free_chunk<TNbElements && not_found_first_free);

						// we reached the end of the chunk.
						if(index_in_free_chunk>= TNbElements){
							if(free_chunk->next!=0){
								free_chunk= free_chunk->next;
								index_in_free_chunk=-1;
							}
							else {

								removeEmptyChunks();
								return;
							}
							/* it is then compacted. Nothing else to do */
						}
					}
					// we relocate the currennt location
					current_chunk = free_chunk;
					index_in_chunk = index_in_free_chunk;
				}
				else if(current_chunk->buf[index_in_chunk]+1<=TNbElements-index_in_free_chunk)
				{
					TCellID* free_elt = &(free_chunk->buf[index_in_free_chunk]);

					int size_of_elt = current_chunk->buf[index_in_chunk]+1;
					// we move the elt to the free space
					memmove(&free_chunk->buf[index_in_free_chunk],
							&current_chunk->buf[index_in_chunk],
							size_of_elt*sizeof(TCellID));

					AMove[&current_chunk->buf[index_in_chunk]] =
										&free_chunk->buf[index_in_free_chunk];

					/* the space where the elt was stored is not free */
					for(int i=0;i<size_of_elt;i++)
						current_chunk->buf[index_in_chunk+i]=NullID;
					// we move the free space
					if(index_in_free_chunk+size_of_elt<TNbElements-1)
						index_in_free_chunk+=size_of_elt;
					else // go to the next chunk
					{
						free_chunk = free_chunk->next;
						index_in_free_chunk=0;
					}
				}
			}
		}


		// we reached the end of the current chunk
		current_chunk = current_chunk->next;
		index_in_chunk=-1;
	};
	removeEmptyChunks();
}
/*----------------------------------------------------------------------------*/
template<int TNbElements>
void LIDVectorAllocator<TNbElements>::
removeEmptyChunks()
{
	/* Traversal of the chunk list for removing empty chunk */
	MemoryChunk* c = first_chunk_;
	MemoryChunk* prev_c = 0;
	while(c!=0)
	{
		if (isAnEmptyChunk(c)){
			MemoryChunk *next_c = c->next;
			delete c;//free(c);
			if(prev_c)
				prev_c->next = next_c;

			c =next_c;
		}
		else{
			prev_c = c;
			c=c->next;
		}
	}
	initFreeElts();
}
/*----------------------------------------------------------------------------*/
template<int TNbElements>
bool LIDVectorAllocator<TNbElements>::
isAnEmptyChunk(MemoryChunk* AChunk) const
{
	/* Traversal of the chunk list for removing empty chunk */
	for(int index = 0;index<TNbElements;index++){
		if (AChunk->buf[index]!=NullID)
			return false;
	}
	return true;
}
/*----------------------------------------------------------------------------*/
template<int TNbElements>
void LIDVectorAllocator<TNbElements>::
initFreeElts()
{
	if(first_chunk_==0)
		return;

	clearFreeElts();
	/* we traverse all the chunks and we refill the free elements structure */
	MemoryChunk* current_chunk = first_chunk_;
	int first_free=0;
	int nb_free=0;
	bool is_first_free = true;
	while(current_chunk){
		TCellID* buffer = current_chunk->buf;
		for(int i =0; i<TNbElements;i++){
			if(buffer[i]==NullID){
				if (is_first_free){
					is_first_free=false;
					first_free = i;
					nb_free = 1;
				}
				else
					nb_free++;
			}
			else if(nb_free!=0){
				newFreeElt(&buffer[first_free],nb_free);
				nb_free=0;
				first_free=0;
				is_first_free = true;
			}
		}
		/* we check if the end of the current chunk is a series of free items*/
		if(nb_free!=0){
			newFreeElt(&buffer[first_free],nb_free);
			nb_free=0;
			first_free=0;
			is_first_free = true;
		}
		current_chunk = current_chunk->next;
	};

}
/*----------------------------------------------------------------------------*/
template<int TNbElements>
void LIDVectorAllocator<TNbElements>::
newFreeElt(TCellID* ALoc, const TInt& ASize){

	FreeElementForLIDVectorAllocator* new_list;

	/* the element which must be made free is simply added in the
	 * list of free elements of size size*/
	if(free_elements_.find(ASize)==free_elements_.end())
	{
		/* We create a new list of free elements whose size is
		 * equal to size */
		new_list = new FreeElementForLIDVectorAllocator(ALoc,ASize,0);
	}
	else
	{
		/* we add a new free element to the list of free elements whose
		 * size is equal to remainder*/
		 new_list = new FreeElementForLIDVectorAllocator(
				 ALoc,ASize,free_elements_[ASize]);
	}

	free_elements_.insert(
		std::pair<TInt,FreeElementForLIDVectorAllocator*>(ASize,new_list));

}
/*----------------------------------------------------------------------------*/
template<typename T, int TNbElements>
ChunkCollection<T,TNbElements>::
ChunkCollection()
: free_elements_(0),first_chunk_(0),last_chunk_(0),nb_stored_elements_(0)
{}
/*----------------------------------------------------------------------------*/
template<typename T, int TNbElements>
ChunkCollection<T,TNbElements>::
~ChunkCollection()
{
	/* we must free the different block */
	//if(first_chunk_)
		removeAllChunks();
}
/*----------------------------------------------------------------------------*/
template<typename T, int TNbElements>
void* ChunkCollection<T,TNbElements>::
allocate(const int& ASize)
{
#ifdef DEBUG_VERBOSE
	std::cout<<"Chunk allocate"<<std::endl;
#endif //DEBUG_VERBOSE

	/* if there is no free elements, a new chunk is added */

	if(!free_elements_)
		addChunk();

	/* the first element of the free list is returned */
	void* result = free_elements_;
	/* we remove this element from the free list */
	free_elements_ = *reinterpret_cast<char**>(free_elements_+sizeof(char*));
	/* we keep in mind that a new atomic memory block is used now */
	++(nb_stored_elements_);

	return result;
}
/*----------------------------------------------------------------------------*/
template<typename T, int TNbElements>
void ChunkCollection<T,TNbElements>::
deallocate(void* APtr)
{
	/* the element which must be made free is simply added in the list of free
	 * elements */

	/* we change the rightest bit to 1 to express this element is free now */
	*reinterpret_cast<char**>(APtr) = (char*)1;

	*reinterpret_cast<char**>((char*)APtr + sizeof(char*)) = free_elements_;

	free_elements_ = reinterpret_cast<char*>(APtr);



	/* if all the elements are free, the memory chunks are free too */
    if(!(--nb_stored_elements_))
    	removeAllChunks();
}
/*----------------------------------------------------------------------------*/
template<typename T, int TNbElements>
void ChunkCollection<T,TNbElements>::
addChunk()
{
#ifdef DEBUG_VERBOSE
	std::cout<<"Add chunks"<<std::endl;
#endif //DEBUG_VERBOSE

	/* Allocation of a new memory chunk */
	MemoryChunk* new_chunk = (MemoryChunk*) malloc(sizeof(MemoryChunk));
	if(!new_chunk)
		throw std::bad_alloc();

  	/* the new chunk is added at the end of the chunk list to preserve order*/
	if(first_chunk_)
	{
		last_chunk_->next = new_chunk;
		new_chunk->prev = last_chunk_;
		new_chunk->next=0;
		last_chunk_ = new_chunk;

	}
	else
	{
		new_chunk->next = 0;
		new_chunk->prev=0;
		first_chunk_ = new_chunk;
		last_chunk_ = new_chunk;
	}

	/* All the elements of this new chunk are added to free elements */
	free_elements_ = new_chunk->buf;

	for(size_t i=0; i<TNbElements; i++){
		char* current_elt = new_chunk->buf+i*sizeof(T);
		char* next_elt = new_chunk->buf+(i+1)*sizeof(T);

		*reinterpret_cast<char**>(current_elt) = (char*)1;

		if(i<TNbElements-1)
			*reinterpret_cast<char**>(current_elt + sizeof(char*)) = next_elt;
		else
			*reinterpret_cast<char**>(current_elt + sizeof(char*)) = 0;
	}
}
/*----------------------------------------------------------------------------*/
template<typename T, int TNbElements>
void ChunkCollection<T,TNbElements>::
removeAllChunks()
{
	/* Traversal of the chunk list for making every chunk free */
	if(first_chunk_!=0)
		do{
			MemoryChunk* current_chunk = first_chunk_->next;
			free(first_chunk_);
			first_chunk_ = current_chunk;
		} while(first_chunk_);

	/* There is no more chunk, thus no more free elements too */
	free_elements_ = 0;
	nb_stored_elements_=0;
}
/*----------------------------------------------------------------------------*/
template<typename T, int TNbElements>
void ChunkCollection<T,TNbElements>::
removeEmptyChunks()
{
	/* Traversal of the chunk list for removing empty chunk */
	MemoryChunk* c = first_chunk_;
	while(c!=0)
	{
		if (isAnEmptyChunk(c)){
			MemoryChunk *prev_c, *next_c;
			prev_c = c->prev;
			next_c = c->next;
			if(prev_c)
				prev_c->next = next_c;
			if(next_c)
				next_c->prev = prev_c;
			free(c);
			c = next_c;
		}
		else
			c=c->next;
	}
	initFreeElts();
}
/*----------------------------------------------------------------------------*/
template<typename T, int TNbElements>
bool ChunkCollection<T,TNbElements>::
isAnEmptyChunk(MemoryChunk* AChunk) const
{
	/* Traversal of the chunk list for removing empty chunk */
	for(int index = 0;index<TNbElements;index++){
		char* c=0;
		c = *reinterpret_cast<char**>(AChunk->buf+index*sizeof(T));
		if (c != (char*)((unsigned long)c|1))
			return false;
	}
	return true;

}
/*----------------------------------------------------------------------------*/
template<typename T, int TNbElements>
void ChunkCollection<T,TNbElements>::
initFreeElts()
{
	free_elements_=0;

	MemoryChunk *ch = first_chunk_;

	if(ch==0)
		return;

	while(ch==0){
		for(int index = 0;index<TNbElements;index++){
			char* c= *reinterpret_cast<char**>(ch->buf+index*sizeof(T));
			//look for free items
			if (c == (char*)((unsigned long)c|1)){
				*reinterpret_cast<char**>(ch->buf+index*sizeof(T)) = free_elements_;
				free_elements_ = ch->buf+index*sizeof(T);
			}
		}
	}
}
/*----------------------------------------------------------------------------*/
template<typename T, int TNbElements>
void ChunkCollection<T,TNbElements>::
print()
{
	/* Traversal of the chunk list for making every chunk free */
	std::cout<<"---------------"<<std::endl;
	MemoryChunk* current_chunk = first_chunk_;
	if(current_chunk==0){
		std::cout<<"- empty -"<<std::endl;
		return;
	}
	int chunk_index=0;
	do{
		char* buffer = current_chunk->buf;
		std::cout<<"Chunk "<<++chunk_index<<": ";
		for(int i =0; i<TNbElements;i++){
			std::cout<<(long)(*reinterpret_cast<char**>(buffer + i*sizeof(T)));
			if(i!=TNbElements-1)
				std::cout<<" - ";
		}

			std::cout<<std::endl;
		current_chunk = current_chunk->next;
	} while(current_chunk);
}
/*----------------------------------------------------------------------------*/
template<typename T, int TNbElements>
int ChunkCollection<T,TNbElements>::
getNbChunks() const
{
	int nb=0;
	MemoryChunk* current_chunk = first_chunk_;
	while(current_chunk){
		nb++;
		current_chunk = current_chunk->next;
	}

	return nb;
}
/*----------------------------------------------------------------------------*/
template<typename T, int TNbElements>
void ChunkCollection<T,TNbElements>::
compact()
{
	if(first_chunk_==0)
		return;
	MemoryChunk *current_chunk = first_chunk_, *free_chunk = first_chunk_;
	int index_in_chunk=0, index_in_free_chunk=0;

	// first, we have to find the first free item
	bool not_found_first_free=true;
	int numbering_free, numbering_current;
	numbering_current=numbering_free=1;
	while(not_found_first_free){
		//traversal of current chunk items
		char* c=0;
		while (index_in_chunk!=TNbElements && not_found_first_free){

			c = *reinterpret_cast<char**>(current_chunk->buf+index_in_chunk*sizeof(T));

			if (c == (char*)((unsigned long)c|1)){
				// it is a free item
				free_chunk = current_chunk;
				index_in_free_chunk = index_in_chunk;
				not_found_first_free=false;
			}
			else{
				index_in_chunk++;
			}
		}
		//while (index_in_chunk!=TNbElements && not_found_first_free);

		if(!not_found_first_free){
			break;
		}

		// we reached the end of the chunk
		if(index_in_chunk== TNbElements && current_chunk->next!=0){
			current_chunk= current_chunk->next;
			index_in_chunk=0;
			index_in_free_chunk=0;
			numbering_current++;
			numbering_free++;
		}
		else{
			// we reached the end of the collection without finding a free element
			// It happens when the collection is compact and stores a multiple
			// of TNbElements elements.
			return;
		}

	};

	//traversal of the chunks to find used items and compact them
	while(current_chunk){
		//traversal of current chunk items
		void* c;
		do{
			if(index_in_chunk<TNbElements-1)
				index_in_chunk++;
			c = *reinterpret_cast<char**>(current_chunk->buf+index_in_chunk*sizeof(T));

			if (c != (char*)((unsigned long)c|1)){

				char* free_elt = free_chunk->buf+index_in_free_chunk*sizeof(T);

				// we move the contain of c to ptr free
				memmove(free_chunk->buf+index_in_free_chunk*sizeof(T),
						current_chunk->buf+index_in_chunk*sizeof(T),sizeof(T));
//				memcpy(	free_chunk->buf+index_in_free_chunk*sizeof(T),
//						current_chunk->buf+index_in_chunk*sizeof(T),sizeof(T));

				/* we change the rightest bit to 1 to express this element is free now */
				*reinterpret_cast<char**>(current_chunk->buf+index_in_chunk*sizeof(T)) = (char*)1;
				*reinterpret_cast<char**>(current_chunk->buf+index_in_chunk*sizeof(T) + sizeof(char*)) = free_elements_;
				free_elements_ = current_chunk->buf+index_in_chunk*sizeof(T);

				// we move ptr_free to the next free space
				if(index_in_free_chunk<TNbElements-1)
					index_in_free_chunk++;
				else // go to the next chunk
				{
					free_chunk = free_chunk->next;
					index_in_free_chunk=0;
					numbering_free++;
				}
			}
		}
		while (index_in_chunk!=TNbElements-1);

		// we reached the end of the current chunk
		current_chunk = current_chunk->next;
		index_in_chunk=-1;
		numbering_current++;
	};

	removeEmptyChunks();
}
/*----------------------------------------------------------------------------*/
template<typename T, int TNbElements>
void ChunkCollection<T,TNbElements>::
serialize(std::ostream& AStr)
{
	MemoryChunk* current_chunk = first_chunk_;
	int nbChunks = getNbChunks();
	AStr.write((char*)&nbChunks,sizeof(int));
	while(current_chunk){
		AStr.write((char*)&current_chunk->buf,TNbElements*sizeof(T));
		current_chunk = current_chunk->next;
	} ;

	AStr.write((char*)&nb_stored_elements_,sizeof(int));
	/* as being a pointer, the free elt is not serialized */
}
/*----------------------------------------------------------------------------*/
template<typename T, int TNbElements>
void ChunkCollection<T,TNbElements>::
unserialize(std::istream& AStr)
{
	removeAllChunks();
	int nbChunks;
	AStr.read((char*)&nbChunks,sizeof(int));
	for(int i=0;i<nbChunks;i++){
		MemoryChunk* new_chunk = (MemoryChunk*) malloc(sizeof(MemoryChunk));
			if(!new_chunk)
				throw std::bad_alloc();

		AStr.read((char*)&new_chunk->buf,TNbElements*sizeof(T));

		/* the new chunk is added at the end of the chunk list to preserve order*/
		if(first_chunk_)
		{
			last_chunk_->next = new_chunk;
			new_chunk->prev = last_chunk_;
			new_chunk->next=0;
			last_chunk_ = new_chunk;

		}
		else
		{
			new_chunk->next = 0;
			new_chunk->prev=0;
			first_chunk_ = new_chunk;
			last_chunk_ = new_chunk;
		}
	}
	AStr.read((char*)&nb_stored_elements_,sizeof(int));
	/* we have to compute the first free elt*/
	initFreeElts();
}
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_CHUNKMEMORYALLOCATOR_H_ */
/*----------------------------------------------------------------------------*/
