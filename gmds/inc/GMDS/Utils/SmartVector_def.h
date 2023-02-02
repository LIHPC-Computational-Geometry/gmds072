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
/** \file    SmartVector.h
 *  \author  F. LEDOUX
 *  \date    04/06/2009
 */
/*----------------------------------------------------------------------------*/
/* \class SmartVector
 *
 * \brief Provide a container storing a collection of T objects.
 * 		  Each elements has its own ids, the item where it is stored in
 * 		  *this. This container must not be used like a traditional vector. ...
 *
 * 		  WARNING; This container can only stores object whose size is greater
 * 		  than sizeof(int). Moreover, it must be used for large collections of
 * 		  items. It si otherwise too expensive in memory occupation.
 */
/*----------------------------------------------------------------------------*/
template<typename T>  class EXPORT_GMDS SmartVector{

public:

	class EXPORT_GMDS Iterator {
	public:

		friend class SmartVector<T>;

		Iterator(const SmartVector<T>* AContainer){
			container_=AContainer;
			current_=0;
			// we get the location of the first real item
			while(current_!=container_->top_ && container_->mark_[current_]==0)
				current_++;
		}
		Iterator(const Iterator& AIt){
			container_ =AIt.container_;
			current_=AIt.current_;
		}
		Iterator(){container_=0; current_=0;}

		virtual void reinit(){current_=0;}
		virtual T currentItem() const {	return container_->vec_[current_];}
		virtual void next(){// move to the first real successor of index_
			do{
				current_++;
			}
			while(current_!=container_->top_ && container_->mark_[current_]==false);
		}
		virtual bool isDone() const {return current_==container_->top_;}

		virtual void operator()(){;}
		virtual T operator*(){	return container_->vec_[current_]; }

		virtual void operator++(){
			// move to the first real successor of index_
			do{
				current_++;
			}while(current_!=container_->top_ && container_->mark_[current_]==0);
		}

	private:

		const SmartVector<T>* container_;
		TInt current_;
	};


	typedef  Iterator iterator;

	/*------------------------------------------------------------------------*/
	/** \brief  Default Constructor
	 */
	SmartVector(const int capacity=GChunkSize);

	/*------------------------------------------------------------------------*/
        /** \brief  Default Constructor
         */
        SmartVector(std::vector<bool> AMarks);

	/*------------------------------------------------------------------------*/
	/** \brief  Copy constructor. Note the algorithm is linear in its size.
	 */
	SmartVector(const SmartVector<T>& vec);

	/*------------------------------------------------------------------------*/
	/** \brief  Destructor. Note the algorithm is linear in the number of holes
	 * 			in the container.
	 */
	~SmartVector();

	/*------------------------------------------------------------------------*/
	/** \brief  Gives the number of elements stored in the container.
	 */
	TInt size() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Gives the next right index after the last used item.
	 */
	TInt top() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Gives the size of the container, i.e, the effective number of
	 * 			items that are available in the container. It corresponds to
	 * 			the number of stored elements + the number of free spaces.
	 */
	TInt capacity() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Indicates if the container is empty.
	 */
	bool empty() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Clear the container and resizes it to 2.
	 */
	void clear();
	/*------------------------------------------------------------------------*/
	/** \brief  resizes the container to ASize. size and top are  both
	 * 			equal to ASize. If ASize>current(size) it is just extended,
	 * 			new items are marked free. Otherwise, all the items between
	 * 			ASize and current(size) are lost.
	 */
	void resize(const TInt& ASize);

	/*------------------------------------------------------------------------*/
	/** \brief  make the container full of elements putting all the marks to 1
	 * 		    and top_ equals to size_;
	 */
	void fillAll();

	/*------------------------------------------------------------------------*/
	/** \brief  This method is necessary when you want to regularize the
	 * 			container after several assignements. Indeed assignement method
	 * 			has the advantage to check nothing before insertion. The problem
	 * 			is then that the container is no more coherent. This method fix
	 * 			the container.
	 */
	void update();

	/*------------------------------------------------------------------------*/
	/** \brief  Give the value stored in the AIndex item but do not allow to
	 * 			modify it.
	 *
	 * 			Warning, in release mode, no test is performed to ensure the
	 * 			choice of AIndex.
	 */
	inline T const& operator[](const TInt& AIndex) const{
#ifdef __DEBUG__
		if(AIndex>=capacity() || mark_[AIndex]==0)
			throw GMDSException("Bad index in SmartVector<T>::operator[]");
#endif //__DEBUG__

		return vec_[AIndex];

	}
	/*------------------------------------------------------------------------*/
	/** \brief  Give the value stored in the AIndex item and allow to
	 * 			modify it.
	 *
	 * 			Warning, in release mode, no test is performed to ensure the
	 * 			choice of AIndex. Moreover, in order to improve performances,
	 * 			the container structure is partially corrupted. Thus an update
	 * 			is necessary at the end (especially for the iterators)
	 */
	inline T& operator[](const TInt& AIndex){
#ifdef __DEBUG__
		if(AIndex>=capacity() || mark_[AIndex]==0)
			throw GMDSException("Bad index in SmartVector<T>::operator[]");
#endif //__DEBUG__
		if(!mark_[AIndex])
		{
			size_++;
			mark_[AIndex]=1;

		}
		return vec_[AIndex];

	}


	/*------------------------------------------------------------------------*/
	/** \brief  Provides the next index where an item can be added. It also
				reserve the necessary space for it and update some internal
				data. It must be used with the assign method.
	 */

	TInt  selectNewIndex();

	/*------------------------------------------------------------------------*/
	/** \brief  Assign AElt to index AIndex. If this index is still occupied,
	 * 			it contents is replaced. If AIndex is out of the vector bondary
	 * 			nothing is specified to avoid the insertion.
	 */
	void assign(T AElt, const TInt& AIndex);

	/*------------------------------------------------------------------------*/
	/** \brief  Add AElt. The container is responsible of finding the index
	 * 			where to add the element. This operation is SAFER than assign.
	 */
	void add(T AElt);

	/*------------------------------------------------------------------------*/
	/** \brief  find if an AElt is inside the vector. It returns true if it is,
	 * 			false otherwise.
	 */
	bool find(T AElt) ;//const;

	/*------------------------------------------------------------------------*/
	/** \brief  Remove the element stored in AIndex.
	 *
	 * 			Warning - In release mode no test is performed on AIndex. In
	 * 			debug mode, AIndex must be in [0, size[
	 */
	void remove(const TInt& AIndex);

	/*------------------------------------------------------------------------*/
	/** \brief  Remove the element AElt
	 */
	void removeElement(const T& AElt);

	/*------------------------------------------------------------------------*/
	/** \brief  Remove the element in position AIterator
	 */
	void remove(Iterator& AIterator);

	/*------------------------------------------------------------------------*/
	/** \brief  Indicates if the AIndex item is available.
	 */
	bool isAvailable(const TInt& AIndex) const;
	/*------------------------------------------------------------------------*/
	/** \brief  Indicates if the AIndex item is out of the container.
	 */
	bool isOutOfContainer(const TInt& AIndex) const;

	/*------------------------------------------------------------------------*/
	/** \brief  Display the content of the container
	 */
	void display();


	void mark(const std::vector<int>& ref);

	/*------------------------------------------------------------------------*/
	/** \brief  Provides an iterator on the first element of the container.
	 */
	iterator begin(){
		return Iterator(this);
	}

	/*------------------------------------------------------------------------*/
	/** \brief  Provides an iterator on the first element of the container.
	 */
//	const_handle_iterator begin() const;

	/*------------------------------------------------------------------------*/
	/** \brief  Compact the container by removing all free items.
	 */
	void compact();
	/*------------------------------------------------------------------------*/
	/** \brief  Compact the container by removing all free items and provide
	 * 			the movement of elements in AMove. Elt in location i is moved
	 * 			to location AMove[i].
	 *
	 */
	void compactWithMemory(std::vector<int>& AMove);

	/*------------------------------------------------------------------------*/
	/** \brief Serialize the variable into stream str. Warning this method does
	 * 		   not support typename T where pointers would be present.
	 */
	inline void serialize(std::ostream& stream);

	/*------------------------------------------------------------------------*/
	/** \brief Unserialize the variable from stream str. Warning this method
	 * 		   does not support typename T where pointers would be present.
	 */
	inline void unserialize(std::istream& stream);

	/*------------------------------------------------------------------------*/
        /** \brief  Return the vector of the boolean mark.
         *
         * \return the vector of the boolean mark
         */
        std::vector<bool> getMarks() {
                return mark_;
        }

protected:
	/* top_ is the index of the the right next item after the last used)*/
	TInt top_;
	TInt size_;
	std::vector<TInt> free_stack_;
	std::vector<T> vec_;
	/* mark_[i] = 0 indicates that vec_[i] is free, otherwise it is not
	 * available */
	std::vector<bool> mark_;
};
/*----------------------------------------------------------------------------*/
