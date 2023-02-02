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
 * SmartBitVector.h
 *
 *  Created on: 6 févr. 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_SMARTBITVECTOR_H_
#define GMDS_SMARTBITVECTOR_H_
/*----------------------------------------------------------------------------*/

#include <GMDS/Utils/CommonTypes.h>
/*----------------------------------------------------------------------------*/
#include <vector>
/*----------------------------------------------------------------------------*/
namespace gmds {
  /*----------------------------------------------------------------------------*/
  class EXPORT_GMDS  SmartBitVector {

  public:

    class EXPORT_GMDS  Iterator{
    public:

      friend class SmartBitVector;

      Iterator(const SmartBitVector* AContainer){
	m_container=AContainer;
	m_current=0;
	// we get the location of the first real item
	while(m_current!=m_container->m_top && m_container->m_bits[m_current]==0)
	  m_current++;
      }
      Iterator(const Iterator& AIt){
	m_container =AIt.m_container;
	m_current=AIt.m_current;
      }
      Iterator(){m_container=0; m_current=0;}

      void reinit(){m_current=0;}
      TInt value() const {return m_current;}
      void next(){// move to the first real successor of index_
	do{
	  m_current++;
	}while(m_current!=m_container->m_top && m_container->m_bits[m_current]==0);
      }
      bool isDone() const {return m_current==m_container->m_top;}

    private:

      const SmartBitVector* m_container;
      TInt m_current;
    };

    typedef  Iterator iterator;



    /*------------------------------------------------------------------------*/
    /** \brief  Default Constructor
     */
    SmartBitVector(const int capacity=GChunkSize);

    /*------------------------------------------------------------------------*/
    /** \brief  Copy constructor. Note the algorithm is linear in its size.
     */
    SmartBitVector(const SmartBitVector& vec);

    /*------------------------------------------------------------------------*/
    /** \brief  Destructor. Note the algorithm is linear in the number of holes
     * 			in the container.
     */
    ~SmartBitVector();

    /*------------------------------------------------------------------------*/
    /** \brief  Gives the number of elements stored in the container.
     */
    TInt size() const;

    /*------------------------------------------------------------------------*/
    /** \brief  Gives the next valid index after the last used item.
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
    bool  operator[](const TInt& AIndex) const;
    /*------------------------------------------------------------------------*/
    /** \brief  Put a new bit from 0 to 1.
     *
     *  @result the index of this new bit
     */
    TInt  selectNewBit();

    /*------------------------------------------------------------------------*/
    /** \brief  Unselect the bit stored in AIndex.
     *
     * 			Warning - In release mode no test is performed on AIndex. In
     * 			debug mode, AIndex must be in [0, size[
     */
    void unselect(const TInt& AIndex);

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


    /*------------------------------------------------------------------------*/
    /** \brief  Provides an iterator on the first element of the container.
     */
    iterator begin(){
      return Iterator(this);
    }
    /*------------------------------------------------------------------------*/
    /** \brief  Compact the container by removing all free items and provide
     * 			the movement of elements in AMove. Elt in location i is moved
     * 			to location AMove[i].
     *
     */
    void compact(std::vector<TInt>& AMove);

    /*------------------------------------------------------------------------*/
    /** \brief Serialize the variable into stream str. Warning this method does
     * 		   not support typename T where pointers would be present.
     */
    void serialize(std::ostream& stream);

    /*------------------------------------------------------------------------*/
    /** \brief Unserialize the variable from stream str. Warning this method
     * 		   does not support typename T where pointers would be present.
     */
    void unserialize(std::istream& stream);

    /*------------------------------------------------------------------------*/
    /** \brief  Provides the next index where an item can be added. It also
	reserve the necessary space for it and update some internal
	data. It must be used with the assign method.
    */
    inline TInt getFreeIndex(){
      TInt free=0;
      TInt cap=m_capacity;

      if(m_free_stack.empty()){
	free= m_top;
	m_top++;
	if(m_top>=cap){
	  //				m_bits.reserve((cap+1)*2);
	  m_bits.resize((cap+1)*2,false);
	  m_capacity = (cap+1)*2;
	}
      }
      else{
	free = m_free_stack.back();
	m_free_stack.pop_back();
      }

      return free;
    }

    /*------------------------------------------------------------------------*/
    /** \brief  Assign the bit AIndex to true. If AIndex is out of the vector
     * 			bounds nothing is specified to avoid the insertion.
     */
    inline void assign(const TInt& AIndex){
      if(m_bits[AIndex]==0) {
	m_size++;
	m_bits[AIndex]=1;
      }

    }

    /*------------------------------------------------------------------------*/
    /** \brief  Return the vector of the boolean mark.
     *
     * \return the vector of the boolean mark 
     */
    std::vector<bool> getBits() {
      return m_bits;
    }

  protected:
    /* top_ is the index of the the right next item after the last used)*/
    TInt m_top;
    TInt m_size;
    TInt m_capacity;
    std::vector<TInt> m_free_stack;

    std::vector<bool> m_bits;

  };
  /*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_SMARTBITVECTOR_H_ */
/*----------------------------------------------------------------------------*/
