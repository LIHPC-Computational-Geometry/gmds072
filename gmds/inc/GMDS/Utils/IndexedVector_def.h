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
/** \file    IndexedVector_def.h
 *  \author  F. LEDOUX
 *  \date    02/06/2014
 */
/*----------------------------------------------------------------------------*/
/* \class IndexedVector
 *
 * \brief Provide a container storing a collection of T objects.
 */
/*----------------------------------------------------------------------------*/
template<typename T> class IndexedVector{

public:
	/*------------------------------------------------------------------------*/
	/** \brief  Default Constructor
	 */
	IndexedVector(const int capacity=2);

	/*------------------------------------------------------------------------*/
	/** \brief  Copy constructor. Note the algorithm is linear in its size.
	 */
	IndexedVector(const IndexedVector<T>& vec);

	/*------------------------------------------------------------------------*/
	/** \brief  Destructor. Note the algorithm is linear in the number of holes
	 * 			in the container.
	 */
	~IndexedVector();
	/*------------------------------------------------------------------------*/
	/** \brief  Gives the size of the container, i.e, the effective number of
	 * 			items that are available in the container. It corresponds to
	 * 			the number of stored elements + the number of free spaces.
	 */
	TInt capacity() const;
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
	/** \brief  Give the value stored in the AIndex item but do not allow to
	 * 			modify it.
	 *
	 * 			Warning, in release mode, no test is performed to ensure the
	 * 			choice of AIndex.
	 */
	 T const& operator[](const TInt& AIndex) const;
	/*------------------------------------------------------------------------*/
	/** \brief  Give the value stored in the AIndex item and allow to
	 * 			modify it.
	 *
	 * 			Warning, in release mode, no test is performed to ensure the
	 * 			choice of AIndex. Moreover, in order to improve performances,
	 * 			the container structure is partially corrupted. Thus an update
	 * 			is necessary at the end (especially for the iterators)
	 */
	 T& operator[](const TInt& AIndex);

	/*------------------------------------------------------------------------*/
	/** \brief  Assign AElt to index AIndex. If this index is still occupied,
	 * 			it contents is replaced. If AIndex is out of the vector
	 * 			limits, nothing is specified to avoid the insertion.
	 */
	void assign(T& AElt, const TInt& AIndex);


	/*------------------------------------------------------------------------*/
	/** \brief Serialize
	 */
	void serialize(std::ostream& stream);

	/*------------------------------------------------------------------------*/
	/** \brief Unserialize
	 */
	void unserialize(std::istream& stream);

protected:
	std::vector<T> m_vec;
};
/*----------------------------------------------------------------------------*/
