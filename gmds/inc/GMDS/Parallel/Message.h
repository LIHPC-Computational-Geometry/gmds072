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
 * DistributedCellData.h
 *
 *  Created on: 10 juil. 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MESSAGE_H_
#define GMDS_MESSAGE_H_
/*----------------------------------------------------------------------------*/
#ifdef GMDS_PARALLEL
/*----------------------------------------------------------------------------*/
#include <iostream>
/*----------------------------------------------------------------------------*/
namespace gmds {
    /*------------------------------------------------------------------------*/
    /** \class Message
     *
     * \brief Provide an object structure to build message that will be sent
     *        using the DistributedManager. Note that the message does not
     *        know its inner structure. If you add 3 double, then 1 string 
     *        into it, the receiver must know this structure in order to read 
     *        it
     */
class Message
{

public:
    /*------------------------------------------------------------------------*/
    /** \brief  Constructor.
     *
     * \param[in] the message size.
     */
    Message(const int ASize=0);
    
    /*------------------------------------------------------------------------*/
    /** \brief  Acces to the message size.
     */
    int size() const;
    
    /*------------------------------------------------------------------------*/
    /** \brief  Acces to the message contents.
     */
    char* contents();
    /*------------------------------------------------------------------------*/
    /** \brief  Resize the message to be equal to  \p ASize
     *
     * \param[in] ASize the new message size
     */
    void reserve(const int ASize);
    
    /*------------------------------------------------------------------------*/
    /** \brief  Add a vector of T-type elements at the end of the message. If 
     *          their is not enough room in the message, its size is double on
     *          the fly.
     *
     * \param[in] AElts the elements to add
     */
    template<class T> void push_back(std::vector<T>& AElts);
    template<class T> operator<<(const std::vector<T>& AElt);
    
    /*------------------------------------------------------------------------*/
    /** \brief  Add a T-type element at the end of the message. If
     *          their is not enough room in the message, its size is double on
     *          the fly.
     *
     * \param[in] AElt the single element to add
     */
    template<class T> void push_back(T& AElt);
    template<class T> operator<<(const T& AElt);
    

    template<class T> void pop_front(T& AElts);
    template<class T> operator>>(const T& AElt);
    template<class T> void pop_front(std::vector<T>& AElts);
    template<class T> operator>>(const std::vector<T>& AElt);
private:

    /** the contents of the message we handle*/
    char* m_contents;
    /** message size */
    int m_size;
    
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif //GMDS_PARALLEL
/*----------------------------------------------------------------------------*/
#endif /* GMDS_DISTRIBUTEDCELLDATA_H_ */
/*----------------------------------------------------------------------------*/
