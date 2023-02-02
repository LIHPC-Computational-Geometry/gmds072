/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux (2015)
 *
 * franck.ledoux@cea.fr
 *
 * The FRAME software is a computer program whose purpose is to provide a set
 * of algorithms to build 2D and 3D meshes using frame field concept. The main
 * focus of these algorithms is quadrilateral and hexahedral meshing.
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
 * PatchVertex2D.h
 *
 *  Created on: sept. 18 2015
 *      Author: Franck Ledoux
 */
/*----------------------------------------------------------------------------*/
#ifndef PATCH_VERTEX_2D_H_
#define PATCH_VERTEX_2D_H_

/*----------------------------------------------------------------------------*/
// GMDS Headers
#include <GMDS/IG/IG.h>
/*----------------------------------------------------------------------------*/
// FRAME Headers
/*----------------------------------------------------------------------------*/
class EXPORT_GMDS PatchVertex2D {

public:
    
    
    /* \enum vertices of patchs can be regular (jointed to 4 patches or 3 if
     *       on on the boundary), extraordinary (if corresponding to a field
     *       singularity or a non-convex boundary point) or t-joint if it is
     *       a non conform point between patchs.
     */
    typedef enum {REGULAR, EXTRAORDINARY, TJOINT} Type;
    
    /*------------------------------------------------------------------------*/
    /* \brief Constructor
     *
     * \param APnt  the vertex location
     * \param AType the type of patch vertex
     */
    PatchVertex2D(const gmds::math::Point& APnt, const Type AType);
    
    /*------------------------------------------------------------------------*/
    /* \brief Accessor on the location
     *
     * \return the vertex location
     */
    gmds::math::Point& location();
    
    /*------------------------------------------------------------------------*/
    /* \brief const accessor on the location
     *
     * \return the vertex location
     */
    gmds::math::Point location() const;
    
    /*------------------------------------------------------------------------*/
    /* \brief Accessor on the vertex type
     *
     * \return the vertex type
     */
    Type type() const;
    
private:
    
    /** vertex location */
    gmds::math::Point m_location;
    /** vertex type */
    Type m_type;
};
/*----------------------------------------------------------------------------*/
#endif /* PATCH_VERTEX_2D_H_ */
/*----------------------------------------------------------------------------*/
