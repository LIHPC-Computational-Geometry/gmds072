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
 * PatchEdge2D.h
 *
 *  Created on: sept. 18 2015
 *      Author: Franck Ledoux
 */
/*----------------------------------------------------------------------------*/
#ifndef PATCH_EDGE_2D_H_
#define PATCH_EDGE_2D_H_
/*----------------------------------------------------------------------------*/
// STL Headers
#include <vector>
/*----------------------------------------------------------------------------*/
// GMDS Headers
#include <GMDS/IG/IG.h>
/*----------------------------------------------------------------------------*/
// FRAME Headers
#include "PatchVertex2D.h"
/*----------------------------------------------------------------------------*/
class EXPORT_GMDS PatchEdge2D {
    
public:
    
    /*------------------------------------------------------------------------*/
    /* \brief Constructor
     *
     * \param AV1 A first patch vertex
     * \param AV2 A second patch vertex
     */
    PatchEdge2D(PatchVertex2D* AV1, PatchVertex2D* AV2);
    
    /*------------------------------------------------------------------------*/
    /* \brief Accessor on  the first end vertex
     *
     * \return the first end vertex
     */
    PatchVertex2D* firstVertex();
    /*------------------------------------------------------------------------*/
    /* \brief Accessor on the second end vertex
     *
     * \return the second end vertex
     */
    PatchVertex2D* secondVertex();
    
    /*------------------------------------------------------------------------*/
    /* \brief Modificator of the first end vertex
     *
     * \param a new vertex
     */
    void setFirstVertex(PatchVertex2D* AV);
    
    /*------------------------------------------------------------------------*/
    /* \brief Modificator of the second end vertex
     *
     * \param a new vertex
     */
    void setSecondVertex(PatchVertex2D* AV);
    
    
    /*------------------------------------------------------------------------*/
    /* \brief Accessor on the discretization of the edge
     *
     * \return a vector of points
     */
    std::vector<gmds::math::Point>& discretization();
    
    /*------------------------------------------------------------------------*/
    /* \brief set the discretization of the edge
     *
     * \param APoints a vector of ordered points
     */
    void setDiscretization(std::vector<gmds::math::Point>& APnts);
    
    /*------------------------------------------------------------------------*/
    /* \brief Accessor on the ids corresponding to the traversed mesh
     *        faces
     *
     * \return a vector of ids
     */
    std::vector<gmds::TCellID>& crossedMeshFaces();
    
    /*------------------------------------------------------------------------*/
    /* \brief set the collection of ids corresponding to the traversed mesh 
     *        faces
     *
     * \param AIDs a vector of face ids
     */
    void setCrossedMeshFaces(std::vector<gmds::TCellID>& AIDs);
    
    
private:
    
    /** first end vertex */
    PatchVertex2D* m_first_vertex;
    /** firs end vertex */
    PatchVertex2D* m_second_vertex;
    
    /** discretization point if the edge is not a segment */
    std::vector<gmds::math::Point> m_discretization;
    
    /** ids of mesh faces crossed by the edge */
    std::vector<gmds::TCellID> m_crossed_faces;
 };
/*----------------------------------------------------------------------------*/
#endif /* PATCH_EDGE_2D_H_ */
/*----------------------------------------------------------------------------*/
