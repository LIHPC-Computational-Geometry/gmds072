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
 * PatchComplex2D.h
 *
 *  Created on: sept. 11, 2015
 *      Author: Franck Ledoux
 */
/*----------------------------------------------------------------------------*/
#ifndef PATCH_COMPLEX_2D_H_
#define PATCH_COMPLEX_2D_H_
/*----------------------------------------------------------------------------*/
// STL Files headers
#include <vector>
#include <map>
/*----------------------------------------------------------------------------*/
// GMDS Headers
#include <GMDS/Math/Point.h>
#include <GMDS/IG/IGMesh.h>
/*----------------------------------------------------------------------------*/
// Patching component headers
#include <Patch2D.h>
#include <PatchEdge2D.h>
#include <PatchVertex2D.h>
/*----------------------------------------------------------------------------*/
/* \class This class describes a complex of 2D patchs. It is made of 2D patchs
 *        that are connected throught a background triangular mesh.
 */

class EXPORT_GMDS PatchComplex2D
{
public:
    
    /*------------------------------------------------------------------------*/
    /* \brief Constructor from a background mesh we work on
     *
     * \param AMesh the background mesh that will structure the complex
     */
    PatchComplex2D(gmds::IGMesh* AMesh);
    /*------------------------------------------------------------------------*/
    /* \brief Destructor
     */
    virtual ~PatchComplex2D();
    
    /*------------------------------------------------------------------------*/
    /* \brief Gives access to the underlying mesh
     */
    gmds::IGMesh* mesh();
    
    /*------------------------------------------------------------------------*/
    /* \brief Create a vertex from a location and a type
     *
     * \param ALoc  the vertex location
     * \param AType the vertex type
     *
     * \return a vertex address
     */
    PatchVertex2D* newVertex(const gmds::math::Point& ALocation,
                             const PatchVertex2D::Type& AType =
                             PatchVertex2D::REGULAR);
    
    /*------------------------------------------------------------------------*/
    /* \brief Create a edge from two vertices from a location and a type
     *
     * \param AV1 a first vertex
     * \param AV2 a second vertex
     *
     * \return an edge address
     */
    PatchEdge2D* newEdge(PatchVertex2D* AV1, PatchVertex2D* AV2);

    /*------------------------------------------------------------------------*/
    /* \brief Create a patch from a set of edges
     *
     * \param AEdges an ordered collection of edges
     *
     * \return a patch address
     */
    Patch2D* newPatch(std::vector<PatchEdge2D*>& AEdges);
    
    
    void write(const std::string& AFileName);

private:
    
    /** Background mesh that structures the complex */
    gmds::IGMesh* m_mesh;
    
    /** the  patches that form the complex */
    std::vector<Patch2D*> m_patchs;
    
    /** the  patche edges that form the complex */
    std::vector<PatchEdge2D*> m_edges;
    
    /** the  patch vertices that form the complex */
    std::vector<PatchVertex2D*> m_vertices;
    
};
/*----------------------------------------------------------------------------*/

#endif /* PATCH_COMPLEX_2D_H_ */
