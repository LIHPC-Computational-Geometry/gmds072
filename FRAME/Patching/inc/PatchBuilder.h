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
 * PatchBuilder.h
 *
 *  Created on: sept. 18 2015
 *      Author: Franck Ledoux
 */
/*----------------------------------------------------------------------------*/
#ifndef PATCH_BUILDER_H_
#define PATCH_BUILDER_H_
/*----------------------------------------------------------------------------*/
// GMDS Headers
#include <GMDS/IG/IGMesh.h>
#include <GMDS/Math/Cross2D.h>
#include <GMDS/Math/Segment.h>
/*----------------------------------------------------------------------------*/
// FRAME Headers
#include "PatchComplex2D.h"
/*----------------------------------------------------------------------------*/
/* \class PatchBuilder
 *
 * \brief this class provides services to build patchs
 */
class EXPORT_GMDS PatchBuilder
{
public:
    /*------------------------------------------------------------------------*/
    /* \brief Constructor starting from a mesh having a cross field defined 
     *        onto it.
     *
     * \param AComplex
     * \param AMesh the background mesh we work on
     * \param AField the cross field associated to AMesh
     */
    PatchBuilder(PatchComplex2D* AComplex,
                 gmds::IGMesh*                        AMesh,
                 gmds::Variable<gmds::math::Cross2D>* AField);
    
    
    /*------------------------------------------------------------------------*/
    /** \brief Initialization of the boolean marks used for the algorithm
     *  \param AMarkNodePnt mark for nodes classified on points
     *  \param AMarkNodeCrv mark for nodes classified on curves
     *  \param AMarkEdgeCrv mark for edges classified on curves
     */
    void initMarks(const int AMarkNodePnt, const int AMarkNodeCrv,
                   const int AMarkEdgeCrv);
    
    /*------------------------------------------------------------------------*/
    /* \brief Create a patch structure
     */
    void execute();
    
    /*------------------------------------------------------------------------*/
    /* \brief Build a single patch surrounding AOrigin. This patch covers the
     *        mesh at a maximum distance of ADistance from AOrigin
     *
     * \param AOrigin a mesh node we start from
     * \param ADistance the distance from AOrigin the patch must go
     */
    void createPatch(gmds::Node&   AOrigin,
                     gmds::TCoord ADistance);

private:
    void createPatchFromAGeometricPoint(gmds::Node&   AOrigin,
                                        gmds::TCoord ADistance);
    void createPatchFromAGeometricCurve(gmds::Node&   AOrigin,
                                        gmds::TCoord ADistance);
    void createPatchFromAnInnerPoint(gmds::Node&   AOrigin,
                                     gmds::TCoord ADistance);

    void getPoint(const gmds::math::Point&  AFromPnt,
                  const gmds::math::Vector& AFromDir,
                  const int&                AFromCellDim,
                  const gmds::TCellID&      AFromCellID,
                  const double              ADistanceMax,
                  gmds::math::Point&        AToPnt,
                  gmds::math::Vector&       AToDir,
                  int&                      AToCellDim,
                  gmds::TCellID&            AToCellID,
                  int&                      AEndOnBnd);
private:
    
    /** The patch complex we work with*/
    PatchComplex2D* m_complex;
    /** Mesh we start from */
    gmds::IGMesh* m_mesh;

    /* Cross field we start from*/
    gmds::Variable<gmds::math::Cross2D>* m_field;
    gmds::Variable<gmds::TCoord>* m_distance_field;
    
    /** mark for nodes classified on geometric points */
    int m_mark_nodes_on_point;
    /** mark for nodes classified on geometric curves */
    int m_mark_nodes_on_curve;
    /** mark for edges classified on geometric curves */
    int m_mark_edges_on_curve;
        
    /** mark for the faces we work on*/
    int m_mark_face;

};
/*----------------------------------------------------------------------------*/
#endif /* PATCH_BUILDER_H_ */
/*----------------------------------------------------------------------------*/

