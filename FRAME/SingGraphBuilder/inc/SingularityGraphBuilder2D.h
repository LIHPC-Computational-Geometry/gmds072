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
 * SingularityGraphBuilder2D.h
 *
 *  Created on: April 10 2014
 *      Author: F. Ledoux
 */
/*----------------------------------------------------------------------------*/
#ifndef SINGULARITYGRAPHBUILDER_2D_H_
#define SINGULARITYGRAPHBUILDER_2D_H_
/*----------------------------------------------------------------------------*/
#include <cstdlib>
#include <iostream>
#include <string>
#include <list>
#include <map>
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IGMesh.h>
#include <GMDS/Math/Quaternion.h>
/*----------------------------------------------------------------------------*/
#include "SingularityGraph.h"
/*----------------------------------------------------------------------------*/
#include <Tools.h>
/*----------------------------------------------------------------------------*/
/** \brief Class providing an algorithm to build a 2D singularity graph from
 *         a triangular mesh and a 2D cross field defined on this mesh
 */
class EXPORT_GMDS SingularityGraphBuilder2D
{
public:
    /*------------------------------------------------------------------------*/
    /** \brief Constructor.
     *
     * \param AMesh the mesh where we work on
     * \param AField the cross field associated to AMesh
     * \param ATolerance the tolerance used to connect singularity lines and
     *                   points. It is a % of the boundary box diagonal. It
     *                   must be set in [0.01,0.1]
     * \param ABuildGeomSing flag to build the geometric singularity points.
     */
    SingularityGraphBuilder2D(gmds::IGMesh* AMesh,
                              gmds::Variable<gmds::math::Cross2D>* AField,
                              const double ATolerance = 0.05,
                              const bool ABuildGeomSing = true);
    
    /*------------------------------------------------------------------------*/
    /** \brief Execution of the algorithm
     */
    void execute();
    
    /*------------------------------------------------------------------------*/
    /** \brief Function to give the directory where we want to put the output
     *		   files (vtk files)
     */
    void setDebugPrefix(const std::string& ADirName) {
        m_output_directory_name = ADirName;
    }
    
    
    /*------------------------------------------------------------------------*/
    /** \brief Initialization of the boolean marks used for the algorithm
     * \param AMarkNodePnt mark for nodes classified on points
     * \param AMarkNodeCrv mark for nodes classified on curves
     * \param AMarkEdgeCrv mark for edges classified on curves
     */
    void initMarks(const int AMarkNodePnt, const int AMarkNodeCrv,
                   const int AMarkEdgeCrv);
    
    
protected:
    void colorFaces(const int markF, const int markE);
    
    /*------------------------------------------------------------------------*/
    /** \brief  Detect the triangles that contains singularity points
     */
    void detectSingularTriangles();
    
    /*------------------------------------------------------------------------*/
    /** \brief  For a a face AFace, which is singular, this method computes the
     *          geometric location of the singularity point in AFace
     *
     * \param[IN ] AFace    a singular face
     * \param[OUT] APosSing the singularity point location
     */
    void computeSingPointInfo( gmds::Face& AFace,
                              gmds::math::Point& APosSing);
    
    /*------------------------------------------------------------------------*/
    /** \brief  For a a face AFace, which is singular, this method computes the
     *          geometric location of the singularity point in AFace
     */
    void createSingPointAndSlots(gmds::Face& AFace);
    
    /*------------------------------------------------------------------------*/
    /** \brief Add the boundary nodes and edges into the singularity graph.
     *
     * \param ABuildGeomSlots indicates if we want to build geometric slots
     *        (default is false). If true, it means the resulting domain
     *        partitioning will be made of 4-sided patches only.
     */
    void addGeometryToSingularityGraph(const bool ABuildGeomSlots = false);
    
    void writeOutput(const std::string& AFileName);
    void writeOutputSingle(const std::string& AFileName);
    /*------------------------------------------------------------------------*/
    /** \brief Creation of singularity points
     */
    //  void createSingularityPoints();
    
    /*------------------------------------------------------------------------*/
    /** \brief Creation of singularity lines
     */
    void createSingularityLines();
    /*------------------------------------------------------------------------*/
    /** \brief Compute the singularity line starting from point ASingPnt in the
     *         direction of ASingSlot
     *
     * \param ASingPnt  the singularity point we start from
     * \param ASingSlot the associated slot we really consider
     */
    void computeSingularityLine(SingularityPoint* ASingPnt,
                                SingularityPoint::Slot* ASingSlot);
    /*------------------------------------------------------------------------*/
    /** \brief Compute the singularity line starting from point ASingPnt in the
     *         direction of ASingSlot
     *
     * \param[IN ] AFromPnt      the singularity point we start from
     * \param[IN ] AFromSlot     the associated slot we really consider
     * \param[OUT] AToSingPnt    the singularity point we arrive at, if it exists
     * \param[OUT] AToSlot       the associated slot
     * \param[OUT] AToPnt        the location we arrive at the end
     * \param[OUT] AToDir        the last direction to arrive at the end
     * \param[OUT] APoints       a discretization of the stream line
     * \param[OUT] ATriangles    the list of traversed triangles
     * \param[OUT] AToCellDim    the dimension of the cell we finish on (0,1,2)
     * \param[OUT] AToCellID     the id of the cell we finish on
     * \param[OUT] AEndOnBnd     indicates if we finish on the boundary (true)
     * \param[OUT] AToSlotIsFree indicates if we finish onto a free slot (true)
     * \param[OUT] APntToCreate  indicates if we must create the end point (true)
     
     */
    void computeStreamLine(SingularityPoint*               AFromPnt,
                           SingularityPoint::Slot*         AFromSlot,
                           SingularityPoint*&              AToSingPnt,
                           SingularityPoint::Slot*&        AToSlot,
                           gmds::math::Point&              AToPnt,
                           gmds::math::Vector&             AToDir,
                           std::vector<gmds::math::Point>& APoints,
                           std::vector<gmds::TCellID>&     ATriangles,
                           int&                            AToCellDim,
                           gmds::TCellID&                  AToCellID,
                           bool&                           AEndOnBnd,
                           bool&                           AToSlotIsFree,
                           bool&                           APntToCreate);
    
    /*------------------------------------------------------------------------*/
    /** \brief Improve a singularity line to avoid side effects due to the way
     *         singularity points are considered.
     *
     * \param ALine     the singularity line we work on
     * \param AFromPnt  the singularity point we start from
     * \param AFromSlot the associated slot we really consider
     * \param AToPnt    the singularity point we go to
     * \param AToSlot   the associated slot we really consider
     */
    void backtrackSingularityLine(SurfaceSingularityLine* ALine,
                                  SingularityPoint*       AFromPnt,
                                  SingularityPoint::Slot* AFromSlot,
                                  SingularityPoint*       AToPnt,
                                  SingularityPoint::Slot* AToSlot);
    
    /*------------------------------------------------------------------------*/
    /** \brief Initialize the confusing ball for singularity APnt
     *
     * \param APnt the singularity point we work on
     */
    void initConfusingBalls(SingularityPoint* APnt);
    
    /*------------------------------------------------------------------------*/
    /** \brief Creates a geometric singularity point.
     *
     * \param [IN]  AInPnt   the point where we are located
     * \param [IN]  AInVec   the geometric direction we arrive
     * \param [IN]  ACellDim the dimension of the cell we are located on
     * \param [IN]  ACellID  the id of the cell we are located on
     * \param [OUT] APnt the created singularity point
     * \param [OUT] ASlot the associated slot directed towards AInVec
     *
     */
    void createGeometricSingularityPoint(const gmds::math::Point&  AInPnt,
                                         const gmds::math::Vector& AInVec,
                                         const int                 ACellDim,
                                         const gmds::TCellID       ACellID,
                                         SingularityPoint*&         APnt,
                                         SingularityPoint::Slot*&   ASlot);
    
    /*------------------------------------------------------------------------*/
    /** \brief Returs the singularity lines that go througt the face AFace
     *
     * \param AFace a mesh face
     *
     * \return the singularity lines going through AFace
     */
    std::vector<SurfaceSingularityLine*> getSingularityLinesIn(const gmds::Face& AFace);
    
    /*------------------------------------------------------------------------*/
    /** \brief This operation detects line intersection in steady areas and
     *         create singularity points
     */
    void detectLineIntersections();
    
    /*------------------------------------------------------------------------*/
    /** \brief This operation create intersection between 2 lines in the
     *         vicinity of a mesh face
     */
    void createLineIntersection(SurfaceSingularityLine *ALine1,
                                SurfaceSingularityLine *ALine2,
                                gmds::Face& AFace,
                                std::vector<gmds::math::Point>& AAddedPoints);
    void createLineIntersection(std::vector<SurfaceSingularityLine*>& ALines,
                                gmds::Face& AFace,
                                std::vector<gmds::math::Point>& AAddedPoints);
    
    void writeSingularityPointsAndSlots();
private:
    
    /** Mesh we start from */
    gmds::IGMesh* m_mesh;
    /* Cross field we start from*/
    gmds::Variable<gmds::math::Cross2D>* m_field;
    
    /** DEBUG variable which stores the index of every single triangle of m_mesh */
    gmds::Variable<int>* m_index;
    
    Tools m_tool;
    /**directory where debug files will be written*/
    std::string m_output_directory_name;
    
    /** flag that indicates if geometric singularity point must be built */
    bool m_build_geometric_singularities;
    /** mark for nodes classified on geometric points */
    int m_mark_nodes_on_point;
    /** mark for nodes classified on geometric curves */
    int m_mark_nodes_on_curve;
    /** mark for edges classified on geometric curves */
    int m_mark_edges_on_curve;
    
    /**the obtained singularity graph*/
    SingularityGraph m_graph;
    /** technical container, which is used to stor all the free slots during the
     *  algorithm */
    std::list<SingularityPoint::Slot*> m_free_slots;
    
    /** radius of the bounding box containing m_mesh. It is computed at the
     initial step*/
    double m_mesh_radius;
    /** confusing distance used for connecting graph lines and points */
    double m_confusing_distance;
    
    /** list of faces containing a 3-valent singularity point */
    std::list<gmds::Face> m_singularities_3;
    /** list of faces containing a 5-valent singularity point */
    std::list<gmds::Face> m_singularities_5;
    
    /** mark for faces containing a singularity point of the cross field*/
    int m_mark_faces_with_sing_point;
    /** mark for faces traversed by a singularity line of the cross field*/
    int m_mark_faces_with_sing_line;
    
    /** map used by the confusing ball. For faces, edges and nodes that are not
     located into a singularity, maps return value 0. */
    std::map<gmds::TCellID, SingularityPoint*> m_faces_to_singularity_on_surf;
    std::map<gmds::TCellID, SingularityPoint*> m_edges_to_singularity_on_surf;
    std::map<gmds::TCellID, SingularityPoint*> m_nodes_to_singularity_on_surf;
    
    /*std::list<gmds::TCellID> m_2SingTetIDsAlongSurf;
     std::map<gmds::TCellID, int> m_region_singularity_type;
     std::map<gmds::TCellID, SurfaceSingularityPoint*> m_face_to_singularity_on_surf;
     
     //boolean marks for executing the algorithm
     
     int m_markNodesOnVert;
     int m_markNodesOnCurv;
     int m_markNodesOnSurf;
     int m_markAloneNodes;
     
     int m_markEdgesOnSurf;
     int m_markEdgesOnCurv;
     
     int m_markFacesOnSurf;
     
     int m_markClusterSingDone;
     int m_markBordVolSingForFaces;
     int m_markBordVolSingForEdges;
     int m_markBordVolSingForNodes;
     
     int m_markBordSurfSingForFaces;//faces containing sing. are marked with it
     int m_markBordSurfSingForEdges;//edges containing sing. are marked with it
     int m_markBordSurfSingForNodes;//nodes containing sing. are marked with it
     */
};
/*----------------------------------------------------------------------------*/
#endif /* SINGULARITYGRAPHBUILDER_H_ */
/*----------------------------------------------------------------------------*/
