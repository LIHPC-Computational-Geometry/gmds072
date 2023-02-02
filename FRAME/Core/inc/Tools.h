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
 * ToolKit.h
 *
 *  Created on: sept. 18, 2015
 *      Author: Franck Ledoux
 */
/*----------------------------------------------------------------------------*/
#ifndef FRAME_TOOLKIT_H_
#define FRAME_TOOLKIT_H_
/*----------------------------------------------------------------------------*/
// STL Headers
/*----------------------------------------------------------------------------*/
// GMDS Headers
#include <GMDS/IG/IG.h>
#include <GMDS/Math/Point.h>
#include <GMDS/Math/Vector.h>
#include <GMDS/Math/Cross2D.h>
#include <GMDS/Math/Chart.h>
#include <GMDS/Math/Quaternion.h>
#include <GMDS/Math/Triangle.h>
#include <GMDS/Math/Segment.h>
#include <GMDS/Math/AxisAngleRotation.h>
/*----------------------------------------------------------------------------*/
// Predicates_psm File Headers
#include <Predicates_psm.h>
/*----------------------------------------------------------------------------*/
/* \class This class gathers elementary algorithms used for different purposes
 *        by the main algorithms provided in FRAME.
 */
class EXPORT_GMDS Tools
{
public:
    
    struct PointVolumetricData{
        gmds::math::Point pnt;      // Point we are located at
        gmds::math::Vector3d dir;   // Direction we propage to
        gmds::Region tet;           // tetra where pnt is located in
        
        PointVolumetricData(const gmds::math::Point& AP,
                            const gmds::math::Vector3d& AV,
                            const gmds::Region& AR)
        :pnt(AP),dir(AV),tet(AR){;}
    };
    struct PointSurfacicData{
        gmds::math::Point pnt;    // Point we are located at
        gmds::math::Vector3d dir; // Direction we propage to
        gmds::Face tri;           // triangle where pnt is located in
        
        PointSurfacicData(const gmds::math::Point& AP,
                            const gmds::math::Vector3d& AV,
                            const gmds::Face& AF)
        :pnt(AP),dir(AV),tri(AF){;}
    };
    /*------------------------------------------------------------------------*/
    /** \brief Constructor.
     *
     * \param AMesh the mesh where we work on
     * \param AField the cross field associated to AMesh
     */
    Tools(gmds::IGMesh* AMesh,
          gmds::Variable<gmds::math::Cross2D>* AField,
          gmds::Variable<gmds::math::AxisAngleRotation>* ARotField=0);

    /*------------------------------------------------------------------------*/
    /** \brief  Compute where the frame field defined by m_axis_angle will
     *          transport \p AData.pnt considering data in \p AData and the max
     *          distance is \p AMaxDist.
     *
     * \details This algorithm uses Heun's scheme + a fuzzy approch where we 
     *          always go in/out of simplex along faces
     *
     * \param[in] AData     point data in the flow
     * \param[in] AMaxDist  Maximum distance allowed for the displacement
     *
     * \param[out] APnt the reached location.
     * 
     * \return true if APnt has been computed, false if we reach a tet 
     *              being FF singular
     */
    bool followFlow(const PointVolumetricData& AData,
                    const double               AMaxDist,
                    gmds::math::Point&         APnt);
    
    bool followFlow(const PointSurfacicData& AData,
                    const double             AMaxDist,
                    const int                AMarkEdgeOnCurve,
                    gmds::math::Point&       APnt);
    
    bool isFFSingular(const gmds::Region& AR);
    bool isFFSingular(const gmds::Face& AF);
    gmds::math::Chart::Mapping getRij(const gmds::TCellID AFrom,
                                      const gmds::TCellID ATo) const;
    
    /*------------------------------------------------------------------------*/
    /** \brief Gives the out point when we go through a triangle following the
     *         cross field
     *
     * \param AFace       the face we work on
     * \param AInPnt      the geometric point we start from
     * \param AInVec      the geometric direction to follow
     * \param AInCellDim  the dimension of the cell start_pnt is located
     * \param AInCellID   the id of the cell start_pnt is located on
     * \param AOutPnt     the geometric point we go out
     * \param AOutVec     the geometric direction to follow after
     * \param AOutCellDim the dimension of the cell we go out
     * \param AOutCellID  the id of the cell we go out
     */
    void	traverseTriangle(const gmds::Face&         AFace,
                             const gmds::math::Point&  AInPnt,
                             const gmds::math::Vector& AInVec,
                             const int                 AInCellDIm,
                             const gmds::TCellID       AInCellID,
                             gmds::math::Point&        AOutPnt,
                             gmds::math::Vector&       AOutVec,
                             int&                      AOutCellDIm,
                             gmds::TCellID&            AOutCellID);
    
    /*------------------------------------------------------------------------*/
    /** \brief Gives the out point when we go through a triangle following the
     *         cross field and starting from a node
     *
     * \param AFace       the face we work on
     * \param ANode       the node we come from
     * \param AInPnt      the geometric point we start from
     * \param AInVec      the geometric direction to follow
     * \param AOutPnt     the geometric point we go out
     * \param AOutVec     the geometric direction to follow after
     * \param AOutCellDim the dimension of the cell we go out
     * \param AOutCellID  the id of the cell we go out
     */
    void	traverseTriangle(const gmds::Face&         AFace,
                             const gmds::Node&         ANode,
                             const gmds::math::Point&  AInPnt,
                             const gmds::math::Vector& AInVec,
                             gmds::math::Point&        AOutPnt,
                             gmds::math::Vector&       AOutVec,
                             int&                      AOutCellDIm,
                             gmds::TCellID&            AOutCellID);
    
    /*------------------------------------------------------------------------*/
    /** \brief Gives the out point when we go through a triangle following the
     *         cross field and starting from an edge
     *
     * \param AFace       the face we work on
     * \param AEdge       the edge we come from
     * \param AInPnt      the geometric point we start from
     * \param AInVec      the geometric direction to follow
     * \param AOutPnt     the geometric point we go out
     * \param AOutVec     the geometric direction to follow after
     * \param AOutCellDim the dimension of the cell we go out
     * \param AOutCellID  the id of the cell we go out
     */
    void	traverseTriangle(const gmds::Face&         AFace,
                             const gmds::Edge&         AEdge,
                             const gmds::math::Point&  AInPnt,
                             const gmds::math::Vector& AInVec,
                             gmds::math::Point&        AOutPnt,
                             gmds::math::Vector&       AOutVec,
                             int&                      AOutCellDIm,
                             gmds::TCellID&            AOutCellID);
    
    /*------------------------------------------------------------------------*/
    /** \brief Performs the Heun's algorithm into a triangle to compute the
     *                  out point when we go through this triangle and we arrive
     *                  from an edge. Note that the algorithm is numerically
     *                  corrected to avoid to go back in the face we come from.
     *
     * \param AInEdge   the edge we come from
     * \param AInPnt    the geometric point we start from
     * \param AInVec    the geometric direction to follow
     * \param AOppNode  the node opposite to AInEdge
     * \param AInNode1  the first node of AInEdge
     * \param AInNode2  the second node of AInEdge
     * \param AOutEdge1 another edge of the face we work on
     * \param AOutEdge2 a second another edge of the face we work on
     * \param AOutPnt   the geometric point we go out
     * \param AOutVec   the geometric direction to follow after
     *
     * \return an index indicating which is the cell we go out through:
     *         - 1 for AOppNode
     *         - 2 for AInNode1
     *         - 3 for AInNode2
     *         - 4 for AOutEdge1
     *         - 5 for AOutEdge2
     *
     */
    int heunsComputation(const gmds::Edge&         AInEdge,
                         const gmds::math::Point&  AInPnt,
                         const gmds::math::Vector& AInVec,
                         const gmds::Node&         AOppNode,
                         const gmds::Node&         AInNode1,
                         const gmds::Node&         AInNode2,
                         const gmds::Edge&         AOutEdge1,
                         const gmds::Edge&         AOutEdge2,
                         gmds::math::Point&        AOutPnt,
                         gmds::math::Vector&       AOutVec);
    
    /*------------------------------------------------------------------------*/
    /** \brief Performs the Heun's algorithm into a triangle to compute the
     *                  out point when we go through this triangle and we arrive
     *                  from a node. Note that the algorithm is numerically
     *                  corrected to avoid to go out from the triangle we get
     *                  in.
     *
     * \param AInNode   the node we come from
     * \param AInPnt    the geometric point we start from
     * \param AInVec    the geometric direction to follow
     * \param AOppNode1 the first node of AOppEdge
     * \param AOppNode2 the second node of AOppEdge
     * \param AOppEdge  the edge opposite to AInNode
     * \param AOutPnt   the geometric point we go out
     * \param AOutVec   the geometric direction to follow after
     *
     * \return an index indicating which is the cell we go out through:
     *         - 0 if it does not intersect
     *         - 1 for AOppNode1
     *         - 2 for AOppNode2
     *         - 3 for AOppEdge
     *
     */
    int heunsComputation(const gmds::Node&         AInNode,
                         const gmds::math::Point&  AInPnt,
                         const gmds::math::Vector& AInVec,
                         const gmds::Node&         AOppNode1,
                         const gmds::Node&         AOppNode2,
                         const gmds::Edge&         AOppEdge,
                         gmds::math::Point&        AOutPnt,		
                         gmds::math::Vector&       AOutVec);
    
    /*------------------------------------------------------------------------*/
    /** \brief Let a ray r=[AInPnt, AInVec) interesecting an edge AEdge, this
     *         method returns the point of intersection P between r and AEdge,
     *         and the output vector at P, that respects the underlying cross
     *         field
     *
     * \param AEdge     the edge we intersect
     * \param AInPnt    the geometric point we start from
     * \param AInVec    the geometric direction to follow
     * \param AOutPnt   the geometric point we go out
     * \param AOutVec   the geometric direction to follow after
     *
     * \return true if the ray intersects the edge, false otherwise
     */
    bool computeOutVectorFromRayAndEdge(const gmds::Edge&         AEdge,
                                        const gmds::math::Point&  AInPnt,
                                        const gmds::math::Vector& AInVec,
                                        gmds::math::Point&        AOutPnt,
                                        gmds::math::Vector&       AOutVec);
    
    /*------------------------------------------------------------------------*/
    /** \brief Arriving at node ANode with direction AInVec, it returns the
     *         best fit direction prescribed by the cross field at ANode.
     *
     * \param ANode   the node we consider
     * \param AInVec    the geometric direction we arrive
     * \param AOutVec   the geometric direction to follow after
     *
     * \return true if the ray intersects the edge, false otherwise
     */
    void computeOutVectorAtPoint(const gmds::Node&         ANode,
                                 const gmds::math::Vector& AInVec,
                                 gmds::math::Vector&       AOutVec); 
    
    
    /*------------------------------------------------------------------------*/
    /** \brief Compute the next cell we will got through. It is a face (dim=2),
     *         or an edge (dim=1).
     *
     * \param[IN]  AFromPnt     the geometric point we start from
     * \param[IN]  AFromVec     the geometric direction we go along
     * \param[IN]  AFromCellDim the dim. of the mesh cell we start from (0 or 1)
     * \param[IN]  AFromCellID  the id of the mesh cell we start from
     * \param[OUT] AToCellDim   the dim. of the mesh cell we start from (0 or 1)
     * \param[OUT] AToCellID    the id of the mesh cell we start from
     */
    void findNextCell(const gmds::math::Point&  AFromPnt,
                      const gmds::math::Vector& AFromVec,
                      const int AFromCellDim,
                      const gmds::TCellID AFromCellID,
                      int& AToCellDim,
                      gmds::TCellID& AToCellID);
    
    /*------------------------------------------------------------------------*/
    /** \brief Compute the next cell we will got through. It is a face (dim=2),
     *         or an edge (dim=1).
     *
     * \param[IN]  AFromPnt   the geometric point we start from
     * \param[IN]  AFromVec   the geometric direction we go along
     * \param[IN]  AFromNode  the node we come from
     * \param[OUT] AToCellDim the dim. of the mesh cell we start from (0 or 1)
     * \param[OUT] AToCellID  the id of the mesh cell we start from
     */
    void findNextCell(const gmds::math::Point&  AFromPnt,
                      const gmds::math::Vector& AFromVec,
                      const gmds::Node& AFromNode,
                      int& AToCellDim,
                      gmds::TCellID& AToCellID);
    /*------------------------------------------------------------------------*/
    /** \brief Compute the next cell we will got through. It is a face (dim=2),
     *         or an edge (dim=1).
     *
     * \param[IN]  AFromPnt   the geometric point we start from
     * \param[IN]  AFromVec   the geometric direction we go along
     * \param[IN]  AFromEdge  the edge we come from
     * \param[OUT] AToCellDim the dim. of the mesh cell we start from (0 or 1)
     * \param[OUT] AToCellID  the id of the mesh cell we start from
     */
    void findNextCell(const gmds::math::Point&  AFromPnt,
                      const gmds::math::Vector& AFromVec,
                      const gmds::Edge& AFromEdge,
                      int& AToCellDim,
                      gmds::TCellID& AToCellID);
    
    /*------------------------------------------------------------------------*/
    /** \brief Indicate if the ray starting from AFromNode and following AVec
     *         is aligned with AEdge.
     *
     * \param AVec      a direction modelized by a vector
     * \param AFromNode the node we start from
     * \param AEdge     an edge incident to AFromNode (otherwise returns false)
     *
     * \return a boolean indicating if the AVec is aligned with AEdge
     */
    bool isAlong(const gmds::math::Vector& AVec,
                 const gmds::Node& AFromNode,
                 gmds::Edge& AEdge);
    
    /*------------------------------------------------------------------------*/
    /** \brief Indicate if the ray starting from APnt and following AVec
     *         is going into the face AFace.
     *
     * \param APnt      the starting point
     * \param AVec      a direction modelized by a vector
     * \param AFromNode the node of AFace, that is located at APnt
     * \param AFace     the face we want to check
     *
     * \return a boolean
     */
    bool isGoingInto(const gmds::math::Point& APnt,
                     const gmds::math::Vector& AVec,
                     const gmds::Node& AFromNode,
                     const gmds::Face& AFace);
    /*------------------------------------------------------------------------*/
    /** \brief Indicate if the ray starting from APnt and following AVec
     *         is going into the face AFace.
     *
     * \param APnt      the starting point
     * \param AVec      a direction modelized by a vector
     * \param AFromEdge the edge of AFace, that contains APnt
     * \param AFace     the face we want to check
     *
     * \return a boolean
     */
    bool isGoingInto(const gmds::math::Point& APnt,
                     const gmds::math::Vector& AVec,
                     const gmds::Edge& AFromEdge,
                     const gmds::Face& AFace);
    
    gmds::math::Chart computeChartIn(const gmds::math::Point& APnt,
                                     const gmds::Face& AFace);
    
    void computeFuzzyHeuns(const gmds::math::Point&                 AFromPnt,
                           const gmds::math::Vector3d&              AFromDir,
                           const std::vector<gmds::Face>&           AFaces,
                           const std::vector<gmds::math::Triangle>& ATri,
                           gmds::math::Point&                       AToPnt,
                           gmds::math::Vector3d&                    AToDir,
                           int&                                     AToFaceId);
    void computeFuzzyHeuns(const gmds::math::Point&                 AFromPnt,
                           const gmds::math::Vector3d&              ADirPnt,
                           const std::vector<gmds::Edge>&           AFaces,
                           const std::vector<gmds::math::Segment>&  ATri,
                           gmds::math::Point&                       AToPnt,
                           gmds::math::Vector3d&                    AToDir,
                           int&                                     AToFaceId);
    
    /*------------------------------------------------------------------------*/
    /** \brief Indicates if a point \p APnt projected onto the plane defining
     *         \p ATri belong to \p ATri or not. Boolean \p AOnEdge0, \p AOnEdge1
     *         and \p \p AOnEdge2 indicates if the points lies on the edge o
     *         opposite to node 0, 1 and 3 respectively
     *
     * \param[in]  APnt a point
     * \param[in]  ATri a triangle
     * \param[out] AOnEdge0 true if APnt lies on the edge opposite to the node
     *                      0 of \p ATri
     * \param[out] AOnEdge1 true if APnt lies on the edge opposite to the node
     *                      1 of \p ATri
     * \param[out] AOnEdge2 true if APnt lies on the edge opposite to the node
     *                      2 of \p ATri
     *
     * \return true if \p APnt is in \p ATri, false otherwise
     */
    
    bool isIn(const gmds::math::Point& APnt,
              const gmds::Face& ATri,
              bool& AOnEdge0, bool& AOnEdge1, bool& AOnEdge2);
    
    char orient3d(const gmds::math::Point& AP0,
                  const gmds::math::Point& AP1,
                  const gmds::math::Point& AP2,
                  const gmds::math::Point& AP3);
    

private:
    
    /** Mesh we start from */
    gmds::IGMesh* m_mesh;
    /* Cross field we start from*/
    gmds::Variable<gmds::math::Cross2D>* m_field;
    
    gmds::Variable<gmds::math::AxisAngleRotation>* m_rot_field;
    
};
/*----------------------------------------------------------------------------*/
#endif /* FRAME_TOOLKIT_H_ */
/*----------------------------------------------------------------------------*/
