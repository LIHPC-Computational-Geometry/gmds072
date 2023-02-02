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
#ifndef SH_TRIANGULAR_SURFACE_MANIPULATOR_H_
#define SH_TRIANGULAR_SURFACE_MANIPULATOR_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <map>

/*----------------------------------------------------------------------------*/
// Predicates_psm File Headers
#include <Predicates_psm.h>

/*----------------------------------------------------------------------------*/
// GMDS File Headers
#include <GMDS/IG/IGMesh.h>
/*----------------------------------------------------------------------------*/
// Tetgen File Headers
#include <tetgen.h>
/*----------------------------------------------------------------------------*/
/** \class  TriangularSurfaceManipulator
 *  \brief  This class gathers all the algorithms impliying some triangular
 *          surfacic mesh modifications in the HexDom mesher.
 */
class EXPORT_GMDS TriangularSurfaceManipulator{
    
public:
    /*------------------------------------------------------------------------*/
    /** \brief Constructor.
     *
     * \param[in] AMesh   the triangular mesh that will be modified
     */
    TriangularSurfaceManipulator(gmds::IGMesh* AMesh);
    /*------------------------------------------------------------------------*/
    /** \brief Destructor.
     */
    virtual ~TriangularSurfaceManipulator();
    
    /*------------------------------------------------------------------------*/
    /** \brief Try to insert  the loop \p ALoop into the background mesh. A 
     *         loop is defined by a set of points and the implicit segments 
     *         that connect a point to the next one (last and first are 
     *         connected too).
     *         For each point classified on a curve, we give its curve id in
     *         \p ACurve. A point classified on a geometric vertex can be 
     *         assigned to several curves.
     *
     *         As a result, it returns the ID of all the faces enclosed by the
     *         loop in \p AEnclosedFaces and all the nodes forming the
     *         loop \p ALoopNodes (in an ordered fashion
     *
     * \param[in] ALoop a cyclic list of points
     */
    
    void insertLoop(std::vector<gmds::math::Point>& ALoop,
                    std::vector<int> & AClassification,
                    std::vector<int> & AID,
                    std::vector<gmds::TCellID>& AEnclosedFaces,
                    std::vector<gmds::TCellID>& ALoopNodes,
                    std::vector<bool>& ALoopNodesFromPnt);
    
    std::vector<gmds::Node>
    insertFreePoints(std::vector<gmds::math::Point>& APntToInsert,
                     std::vector<gmds::TCellID>& AInFaces);
    
    /*------------------------------------------------------------------------*/
    /** \brief Extract the faces enclosed by the loop made of the ordered list
     *         of nodes \p ALoopNodes. When this function is called, all the 
     *         nodes are linked by edges in the mesh.
     *
     * \param[in] ALoopNodes a cyclic list of nodes
     *
     * \return the ids of the enclosed faces.
     */

    std::vector<gmds::TCellID>
    extractEnclosedFaces(const std::vector<gmds::Node>& ALoopNodes);
    std::vector<gmds::TCellID>
    extractEnclosedFaces(const std::vector<gmds::TCellID>& ALoopNodes);
    
    /*------------------------------------------------------------------------*/
    /** \brief Find the face containing the point \p APnt. This inclusion
     *         is approximated due to the fact that APnt does not necessary
     *         lies upon the surface mesh
     *
     *         This procedure is not optimized at all
     *
     * \param[in]  APnt a point we look the face containing it
     * \param[out] AOutNode the node created or found in the mesh
     *
     * \return the face that contains APnt.
     */
    std::vector<gmds::Face> findFaces(const gmds::math::Point& APnt,
                                      gmds::Node& AOutNode);
    
    std::vector<gmds::Face> findFaces(const gmds::math::Point& APnt,
                                      const std::vector<gmds::Face>& AFaces,
                                      gmds::Node& AOutNode);
    /*------------------------------------------------------------------------*/
    /** \brief Walk from \p AFrom to \p ATo starting from the face \p AF which
     *         must contain \p AFrom. Point generated are put into \p APoints
     *         as well as the faces containing each of those points which are
     *         in \p AData. Note that datas relative to \p AFrom must be in 
     *         \APoints and \p AData, while those relative to ATo will be added
     *         by this process.
     *
     * \p ATo can be modified if he lies very closely to a vertex or edge of the
     *        triangular mesh.
     *
     * \param[in]     AFrom      the point we start from
     * \param[in]     ATo        the point we go to
     * \param[in]     AFromFaces Faces containing AFrom
     * \param[in/out] APoints    the complete set of generated points, in fact 
     *                           nodes are generated or listed
     * \param[in/out] AData      the faces containing each generated point
     * \param[in]     AIsLast    when we insert the last segment, a specific behaviour
     *                           is done for the last point (same as the first one)
     *
     * \return the face that contains APnt.
     */
    void walk(const gmds::math::Point& AFrom,
              gmds::math::Point& ATo,
              const std::vector<gmds::Face>& AFromFaces,
              std::vector<gmds::Node>& ANewNodes,
              std::vector<bool>& ANewFromPnt,
              std::map<gmds::TCellID, std::vector<int> >& AData,
              std::vector<gmds::Face>& AOutFaces,
              const bool AIsLast);
    
    void walkOnCurve(const gmds::math::Point& AFrom,
                     gmds::math::Point& ATo,
                     const int AFromClassif,
                     const int AToClassif,
                     const int AFromGeomID,
                     const int AToGeomID,
                     const std::vector<gmds::Face>& AFromFaces,
                     std::vector<gmds::Node>& ANewNodes,
                     std::vector<bool>& ANewFromPnt,

                     std::map<gmds::TCellID, std::vector<int> >& AData,
                     std::vector<gmds::Face>& AOutFaces,
                     const bool AIsLast);

    std::vector<gmds::Node> getNodes(const gmds::Node& AFrom,
                                     const int ASurfColor0,
                                     const int ASurfColor1);

    void findOutFace(const gmds::Node& AFrom,
                     const gmds::math::Point& ATo,
                     const gmds::Face& ACurFace,
                     gmds::Face& ANextFace);
    /*------------------------------------------------------------------------*/
    /** \brief Split the face \p AF by inserting into it all the nodes of
     *         \p AN
     *
     * \param[in] AF a face to split
     * \param[in] AN nodes to be inserted
     *
     * \return true if the face is effectively splitted, false otherwise
     */
    bool split(gmds::Face& AF, std::vector<gmds::Node>& AN);
    
    void triangulate(gmds::Face& AFace,
                     std::vector<gmds::math::Point>& APoints);
    
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
    
    bool isIn(const gmds::math::Point& APnt,
              const gmds::Face& ATri,
              std::vector<gmds::TCoord>& ACoord);
    
    
    void simplify(std::vector<gmds::Node>& AToKeep);
    /*------------------------------------------------------------------------*/
    /** \brief Return the faces sharing edge [\p AN1, \p AN2].
     * \param[in]  AN1 a first node
     * \param[in]  AN2 a second node
     *
     * \return all the faces sharing edge [\p AN1, \p AN2]
     */
    std::vector<gmds::Face> getFaces(const gmds::Node& AN1,
                                     const gmds::Node& AN2);
    

   
    char orient3d(const gmds::math::Point& AP0,
                  const gmds::math::Point& AP1,
                  const gmds::math::Point& AP2,
                  const gmds::math::Point& AP3);
    
    std::vector<gmds::Face>
    getFaces(const gmds::Face& ATri,
             const bool AOnEdge0,
             const bool AOnEdge1,
             const bool AOnEdge2);
    
    void writePoints(std::vector<gmds::math::Point>& ALoop);
    
    void buildEdge(gmds::Node& ANI, gmds::Node& ANJ);
    /*------------------------------------------------------------------------*/
    /** \brief Give the face adjacent to \p AN such that ray coming from \p AN
     *         and going towards \p ATo intersect it. Some coplanarity issues
     *         can occur if the point \p ATo is far from being coplanar with
     *         the faces of \p AN.
     *
     * \param[in]  AN a starting node
     * \param[in]  ATo the point we go towards
     *
     * \return A face adjacent to \p AN
     */

    gmds::Face getOutFace(const gmds::Node& ANI, const gmds::math::Point ATo);
    
    
protected:
    gmds::IGMesh* m_mesh;
    gmds::Variable<gmds::TCellID>* m_cavity;
    /*---------------------------------------------------------------------------*/
    /** Local node numbering so that the triangle formed with:
     *   - local node  i
     *   - m_local_tri_node2edge[i][0]
     *   - m_local_tri_node2edge[i][1]
     *   has the same orientation as the original triangle for any vertex i.
     */
    static const int m_local_tri_node2edge[3][2];
    
};

/*----------------------------------------------------------------------------*/
#endif /* SH_TRIANGULAR_SURFACE_MANIPULATOR_H_ */
/*----------------------------------------------------------------------------*/
