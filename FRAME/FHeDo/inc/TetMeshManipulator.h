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
#ifndef SH_TETMESH_MANIPULATOR_H_
#define SH_TETMESH_MANIPULATOR_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <map>

/*---------------------------------------------------------------------------*/
// Predicates_psm File Headers
#include <Predicates_psm.h>

/*----------------------------------------------------------------------------*/
// GMDS File Headers
#include <GMDS/IG/IGMesh.h>
/*----------------------------------------------------------------------------*/
// Tetgen File Headers
#include <tetgen.h>
/*----------------------------------------------------------------------------*/
namespace fhedo{
/*----------------------------------------------------------------------------*/
/** \class  TetMeshManipulator
 *  \brief  This class gathers all the algorithms impliying some tet mesh
 *          modifications in the HexDom mesher. In practice, tetgen is used
 *          to perform those modifications
 *
 *  The findTet() function is based on the following two references.
 *  The first one randomizes the choice of the next tetrahedron.
 *  The second one uses an inexact locate() function to initialize
 *  the exact one (it is called "structural filtering"). The first
 *  idea is used in both CGAL and tetgen, and the second one is used
 *  in CGAL.
 *  - Walking in a triangulation, O Devillers, S Pion, M Teillaud
 *   17th Annual Symposium on Computational geometry, 106-114
 *  - Stefan Funke , Kurt Mehlhorn and Stefan Naher, "Structural filtering,
 *  a paradigm for efficient and exact geometric programs", 1999
 */
class EXPORT_GMDS TetMeshManipulator{
    
public:
    /*------------------------------------------------------------------------*/
    /** \brief Constructor.
     *
     * \param[in] AMesh   background tetrahedral mesh \p APnts was built
     *                    from
     */
    TetMeshManipulator(gmds::IGMesh* AMesh);
    /*------------------------------------------------------------------------*/
    /** \brief Destructor.
     */
    virtual ~TetMeshManipulator();
    
    /*------------------------------------------------------------------------*/
    /** \brief Try to insert  the loop \p ALoop into the background mesh. A 
     *         loop is defined by a set of points and the implicit segments 
     *         that connect a point to the next one (last and first are 
     *         connected too).
     *
     * \param[in] ALoop a cyclic list of points
     */

    void insertLoopOnSurface(gmds::IGMesh* ATriMesh,
                             std::vector<gmds::math::Point>& ALoop);
    
    void insertPoints(gmds::IGMesh* AInMesh,
                      std::vector<gmds::math::Point>& APoints,
                      const std::vector<int>& AClassification);
    /*------------------------------------------------------------------------*/
    /** \brief Try to insert \p APnt into the tet \p ATet
     *
     * \param[in] APnt the point to insert
     * \param[in] ATet the tetrahedral element we start the process from
     * \param[in] AC       the presumed classificatio of \p APnt
     * \param[out] ANewTet  a tet of the new cavity
     *
     * \true if the insertion succeeded, false otherwise
     */

    bool insert(gmds::math::Point &APnt, gmds::Region& ATet, const int AC,
                gmds::Region& ANewTet);

    /*------------------------------------------------------------------------*/
    /** \brief Insert new points \p APoints into the tetrahedral mesh
     *         \p AMesh to get the output mesh \p AOutMesh
     *
     * \param[in] AInMesh   the tet mesh we want to work on
     * \param[in] APoints   the points to be inserted
     * \param[out] AOutMesh the resulting mesh
     */
    void insertPointsTetgen(gmds::IGMesh* AInMesh,
                            std::vector<gmds::math::Point>& APoints,
                            gmds::IGMesh* AOutMesh);
    /*------------------------------------------------------------------------*/
    /** \brief Insert new points \p APoints into the tetrahedral mesh
     *         \p AMesh to get the output mesh \p AOutMesh
     *
     * \param[in] AInMesh   the tet mesh we want to work on
     * \param[in] APoints   the points to be inserted
     * \param[out] AOutMesh the resulting mesh
     */
    void tetrahedralizeTetgen(gmds::IGMesh* ABndMesh,
                              gmds::IGMesh* ATetMesh);
    /*------------------------------------------------------------------------*/
    /** \brief Check if m_mesh is Delaunay admissible. Warning this function
     *         is not optmized.
     */

    bool isDelaunay();
    
protected:
    
    
    /*------------------------------------------------------------------------*/
    struct facet{
        gmds::TCellID n[3];
        
        facet(gmds::TCellID AI = gmds::NullID,
              gmds::TCellID AJ = gmds::NullID,
              gmds::TCellID AK = gmds::NullID) {
            n[0]=AI;
            n[1]=AJ;
            n[2]=AK;
        }
        
        gmds::TCellID  operator() (const int AI) const {return n[AI];}
        gmds::TCellID& operator() (const int AI) {return n[AI];}
    };
    
    /*------------------------------------------------------------------------*/
    /** \brief Finds the tetrahedron that contains a point.
     * \details If the point is on a face, edge or vertex,
     *  the function returns one of the tetrahedra incident
     *  to that face, edge or vertex.
     * \param[in] p a pointer to the coordinates of the point
     * \param[out] orient a pointer to an array of four Sign%s
     *  or nil. If non-nil, returns the orientation with respect
     *  to the four facets of the tetrahedron that contains \p p.
     * \return the index of a tetrahedron that contains \p p.
     *  If the point is outside the convex hull of
     *  the inserted so-far points, then the returned tetrahedron
     *  is a virtual one (first vertex is the "vertex at infinity"
     *  of index -1) or NO_TETRAHEDRON if the virtual tetrahedra
     *  were previously removed.
     */
    gmds::Region findTet(gmds::math::Point& APnt,
                         const gmds::Region& ASeed,
                         const int AC);
    
    /** Returns true if the point \p AP is contained in region 
     *  \p AR knowing that \p AP is classified as \p AC.
     */
    bool contains(gmds::math::Point& AP,
                         const gmds::Region& AR,
                         const int AC);
    
    /*------------------------------------------------------------------------*/
    /** \brief Init the cavity of a \p AP knowing it is classified on a curve 
     *         (\p AC=1), a surface (\p AC=2) or inside the
     *         volume (\p AC=3) and that it is "in" \p AR.
     *
     *         Points on vertices are not put since they already exist in the 
     *         mesh so AC=0 cannot occur.
     *
     * \param[in]  AP the point we start from
     * \param[in]  AC the point classification
     * \param[in]  AR the region that contains AP
     * 
     * \return the initial cavity before expansion
     */
    std::vector<gmds::Region> initCavity(gmds::math::Point& APnt,
                                         const int AC,
                                         const gmds::Region& AR);
    
    std::vector<gmds::Region> initFuzzyCavity(gmds::math::Point& APnt,
                                         const int AC,
                                         const gmds::Region& AR);
    
    /*------------------------------------------------------------------------*/
    /** \brief Starting from \p AR and three of its node \p AN0, \p AN1 and
     *         \p AN2 it returns the tetrahedral element sharing \p ANodes with
     *         \p AR. If there is no tet sharing those nodes with  \p AR, the
     *         NullID region is returned (boundary case).
     *
     * \param[in] AR a tetrahedral region
     * \param[in] AN0 a first node of AR
     * \param[in] AN1 a second node of AR
     * \param[in] AN2 a third node of AR
     *
     * \return the region sharing \p AN0, \p AN1 and \p AN2 with \p AR
     */
    
    gmds::Region getOppositeRegion(const gmds::Region& AR,
                                   const gmds::Node&   AN0,
                                   const gmds::Node&   AN1,
                                   const gmds::Node&   AN2);
    /*------------------------------------------------------------------------*/
    /** \brief Reorient each tet of the mesh in order to get a local numbering
     *         such that facet normal are going outside
     */
    void correctTetOrientation();
    /*------------------------------------------------------------------------*/
    /** \brief Build outside tets which are made of a boundary facet +
     *         an infinity node (denoted by gmds::InfinityID
     */
    void buildOutTets();
    
    /*------------------------------------------------------------------------*/
    /** \brief Delete all the tet having at least one infinity node (denoted
     *         by gmds::InfinityID)
     */
    void deleteOutTets();
    
    
    /*------------------------------------------------------------------------*/
    /** \brief Returns if \p ATet is a real tetrahedron or a ghost one
     *         (connected to an infinityPoint)
     */
    bool isInfinityTet(const gmds::Region& ATet) const;
    /*------------------------------------------------------------------------*/
    /** \brief Returns the real facet of \p ATet assuming \p ATet is a ghost
     *         tetrahedral element. The obtained facet is well-oriented.
     */
    facet realFacet(const gmds::Region& ATet);
    
    /*------------------------------------------------------------------------*/
    /** \brief Compute the cavity made of tetrahedra that are in conflict
     *  with \p APnt.
     *
     * \param[in]     APnt the point to be inserted
     * \param[in/out] ACav the set of regions defining the cavity of APnt
     * \param[in]     ACavBnd elements of a \p ACav whose one facet is in the
     *                boundary of ACav
     * \param[in]     ACavID cavity id
     */
    void expandCavity(const gmds::math::Point &APnt,
                      std::vector<gmds::Region>& ACav,
                      std::vector<gmds::Region>& ACavBnd,
                      const int ACavID);
    void expandCavityND(const gmds::math::Point &APnt,
                        const int AClass,
                        std::vector<gmds::Region>& ACav,
                        std::vector<gmds::Region>& ACavBnd,
                        const int ACavID);
    
    bool isInConflict(const gmds::math::Point& APnt,
                      gmds::Region& AR);
    /*------------------------------------------------------------------------*/
    /** \brief The cavity computation relies on floating point arithmetic. This
     *         procedure is here to correct the cavity (by expanding it) in case
     *         it is not star-shaped from \p APnt.
     *
     * \param[in]     APnt the point to be inserted
     * \param[in/out] ACav the set of regions defining the cavity of APnt
     */
    void checkAndCorrectCavity(const gmds::math::Point &APnt,
                       std::vector<gmds::Region>& ACav);
        gmds::Region seed();
    
    int random(const int& AMax) const;
    bool isIn(const gmds::math::Point APnt,
              const gmds::Region& ATet,
              gmds::math::Vector4d& ACoord) const;

    
    bool stellateCavity(const gmds::math::Point& APnt,
                        std::vector<gmds::Region>& ACavBnd,
                        const int ACavID,
                        gmds::Region& ANewTet);
    
    
    double quality(const gmds::math::Point& AN1, const gmds::math::Point& AN2,
                   const gmds::math::Point& AN3, const gmds::math::Point& AN4);
    
    /*------------------------------------------------------------------------*/
    /** \brief Remove the edge [\p AN1, \p AN2] by merging \p AN1 onto \p AN2
     *
     * \param[in]  AN1 a first node mesh
     * \param[out] AN2 a second node mesh
     */
    void removeNodeViaEdgeContraction(gmds::Node& AN1, gmds::Node& AN2);
    
    std::vector<gmds::Region> getRegions(const gmds::Node& AN1,
                                         const gmds::Node& AN2);
    std::vector<gmds::Region> getRegions(const gmds::Node& AN1,
                                         const gmds::Node& AN2,
                                         const gmds::Node& AN3);
    /*------------------------------------------------------------------------*/
    /** \brief Get the regions sharing a face with \p AR
     * \param[in]  AR a region
     */    std::vector<gmds::Region> getRegions(const gmds::Region& AR);
    
    /*------------------------------------------------------------------------*/
    /** \brief load the GMDS Mesh we want to work on.
     *
     * \param[in] AMesh the tet mesh we want to work on
     */
    void loadTetMesh(gmds::IGMesh* AMesh,
                     tetgenio& ATetMesh,
                     std::map<gmds::TCellID,int>& AGMDS2Tet,
                     std::map<int,gmds::TCellID>& ATet2GMDS);
    
    
    /*------------------------------------------------------------------------*/
    /** \brief load the GMDS Mesh we want to work on.
     *
     * \param[in] AMesh the triangular surfacic mesh we want to work on
     */
    void loadTriangularMesh(gmds::IGMesh* AMesh,
                            tetgenio& ATetMesh,
                            std::map<gmds::TCellID,int>& AGMDS2Tet,
                            std::map<int,gmds::TCellID>& ATet2GMDS);
    
    
    /*------------------------------------------------------------------------*/
    /** \brief load a set of points into a tetgenio structure
     *
     * \param[in] AMesh the tet mesh we want to work on
     */
    void loadNodes(std::vector<gmds::math::Point>& APnts,
                   tetgenio& ATetPnts);
    
    bool isIn(const gmds::math::Point& AP,
              const gmds::Region& AR,
              std::vector<char>& AS) const;

    
    /*------------------------------------------------------------------------*/
    /** \brief Check if \p AP1 and \p AP2 are the same point
     *
     * \param[in] AP1 first point
     * \param[in] AP2 second point
     *
     * \return true if \p AP1 and \p AP2 have  same floating-point coordinate,
     *         false otherwise
     */
    bool same(gmds::math::Point& AP1, gmds::math::Point& AP2);
    
    /*------------------------------------------------------------------------*/
    /** \brief Check 3D point colinearity using PCK Predicates
     *
     * \param[in] AP1 first point
     * \param[in] AP2 second point
     * \param[in] AP3 third point
     *
     * \return true if \p AP1, \p AP2 and \p AP3 are colinear,  false otherwise
     */
    bool colinear(gmds::math::Point& AP1,
                  gmds::math::Point& AP2,
                  gmds::math::Point& AP3);
    char orient3d(const gmds::math::Point& AP0,
                  const gmds::math::Point& AP1,
                  const gmds::math::Point& AP2,
                  const gmds::math::Point& AP3);
    
    bool computeProj(const gmds::math::Point& AP,
                     const gmds::math::Point& AT1,
                     const gmds::math::Point& AT2,
                     const gmds::math::Point& AT3,
                     gmds::math::Point& ANewP);
    
    char orient3d(const gmds::math::Point& AP,
                       const gmds::Node& AN1,
                       const gmds::Node& AN2,
                       const gmds::Node& AN3);
    
    void writeTetMesh();
    void writeCavity(std::vector<gmds::Region>& ACav,
                     std::vector<gmds::Region>& ACavBnd);

protected:
    gmds::IGMesh* m_mesh;
    gmds::Variable<gmds::TCellID>* m_cavity;
    static const int m_local_tet_node2facet[4][3];
    
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* SH_TETMESH_MANIPULATOR_H_ */
/*----------------------------------------------------------------------------*/
