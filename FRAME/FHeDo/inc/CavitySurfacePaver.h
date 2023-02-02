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
#ifndef SH_CAVITY_SURFACE_PAVER_H_
#define SH_CAVITY_SURFACE_PAVER_H_

/*----------------------------------------------------------------------------*/
// STL File Headers
#include <vector>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
#include <GMDS/IG/IGMesh.h>
#include <GMDS/Math/Triangle.h>
/*----------------------------------------------------------------------------*/
/** \struct CavitySurfacePaver
 *  \brief This structure defines the boundary of a cavity surface
 */
struct EXPORT_GMDS BndPoint{
    enum vertex_type{
        END,
        SIDE,
        CORNER,
        REVERSAL,
        UNDEFINED
    };
    
    BndPoint(const gmds::Node& AN, const gmds::math::Vector3d& AV);
    gmds::Node  node;
    /** Vector at node going into the domain to be meshed*/
    gmds::math::Vector3d inward;
    vertex_type type;
};

/*----------------------------------------------------------------------------*/
/** \struct CavitySurfacePatch
 *  \brief Local definition of a surface as a set of triangles to perform 
 *         geometric projection
 */
struct EXPORT_GMDS CavitySurfacePatch{
    
    /*------------------------------------------------------------------------*/
    /** \brief Default constructor.
     *
     
     * \param AT a set of triangles defining the patch
     * \param AN output normal in each triangle
     */
    CavitySurfacePatch(const std::vector<gmds::math::Triangle>& AT,
                       const std::vector<gmds::math::Vector3d>& AN);

    /*------------------------------------------------------------------------*/
    /** \brief compute the projected point of \p AP onto *this
     *
     * \param[in] AP the point to be projected
     *
     * /return the projected point
     */
    gmds::math::Point project(gmds::math::Point& AP) const;

    gmds::math::Vector3d normal( gmds::math::Point& AP) const;
    void write();

    std::vector<gmds::math::Triangle> triangles;

    std::vector<gmds::math::Vector3d> normals;
};
/*----------------------------------------------------------------------------*/
/** \struct CavitySurfacePaver
 *  \brief This structure defines the boundary of a cavity surface
 */
struct EXPORT_GMDS CavitySurfaceBnd{
    
    /*------------------------------------------------------------------------*/
    /** \brief Default constructor.
     */
    CavitySurfaceBnd(const std::vector<gmds::Node>& ABnd,
                     CavitySurfacePatch* APatch);
    
    
    void setPatch(CavitySurfacePatch* APatch){
        patch=APatch;
    }
    /*------------------------------------------------------------------------*/
    /** \brief Returns the element before the one of index \p AIndex-1. It takes
     *         the cyclic behaviour of the boundary into account
     *
     * \param[in] AIndex the index of the current point
     *
     * \return the point before the one with index \p AIndex
     */
    BndPoint& prev(const int AIndex);
    
    /*------------------------------------------------------------------------*/
    /** \brief Returns the element after the one of index \p AIndex-1. It takes
     *         the cyclic behaviour of the boundary into account
     *
     * \param[in] AIndex the index of the current point
     *
     * \return the point after the one with index \p AIndex
     */
    BndPoint& next(const int AIndex);
    
    /*------------------------------------------------------------------------*/
    /** \brief Returns the element of index \p AI
     *
     * \param[in] AI the index of the point we want to access to
     *
     * \return the point of index \p AIndex
     */
    BndPoint& operator[](const int AIndex);
    
    /*------------------------------------------------------------------------*/
    /** \brief Flag each node as end, side, corner or reversal
     */
    void initTypes();
    
    
    /*------------------------------------------------------------------------*/
    /** \brief Compute the type of each point between indices \p AI and  /p AJ
     *         both included
     *
     * \param[in] AI index of the first point to update
     * \param[in] AJ index of the last  point to update
     */
    void updateTypes(const int& AI, const int& AJ);
    
    int size() const;
    
    void clear();
    /*------------------------------------------------------------------------*/
    /** \brief compute the angle between \p AV1 and  \p AV2
     *
     * /return a value between -180 and 180 degree
     */
    static double angle(const gmds::math::Vector3d& AV1,
                        const gmds::math::Vector3d& AV2,
                        const gmds::math::Vector3d& ANormal);
    /*------------------------------------------------------------------------*/
    /** \brief compute the angle in point of index AI
     */
    double angle(const int AI);
    /*------------------------------------------------------------------------*/
    /** \brief We replace all the points in [\p AI, \p AJ] by \p ANodes. Points
     *         \p AI and \p AJ are included in the replaced points. This 
     *         modification leads to an angle and type computation for some 
     *         boundary points.
     *
     * \param[in] AI   index of the first point to replace
     * \param[in] AJ   index of the last  point to replace
     * \param[in] AN   nodes to be inserted in replacement
     * \param[in] ADir inward vectors computed
     */

    void replace(const int AI, const int AJ, std::vector<gmds::Node>& AN,
                 std::vector<gmds::math::Vector3d>& ADir);
    
    /** the boundary of a cavity must form a loop, so  a closed ordered list
     *  of nodes */
    std::vector<BndPoint> points;
    
    CavitySurfacePatch* patch;
};

/*----------------------------------------------------------------------------*/
/** \class CavitySurfacePaver
 *  \brief This class provides an algorithm to pave a cavity boundary on a
 *         surface. This class is a simpler paver in the meaning that it only
 *         fills surface hole at the end of the hex-dominante algorithm
 */
class EXPORT_GMDS CavitySurfacePaver{
    
public:
    /*------------------------------------------------------------------------*/
    /** \brief Constructor.
     *
     * \param[in] AMesh the mesh we want to add the element in
     * \param[in] ABnd an ordered list of nodes such that the surface to pave is
     *            on the left of each edge; Nodes of \p ABnd must be in \p AMesh
     
     * \param[in] ATri discrete definition of the geometrical patch we have to 
     *                 mesh
     */
    CavitySurfacePaver(gmds::IGMesh*                            AMesh,
                       const std::vector<gmds::Node>&           ABnd,
                       const std::vector<gmds::math::Triangle>& ATri,
                       const std::vector<gmds::math::Vector3d>& AN);
    
    /*------------------------------------------------------------------------*/
    /** \brief Function to be called for generating the mesh
     */
    void execute();
    
    void write();
protected:
    
    /*------------------------------------------------------------------------*/
    /** \brief Finalize the mesh when the boundary has 6 nodes or less
     */
    void finalize();
    
    /*------------------------------------------------------------------------*/
    /** \brief Finalize the mesh when the boundary has 5 nodes exactly
     */
    void finalize5();
    /*------------------------------------------------------------------------*/
    /** \brief Finalize the mesh when the boundary has 6 nodes exactly
     */
    void finalize6();
    
    /*------------------------------------------------------------------------*/
    /** \brief Check wich 6 patterns we encounter
     */
    bool checkPattern6(const BndPoint::vertex_type APattern[],
                       int& AStartingIndex);
    
    /*------------------------------------------------------------------------*/
    /** \brief Check if we have a 4 END corner-shape. If yes, it returns the
     *         indices of the 4 END corners in  \p AIndex
     *
     * \param[out] AIndex the index of the the 4 corners if found
     *
     * \return true if we encounter a 4-sided configuration, false otherwise
     *
     */
    bool check4Sided(int (&AIndex)[4]);
    /*------------------------------------------------------------------------*/
    /** \brief if we have a grid structure, that is a quad patch with the
     *         same number of edges on opposite patch edges, this algorithm
     *         directly mesh the patch.
     *
     * \return true if it is done, false otherwise
     *
     */
    bool meshGridPatch();
    
    void mesh1RowPatch(const int& AI, const int& AJ);
    
    void meshTFIPatch(const int AIndex[4]);
    void findSmallestRowToInsert(int& AI, int& AJ);

    std::vector<gmds::Node> insertRow(const int& AI, const int& AJ);

    /*------------------------------------------------------------------------*/
    /** \brief When a row is inserted, free nodes, that is those inserted by
     *         the paver are smoothed by a simple Laplacian algorithm
     */
    void smooth(std::vector<gmds::Node>& ANodes);
    
    /*------------------------------------------------------------------------*/
    /** \brief compute the quad quality. It is included in [0,1] and based on
     *         the angle. Quality angle is based onto the difference to 90
     *         degrees, and this function returns the worst corner angle.
     */
    double quality(const gmds::Node& AN1, const gmds::Node& AN2,
                   const gmds::Node& AN3, const gmds::Node& AN4);
    /*------------------------------------------------------------------------*/
    /** \brief add a quad and performs the right connectivity updates
     */
    gmds::Face addQuad(gmds::Node& AN1, gmds::Node& AN2,
                       gmds::Node& AN3, gmds::Node& AN4);
    
    /*------------------------------------------------------------------------*/
    /** \brief add a quad and performs the right connectivity updates
     */
    gmds::Face addQuad(gmds::TCellID& AN1, gmds::TCellID& AN2,
                       gmds::TCellID& AN3, gmds::TCellID& AN4);
    /*------------------------------------------------------------------------*/
    /** \brief add a triangle and performs the right connectivity updates
     */
    gmds::Face addTriangle(gmds::Node& AN1,
                           gmds::Node& AN2,
                           gmds::Node& AN3);
    /*------------------------------------------------------------------------*/
    /** \brief add a quad and performs the right connectivity updates
     */
    gmds::Node addNode(gmds::math::Point& AP);
    
    /*------------------------------------------------------------------------*/
    /** \brief Create a row side point from a SIDE boundary point of index \p AI
     */
    gmds::math::Point createRowSideNode(const int& AI);
    /*------------------------------------------------------------------------*/
    /** \brief Create 3 points from a CORNER boundary point of index \p AI
     */
    std::vector<gmds::math::Point> createRowCornerNode(const int& AI);
    /*------------------------------------------------------------------------*/
    /** \brief Create 5 points from a REVERSAL boundary point of index \p AI
     */
    std::vector<gmds::math::Point> createRowReversalNode(const int& AI);
 
    gmds::IGMesh* m_mesh;
    
    CavitySurfacePatch m_patch;

    CavitySurfaceBnd m_bnd;
    /*We keep the list of the initial nodes that make the boundary loop*/
    std::vector<gmds::Node> m_initial_nodes;
    
    /*All inserted points are kept in mind*/
    std::vector<gmds::Node> m_free_nodes;
    
    enum node_type{
        FIXED,
        FLOATING,
        FRONT,
    };
    /** for each node handled by the advancing front algorithm, we know
     * if it is fixed (initial boundary node), on the front (ie will be
     *  used to insert a new face, or floating, that is no more in the 
     *  front */
    std::map<gmds::TCellID, node_type> m_node_type;
    std::vector<gmds::Face> m_faces;

};

/*----------------------------------------------------------------------------*/
#endif /* SH_CAVITY_SURFACE_PAVER_H_ */
/*----------------------------------------------------------------------------*/
