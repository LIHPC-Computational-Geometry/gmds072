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
#ifndef SH_HEX_DOM_MESH_GENERATOR_H_
#define SH_HEX_DOM_MESH_GENERATOR_H_

/*----------------------------------------------------------------------------*/
// STL File Headers
#include <vector>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
#include <GMDS/IG/IGMesh.h>
#include <GMDS/Math/Chart.h>
#include <GMDS/Math/AxisAngleRotation.h>
#include <GMDS/Utils/OrientedGraph.h>
/*----------------------------------------------------------------------------*/
#include "Params.h"
/*----------------------------------------------------------------------------*/
namespace fhedo{
/*----------------------------------------------------------------------------*/
/** \class  HDMeshGenerator
 *  \brief  This class provides an algorithm to generate hex-dominant meshes
 *          from a set of points obtained via the PointGenerator class.
 */
class EXPORT_GMDS HDMeshGenerator{
    
public:
    /*------------------------------------------------------------------------*/
    /** \brief Constructor.
     *
     * \param[in] AMesh   background tetrahedral mesh \p APnts was built
     *                    from
     * \param[in] ARot    frame field defined in the nodes of \p AMesh
     * \param[in] AMarkFaceOnSurface mark for faces classified on surfaces
     * \param[in] AMarkEdgeOnCurve   mark for edges classified on curves
     * \param[in] AMarkEdgeOnSurface mark for edges classified on surfaces
     * \param[in] AMarkNodeOnPoint   mark for nodes classified on points
     * \param[in] AMarkNodeOnCurve   mark for nodes classified on curves
     * \param[in] AMarkNodeOnSurface mark for nodes classified on surfaces
     * \param[in] APnts   the points we build the hex-dom mesh from
     * \param[in] ACharts the charts associated to each point in \p APnts
     * \param[in] ATypes  a flag indicating if the corresponding point was
     *                    extracted from a stable tet (0), a PGP sing tet (1)
     *                    or a FF sing tet(2)
     * \param[in] AClass  a flag indicating if the corresponding point is
     *                    classified onto a point (0), a curve (1), a surf(2)
     *                    or a volume (3)
     * \param[in] ACurv   a flag indicating the curve number it is classified
     *                    with 0 value if it is in the volume, on a surface or a
     *                    point
     
     * \param[in] ASurf   a flag indicating the surface number it is classified
     *                    with 0 value if it is in the volume, on a curve or a
     *                    point
     * \param[in] ANormal normal for each boundary node
     */
    HDMeshGenerator(gmds::IGMesh* AMesh,
                    const ParamsGlobal& AParamGL,
                    const ParamsMark& AMarks,
                    const ParamsHexDom& AParamHD,
                    const std::vector<gmds::math::Point>& APnts,
                    const std::vector<gmds::math::Chart>& ACharts,
                    const std::vector<gmds::Cell::Data>&  AData,
                    const std::vector<int>&               ATypes,
                    const std::vector<int>&               AClass,
                    const std::vector<int>&               ACurv,
                    const std::vector<int>&               ASurf,
                    const std::map<int, gmds::math::Vector3d>& ANormal,
                    const double ASpacing);
    
    /*------------------------------------------------------------------------*/
    /** \brief Function to be called for generating the mesh
     */
    void execute();
    
    /*------------------------------------------------------------------------*/
    /** \brief Access to the generated mesh
     * 
     * \return the points generated by the algorithm
     */
    const gmds::IGMesh& mesh() const {
        return m_hexdom;
    }

protected:
    
    enum EPointClassification{
        ON_VERTEX =0,
        ON_CURVE  =1,
        ON_SURFACE=2,
        IN_VOLUME =3
    };
    
    enum EPointType{
        REGULAR    =0,
        PARAM_SING =1,
        FRAME_SING =2
    };
    
    /*------------------------------------------------------------------------*/
    /** \struct HexCorner
     * \brief Local Structure to store hex corner data for each input point. A
     *        corner data is made of the index of the current point (p), the
     *        index of the 3 points helping to form the corner (adj). Each 
     *        "edge" [p,adj[i]] is geometrically defined by vec[i].
     *
     *        The free field indicates if the corner is already used for a hex.
     *        The index field give the position of *this in the list of corners
     *        defined in point this->p
     *
     */
    struct HexCorner {
        int p;
        int index;
        bool free;
        int adj[3];
        gmds::math::Vector3d vec[3];
    };
    /*------------------------------------------------------------------------*/
    /** \struct OrientedEdge
     * \brief Local Structure to store an oriented edge made of the index of 
     *        the first and second points and the chart data used in the first
     *        point to select this OrientedEdge. A chart data (-1,-1) means
     *        that the selected end point finally do not correspond to a chart
     *        alignment (will happen for boundary points for instance)
     */
    struct OrientedEdge {
        int first;
        int second;
        int axis;
        int dir;
        OrientedEdge(const int AF=-1, const int AS=-1,
                     const int AAxis=-1, const int ADir=-1):
        first(AF),second(AS),axis(AAxis), dir(ADir){;}
        
        bool operator==(const OrientedEdge& AE)const{
            return (first==AE.first) && (second==AE.second);
        }
        void setInvalid(){
            first=-1; second=-1;
        }
        bool isInvalid(){
            return (first==-1 && second==-1);
            
        }
        
        bool isInverse(OrientedEdge& AE){
            return (AE.first==second && AE.second==first);
        }
        bool isEqual(OrientedEdge& AE){
            return (AE.first==first && AE.second==second);
        }
        
        /*--------------------------------------------------------------------*/
        /** \brief Weak equality (equal or inv) between  *this and \p AE
         *
         * \param[in] AE an oriented edge bo be compared with
         *
         * \return true if \p AE is the same or the inverse edge of *this, 
         *         false otherwise.
         */

        bool operator==(OrientedEdge& AE){
            return isEqual(AE) || isInverse(AE);
        }
        
    };
    
    /*------------------------------------------------------------------------*/
    /** \brief Build for each point a list of close points to be compared with
     *         during edge creation.
     */
    void createDistanceFilter();
    
    
    /*------------------------------------------------------------------------*/
    /** \brief Associate each point to a mesh cell
     */
    void computeMeshAssociation();
    

    /*------------------------------------------------------------------------*/
    /** \brief Build oriented edges from each point to its neighbor points in a
     *         local point of view. \p AEdges[i] contains all the oriented edges
     *         coming from point i.
     *
     * \param[out] AEdges the list of oriented edges that will be created.
     */
    void buildOrientedEdges(std::vector<std::vector<OrientedEdge> >& AEdges);
    
    /*------------------------------------------------------------------------*/
    /** \brief Build oriented edges for the point \p APntId which is located
     *         into the volume. The input point is necessary REGULAR.
     *
     * \param[in]  APntID index of the point we work on
     * \param[in]  AVec   the [3][2] vectors of the frame in point APntID
     * \param[out] APnt   the [3][2] potential points location found
     * \param[out] AIndex the [3][2] potential point  index found
     * \param[out] AFound the [3][2] Indicate if a candidate is found or not
     */
    void buildOrientedEdgesInVolume(const int                  APntID,
                                    const gmds::math::Vector3d AVec[][2],
                                    gmds::math::Point          APnt[][2],
                                    int                        AIndex[][2],
                                    bool                       AFound[][2]);
    
    /*------------------------------------------------------------------------*/
    /** \brief Build oriented edges for the point \p APntId which is located
     *         onto a surface. The input point is necessary REGULAR.
     *
     *
     * \param[in]  APntID index of the point we work on
     * \param[in]  AVec   the [3][2] vectors of the frame in point APntID
     * \param[out] APnt   the [3][2] potential points location found
     * \param[out] AIndex the [3][2] potential point  index found
     * \param[out] AFound the [3][2] Indicate if a candidate is found or not
     */
    void buildOrientedEdgesOnSurface(const int                  APntID,
                                     const gmds::math::Vector3d AVec[][2],
                                     gmds::math::Point          APnt[][2],
                                     int                        AIndex[][2],
                                     bool                       AFound[][2]);

    /*------------------------------------------------------------------------*/
    /** \brief Build oriented edges for the point \p APntId which is located
     *         on a curve. The input point is necessary REGULAR.
     *
     * \param[in]  APntID index of the point we work on
     * \param[in]  AVec   the [3][2] vectors of the frame in point APntID
     * \param[out] APnt   the [3][2] potential points location found
     * \param[out] AIndex the [3][2] potential point  index found
     * \param[out] AFound the [3][2] Indicate if a candidate is found or not
     */
    void buildOrientedEdgesOnCurve(const int                  APntID,
                                   const gmds::math::Vector3d AVec[][2],
                                   gmds::math::Point          APnt[][2],
                                   int                        AIndex[][2],
                                   bool                       AFound[][2]);

    
    /*------------------------------------------------------------------------*/
    /** \brief Select among different point which is the bes aligned and/or
     *         at a good distance to build an oriented edge
     *
     * \param[in] AID   the candidate point we work
     * \param[in] ADot  the dot product of the ref pnt and the candidate pnt
     * \param[in] ADist the distance between the ref pnt and the candidate pnt
     *
     * \return    the index of the "best" candidate
     */

    int filterPointsForBuildingOrientedEdge(const std::vector<int>& AID,
                                            const std::vector<double>& ADot,
                                            const std::vector<double>& ADist);
    /*------------------------------------------------------------------------*/
    /** \brief Correct oriented edges to get a onsistent structure
     *
     * \param[in]  AInEdges  the list of oriented edges that we start from.
     * \param[out] AOutEdges the list of oriented edges that will be created.
     */
    void buildEdges(std::vector<std::vector<OrientedEdge> >& AInEdges,
                    std::vector<std::vector<OrientedEdge> >& AOutEdges);
    /*------------------------------------------------------------------------*/
    /** \brief ...
     */
    bool computeVolumePointFrom(const int APntIndex,
                                const gmds::math::Vector3d& AV,
                                gmds::math::Point& APnt);
    /*------------------------------------------------------------------------*/
    /** \brief ...
     */
    bool computeSurfacePointFrom(const int APntIndex,
                                 const gmds::math::Vector3d& AV,
                                 gmds::math::Point& APnt);

    /*------------------------------------------------------------------------*/
    /** \brief Check if an oriented edge going from \p AFrom to \p ATo exists
     *         in \p AEdgeSet.
     *
     * \param[in]  AFrom    A first point index
     * \param[in]  ATo      A second point index
     * \param[in]  AEdgeSet A set of oriented edges
     * \param[out] AOutEdge The edge from \p AFrom to \p ATo if it exists
     *
     * \return true if an edge going from \p AFrom to \p ATo exists in 
     *         \p AEdgeSet, false otherwise.
     */
    bool isIn(const int AFrom, const int ATo,
              const std::vector<OrientedEdge>& AEdgeSet,
              OrientedEdge& AOutEdge) const;
    
    /*------------------------------------------------------------------------*/
    /** \brief Fill the m_hc_mapping attribute for each stable points. Such a
     *         point should be the corner of 1, 2, 4 or 8 hexahedral elts
     *         depending if it is classified on a geometric vertex, curve,
     *         surface or volume.
     */
    void buildHexCorners(std::vector<std::vector<OrientedEdge> >& AEdges);
    
    /*------------------------------------------------------------------------*/
    /** \brief Build hexahedral elements
     */
    void buildHexahedral();
    /*------------------------------------------------------------------------*/
    /** \brief Build remaining cavities
     */
    void buildCavities();
    
    /*------------------------------------------------------------------------*/
    /** \brief Build 3-sided base prism elements
     */
    void buildPrism3(const std::vector<std::vector<OrientedEdge> >& AEdges);
    
    /*------------------------------------------------------------------------*/
    /** \brief Find among the free corners of \p AI and \p AJ an end point 
     *         which is different of \p AFrom but both in a free corner of 
     *         \p AI and \p AJ
     *
     * \param[in] AFrom point that must be appeared in \p AI and \p AJ corners
     * \param[in] AI    an end point
     * \param[in] AJ    another end point

     *
     * \return the id of common points
     */
    std::vector<int> findCommonPoints(const int AFrom,
                                      const int AI,
                                      const int AJ);
    /*------------------------------------------------------------------------*/
    /** \brief Starting from three corner points, it returns the free corner of
     *         the point that apperas in all of them and which is pointed to all
     *         of their origins
     *
     * \param[in]  ACorner1   first corner
     * \param[in]  ACorner2   second corner
     * \param[in]  ACorner3   third corner
     * \param[out] ACornerOut resulting corner
     *
     * \return true if we find the out corner, false otherwise
     */
    bool findCommmonLastCorner(const HexCorner& ACorner1,
                               const HexCorner& ACorner2,
                               const HexCorner& ACorner3,
                               HexCorner&       ACornerOut);
    /*------------------------------------------------------------------------*/
    /** \brief Looking at all the corners defined for origin point \p AOrigin,
     *         it returs the corner corresponding to (\p AOrigin, \p AI, \p AJ,
     *         \p AK) if it exist
     * \param[in]  AOrigin origin point
     * \param[in]  AI      an end point
     * \param[in]  AJ      another end point
     * \param[in]  AK      another end point
     * \param[out] AOut    the found hex corner if it exists
     *
     * \return true if we find the out corner, false otherwise
     */
    bool getCorner(const int AOrigin, const int AI, const int AJ,
                   const int AK, HexCorner& AOut);
    
    /*------------------------------------------------------------------------*/
    /** \brief Check if \p AC is a corner corresponding to (\p AOrigin, \p AI, 
     *         \p AJ, \p AK)
     *
     * \param[in] AC      the found hex corner if it exists
     * \param[in] AOrigin origin point
     * \param[in] AI      an end point
     * \param[in] AJ      another end point
     * \param[in] AK      another end point
     *
     * \return true if they correspond, false otherwise
     */
    bool isCorner(const HexCorner& AC, const int AOrigin,
                  const int AI, const int AJ, const int AK);
    
    /*------------------------------------------------------------------------*/
    /** \brief Find among the corners in \p AIn those having \p AFrom as an
     *         adjacent point.
     *
     * \param[in]  AIn   the set of corners to be parsed
     * \param[in]  AFrom the point they must match
     * \param[out] AOut  corners of \p AIn having \p AFrom as end point
     *
     * \return the number of found compatible corners
     */
    int findFreeCorners(const std::vector<HexCorner>& AIn,
                        const int AFrom,
                        std::vector<HexCorner>& AOut);
    
    /*------------------------------------------------------------------------*/
    /** \brief Return the free corners of \p AOrigin with \p AI and \p AI as
     *         end points
     *
     * \param[in] AOrigin corner origin points
     * \param[in] AI      an end point
     * \param[in] AJ      another end point
     *
     * \return free corners of \p AOrigin with \p AI and \p AI as end points
     */
    std::vector<HDMeshGenerator::HexCorner>
    findFreeCorners(const std::vector<int>& AOrigin,
                    const int AI, const int AJ);
    
    /*------------------------------------------------------------------------*/
    /** \brief Add a corner data to the point with index \p AIndex
     *
     * \param[in] AIndex  the point index we add a corner to
     * \param[in] AIndex1 the point index we go to in direction 1
     * \param[in] AV1     the vector in direction 1
     * \param[in] AIndex2 the point index we go to in direction 2
     * \param[in] AV2     the vector in direction 2
     * \param[in] AIndex3 the point index we go to in direction 3
     * \param[in] AV3     the vector in direction 3
     */
    void addCorner(const int AIndex,
                   const int AIndex1, const gmds::math::Vector3d& AV1,
                   const int AIndex2, const gmds::math::Vector3d& AV2,
                   const int AIndex3, const gmds::math::Vector3d& AV3);
    
    /*------------------------------------------------------------------------*/
    /** \brief Build all the corners associated to \p AOrigin by considering
     *         existing edges \p AEdges. This construction is based on geom.
     *         criteria
     *
     * \param[in] AOrigin origin point
     * \param[in] AEdges  existing edges connected to \p AOrigin
     */
    void buildCornersAsSolidAngles(const int AOrigin,
                                   std::vector<OrientedEdge>& AEdges);
    
    
    /*------------------------------------------------------------------------*/
    /** \brief Compute the mesh association for point AID. This computation is
     *         performed knowing how the point is classified (0,1,2 or 3) and
     *         the id of the curve
     *
     * \param[in] AID point id
     */
    void computeMeshAssociation(const int AID);

    /*------------------------------------------------------------------------*/
    /** \brief Return the region information (dim=3, id) of the region that
     *         contains \p APnt.
     * \param[in] APnt   the point we start from
     * \param[in] ATetID an initial tetrahedron we start from
     *
     * \return the info (3, id) that characterizes the tet containing APnt
     */
    gmds::Cell::Data getRegionContaining(const gmds::math::Point& APnt,
                                         const gmds::TCellID ATetID);
    /*------------------------------------------------------------------------*/
    /** \brief Return the face information (dim=2, id) of the boundary face that
     *         contains \p APnt.
     * \param[in] APnt    the point we start from
     * \param[in] ATetID  an initial tetrahedron we start from
     * \param[in] ASurfID the id of the surface \p APnt is supposed to be on
     *
     * \return the info (2, id) that characterizes the face containing APnt
     */
    gmds::Cell::Data getBoundaryFaceContaining(gmds::math::Point& APnt,
                                               const gmds::TCellID ATetID,
                                               const double ADistance,
                                               const int ASurfID);
    
    /*------------------------------------------------------------------------*/
    /** \brief Return the edge information (dim=2, id) of the boundary edge that
     *         contains \p APnt.
     * \param[in] APnt    the point we start from
     * \param[in] ATetID  an initial tetrahedron we start from
     * \param[in] ASurfID the id of the curve \p APnt is supposed to be on
     *
     * \return the info (1, id) that characterizes the edge containing APnt
     */
    gmds::Cell::Data getBoundaryEdgeContaining( gmds::math::Point& APnt,
                                               const gmds::TCellID ATetID,
                                               const double ADistance,
                                               const int ACurvID);
    /*------------------------------------------------------------------------*/
    /** \brief This function give the list of region ids such that generated
     *          points associated to those regions could be at a distance
     *          less or equal to \p AEps from  \p APnt
     *
     * \param[in] AP    the point we start from
     * \param[in] AT    the tetrahedron which should containt AP
     * \param[in] AEps  the distance we look for
     *
     * \return the ids of the tetrahedra that contains points at the specified
     *         distance
     */
    
    std::set<gmds::TCellID>
    getCloseRegionsFrom(const gmds::math::Point& AFromPnt,
                        const gmds::Region& AFromTet,
                        const double AEpsilon);
   bool getOppositeRegion(const gmds::TCellID  ANodeID,
                           const gmds::Region&  AR,
                           gmds::Region&        AOut);
    bool getOppositeFace(const gmds::TCellID  ANodeID,
                         const gmds::Region&  AR,
                         gmds::Face&        AOut);
    bool getOppositeBndFace(const gmds::TCellID ANodeID,
                            const gmds::Face&   AR,
                            gmds::Face&         AOut);
    void updateTetMesh();
    /*------------------------------------------------------------------------*/
    /** \brief Return the face of \p AR sharing nodes \p ANI and \p ANJ
     *         with face \p AFrom
     *
     * \param[in] AFrom  a face of \p AR whose \p ANI and \p ANJ are two nodes
     * \param[in] AR     a region
     * \param[in] ANI    a node id
     * \param[in] ANJ    a node id
     *  
     * \return the expected face
     */

    gmds::Face getFace(const gmds::Face& AFrom,
                       const gmds::Region& AR,
                       const gmds::TCellID ANI,
                       const gmds::TCellID ANJ);
    
    
    bool extractBndLoops(std::vector<gmds::Face>& AFaces,
                         std::vector<std::vector<gmds::Node> > & ALoops);
    
    /*------------------------------------------------------------------------*/
    /** \brief Starting from \p ALoops, i.e. nodes defined loops of cavity in 
     *         the hex_dom mesh, this function attempts to close all the cavities
     *         by quad and triangles on surfaces. Each patch of faces is in \p
     *         return in \p AFaces.
     *
     * \param[in ] ALoops cavity loops
     *
     */
    bool imprintLoop(std::vector<std::vector<gmds::Node> >& ALoops);
    
    
    void meshSurfacePatch(std::vector<gmds::Face>& ASurfFaceRep,
                          std::vector<gmds::math::Vector>& ASurfNormalRep,
                          std::vector<gmds::math::Point>& ABndPoints,
                          gmds::IGMesh& APatch);

    void extractFinalCavityBoundaries(std::vector<std::vector<gmds::Face> >& AFinalCavities);

    void closeCavity(std::vector<gmds::Face>& AC);
    
    /*------------------------------------------------------------------------*/
    /** \brief Use a tet meshing algorithm to fill in the cavity made of
     *         faces \p AC
     *
     * \param[in] AC Faces defining the boundary of a closed cavity
     */
    void buildTetMesh(std::vector<gmds::Face>& AC);

    
    /*------------------------------------------------------------------------*/
    /** \brief Use the whisker weaving algorithm to mesh the cavity made of 
     *         faces \p AC
     *
     * \param[in] AC Faces defining the boundary of a closed cavity
     */
    void applyWhiskerWeaving(std::vector<gmds::Face>& AC);

    
    gmds::math::Vector3d getOutputNormal(gmds::Face& AFace, gmds::Region& ARegion);
    gmds::math::Vector3d getOutputNormalOfABoundaryFace(gmds::Face& AFace);
    
    /*------------------------------------------------------------------------*/
    /** \brief Starting from the boundary faces of all cavities in \p AInFaces,
     *         it splits them into separate cavities stored in \p ACavFAces
     *
     * \param[in]  AInFaces  cavity boundary faces
     * \param[out] ACavFaces faces splitted by cavity ids
     */
    void flagCavityFaces(std::vector<gmds::Face>& AInFaces,
                         std::vector<std::vector<gmds::Face> >& ACavFaces);
    
    /*------------------------------------------------------------------------*/
    /** \brief Returns the id of the first item of \AV being false
     *
     * \param[in]  AV a vector of boolean values
     *
     * \return the id of the first false item, gmds::NullID otherwise
     */
    gmds::TCellID foundFreeIndex(const std::map<gmds::TCellID, bool>& AVec) const;
    
    /*------------------------------------------------------------------------*/
    /** \brief Check if \p AI is in \p AV
     *
     * \param[in] AI an id
     * \param[in] AV a vector of ids
     *
     * \return true if \p AI is in \p AV, false otherwise
     */
    bool isIn(const gmds::TCellID AI, const std::vector<gmds::TCellID>& AV) const;
    
    std::vector<gmds::TCellID> getFaces(const gmds::Node& ANI, const gmds::Node& ANJ) const;
    
    //The point is moved during this process
    gmds::Face closestFace(gmds::math::Point& AP, const int ASurfID);
    gmds::Edge closestEdge(gmds::math::Point& AP, const int ACurvID);
    
    
    void writeInput();
    void writeOutput();
    void writeEdges(std::vector<std::vector<OrientedEdge> >& AInEdges,
                    const std::string& AFileName);
    void writeLoops(std::vector<std::vector<gmds::Node> >& ALoops);
protected:
    /** Mesh we start from */
    gmds::IGMesh* m_mesh;
    /** FHeDo global parameters*/
    ParamsGlobal m_param_gl;
    
    /** Hex Dominant parameters*/
    ParamsHexDom m_param_hexdom;
    
    /** All boolean marks*/
    ParamsMark m_bm;
    /** rotational field defined on vertices of m_mesh*/
    gmds::Variable<gmds::math::AxisAngleRotation>* m_rot_field;

    
    /** points used to build the mesh*/
    std::vector<gmds::math::Point> m_pnt;
    /** rotation associated to each point used to build the mesh*/
    const std::vector<gmds::math::Chart> m_chart;
    std::vector<gmds::Cell::Data> m_mesh_data;
    /** type associated to each point used to build the mesh*/
    const std::vector<int>               m_type;
    /** geometric classification of each point used to build the mesh*/
    const std::vector<int>               m_classification;
    /** Curve numbering of each point used to build the mesh*/
    const std::vector<int>               m_curve;
    /** surface numbering of each point used to build the mesh*/
    const std::vector<int>               m_surface;
    /** Defines a normal constraint alignment for boundary nodes
     * we work on */
     std::map<int, gmds::math::Vector3d> m_normal;
    /** initial required spacing between points*/
    double m_spacing;
    
    /* vector given for each point if it is already used as a hex corner*/
    std::vector<int>               m_used;
    /** the generated hex-dom mesh we work on */
    gmds::IGMesh                   m_hexdom;
    
    
    double m_dot_tolerance;
    

    /** mapping from each point to the number of hex corner it must be in*/
    std::map<int, std::vector<HexCorner> > m_hc_mapping;
    std::map<int, gmds::Node> m_node_mapping;
    
    std::map<int,std::vector<int> > m_filter;
};
    /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* SH_HEX_DOM_MESH_GENERATOR_H_ */
/*----------------------------------------------------------------------------*/
