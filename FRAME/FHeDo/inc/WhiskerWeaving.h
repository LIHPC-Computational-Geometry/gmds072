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
#ifndef FHEDO_WHISKER_WEAVING_H_
#define FHEDO_WHISKER_WEAVING_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <map>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
#include <GMDS/IG/IGMesh.h>
#include <GMDS/Math/VectorND.h>
/*----------------------------------------------------------------------------*/
namespace fhedo{
    
    /*------------------------------------------------------------------------*/
    /** \class WhiskerWeaving
     *  \brief Implementation of a geometrical whisker weaving in the context
     *         of hex-dominant mesh generation
     */
    class EXPORT_GMDS WhiskerWeaving{
    public:
        /*-------------------------------------------------------------------*/
        /** \brief Constructor
         */
        WhiskerWeaving();
        /*-------------------------------------------------------------------*/
        /** \brief Destructor
         */
        virtual ~WhiskerWeaving();
        
        /*-------------------------------------------------------------------*/
        /** \brief Execute the algorithm starting from a boundary mesh \p AM
         *         knowing that the inner side of the volume to be meshed is
         *         given by the direction \p ADir in face of id \p AF
         *
         * \param[in] AM a boundary mesh we start from
         * \param[in] AF the face the inward direction is defined on
         * \param[in] ADir inward direction at face \p AF
         *
         * \return an algorithm status indicating if the algorithm proceeded
         *         as expected.
         */
        int execute(gmds::IGMesh* AM,
                    gmds::TCellID AF,
                    gmds::math::Vector3d& ADir);
        
        /*-------------------------------------------------------------------*/
        /** \briefGive the main name of debug files
         * \param[in] AName debug file prefix
         */
        void setDebugFile(const std::string& AName);
        
        
        struct Loop{
            std::vector<gmds::TCellID> faces;
            //Going from face[0] to face[1] gives the orientation (left/right) using also the inward normal
            bool self_intersect;
            
            std::vector<gmds::TCellID> right_nodes[2];
            std::vector<gmds::TCellID> left_nodes[2];
            std::vector<gmds::TCellID> right_order_nodes;
            std::vector<gmds::TCellID> left_order_nodes;
            /*nb sharp edges traversed by face*/
            int flat_index;
            /*nb sharp CONVEX edges on the left*/
            int left_sharp;
            /*nb sharp CONVEX edges on the right*/
            int right_sharp;
            /*nb faces enclosed on the left*/
            int left_nb_faces;
            /*nb faces enclosed on the right*/
            int right_nb_triangles;
            /*nb triangles enclosed on the left*/
            int left_nb_triangles;
            /*nb triangles enclosed on the right*/
            int right_nb_faces;

            /*max_depth on the left*/
            int left_depth;
            /*max_depth on the right*/
            int right_depth;
            
            void init(){
                self_intersect    = false;
                flat_index        = 0;
                left_sharp        = 0;
                right_sharp       = 0;
                left_nb_faces     = 0;
                right_nb_faces    = 0;
                left_nb_triangles = 0;
                right_nb_triangles= 0;
                left_depth        = 0;
                right_depth       = 0;
            }
        };
    private:

        /*-------------------------------------------------------------------*/
        /** \brief Create the boundary mesh we will work on from input data
         *
         * \param[in] AM a boundary mesh we start from
         * \param[in] AF the face the inward direction is defined on
         * \param[in] ADir inward direction at face \p AF
         */
        void prepareMesh(gmds::IGMesh* AM,
                         gmds::TCellID AF,
                         gmds::math::Vector3d& ADir);
        /*-------------------------------------------------------------------*/
        /** \brief build loops that can be contracted on the current front
         */
        std::vector<Loop> buildLoops();
        /*-------------------------------------------------------------------*/
        /** \brief Compute the data associate to loop \p AL
         *
         * \param[in] AL a loop
         */
        void computeData(Loop& AL);
        
        /*-------------------------------------------------------------------*/
        /** \brief Select the loop to shrink in \p ALoops
         *
         * \param[in]  ALoops   a set of loops
         * \param[out] AL       selected loop
         * \param[out] AToRight if yes shrink to right, otherwise to left
         *
         * \return true if a chord is selected, false otherwise
         */
        bool select(const std::vector<Loop>& ALoops,
                    Loop& ALoop, bool& AToRight);
        /*-------------------------------------------------------------------*/
        /** \brief Returns the ids of the front faces enclosed by \p AL in 
         *         direction \p ADir
         *
         * \param[in]  AL   loop
         * \param[in]  ADir direction to go
         * \param[out] AF   enclosed faces
         */
        void getEnclosedFaces(const Loop& ALoop, const bool AToRight,
                              std::vector<gmds::TCellID>& AF);
        
        /*-------------------------------------------------------------------*/
        /** \brief Shrink the loop \p AL. If it fails to do it, nothing is 
         *         done, and it returns false
         *
         * \param[in] AL a loop
         *
         * \return true if the shrink process was done, false otherwise
         */
        bool shrink(const Loop& AL, const bool AToRight);
        
        /*-------------------------------------------------------------------*/
        /** \brief Compute the side info of a loop \p ALoop where all faces
         *         are marked with \p AMark.
         *
         * \param[in] AL a loop
         * \param[in] AM a mark that distinguish all faces of \p AL
         * \param[in] AS the side we compute info (true - right, letf - false)
         */
        void computeLoopSideInfo(Loop& AL, const int AM, const bool side);
        
        /*-------------------------------------------------------------------*/
        /** \brief Returns the front face sharing node \p AN1 and \p AN2
         *         with \p AF
         *
         * \param[in] AF A face
         * \param[in] AN1 one node of AF
         * \param[in] AN2 another node of AF
         *
         * \return the id of the face sharing \p AN1 and \p AN2 with \p AF.
         */
        gmds::TCellID getFace(const gmds::Face& AF,
                              const gmds::TCellID AN1,
                              const gmds::TCellID AN2);
        
        bool isInPolygon(const gmds::TCellID AF,
                         const std::vector<gmds::math::Point>& APolyPnts,
                         const gmds::math::Plane& APolyPl);
        
        /*-------------------------------------------------------------------*/
        /** \brief Check if a face of \p AF has exactly the nodes of \p AN.
         *
         * \param[in] AF A set of faces
         * \param[in] AN A set of nodes
         *
         * \return true if a face is based on nodes \p AN, false otherwise
         */
        bool existFace(const std::vector<gmds::TCellID>& AF,
                       const std::vector<gmds::TCellID>& AN);
        /*-------------------------------------------------------------------*/
        /** \brief Returns all the faces sharing an edge with \p AF. If a
         *         same face share two edges with \p AF, it will appear twice
         *         in the result.
         *
         * \param[in] AF A face
         *
         * \return the id of the face sharing an edge with \p AF.
         */
        std::vector<gmds::TCellID> getAdjFaces(const gmds::Face& AF);
        /*-------------------------------------------------------------------*/
        /** \brief Returns the two nodes of the quad faces \p AF that are not
         *         \p AN1 and \p AN2
         *
         * \param[in]  AF A face
         * \param[in]  AN1 one node of AF
         * \param[in]  AN2 another node of AF
         * \param[out] AN3 other node of AF
         * \param[out] AN4 other node of AF
         */
        void getOtherNodesFace(const gmds::Face& AF,
                               const gmds::TCellID AN1,
                               const gmds::TCellID AN2,
                               gmds::TCellID& AN3,
                               gmds::TCellID& AN4);
        /*-------------------------------------------------------------------*/
        /** \brief Clean the front after a loop shrink. Remove adjacent faces
         *         sharing 2 edges.
         */

        void cleanFront();
        /*-------------------------------------------------------------------*/
        /** \brief Reverse local node numerotation of \p AF if it does not go
         *         from \p A From to \p ATo, two ids of nodes of \p AF.
         *
         * \param[in] AF A face
         * \param[in] AFrom one node of AF
         * \param[in] ATo another node of AF
         */
        void reverse(gmds::Face& AF,
                     gmds::TCellID AFrom,
                     gmds::TCellID ATo);
        
        /*-------------------------------------------------------------------*/
        /** \brief Indicates if loops \p AL1 and \p AL2 are the same
         *
         * \param[in] AL1 a first loop
         * \param[in] AL2 a second loop
         *
         * \return true if AL1==AL2, false otherwise
         */

        bool same(const Loop& AL1, const Loop& AL2);
        /*-------------------------------------------------------------------*/
        /** \brief Indicates if angle \p AAngle is sharp for this algorithm
         *
         * \param[in] An angle in degree
         *
         * \return true if sharp, false otherwise
         */
        bool isSharpAngle(const double& AAngle);
        /*-------------------------------------------------------------------*/
        /** \brief Indicates if angle \p AAngle is flat for this algorithm
         *
         * \param[in] An angle in degree
         *
         * \return true if flat, false otherwise
         */
        bool isFlatAngle(const double& AAngle);
        /*-------------------------------------------------------------------*/
        /** \brief Indicates if angle \p AAngle is convex for this algorithm
         *
         * \param[in] An angle in degree
         *
         * \return true if convex, false otherwise
         */
        bool isConvexAngle(const double& AAngle);
        /*-------------------------------------------------------------------*/
        /** \brief Indicates if angle \p AAngle is concave for this algorithm
         *
         * \param[in] An angle in degree
         *
         * \return true if convex, false otherwise
         */
        bool isConcaveAngle(const double& AAngle);

        /*-------------------------------------------------------------------*/
        /** \brief Return the oriented dihedral angle between \p AF1 and \p AF2
         *         in degree
         *
         * \param[in] AF1 A face
         * \param[in] AF2 Another face
         *
         * \return the oriented angle between \p AF1 and \p AF2
         */
        double dihedralAngle(const gmds::Face& AF1, const gmds::Face& AF2);
        
        /*-------------------------------------------------------------------*/
        /** \brief Give the ids of the nodes shared by \p AF1 and \p AF2
         *
         * \param[in] AF1 A face
         * \param[in] AF2 Another face
         * \param[out] AN1 first node id
         * \param[out] AN2 second node id
         */
        void commonNodes(const gmds::Face& AF1, const gmds::Face& AF2,
                         gmds::TCellID& AN1, gmds::TCellID& AN2);
        /*-------------------------------------------------------------------*/
        /** \brief Starting from quad \p AQ through adjacent edge 
         *         [\p AN1, \p AN2], it build a loop on the front set of faces
         *
         * \param[in] AQ A quad face
         * \param[in] AN1 first node id
         * \param[in] AN2 second node id
         * \param[in/out] AC1 color for faces trav. through their first edge
         * \param[in/out] AC2 color for faces trav. through their first edge
         *
         * \return true if the loop has been build (come back on the first face)
         */

        bool buildLoop(const gmds::Face& AQ,
                       const gmds::TCellID& AN1,
                       const gmds::TCellID& AN2,
                       std::map<gmds::TCellID,int>& AC1,
                       std::map<gmds::TCellID,int>& AC2,
                       Loop& ALoop);
        
        void writeDebugMesh();
    private:
        /** the mesh we work on*/
        gmds::IGMesh m_mesh;
        
        /** map from ids of the bnd nodes to ids of inner volume nodes */
        std::map<gmds::TCellID, gmds::TCellID> m_from_bnd_nodes;
        /** map from ids of the inner volume nodes to ids the bnd nodes */
        std::map<gmds::TCellID, gmds::TCellID> m_to_bnd_nodes;
        /** map from ids of the bnd faces to ids of inner volume faces */
        std::map<gmds::TCellID, gmds::TCellID> m_from_bnd_faces;
        /** map from ids of the inner volume faces to ids the bnd faces */
        std::map<gmds::TCellID, gmds::TCellID> m_to_bnd_faces;
        /** debug file prefix*/
        std::string m_debug_file;
        
        /** mark for faces in the current front */
        int m_front_mark;
        
        /** variable attached to faces to give the distance of a face to
         *  the initial boundary, used to give priority during loop 
         *  selection*/
        gmds::Variable<int>* m_depth;
        
        
    };
    
    /*------------------------------------------------------------------------*/
}
std::ostream& operator<<(std::ostream& AO,const fhedo::WhiskerWeaving::Loop& AL);
/*----------------------------------------------------------------------------*/
#endif /* FHEDO_WHISKER_WEAVING_H_ */
/*----------------------------------------------------------------------------*/
