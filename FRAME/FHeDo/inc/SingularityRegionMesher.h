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
#ifndef SH_SINGULARITY_REGION_MESHER_H_
#define SH_SINGULARITY_REGION_MESHER_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
/*----------------------------------------------------------------------------*/
// GMDS File Headers
#include <GMDS/IG/IGMesh.h>
#include <GMDS/Math/Chart.h>
/*---------------------------------------------------------------------------*/
// SingularityGraph File Headers
#include <SingularityGraph.h>
/*----------------------------------------------------------------------------*/
// FHeDO File Headers
#include "Params.h"
/*----------------------------------------------------------------------------*/
namespace fhedo{
/*----------------------------------------------------------------------------*/
/** \class  SingularityRegionMesher
 *  \brief  Starts from a tet mesh with a frame field, this algorithm extracts
 *          3D singularity lines of the frame field, build hexahedral regions
 *          around them, and finally modify the initial mesh to be the same
 *          but without new hexahedral regions.
 */
class EXPORT_GMDS SingularityRegionMesher{
    
public:
    /*------------------------------------------------------------------------*/
    /** \brief Constructor.
     *
     * \param[in] AMesh background tetrahedral mesh we start from
     */
    SingularityRegionMesher(gmds::IGMesh* AMesh,
                            ParamsGlobal& AGlobal,
                            double AEdgeSize,
                            ParamsMark& AMarks);
    /*------------------------------------------------------------------------*/
    /** \brief Destructor.
     */
    virtual ~SingularityRegionMesher();
    
    /*------------------------------------------------------------------------*/
    /** \brief Execute the algorithm
     */

    void execute();

protected:
    /*------------------------------------------------------------------------*/
    /** \brief Detects the singular cells of the frame field
     */
    void detectSingularCells();
    /*------------------------------------------------------------------------*/
    /** \brief initialize the used boolean marks specific to this algorithm
     */
    void initializeBooleanMarks();
    
    /*------------------------------------------------------------------------*/
    /** \brief clean the used boolean marks specific to this algorithm
     */
    void finalizeBooleanMarks();

    void createVolumeSingularityLines();
    void createVolumeSingularityPoints();
    void createOneVolumeSingularityPoint(std::vector<gmds::Region>& ACluster);
    void createBoundarySingularityPoints();
    
    std::vector<gmds::math::Vector3d>
    defineBoundarySlotsViaAngles(const std::vector<gmds::Node>& ANodes,
                                 const gmds::math::Point&       ASingLoc,
                                 const gmds::math::Vector&      ANormal);
    std::vector<gmds::math::Vector3d>
    defineBoundarySlotsViaVectors(const std::vector<gmds::Node>& ANodes,
                                  const gmds::math::Point&       ASingLoc,
                                  const gmds::math::Vector&      ANormal);
    void createOneSingularityLineFrom(SingularityPoint* ASingPoint,
                                      SingularityPoint::Slot* ASlot,
                                      const int AMarkFaceUsedForSep);
    
    gmds::math::Vector getOutputNormal(const gmds::Face& AF,
                                       const gmds::Region& AR);
    gmds::math::Vector getInputNormal(const gmds::Face& AF,
                                      const gmds::Region& AR);
    
    void  writeOutput(const std::string& AFileName);
    void  writeOutputSingle(const std::string& AFileName);
    void createSurfaceQuadPatches();
    void createLineTubes();
protected:
    /** the mesh we work on*/
    gmds::IGMesh* m_mesh;
    /** fhedo global parameters*/
    ParamsGlobal m_params_gl;
    
    /** expected target mesh size*/
    double m_edge_size;
    
    /** rotation field defined on the vertices of m_mesh*/
    gmds::Variable<gmds::math::AxisAngleRotation>* m_rot;

    /** set of marks for boundary cells (nodes, edges, faces) */
    ParamsMark m_marks;
    
    
    gmds::IGMesh* m_line_mesh;
    int m_mark_line_sing;
    int m_mark_volume_pnt_sing;
    int m_mark_face_sing;
    int m_mark_edge_sing;
    int m_mark_node_sing;
    SingularityGraph m_graph;
    
    /** store the ids of singular region*/
    std::vector<gmds::TCellID> m_sing_tet;
    /** store for each  singular tet, the number of singular faces it has*/
    std::map<gmds::TCellID, int> m_sing_tet_type;
    /** store the ids of singular boundary faces*/
    std::vector<gmds::TCellID> m_sing_bnd_tri;
    /** for each sing point, we store the quad patch created for it*/
    std::map<void*,std::vector<gmds::TCellID> > m_sing_pnt_to_patch_circle;
    std::map<void*,gmds::TCellID > m_sing_pnt_to_patch_center;
};
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* SH_TETMESH_MANIPULATOR_H_ */
/*----------------------------------------------------------------------------*/
