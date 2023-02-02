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
#ifndef SH_FIELD_GENERATOR_H_
#define SH_FIELD_GENERATOR_H_

/*----------------------------------------------------------------------------*/
// STL File Headers
#include <vector>
#include <map>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
#include <GMDS/IG/IGMesh.h>
#include <GMDS/Math/Chart.h>
#include <GMDS/Math/SHarmonicL4.h>
/*---------------------------------------------------------------------------*/
// Frame File Headers
#include "FieldSolverStrategyItf.h"
#include "Params.h"
/*---------------------------------------------------------------------------*/
namespace fhedo{
/*----------------------------------------------------------------------------*/
/** \class  FieldGenerator
 *  \brief  Computes a frame field onto a 3D mesh following a direct approach
 *          where a linear system is solved and frames are represented by
 *          Spherical Harmonics (SH)
 */
class EXPORT_GMDS FieldGenerator{
    
public:
    /*------------------------------------------------------------------------*/
    /** \brief Constructor.
     *
     * \param[in] ASolver the way the linear solver will be build and solve
     * \param[in] AMesh the mesh where we work on
     * \param[in] MarkNodeOnSurface mark for nodes classified on a surface
     * \param[in] MarkNodeOnCurve   mark for nodes classified on a curve
     * \param[in] MarkNodeOnPoint   mark for nodes classified on a point
     * \param[in] MarkEdgeOnSurface mark for edges classified on a surface
     * \param[in] MarkEdgeOnCurve   mark for edges classified on a curve
     * \param[in] MarkFaceOnSurface mark for faces classified on a surface
     * \param[in] ANbMaxIterationSmoothing the max number of iteration during
     *            the smoothing process
     */
    FieldGenerator(FieldSolverStrategyItf* ASolver,
                   gmds::IGMesh* AMesh,
                   const ParamsGlobal& AGParam,
                   const ParamsFrameField& AParam,
                   const ParamsMark& ABM);
    
    
    /*------------------------------------------------------------------------*/
    /** \brief Function to be called for generating the frame field
     */
    void execute();
    
    
    /*------------------------------------------------------------------------*/
    /** \brief Function returning the computed field energy as defined in
     *         RAY, SOKOLOF paper "On Smooth Frame Field"
     */
    double computeFieldEnergy();
    
    /*------------------------------------------------------------------------*/
    /** \brief Returns the chart field build by this algorithm with 1 chart
     *         per mesh node.
     */
    gmds::Variable<gmds::math::Chart>* chartField(){
        return m_chart_field;
    }
protected:
    
    
    std::vector<gmds::Edge> getEdgesOnCurve(gmds::Node& ANode) const;
    std::vector<gmds::Face> getFacesOnSurface(gmds::Edge& AEdge) const;
    gmds::Node  getNeighboorOn(gmds::Node& AN, gmds::Edge& AE) const;
    
    
    
    virtual void buildFrameField() {;}
    virtual void computeDistanceFields() {;}
    virtual void initQuaternionsBySnapping() {;}
    
    virtual void smoothAll(){;}
    /*------------------------------------------------------------------------*/
    /** \brief Sort nodes to gather boundary nodes first then inner nodes
     */
    void sortNodes();
    
    /*------------------------------------------------------------------------*/
    /** \brief Compute and store the normal to each boundary node. A single
     *         normal is assigned to nodes classified on surfaces, while a
     *         complete frame is computed for nodes classified on curves. Nodes
     *         classified on vertices are not taken into account.
     */
    void computeBoundaryData();
    
    void initHarmonicOnPoint  (gmds::Node& ANode);
    void initHarmonicOnCurve  (gmds::Node& ANode);
    void initHarmonicOnSurface(gmds::Node& ANode);
    /*------------------------------------------------------------------------*/
    /** \brief Build the system to solve for iteration \p AI. Default value for
     *         \p AI is 0 meangin the intilialization stage. Previous solution
     *         is also provided if necessary.
     *
     * \param[in] APrevSolution a representation of the previous solution
     * \param[in] AI the iteration number for assembling the system
     */
    void buildSystem(const int AI=0);
    
  
    void writeSolution();
    gmds::math::Vector9d extractSolutionFromX(const gmds::Node& AN);

    
    /*------------------------------------------------------------------------*/
    /** \brief Initialize the boolean marks
     */
    void initMarks();
    
    /*------------------------------------------------------------------------*/
    /** \brief Clean the boolean marks
     */
    void cleanMarks();
    /*------------------------------------------------------------------------*/
    /** \brief Mark the boundary cells
     */
    void markBoundaryCells();

protected:
    
    /** the object used to solve the linear system */
    FieldSolverStrategyItf* m_solver;
    /** the mesh we work on */
    gmds::IGMesh* m_mesh;
    
    ParamsGlobal m_global_params;
    ParamsFrameField m_params;
    ParamsMark m_bm;
    
    /*** Node mark for all the nodes classified on points and curves*/
    int m_markNodeFixed;
    /*** Node mark for all the nodes on the boundary*/
    int m_markNodeOnBnd;
    /*** Edge mark for all the edges used in the minimization pb*/
    int m_markEdgeSmooth;
    
    gmds::Variable<gmds::math::SHarmonicL4>* m_harmonic_field;
    gmds::Variable<gmds::math::Chart>*       m_chart_field;

    gmds::Variable<int>*       m_ordering;
    std::map<gmds::TCellID, gmds::math::SHarmonicL4> m_H0;
    std::map<gmds::TCellID, gmds::math::SHarmonicL4> m_H4;
    std::map<gmds::TCellID, gmds::math::SHarmonicL4> m_H8;

    
    std::vector<gmds::Node> m_surface_nodes;

    /*** number of inner, boundary nodes */
    int m_nb_inner_nodes; //nodes inside the volume
    int m_nb_fixed_nodes; //nodes on curves and pnts
    int m_nb_surface_nodes; //nodes on surface
    int m_nb_free_nodes; //inner+surface

    // All the edges, but those connecting two fixed nodes
    int m_nb_smooth_edges;
    int m_nb_strict_smooth_edges;
    
    
};
    /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* SH_FIELD_GENERATOR_H_ */
/*----------------------------------------------------------------------------*/
