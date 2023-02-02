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
#ifndef SH_FIELD_SOLVER_STRATEGY_ITF_H_
#define SH_FIELD_SOLVER_STRATEGY_ITF_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <map>
/*---------------------------------------------------------------------------*/
// GMDS File Headers
#include <GMDS/Utils/CommonTypes.h>
#include <GMDS/IG/IGMesh.h>
#include <GMDS/Math/Chart.h>
#include <GMDS/Math/SHarmonicL4.h>
#include <GMDS/Math/AxisAngleRotation.h>
/*----------------------------------------------------------------------------*/
namespace fhedo{
/*---------------------------------------------------------------------------*/
// Frame File Headers
/*----------------------------------------------------------------------------*/
/** \class  FieldSoverStrategyItf
 *  \brief  Interface defining the operations that must be provided to the
 *          FieldGenerator class in order to build a 2D or 3D frame field
 *          (based on Spherical harmonics)
 */
class EXPORT_GMDS FieldSolverStrategyItf{

public:
    
    /*------------------------------------------------------------------------*/
    /** \brief Default constructor
     */
    FieldSolverStrategyItf():
    m_mesh(0),
    m_harmonic_field(0),
    m_chart_field(0),
    m_ordering(0),
    m_H0(0),
    m_H4(0),
    m_H8(0),
    m_surface_nodes(0),
    m_nb_unknowns(0),
    m_markNodeLocked(0),
    m_penalty(10)
    {;}
    
    /*------------------------------------------------------------------------*/
    /** \brief Default destructor
     */
    virtual ~FieldSolverStrategyItf(){};
    /*------------------------------------------------------------------------*/
    /** \brief Gives the mesh to work on
     *
     * \param[in] AMesh The mesh we work on
     */
    void setMesh(gmds::IGMesh* AMesh){m_mesh = AMesh;}
    /*------------------------------------------------------------------------*/
    /** \brief Gives the Sph. Harmonic variable used to store SH in the mesh
     *
     * \param[in] AV the variable where SH are stored and updated
     */
    void setHarmonicVariable(gmds::Variable<gmds::math::SHarmonicL4>* AV){
        m_harmonic_field = AV;
    }
    /*------------------------------------------------------------------------*/
    /** \brief Gives the chart variable used to store orientations in the mesh
     *
     * \param[in] AV the variable where charts are stored and updated
     */
    void setChartVariable(gmds::Variable<gmds::math::Chart>* AV){
        m_chart_field = AV;
    }
    
    void setOrderingVariable(gmds::Variable<int>* AV){
        m_ordering = AV;
    }
    void setH0(std::map<gmds::TCellID, gmds::math::SHarmonicL4>* AMap){
        m_H0 = AMap;
    }
    void setH4(std::map<gmds::TCellID, gmds::math::SHarmonicL4>* AMap){
        m_H4 = AMap;
    }
    void setH8(std::map<gmds::TCellID, gmds::math::SHarmonicL4>* AMap){
        m_H8 = AMap;
    }
    void setNbUnknowns(const int ANb){m_nb_unknowns = ANb;}
    void setNbFreeNodes(const int ANb){m_nb_free_nodes = ANb;}
    void setMarkNodeLocked(const int AMark){m_markNodeLocked = AMark;}
    void setMarkNodeOnSurface(const int AMark){m_markNodeOnSurf = AMark;}
    void setPenaltyTerm(const double AD){m_penalty = AD;}
    void setSurfaceNodes(std::vector<gmds::Node>* AV){m_surface_nodes = AV;}
    /*------------------------------------------------------------------------*/
    /** \brief Function to be called for solving the system
     */
    virtual void setX()=0;
    virtual void addSmoothingTerms()=0;
    virtual void addBoundaryConstraints()=0;
    virtual void addLockedTerms()=0;
    virtual void initializeAssembly()=0;
    virtual void finalizeAssembly()=0;
    virtual void addLocalOptimConstraints()=0;
    virtual void getFeasibleSolution()=0;
    /*------------------------------------------------------------------------*/
    /** \brief Function to be called for solving the system
     */
    virtual void solve()=0;

    virtual void clean()=0;
protected:
    gmds::IGMesh* m_mesh;

    gmds::Variable<gmds::math::SHarmonicL4>* m_harmonic_field;
    gmds::Variable<gmds::math::Chart>* m_chart_field;
    gmds::Variable<int>* m_ordering;
    std::map<gmds::TCellID, gmds::math::SHarmonicL4>* m_H0;
    std::map<gmds::TCellID, gmds::math::SHarmonicL4>* m_H4;
    std::map<gmds::TCellID, gmds::math::SHarmonicL4>* m_H8;
    std::vector<gmds::Node>* m_surface_nodes;
    
    int m_nb_unknowns;
    int m_nb_free_nodes;
    int m_markNodeLocked;
    int m_markNodeOnSurf;
    double m_penalty;
};
    /*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* SH_FIELD_SOLVER_STRATEGY_ITF_H_ */
/*----------------------------------------------------------------------------*/
