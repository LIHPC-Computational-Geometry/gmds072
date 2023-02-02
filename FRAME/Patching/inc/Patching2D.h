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
 * Patching2D.h
 *
 *  Created on: 13 juil. 2014
 *      Author: bibi
 */
/*----------------------------------------------------------------------------*/
#ifndef PATCHING_2D_H_
#define PATCHING_2D_H_
/*----------------------------------------------------------------------------*/
// STL Headers
#include <map>
/*----------------------------------------------------------------------------*/
// GMDS Headers
#include <GMDS/IG/IGMesh.h>
#include <GMDS/Math/Cross2D.h>
#include <GMDS/Math/Segment.h>
/*----------------------------------------------------------------------------*/
// FRAME Headers
#include <PatchComplex2D.h>
/*----------------------------------------------------------------------------*/
class EXPORT_GMDS Patching2D
{
public:
    /*------------------------------------------------------------------------*/
    /* \brief Constructor starting from a mesh having a cross field defined onto
     *        it.
     *
     * \param AMesh the background mesh we work on
     * \param AField the cross field associated to AMesh
     */
    Patching2D(gmds::IGMesh* AMesh,
               gmds::Variable<gmds::math::Cross2D>* AField);
    
    /*------------------------------------------------------------------------*/
    /* \brief Destructor
     */
    virtual ~Patching2D();
    
    /*------------------------------------------------------------------------*/
    /* \brief Main method called to launch the algorithm
     *
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
    
    
private:
    
    /*------------------------------------------------------------------------*/
    /* \brief method used for detecting triangles that contain a field
     *        singularity
     *
     */
    void detectSingularTriangles();
    
    void growPatchs();
    /*------------------------------------------------------------------------*/
    /* \brief Build individual patchs for every triangles that is regular. Built
     *        patches are stored in m_patchs
     *
     */
    void buildRegularPatchs();
    
    /*------------------------------------------------------------------------*/
    /* \brief ...
     *
     */
    void buildLocalPatch(gmds::Face& AFace,
                         const std::vector<gmds::math::Point>& ACorners);
    
    /*------------------------------------------------------------------------*/
    /* \brief Methods used for debug purposes
     *
     */
    void writeOutput(const std::string& AFileName);
    void writeOutputSingle(const std::string& AFileName);

    bool tryToAdd(const gmds::math::Point& AP,
                  std::vector<gmds::math::Point>& AV,
                  const double& ATol) const;
    
    /*------------------------------------------------------------------------*/
    /* \brief  Compute the target size has 1/3 of the smallest edge
     *
     */
    void computeTargetSize();
    

private:
    /** Mesh we start from */
    gmds::IGMesh* m_mesh;
    
    /** the patch complex we build on m_mesh*/
    PatchComplex2D m_complex;
    /** Mesh we create */
    gmds::IGMesh m_quad_mesh;
    /* Cross field we start from*/
    gmds::Variable<gmds::math::Cross2D>* m_field;
    
    
    /** Seeds for creating growing patchs are located at the geometrical
     *  vertices or in the field-singular tirangles. A seed cell is defined
     *  by a pair id+dim*/
    std::vector<gmds::TCellID>  m_seed_ids;
    std::vector<int>            m_seed_dims;
    
    /** DEBUG variable which stores the index of every single triangle of m_mesh */
    gmds::Variable<int>* m_index;
    
    /**directory where debug files will be written*/
    std::string m_output_directory_name;
    
    /** mark for nodes classified on geometric points */
    int m_mark_nodes_on_point;
    /** mark for nodes classified on geometric curves */
    int m_mark_nodes_on_curve;
    /** mark for edges classified on geometric curves */
    int m_mark_edges_on_curve;
    /** mark for faces which contain singularities */
    int m_mark_faces_with_sing;
    
    /** we keep in mind the triangles that contains a singularty point */
    std::vector<gmds::Face> m_singular_triangles;
    
    /** target size*/
    double m_target_size;
    
};
/*----------------------------------------------------------------------------*/
#endif /* PATCHING_2D_H_ */
/*----------------------------------------------------------------------------*/
