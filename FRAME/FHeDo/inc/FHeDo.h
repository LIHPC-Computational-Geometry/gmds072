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
#ifndef FHEDO_H_
#define FHEDO_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <string>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
#include <GMDS/Utils/Parameters.h>
#include <GMDS/IG/IGMesh.h>
/*----------------------------------------------------------------------------*/
// FHedo File Headers
#include "Params.h"
#include "PointGenerator.h"
/*----------------------------------------------------------------------------*/
namespace fhedo{
    /*------------------------------------------------------------------------*/
    /** \class  FHeDo
     *  \brief  Front-end class to call a frame-based hex-dom mesher.
     */
    class EXPORT_GMDS FHeDo{
        
    public:
        /*--------------------------------------------------------------------*/
        /** \brief Constructor.
         *
         * \param[in] ...
         */
        FHeDo();
        /*--------------------------------------------------------------------*/
        /** \brief Destructor.
         */
        ~FHeDo();
        
        /*--------------------------------------------------------------------*/
        /** \brief Update of the algorithm parameters
         *
         * \param[in] AFileName name of the file containing algorithm param.
         *
         * \return true if the parameter file is correct, false otherwise
         */
        bool setParameters(const std::string& AFileName);
        
        /*--------------------------------------------------------------------*/
        /** \brief Function to be called for starting the algorithm
         */
        void execute();
        
    private:
        /*--------------------------------------------------------------------*/
        /** \brief Initialize the parameter structure of this algorithm
         *
         */
        void initParams();
        /*--------------------------------------------------------------------*/
        /** \brief Read the file given in parameters and initialize m_mesh.
         */
        void readFile();
        
        /*--------------------------------------------------------------------*/
        /** \brief Read the file given in parameters and initialize m_mesh.
         */
        void prepareMesh();
        /*--------------------------------------------------------------------*/
        /** \brief Compute the target size of the mesh
         */
        void computeTargetSize();
        
        /*--------------------------------------------------------------------*/
        /** \brief Build the boundary informations that are necessary for the
         *         algorithm
         */
        void buildBoundaryData();
        /*--------------------------------------------------------------------*/
        /** \brief Removes the slivers that are made of boundary vertices only
         *         Boundary information are then updated
         *
         *  \return true if slivers are removed, false otherwise
         */
        bool removeBoundarySlivers();
        
        /*--------------------------------------------------------------------*/
        /** \brief Generate the frame field
         */
        void generateFrameField();
        
        /*--------------------------------------------------------------------*/
        /** \brief Generate the point field
         */
        void generatePointField();
        
        
        /*--------------------------------------------------------------------*/
        /** \brief Generate the hex-dom mesh
         */
        void generateHexDomMesh();
        /*--------------------------------------------------------------------*/
        /** \brief Initialize all the marks that are used
         */
        void initMeshMarks();
        /*--------------------------------------------------------------------*/
        /** \brief Clean all the marks that are used
         */
        void cleanMeshMarks();
        /*--------------------------------------------------------------------*/
        /** \brief compute the cotangent weights and assigned them as an edge
         *         variable.
         */
        void computeCotangentWeights();

        /*--------------------------------------------------------------------*/
        /** \brief Function called at the end of the point generation algorithm
         *         to write points into a file for debug/vizualisation purpose
         */
        void writePoints(const std::vector<gmds::math::Point>& APnts);
    private:
        
        gmds::IGMesh m_mesh;
        
        gmds::Parameters m_param;

        /** global parameters */
        ParamsGlobal m_param_gl;
        /** frame field generation parameters */
        ParamsFrameField m_param_ff;
        /** point field generation parameters */
        ParamsPointField m_param_pf;
        /** hex-dom algorithm parameters */
        ParamsHexDom m_param_hd;
        /** Boolean marks */
        ParamsMark m_bm;
        
        /** map that stores normal vector only for boundary nodes */
        std::map<gmds::TCellID, gmds::math::Vector3d> m_bnd_normals;

        PointGenerator* m_pg;
    };
    /*------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* FHEDO_H_ */
/*----------------------------------------------------------------------------*/
