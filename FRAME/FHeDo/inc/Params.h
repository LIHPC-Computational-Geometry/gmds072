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
#ifndef FHEDO_PARAMS_H_
#define FHEDO_PARAMS_H_
/*----------------------------------------------------------------------------*/
// STL File Headers
#include <string>
/*----------------------------------------------------------------------------*/
// GMDS File Headers
#include <GMDS/Utils/Parameters.h>
#include <GMDS/IG/IGMesh.h>
/*----------------------------------------------------------------------------*/
/** \brief This file gathers all the parameter structures used in the FHeDo
 *         component.
 */
/*----------------------------------------------------------------------------*/
namespace fhedo{

    /*------------------------------------------------------------------------*/
    /** \struct ParamsMark
     *  \brief  structure that gathers FHeDo boolean marks used during the
     *          different algorithms
     */
    struct EXPORT_GMDS ParamsMark{

        /** node mark for nodes classified on surfaces */
        int mark_node_on_surf;
        /** node mark for nodes classified on curves */
        int mark_node_on_curv;
        /** node mark for nodes classified on points*/
        int mark_node_on_pnt;
        /** node mark for isolated nodes (connected to no cells) */
        int mark_node_isolated;
        /** edge mark for edges classified on surfaces */
        int mark_edge_on_surf;
        /** edge mark for edges classified on curves */
        int mark_edge_on_curv;
        /** face mark for faces classified on surfaces */
        int mark_face_on_surf;
    };
    /*------------------------------------------------------------------------*/
    /** \struct ParamsGlobal
     *  \brief  structure that gathers FHeDo global parameters
     */
    struct EXPORT_GMDS ParamsGlobal{

        enum AlgoPhaseType{
            FF_GEN,     /** Frame field generation*/
            FF_SMOOTH,  /** Frame field smoothing*/
            PF_GEN,     /** Point field generation*/
            HEX_DOM,    /** Hexahedral  generation*/
            HEX_CAV,    /** Cavities filling*/
        };

        /** last phase we do */
        AlgoPhaseType stop_at;
        /** flag indicating that we want debug files to be generated */
        bool with_debug_files;
        /** Input mesh file name */
        std::string mesh_file;
        /** Output directory */
        std::string output_dir;
    };
    /*------------------------------------------------------------------------*/
    /** \struct ParamsFrameField
     *  \brief  structure that gathers parameters of the Frame Field generation
     */
    struct EXPORT_GMDS ParamsFrameField{
        enum SolverType{
            OPENNL=0,
            EIGEN=1
        };
        enum SmoothingType{
            RAY=0,
            LIU=1
        };

        /** solver used for building the frame field */
        SolverType solver_type;
        /** use cotangent_weight or not*/
        bool with_cotangent_weights;
        /** use smoothing or not*/
        bool with_smoothing;
        /** type of smoothing algorithm*/
        SmoothingType smoothing_algo;
        /** Number of smoothing iterations*/
        int smoothing_nb_iter;
        /** Convergence criteria smoothing steps*/
        double smoothing_epsilon;
        /** with or without tet mesh adaptation during FF*/
        bool with_mesh_adaptation;
        /** with or without premeshing around volume singularity lines*/
        bool premeshing_sing_lines;
        /** generate streamlines after frame field computation*/
        bool generate_streamlines;
    };
    /*------------------------------------------------------------------------*/
    /** \struct ParamsPointField
     *  \brief  structure that gathers parameters of the Point Field generation
     */
    struct EXPORT_GMDS ParamsPointField{
        /** Curl parameter for the PGP algorithm*/
        double curl;
        /** indicates if the user provides a spacing target value */
        bool with_user_spacing;
        /** provided spacing value*/
        double spacing;
    };
    /*------------------------------------------------------------------------*/
    /** \struct ParamsHexDom
     *  \brief  structure that gathers parameters of the Hex-dom generation
     */
    struct EXPORT_GMDS ParamsHexDom{
        /** Gives the method to create edge. With interpolation, Heuns' method
         *  is used, without only vector directions are taken into account*/
        bool with_edge_interpolation;
        /** Cone tolerance for connecting points to form edges*/
        double edge_cone_tolerance;
        /** Remaining surface patch are built with quad (triangles otherwise)*/
        bool with_quad_surface;
        /** We try to add hexahedral elements with an advancing front algo.*/
        bool with_whisker_weaving;
        /** Pyramids are added to connect hex and tet. Otherwise non-conform
         * meshes are generated */
        bool with_pyramids;

    };


    /*------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* FHEDO_PARAMS_H_ */
/*----------------------------------------------------------------------------*/
