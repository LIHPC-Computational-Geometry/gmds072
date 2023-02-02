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
#ifndef SMOOTHING_HLBFGS_H_
#define SMOOTHING_HLBFGS_H_
/*---------------------------------------------------------------------------*/
#include <iostream>
/*---------------------------------------------------------------------------*/
#include <GMDS/IG/IGMesh.h>
#include <GMDS/Math/Quaternion.h>
#include <GMDS/Math/Matrix.h>
/*---------------------------------------------------------------------------*/
class SmoothingHLBFGS
{

 public:

  // the algorithm i can be applied on a mesh where Quaternions are defined
  // for some nodes and outward normal is given for nodes on a boundary 
  // surface
  SmoothingHLBFGS(gmds::IGMesh* AMesh,
		  gmds::Variable<gmds::math::Quaternion>*& AField,
		  gmds::Variable<gmds::math::Vector>*& ANormal);

  /* All the node with AMark must be involved in the smoothing process
   */
  void selectNodes(const int AMark);

  // nodes classified on curves, surfaces and points must be marked
  // the option is to provide this mark or to recompute them
  // Only the first option is provided right now.
  void initBoundaryMarks(const int AMarkPnt, const int AMarkCurve,
			 const int AMarkSurf);

  void execute();
  void initCandidates();
  int  initOptimizationData(double*& eulerAngles);
  void rebuildQuaternions(double*& eulerAngles);

  static void evalF(int N, double* x, double *prev_x, double* func, double* grad);

  static void initMatrix(gmds::Node& ANode,
			 gmds::math::Matrix<3, 3,double>& Mix,
			 gmds::math::Matrix<3, 3,double>& Miy,
			 gmds::math::Matrix<3, 3,double>& Miz,
			 gmds::math::Matrix<3, 3,double>& DMix,
			 gmds::math::Matrix<3, 3,double>& DMiy,
			 gmds::math::Matrix<3, 3,double>& DMiz,
			 std::vector<gmds::TCoord>& trig_cos_vectors,
			 std::vector<gmds::TCoord>& trig_sin_vectors);
  static void newiteration(int iter, int call_iter, double *x, double* f, double *g, double* gnorm);
    
 private:
  void finalSmoothing();
 private:
	    
  static gmds::IGMesh* m_mesh;
  gmds::Variable<gmds::math::Quaternion>* m_frame_field;
  gmds::Variable<gmds::math::Vector>* m_surf_normal;
  //mark to identify the nodes to work on
  int m_mark_candidates;
  // boolean indicating if the candidate mark has been defined
  bool m_candidate_mark_initialized;
  // mark to identify the nodes that are classified on points
  static int m_mark_point;
  // mark to identify the nodes that are classified on curves
  static int m_mark_curve;
  // mark to identify the nodes that are classified on surfaces
  static int m_mark_surface;
  // boolean indicating if the candidate mark has been defined
  bool m_boundary_marks_initialized;
  //we keep in mind the number of call to HLBFGS
  int nb_HLBFGS_calls;
      
  //the ordered list of candidate nodes
  static std::vector<gmds::Node> m_candidates;
  //the list of edges that connect candidates
  // these edges will be used to compute smoothing values between
  // end nodes
  static std::vector<gmds::Edge> m_edges;
  //the index of a candidate node in m_candidates
  static std::map<gmds::TCellID, int> m_candidates_index;
  //We keep in mind a chart for each boundary node we work on
  static std::map<gmds::TCellID, gmds::math::Chart> m_boundary_triad;
};
/*---------------------------------------------------------------------------*/
#endif //SMOOTHING_HLBFGS_H_
/*---------------------------------------------------------------------------*/
