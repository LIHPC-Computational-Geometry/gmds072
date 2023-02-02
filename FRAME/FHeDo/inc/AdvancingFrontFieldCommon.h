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
#ifndef ADVANCING_FRONT_FIELD_GEN_COMMON_H_
#define ADVANCING_FRONT_FIELD_GEN_COMMON_H_
/*----------------------------------------------------------------------------*/
#include <iostream>
#include <set>
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IGMesh.h>
#include <GMDS/Math/Quaternion.h>
/*----------------------------------------------------------------------------*/
class EXPORT_GMDS AdvancingFrontFieldCommon 
{
 public:

  /*------------------------------------------------------------------------*/
  /** \brief Function to be called for executing a frame field generation 
   *		   algorithm. This function follows the template method DP. It
   *		   means that it drives the global execution algorithm, each 
   *		   specific part being delegated to child classes.
   */
  virtual void execute()=0;

  /*------------------------------------------------------------------------*/
  /** \brief Function to give the directory where we want to put the output
   *		   files (vtk files)
   */
  void setDebugDirectory(const std::string& ADirName);
 protected:


  struct Contribution{
    double distance;
    double angle;

    /*-----------------------------------------------------------------*/
    /** \brief Comparison operator < between two contributions
     */
    friend bool operator<(const Contribution& AC1,
			  const Contribution& AC2){
      return(AC1.distance < AC2.distance);
    }
  };

  /*------------------------------------------------------------------------*/
  /** \brief Constructor.
   *
   *  \param AMesh the mesh where we work on
   */
  AdvancingFrontFieldCommon(gmds::IGMesh* AMesh);

  /* mark all the cells on the boundary using m_markXXX attributes*/
  void markBoundaryCells();
  // Compute the quaternions that are put on the node classified on curves
  void initQuaternionOnCurves();
  // Compute the quaternions that are put on the node classified on points
  void initQuaternionOnPoints();

  std::vector<gmds::Edge> getEdgesOnCurve(gmds::Node& ANode) const;
  std::vector<gmds::Face> getFacesOnSurface(gmds::Edge& AEdge) const;
  gmds::Node  getNeighboorOn(gmds::Node& AN, gmds::Edge& AE) const;


  gmds::math::Vector getOutputNormal(gmds::Face& AFace, gmds::Region& ARegion);
  gmds::math::Vector getOutputNormalOfABoundaryFace(gmds::Face& AFace);

  std::vector<gmds::Node> getAdjacentNodes(gmds::Node& ANode, const int AMark);
  std::vector<gmds::Node> getAdjacentNodesByEdge(gmds::Node& ANode, const int AMark);
  /* Computes the smoothing value of a node
   */
  virtual double computeStability(gmds::Node& ANode, const int AMark);
  virtual void updateStabilityInfo(gmds::Node& AN, const int AMark);
  virtual double computeStability(std::set<Contribution>& AContributions);
	
  virtual void writeForDebug(const std::string AFileName="");


  void initNormalOnBoundarySurfaceNodes();
  void extrapolateQuaternion(gmds::Node& ATo, 
			     gmds::Node& AFrom1, 
			     gmds::Node& AFrom2,
			     gmds::Face& AFace);

  void extrapolateQuaternion(gmds::Node& ATo, std::vector<gmds::Node>& AFrom);
  void extrapolateQuaternionOnSurf(gmds::Node& ATo, std::vector<gmds::Node>& AFrom);
  void extrapolateQuaternionInVol(gmds::Node& ATo, std::vector<gmds::Node>& AFrom);

  bool belongToTheSameFace(
			   gmds::Node& AN1,
			   gmds::Node& AN2,
			   gmds::Node& AN3,
			   gmds::Face& AFace);

  //for the advancing front algorithm, nodes on curves and points are 
  //initialized as being alive as they have an already defined
  //quaternion.
  void initAliveNodes();
  void cleanMarks();
  void initMarks();

  inline double quaternionAngle(gmds::Node& AN1, gmds::Node& AN2);

  virtual void buildFrameField() = 0;
  virtual void computeDistanceFields()=0;
  virtual void initQuaternionsBySnapping() = 0;

  virtual void initDistanceBalls() { ; }

  virtual void smoothAll()=0;
    
  // Simplices can be colored using a variable. This behavior is
  // let to children classes.
  virtual void colorSimplices();
    
    
 protected:
  gmds::IGMesh* m_mesh; //background mesh
  std::string m_output_directory_name;

  double m_epsilon_surface;
  gmds::Variable<gmds::math::Quaternion>* m_cross_field;
  gmds::Variable<gmds::math::Vector>* m_surf_normal;

  std::map<gmds::TCellID, std::vector<gmds::TCellID> > m_ball; //all the nodes that are used to compute the smoothing value in m_dist_ball[i]
  std::map<gmds::TCellID, std::map<gmds::TCellID, int> > m_ball_location; //all the nodes that are used to compute the smoothing value in m_dist_ball[i]
  std::map<gmds::TCellID, std::vector<double> > m_distance_ball; //corresponding distance to each node
  std::map<gmds::TCellID, std::vector<double> > m_smooth_ball;  //corresponding smoothing contribution of each node

  std::map<gmds::TCellID, double > m_stability_value;


  int m_markNodeOnBnd;
  int m_markNodeOnSurf;
  int m_markNodeOnCurv;
  int m_markNodeOnPnt;
  int m_markEdgeOnCurv;
  int m_markEdgeOnSurf;
  int m_markFaceOnSurf;
  int m_markIsolated;

  // boolean marks used to perform the advancing front algorithm
  int m_mark_alive;
  int m_mark_narrow;

  std::vector<gmds::Node> m_boundary_nodes; //all the nodes on the boundary (surf, curve, point)
  std::vector<gmds::Node> m_curve_nodes; //all the nodes on a curve or a point
  std::vector<gmds::Node> m_surf_nodes; //all the nodes on a surface (not curve, not point)
};
/*----------------------------------------------------------------------------*/

#endif //ADVANCING_FRONT_FIELD_GEN_COMMON_H_
