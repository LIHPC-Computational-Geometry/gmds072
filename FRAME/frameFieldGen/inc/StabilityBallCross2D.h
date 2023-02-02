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
#ifndef STABILITY_BALL_CROSS_2D_H_
#define STABILITY_BALL_CROSS_2D_H_
/*----------------------------------------------------------------------------*/
#include <map>
#include <list>
#include <set>
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IGMesh.h>
#include <GMDS/Math/Cross2D.h>
/*----------------------------------------------------------------------------*/
/** \class  StabilityBallCross2D
 *  \brief  Propagate a cross field from crosses defined on the boundary 
 *          into the mesh
 */
class EXPORT_GMDS StabilityBallCross2D
{
 public:

  /*------------------------------------------------------------------------*/
  /** \brief Constructor.
   *
   *  \param AMesh     the mesh where we work on
   *  \param AField    the cross field to propagate
   *  \param AFrom     the nodes we start from
   *  \param ATo       the nodes we propagate through
   */
  StabilityBallCross2D(gmds::IGMesh* AMesh,
		       gmds::Variable<gmds::math::Cross2D>* AField,
		       std::vector<gmds::Node>& AFromNodes,
		       std::vector<gmds::Node>& AToNodes);
  
  /*------------------------------------------------------------------------*/
  /** \brief Function to be called for executing the cross field generation 
   *	     algorithm. 
   */
  void execute();
     
 private:  



  /*------------------------------------------------------------------------*/
  /** \brief Local structure to compute and store stability data
   */
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
  /** \brief Compute the distance field that will be used to define the
   *         insertion order of the advancing front process
   */
  void computeDistanceFields();
  /*------------------------------------------------------------------------*/
  /** \brief Initialize the list of alive nodes, i.e. the nodes that define
   *         the propagation front
   */
  void initAliveNodes();
	       
  /*------------------------------------------------------------------------*/
  /** \brief Initialize the narrow band, i.e. the list of nodes adjacent
   *         to at least 2 alive nodes and that are not alive
   */
   void initNarrowBand(std::list<gmds::Node>& ANarrowBand,
		       double&	ASmallestStab,
		       gmds::TCellID& ASmallestStabID);

  /*------------------------------------------------------------------------*/
  /** \brief Return all the nodes sharing a face with ANode
   *  
   * \param ANode the node we want the neighborhood
   *
   * \return the nodes sharing a face with ANode 
   */
   std::vector<gmds::Node> getAdjacentNodes(gmds::Node& ANode) const;

  /*------------------------------------------------------------------------*/
  /** \brief Check if ATo, AFrom1 and AFrom2 belongs to the same face
   *  
   * \param ATo    a node
   * \param AFrom1 a node
   * \param AFrom2 a node
   * \param AFace  the face that contains ATo, AFrom1 and AFrom2 if it exists. 
   *               It is an output param.
   *
   * \return true if ATo, AFrom1 and AFrom2 share a face, false otherwise.
   */
   bool belongToTheSameFace(gmds::Node& ATo   , gmds::Node& AFrom1,
			    gmds::Node& AFrom2, gmds::Face& AFace) const;

  /*------------------------------------------------------------------------*/
  /** \brief Extrapolate a cross value in ATo from the crosses already 
   *         defined in AFrom1 and AFrom2
   *  
   * \param ATo    a node
   * \param AFrom1 a node
   * \param AFrom2 a node
   */
   void extrapolateCross(gmds::Node& ATo, gmds::Node& AFrom1, 
			 gmds::Node& AFrom2);
   /*------------------------------------------------------------------------*/
   /** \brief Extrapolate a cross value in ATo from the crosses already 
    *         defined in AFrom2
    *  
    * \param ATo    a node
    * \param AFrom  a vector of nodes
    */
   void extrapolateCross(gmds::Node& ATo, std::vector<gmds::Node>& AFrom);

   /*------------------------------------------------------------------------*/
   /** \brief Initialize the distance ball for each node
    */
   void initDistanceBalls();

   /*------------------------------------------------------------------------*/
   /** \brief Iterative process to assign a cross to every node of m_mesh
    */
   void propagateCrosses(std::list<gmds::Node>& ANarrowBand,
			 double& ASmallestStab,
			 gmds::TCellID& ASmallestStabID);

   /*------------------------------------------------------------------------*/
   /** \brief Local smoothing restricted to nodes marked with AMark
    *
    * \param AMark a mark on nodes
    */
   void smooth(const int AMark);


   /*------------------------------------------------------------------------*/
   /** \brief Compute the stability data for AN considering the nodes located
    *         in the distance ball of AN and marked with AMark
    *
    * \param ANode the node to work on
    * \param AMark a mark on nodes
    */
   virtual double computeStability(gmds::Node& ANode, const int AMark);

   /*------------------------------------------------------------------------*/
   /** \brief Recompute the stability data for AN considering the nodes located
    *         in the distance ball of AN and marked with AMark
    *
    * \param ANode the node to work on
    * \param AMark a mark on nodes
    */
   virtual void updateStabilityInfo(gmds::Node& AN, const int AMark);

   /*------------------------------------------------------------------------*/
   /** \brief Compute a stability data from a set of contributions
    *
    * \param AContributions set of contribution we want to extract stability
    *                       data from.
    */
   virtual double computeStability(std::set<Contribution>& AContributions);

   void putInBalls(gmds::Node& AN1, gmds::Node& AN2);

 private:
   /** Background triangular mesh */
   gmds::IGMesh* m_mesh; 
  
   /** Cross field that we build on the mesh (node association) */
   gmds::Variable<gmds::math::Cross2D>* m_field;

   /*** set of nodes we start from */
   std::vector<gmds::Node> m_from_nodes;

   /*** set of nodes we must defined a cross on */
   std::vector<gmds::Node> m_to_nodes;

   // distance field that is going to be computed on m_mesh */
   gmds::Variable<double>* m_distance_field; 
	  

   /*** Boolean mark for nodes we start from */
   int m_mark_from;
   /*** Boolean mark for nodes we go through */
   int m_mark_to;
   /*** Boolean mark for nodes where a cross is already defined */
   int m_mark_alive;
   /*** Boolean mark for nodes that are candidates to get a defined cross*/
   int m_mark_narrow;


   /*** All the nodes that are used to compute the smoothing value in m_dist_ball[i] */
   std::map<gmds::TCellID, std::vector<gmds::TCellID> > m_ball; 
   /*** All the nodes that are used to compute the smoothing value in m_dist_ball[i] */
   std::map<gmds::TCellID, std::map<gmds::TCellID, int> > m_ball_location; 
   /*** Corresponding distance to each node */
   std::map<gmds::TCellID, std::vector<double> > m_distance_ball; 
   /*** Corresponding smoothing contribution of each node */
   std::map<gmds::TCellID, std::vector<double> > m_smooth_ball;  

   /*** The stability value computed for each node */
  std::map<gmds::TCellID, double > m_stability_value;
};
/*----------------------------------------------------------------------------*/
#endif //ADVANCING_FRONT_FIELD_GEN_2D_H_
/*----------------------------------------------------------------------------*/
