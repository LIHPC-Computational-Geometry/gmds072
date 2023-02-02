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
/*---------------------------------------------------------------------------*/
#include <iostream>
/*---------------------------------------------------------------------------*/
#include <GMDS/Algo/BoundaryOperator.h>
#include <GMDS/Algo/DistanceFieldBuilder2D.h>
/*---------------------------------------------------------------------------*/
#include <GMDS/IO/VTKWriter.h>
#include <sstream>
/*---------------------------------------------------------------------------*/
#include "StabilityBallCross2D.h"
/*---------------------------------------------------------------------------*/
using namespace gmds;
/*---------------------------------------------------------------------------*/
StabilityBallCross2D::
StabilityBallCross2D(gmds::IGMesh* AMesh,
		      gmds::Variable<gmds::math::Cross2D>* AField,
		      std::vector<gmds::Node>& AFromNodes,
		      std::vector<gmds::Node>& AToNodes)
  : m_mesh(AMesh), m_field(AField), m_from_nodes(AFromNodes),
    m_to_nodes(AToNodes)
{ 
  m_mark_from = m_mesh->getNewMark<Node>();
  m_mark_to   = m_mesh->getNewMark<Node>();
  for(unsigned int i=0; i <m_from_nodes.size(); i++)
    m_mesh->mark(m_from_nodes[i],m_mark_from);
  
  for(unsigned int i=0; i <m_to_nodes.size(); i++)
    m_mesh->mark(m_to_nodes[i],m_mark_to);
  
  m_distance_field = 0;
}
/*---------------------------------------------------------------------------*/
void StabilityBallCross2D::execute()
{
  m_mark_alive = m_mesh->getNewMark<Node>();
  m_mark_narrow = m_mesh->getNewMark<Node>();

  //==================================================================
  //STEP 1- Computes the distance field(s)
  //==================================================================
  std::cout << "======================================"<< std::endl;
  std::cout << "Computation of distance fields" << std::endl;
  computeDistanceFields(); //CHILDREN'S DELEGATED BEHAVIOUR
  
  std::cout << "    DONE" << std::endl;
  

  //==================================================================
  //STEP 2 - buid and compute the distance to the approximated medial
  // axis - NOT DONE NOW
  //==================================================================
  //std::vector<gmds::Node>  nodes_of_approximate_MA =
  //	computeApproximateMedialAxis(surf_nodes);
  //m_medial_axis_distance_field =
  // 	dist_computer.computeDistance(nodes_of_approximate_MA,
  //			surface_nodes_ids);


  //==================================================================
  //STEP 3 - Marks nodes where a quaternion is defined on
  //==================================================================
  std::cout << "======================================"<< std::endl;
  std::cout << "Alive nodes initialization" << std::flush;
  initAliveNodes();
  std::cout << "    DONE" << std::endl;
  //=================================================================
  // Step 0 - We initialize the distance ball only for surf elements
  //==================================================================
  std::cout << std::endl << "2D init distance balls" << std::endl;
  initDistanceBalls();

      //==================================================================
  // Now, we assign a cross to every  nodes !!!
  // We follow the same process than the one used to build the distance
  // field:
  // 1)a list of candidates is computed (a front, or narrow band)
  // 2) for each candidate we compute a stability value
  // 3) We make alive the candidate with the smallest stab. value. Let
  //    C this candidate
  // 4) C is removed from the narrow band and some of its neighbors
  //    are added in the narrow band.
  // Step (3) and (4) are repeated while the narrow band list is 
  // not empty. 
  //==================================================================
  // Step 4 - Initialization of the narrow band list
  //==================================================================
  // A node of the narrow band must be adjacent to 2 alive nodes, 
  // while being not alive.
  std::cout << "======================================"<< std::endl;
  std::cout << "Init Narrow Band" << std::flush;
  std::list<Node> narrow_band;
  // we keep in mind the smallest stability encountered
  TCellID smallest_stab_node_id = NullID;
  double smallest_stab = 10000000;

  initNarrowBand(narrow_band, smallest_stab, smallest_stab_node_id);
  std::cout << "    DONE" << std::endl;


  //======================================================================
  // Step 5 - Advancing front loop
  //======================================================================
  std::cout << "======================================"<< std::endl;
  std::cout << " Advancing-front progress of the narrow band" << std::flush;
  propagateCrosses(narrow_band, smallest_stab, smallest_stab_node_id);
  
  std::cout << "    DONE" << std::endl;
  
  std::cout << "    DONE" << std::endl;

  //==================================================================
  // CLEANING - Boolean marks are cleaned
  //==================================================================
  m_mesh->freeMark<Node>(m_mark_alive);
  m_mesh->freeMark<Node>(m_mark_narrow);

  m_mesh->unmarkAll<Node>(m_mark_from);
  m_mesh->unmarkAll<Node>(m_mark_to  );
  m_mesh->freeMark<Node> (m_mark_from);
  m_mesh->freeMark<Node> (m_mark_to  );
}
/*---------------------------------------------------------------------------*/
void StabilityBallCross2D::initDistanceBalls()
{
  std::vector<Node> all_nodes;
  all_nodes.insert(all_nodes.end(), m_from_nodes.begin(), m_from_nodes.end());
  all_nodes.insert(all_nodes.end(), m_to_nodes.begin()  , m_to_nodes.end()  );

  for (unsigned int i = 0; i < all_nodes.size(); i++) {
    Node ni = all_nodes[i];
    m_ball[ni.getID()].clear();
    m_ball_location[ni.getID()].clear();
    m_distance_ball[ni.getID()].clear();
    m_smooth_ball[ni.getID()].clear();
  }
  return;
  for (unsigned int i = 0; i < all_nodes.size(); i++) {
    Node ni = all_nodes[i];
    math::Point pi = ni.getPoint();
    double di = (*m_distance_field)[ni.getID()];
    /*	if (di>m_max_distance_2D_div2)
	di = m_max_distance_2D_div2;*/
     
    double di2 = di*di;
    int ball_size = 0;
    for (unsigned int j = i + 1; j < all_nodes.size(); j++) {
      //i!=j by definition
      Node nj = all_nodes[j];
      math::Point pj = nj.getPoint();
      double dj = (*m_distance_field)[nj.getID()];
      /*	if (dj>m_max_distance_2D_div2)
		dj = m_max_distance_2D_div2;
      */
      double dj2 = dj*dj;
      math::Vector vij(pi, pj);
      double d = (di2 > dj2) ? di2 : dj2;
      if (vij.norm2() <= d)
	{
	  /* As each point has its own visibility ball, a point P1 can
	   * see P2 while P2 does not see P1.
	   * The current function corrects that by making visibility a
	   * one-to-one mapping
	   */
	  putInBalls(ni, nj);
	  ball_size++;
	}
    }
    // All the nodes can be too far from ni. In this case, we put adjacent
    //nodes in the ball (here adjacent nodes sharing a boundary face
    if (ball_size == 0) {
      std::vector<Node> adj_nodes = getAdjacentNodes(ni);
      for (unsigned int k = 0; k < adj_nodes.size(); k++) {
	putInBalls(ni, adj_nodes[k]);
				
      }
    }
		      
  }
}
/*----------------------------------------------------------------------------*/
void StabilityBallCross2D::putInBalls(Node& AN1, Node& AN2)
{
  math::Point p1 = AN1.getPoint();
  math::Point p2 = AN2.getPoint();
  math::Vector v12(p1, p2);
  //each node is in the ball of the other...
  m_ball[AN1.getID()].push_back(AN2.getID());
  m_ball[AN2.getID()].push_back(AN1.getID());

  m_ball_location[AN1.getID()][AN2.getID()] = m_ball[AN1.getID()].size() - 1;
  //m_ball_location[AN2.getID()][AN1.getID()] = m_ball[AN2.getID()].size() - 1;
  //we keep their distance ...
  double distance = v12.norm();
  m_distance_ball[AN1.getID()].push_back(distance);
  m_distance_ball[AN2.getID()].push_back(distance);

  // Warning, we initialize the contribution of AN1 to AN2 to -10! It is just a
  // way to field the values at the initialization, but this value is
  // ill-defined. It will have to be updated when a quaternion will be assigned
  // to AN1 and AN2.

  m_smooth_ball[AN1.getID()].push_back(-10);
  m_smooth_ball[AN2.getID()].push_back(-10);

}
/*---------------------------------------------------------------------------*/
void StabilityBallCross2D::computeDistanceFields()
{

  gmds::DistanceFieldBuilder2D distanceBuilder(m_mesh);
  std::cout << "Computation of a 2D distance field" << std::endl;
  if (!distanceBuilder.isValid())
    {
      std::cout << "Invalid model for distance computation" << std::endl;
      throw GMDSException("Invalid model for distance computation");
    }
  m_distance_field = distanceBuilder.computeDistance(m_from_nodes, 
						     m_to_nodes,
						     m_mark_to);
}
/*---------------------------------------------------------------------------*/
void StabilityBallCross2D::initAliveNodes()
{
  IGMesh::node_iterator it_node = m_mesh->nodes_begin();
  for (; !it_node.isDone(); it_node.next())
    {
      Node current_node = it_node.value();

      if (m_mesh->isMarked(current_node, m_mark_from))
	m_mesh->mark(current_node, m_mark_alive);
    }

}
/*---------------------------------------------------------------------------*/
void StabilityBallCross2D::
propagateCrosses(std::list<gmds::Node>& ANarrowBand,
		 double& ASmallestStab,
		 TCellID& ASmallestStabID)
{
  //  double step=0;
  while (!ANarrowBand.empty())
    {
      /*  step++;
	  if(step=20){
	  
	  step=0;
	  }*/
      //we take the node with the smallest distance computed in
      Node trial_node = m_mesh->get<Node>(ASmallestStabID);
      //the cross extrapolated in trial_node has been already assigned
      //so, nothing to do...? NO, as trial_node becomes alive, all the nodes
      // in his ball must update their stability info.
      std::vector<TCellID> trial_ball = m_ball[trial_node.getID()];
      for (unsigned int i = 0; i < trial_ball.size(); i++)
	{
	  Node current_node = m_mesh->get<Node>(trial_ball[i]);
	  updateStabilityInfo(current_node, m_mark_alive);
	}

      //the trial node becomes alive ...
      m_mesh->mark(trial_node, m_mark_alive);

      //... and is removed from the narrow band
      m_mesh->unmark(trial_node, m_mark_narrow);
      ANarrowBand.remove(trial_node);

      //adjacent nodes that are not alive or in the narrow band are added
      //in the narrow band if they are now adjacent to 2 alive nodes
      std::vector<Node> trial_adj_nodes =
	getAdjacentNodes(trial_node);

      for (unsigned int j = 0; j < trial_adj_nodes.size(); j++)
	{
	  Node n_j = trial_adj_nodes[j];
	  // we look for nodes that are not yet in the narrow band or alive
	  if (m_mesh->isMarked(n_j, m_mark_narrow) ||
	      m_mesh->isMarked(n_j, m_mark_alive))
	    continue;
	  //we check the mark of the nodes adjacent to n_j
	  std::vector<Node> adj_nodes_j = getAdjacentNodes(n_j);
	  std::vector<Node> ext_nodes;
	  int nb_alive_adj_nodes = 0;
	  for (unsigned int k = 0; k < adj_nodes_j.size(); k++)
	    {
	      Node n_k = adj_nodes_j[k];
	      if (m_mesh->isMarked(n_k, m_mark_alive))
		{
		  nb_alive_adj_nodes++;
		  ext_nodes.push_back(n_k);
		}
	    }

	  if (nb_alive_adj_nodes == 2)
	    {
	      // we have to be sure that n_j, ext_nodes[0], ext_nodes[1]
	      // belong to the same face. Iy yes we can put n_j in the 
	      // narrow band.
	      Face commonFace;
	      if (belongToTheSameFace(n_j, ext_nodes[0], ext_nodes[1], commonFace))
		{
		  ANarrowBand.push_back(n_j);
		  m_mesh->mark(n_j, m_mark_narrow);

		  // We can now extrapolate a quaternion in this node... 
		  extrapolateCross(n_j, ext_nodes[0], ext_nodes[1]);
		  // then compute a stability value in n_j
		}
	    }
	  else if (nb_alive_adj_nodes > 2)
	    {
	      ANarrowBand.push_back(n_j);
	      m_mesh->mark(n_j, m_mark_narrow);

	      // We can now extrapolate a quaternion in this node... 
	      extrapolateCross(n_j, ext_nodes);
	    }
	}//for (unsigned int k = 0; k < adj_nodes_j.size(); k++)

      //now we select the next trial element
      ASmallestStab = 10000000;
      for (std::list<Node>::const_iterator it_nodes = ANarrowBand.begin();
	   it_nodes != ANarrowBand.end(); it_nodes++)
	{
	  Node current_node = *it_nodes;

	  double local_stab = m_stability_value[current_node.getID()];

	  if (local_stab < ASmallestStab)
	    {
	      ASmallestStab = local_stab;
	      ASmallestStabID = current_node.getID();
	    }
	}

    }
}
/*----------------------------------------------------------------------------*/
std::vector<Node> StabilityBallCross2D::
getAdjacentNodes(Node& ANode) const
{
  std::set<Node> adj_nodes;
  std::vector<Face>  adj_faces = ANode.get<Face>();

  for (unsigned int j = 0; j < adj_faces.size(); j++)
    {
      Face f_j = adj_faces[j];
      std::vector<Node> nodes_fj = f_j.get<Node>();
      for (unsigned int j = 0; j < nodes_fj.size(); j++)
	{
	  Node n_j = nodes_fj[j];
	  if (n_j.getID() != ANode.getID())
	    {
	      adj_nodes.insert(n_j);
	    }
	}
    }
  std::vector<Node> adj;
  adj.insert(adj.end(), adj_nodes.begin(), adj_nodes.end());

  return adj;
}
/*---------------------------------------------------------------------------*/
bool StabilityBallCross2D::belongToTheSameFace(gmds::Node& ATo,
						gmds::Node& AFrom1,
						gmds::Node& AFrom2,
						gmds::Face& AFace) const
{
  std::vector<TCellID> common_faces_12;
  std::vector<TCellID> common_faces_123;
  std::vector<TCellID> fs1 = ATo.getIDs<Face>();
  std::vector<TCellID> fs2 = AFrom1.getIDs<Face>();
  for (unsigned int i = 0; i < fs1.size(); i++)
    {
      TCellID id1 = fs1[i];
      bool found = false;
      for (unsigned int j = 0; !found && j < fs2.size(); j++)
	{
	  TCellID id2 = fs2[j];
	  if (id1 == id2)
	    {
	      found = true;
	      common_faces_12.push_back(id1);
	    }
	}
    }
  fs2 = AFrom2.getIDs<Face>();

  for (unsigned int i = 0; i < common_faces_12.size(); i++)
    {
      TCellID id1 = common_faces_12[i];
      bool found = false;
      for (unsigned int j = 0; !found && j < fs2.size(); j++)
	{
	  TCellID id2 = fs2[j];
	  if (id1 == id2)
	    {
	      found = true;
	      common_faces_123.push_back(id1);
	    }
	}
    }
  if (common_faces_123.size() == 1){
    AFace = m_mesh->get<Face>(common_faces_123[0]);
    return true;
  }
  return false;
}
/*----------------------------------------------------------------------------*/
void StabilityBallCross2D::initNarrowBand(std::list<gmds::Node>& ANarrowBand,
					   double& ASmallestStab,
					   TCellID& ASmallestStabID)
{
  IGMesh::node_iterator it_node = m_mesh->nodes_begin();
  for (; !it_node.isDone(); it_node.next())
    {
      Node current_node = it_node.value();
      if (m_mesh->isMarked(current_node, m_mark_alive))
	{
	  //We get adjacent nodes that are on the surface
	  std::vector<Node> adj_nodes = getAdjacentNodes(current_node);

	  for (unsigned int j = 0; j < adj_nodes.size(); j++)
	    {
	      Node n_j = adj_nodes[j];
	      // we look for nodes that are not yet in the narrow band or alive
	      if (m_mesh->isMarked(n_j, m_mark_narrow) ||
		  m_mesh->isMarked(n_j, m_mark_alive))
		continue;

	      std::vector<Node> adj_nodes_j = getAdjacentNodes(n_j);

	      std::vector<Node> ext_nodes;
	      int nb_alive_adj_nodes = 0;
	      for (unsigned int k = 0; k < adj_nodes_j.size(); k++)
		{
		  Node n_k = adj_nodes_j[k];
		  if (m_mesh->isMarked(n_k, m_mark_alive))
		    {
		      nb_alive_adj_nodes++;
		      ext_nodes.push_back(n_k);
		    }
		}
	      if (nb_alive_adj_nodes == 2)
		{
		  // we have to be sure that n_j, ext_nodes[0], ext_nodes[1]
		  // belong to the same face. Iy yes we can put n_j in the 
		  // narrow band.
		  Face commonFace;
		  if (belongToTheSameFace(n_j, ext_nodes[0], ext_nodes[1], commonFace))
		    {
		      ANarrowBand.push_back(n_j);
		      m_mesh->mark(n_j, m_mark_narrow);
		      // We can now extrapolate a quaternion in this node... 
		      extrapolateCross(n_j, ext_nodes[0], ext_nodes[1]);
		      // then compute a stability value in n_j
		      double local_stab = computeStability(n_j, m_mark_alive);
		      m_stability_value[n_j.getID()] = local_stab;
		      if (local_stab < ASmallestStab)
			{
			  ASmallestStab = local_stab;
			  ASmallestStabID = n_j.getID();
			}
		    }

		}//if (nb_alive_adj_nodes == 2)

	    }//for (unsigned int j = 0; j < adj_nodes.size(); j++)

	}//if (m_mesh->isMarked(current_node, AMarkAlive))

    }//for (; !it_node.isDone(); it_node.next())

  std::cout << "Narrow Band Size: " << ANarrowBand.size() << std::endl;
}
/*----------------------------------------------------------------------------*/
void StabilityBallCross2D::smooth(const int AMark)
{

  //	SmoothingHLBFGS smoother(m_mesh, m_cross_field_2D, m_surf_normal);
  /*	FrameFieldLaplacianSmoothing smoother(m_mesh, m_cross_field_2D, m_surf_normal);
	smoother.initBoundaryMarks(m_markNodeOnPnt, m_markNodeOnCurv, m_markNodeOnSurf);
	smoother.selectNodes(AMark);
	smoother.execute();
  */
}
/*----------------------------------------------------------------------------*/
void StabilityBallCross2D::extrapolateCross(gmds::Node& ATo, 
					     gmds::Node& AFrom1, 
					     gmds::Node& AFrom2)
{

  //===================================================================
  //We get quaternions
  //===================================================================
  math::Cross2D c1 = (*m_field)[AFrom1.getID()];
  math::Cross2D c2 = (*m_field)[AFrom2.getID()];

  //===================================================================
  //We get coeff (inverse of distance)
  //===================================================================
  std::vector<TCoord> cs;
  math::Point p0 = ATo.getPoint();
  math::Point p1 = AFrom1.getPoint();
  math::Point p2 = AFrom2.getPoint();

  math::Vector v01(p0, p1);
  math::Vector v02(p0, p2);
  TCoord w1 =1.0/v01.norm();
  TCoord w2 =1.0/v02.norm();

  //===================================================================
  // extrapolated cross
  //===================================================================
  math::Cross2D c  = c1 + c2;

  (*m_field)[ATo.getID()] = c;

}
/*----------------------------------------------------------------------------*/
void StabilityBallCross2D::extrapolateCross(gmds::Node& ATo, 
						std::vector<gmds::Node>& AFrom)
{
  math::Point p0 = ATo.getPoint();
  
  //===================================================================
  //We get crosses and  coeff (inverse of distance)
  //===================================================================
  std::vector<math::Cross2D> from_crosses;
  std::vector<TCoord> from_weights;
  for (unsigned int i = 0; i < AFrom.size(); i++) {
    Node current = AFrom[i];
    math::Cross2D current_c = (*m_field)[current.getID()];
    from_crosses.push_back(current_c);
    //coeef
    math::Point p1 = current.getPoint();
    math::Vector v01(p0, p1);
    from_weights.push_back(1.0/v01.norm());
  }

  //===================================================================
  //2D cross extrapolated
  //===================================================================
  math::Cross2D c_mean = math::Cross2D::mean(from_crosses,from_weights);


  //===================================================================
  //We assign the 2D cross on the mesh node
  //===================================================================
  (*m_field)[ATo.getID()] = c_mean;
}

/*----------------------------------------------------------------------------*/
void StabilityBallCross2D::updateStabilityInfo(Node& ANode, 
						   const int AMark)
{
  double stability = 100000;
  //ids of the node, which are in the distance ball of ANode
  std::vector<TCellID>& ball = m_ball[ANode.getID()];
  std::map<TCellID, int>& ball_location = m_ball_location[ANode.getID()];
  std::vector<double>& distances = m_distance_ball[ANode.getID()];
  std::vector<double>& angles = m_smooth_ball[ANode.getID()];

  //==================================================================
  // Among the nodes in the distance ball, we keep the ones where 
  // a Quaternion is defined (i.e. marked with AMark)
  //==================================================================
  std::vector<TCellID> defined_node_ids;
  std::vector<math::Quaternion> defined_quat;
  for (unsigned int i = 0; i < ball.size(); i++)
    {
      Node ni = m_mesh->get<Node>(ball[i]);
      if (m_mesh->isMarked(ni, AMark))
	{
	  defined_node_ids.push_back(ball[i]);
	}
    }
  //==================================================================
  // Now we order by growing distance the available contributions
  //==================================================================

  std::set<Contribution> ordered_contrib;
  for (unsigned int i = 0; i < defined_node_ids.size(); i++)
    {
      int location_in_ball = ball_location[defined_node_ids[i]];
      Contribution ci;
      ci.distance = distances[location_in_ball];
      ci.angle = angles[location_in_ball];
    }

  //==================================================================
  //now we compute the stability value
  //==================================================================
  computeStability(ordered_contrib);
}

/*----------------------------------------------------------------------------*/
double StabilityBallCross2D::
computeStability(Node& ANode, const int AMark) {
  //ids of the node, which are in the distance ball of ANode
  std::vector<TCellID>& ball      = m_ball[ANode.getID()];
  std::vector<double >& distances = m_distance_ball[ANode.getID()];
    //  std::vector<double >& angles    = m_smooth_ball[ANode.getID()];

  //as we compute the stability of ANode, a quaternion is defined on it
  math::Cross2D cross_ref = (*m_field)[ANode.getID()];
  //==================================================================
  // Among the nodes in the distance ball, we keep the ones where 
  // a quaternion is defined (i.e. marked with AMark)
  //==================================================================
  std::vector<TCellID> defined_node_ids;
  std::vector<math::Cross2D> defined_crosses;
  for (unsigned int i = 0; i < ball.size(); i++) {
    Node ni = m_mesh->get<Node>(ball[i]);
    if (m_mesh->isMarked(ni, AMark)) {
      defined_node_ids.push_back(i);
    }
  }
  //==================================================================
  // Now we order by growing distance the available contributions
  //==================================================================

  std::set<Contribution> ordered_contrib;
  for (unsigned int i = 0; i < defined_node_ids.size(); i++) {
    Contribution ci;
    ci.distance = distances[defined_node_ids[i]];

      // a defined node is alive, so a quaternion is associated to it
    math::Cross2D cross_i = (*m_field)[defined_node_ids[i]];
    ci.angle = cross_ref.angle(cross_i);
    //angles[defined_node_ids[i]];
    ordered_contrib.insert(ci);
  }

  //==================================================================
  //now we compute the stability value
  //==================================================================
  double stab = computeStability(ordered_contrib);
  return stab;
}
/*----------------------------------------------------------------------------*/
double  StabilityBallCross2D::
computeStability(std::set<Contribution>& AContributions)
{
  double res = 0.0;
  double lastDist = 0.0;
  double delta = 0.0;

  std::set<Contribution>::const_iterator it = AContributions.begin();
  const std::set<Contribution>::const_iterator it_end = AContributions.end();
  while (it != it_end) {
    delta = (*it).angle;
    res += delta;// *((*it).distance - lastDist); //NB : on a tj *it.dist >= lastDist
    lastDist = (*it).distance;
    it++;
  }

  // This information should be scaled to the ball size, or in a discrete 
  // manner, to the number of contributors??? Otherwise, we could not compare two
  // stability data. The choice is not so obvious. After all we want to introduce
  // quaternions in thin area first. And the one in the middle object. 

  return res;
}
/*----------------------------------------------------------------------------*/
