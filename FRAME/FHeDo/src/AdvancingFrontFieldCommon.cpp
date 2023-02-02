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
#include <sstream>
#include <set>
/*---------------------------------------------------------------------------*/
#include <GMDS/Algo/BoundaryOperator.h>
#include <GMDS/Algo/DistanceFieldBuilder2D.h>
/*---------------------------------------------------------------------------*/
#include <GMDS/IO/VTKWriter.h>
/*---------------------------------------------------------------------------*/
#include "AdvancingFrontFieldCommon.h"
#include "SmoothingHLBFGS.h"
#include "FrameFieldLaplacianSmoothing.h"
/*---------------------------------------------------------------------------*/
using namespace gmds;
/*---------------------------------------------------------------------------*/
AdvancingFrontFieldCommon::AdvancingFrontFieldCommon(IGMesh* AMesh)
  : m_mesh(AMesh), m_output_directory_name("."), m_epsilon_surface(0.1)
{
  try{
    m_cross_field = m_mesh->getVariable<math::Quaternion>(GMDS_NODE, "quaternion");
  }
  catch (GMDSException& e){
    std::cout << e.what() << std::endl;
    std::cout << "\t -> generation of the quaternion variable on nodes"<<std::endl;
    m_cross_field = m_mesh->newVariable<math::Quaternion>(GMDS_NODE, "quaternion");

  }
  try{
    m_surf_normal = m_mesh->getVariable<math::Vector>(GMDS_NODE, "normal");
  }
  catch (GMDSException& e){
    std::cout << e.what() << std::endl;
    std::cout << "\t -> generation of the normal variable on nodes"<<std::endl;
    m_surf_normal = m_mesh->newVariable<math::Vector>(GMDS_NODE, "normal");

  }
}

/*---------------------------------------------------------------------------*/
void AdvancingFrontFieldCommon::setDebugDirectory(const std::string& ADirName)
{
  m_output_directory_name = ADirName;
}
/*---------------------------------------------------------------------------*/
void AdvancingFrontFieldCommon::colorSimplices()
{
  std::cout<<"No coloring scheme by default!!! This behaviour is delegated to child classes"<<std::endl;
}

/*---------------------------------------------------------------------------*/
void AdvancingFrontFieldCommon::initAliveNodes()
{
  IGMesh::node_iterator it_node = m_mesh->nodes_begin();
  for (; !it_node.isDone(); it_node.next())
    {
      Node current_node = it_node.value();

      //only nodes on curves are considered as being alive. In specific, 
      // nodes on geometric points must be ignored. Their associated 
      // quaternion having any geometric meaning.
      if (m_mesh->isMarked(current_node, m_markNodeOnCurv) ||
	  m_mesh->isMarked(current_node, m_markNodeOnPnt))
	m_mesh->mark(current_node, m_mark_alive);
    }

}
/*---------------------------------------------------------------------------*/
void AdvancingFrontFieldCommon::initNormalOnBoundarySurfaceNodes()
{
  IGMesh::node_iterator it_node = m_mesh->nodes_begin();
  for (; !it_node.isDone(); it_node.next())
    {
      Node current_node = it_node.value();
      if (m_mesh->isMarked(current_node, m_markNodeOnSurf) &&
	  !m_mesh->isMarked(current_node, m_markNodeOnCurv) &&
	  !m_mesh->isMarked(current_node, m_markNodeOnPnt))
	{
	  std::vector<Face> adj_faces = current_node.get<Face>();
	  std::vector<Face> adj_bnd_faces;
	  for (unsigned int i = 0; i < adj_faces.size(); i++)
	    {
	      Face fi = adj_faces[i];
	      if (m_mesh->isMarked(fi, m_markFaceOnSurf))
		adj_bnd_faces.push_back(fi);
	    }
	  math::Vector out_vec;
	  for (unsigned int i = 0; i < adj_bnd_faces.size(); i++)
	    {
	      Face fi = adj_bnd_faces[i];
	      out_vec = out_vec + getOutputNormalOfABoundaryFace(fi);
	    }
	  out_vec = out_vec / adj_bnd_faces.size();
	  (*m_surf_normal)[current_node.getID()] = out_vec;
	}
    }
}
/*---------------------------------------------------------------------------*/
void AdvancingFrontFieldCommon::markBoundaryCells()
{
  BoundaryOperator boundaryOp(m_mesh);
  if (!boundaryOp.isValid())
    {
      std::cout << "Invalid model for boundary operations" << std::endl;
      throw GMDSException("Invalid model for boundary operations");
    }



  boundaryOp.markCellOnGeometry(
				m_markFaceOnSurf, m_markEdgeOnSurf, m_markNodeOnSurf,
				m_markEdgeOnCurv, m_markNodeOnCurv, m_markNodeOnPnt,
				m_markIsolated);

  IGMesh::node_iterator it_node = m_mesh->nodes_begin();
  for (; !it_node.isDone(); it_node.next())
    {
      Node current_node = it_node.value();
      if (m_mesh->isMarked(current_node, m_markNodeOnSurf) ||
	  m_mesh->isMarked(current_node, m_markNodeOnCurv) ||
	  m_mesh->isMarked(current_node, m_markNodeOnPnt))
	{
	  m_mesh->mark(current_node, m_markNodeOnBnd);
	}
    }
}
/*---------------------------------------------------------------------------*/
void AdvancingFrontFieldCommon::initMarks()
{
  m_markNodeOnBnd = m_mesh->getNewMark<Node>();
  m_markNodeOnSurf = m_mesh->getNewMark<Node>();
  m_markNodeOnCurv = m_mesh->getNewMark<Node>();
  m_markNodeOnPnt = m_mesh->getNewMark<Node>();

  m_markEdgeOnSurf = m_mesh->getNewMark<Edge>();
  m_markEdgeOnCurv = m_mesh->getNewMark<Edge>();

  m_markFaceOnSurf = m_mesh->getNewMark<Face>();

  m_markIsolated = m_mesh->getNewMark<Node>();

  m_mark_alive = m_mesh->getNewMark<Node>();
  m_mark_narrow = m_mesh->getNewMark<Node>();
}
/*---------------------------------------------------------------------------*/
void AdvancingFrontFieldCommon::cleanMarks()
{
  m_mesh->unmarkAll<Node>(m_markNodeOnBnd);
  m_mesh->unmarkAll<Node>(m_markNodeOnSurf);
  m_mesh->unmarkAll<Node>(m_markNodeOnCurv);
  m_mesh->unmarkAll<Node>(m_markNodeOnPnt);
  m_mesh->unmarkAll<Node>(m_markIsolated);
  m_mesh->unmarkAll<Node>(m_mark_alive);
  m_mesh->unmarkAll<Node>(m_mark_narrow);

  m_mesh->unmarkAll<Edge>(m_markEdgeOnCurv);

  m_mesh->unmarkAll<Edge>(m_markEdgeOnSurf);
  m_mesh->unmarkAll<Face>(m_markFaceOnSurf);

  m_mesh->freeMark<Node>(m_markNodeOnBnd);
  m_mesh->freeMark<Node>(m_markNodeOnSurf);
  m_mesh->freeMark<Node>(m_markNodeOnCurv);
  m_mesh->freeMark<Node>(m_markNodeOnPnt);
  m_mesh->freeMark<Node>(m_markIsolated);
  m_mesh->freeMark<Node>(m_mark_alive);
  m_mesh->freeMark<Node>(m_mark_narrow);

  m_mesh->freeMark<Edge>(m_markEdgeOnCurv);
  m_mesh->freeMark<Edge>(m_markEdgeOnSurf);

  m_mesh->freeMark<Face>(m_markFaceOnSurf);

}
/*----------------------------------------------------------------------------*/
void AdvancingFrontFieldCommon::initQuaternionOnPoints()
{
    IGMesh::node_iterator it_node = m_mesh->nodes_begin();
    for (; !it_node.isDone(); it_node.next())
    {
        Node current_node = it_node.value();
        
        if (m_mesh->isMarked(current_node, m_markNodeOnPnt))
        {
            std::vector<Node> curve_nodes;
            std::vector<Edge> curve_edges = getEdgesOnCurve(current_node);
            for (unsigned int i = 0; i < curve_edges.size(); i++)
            {
                Edge ei = curve_edges[i];
                curve_nodes.push_back(getNeighboorOn(current_node, ei));
            }
            
            std::vector<math::Quaternion> adj_quat;
            std::vector<TCoord> adj_coef;
            adj_quat.resize(curve_nodes.size());
            adj_coef.resize(curve_nodes.size());
            for (unsigned int i = 0; i < curve_nodes.size(); i++)
            {
                adj_quat.push_back((*m_cross_field)[curve_nodes[i].getID()]);
                adj_coef.push_back(1);
            }
            (*m_cross_field)[current_node.getID()] =
            math::Quaternion::mean(adj_quat, adj_coef);
        }
        
    }
}
/*----------------------------------------------------------------------------*/
void AdvancingFrontFieldCommon::updateStabilityInfo(Node& ANode, const int AMark)
{
  //double stability = 100000;
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
double AdvancingFrontFieldCommon::
quaternionAngle(Node& AN1, Node& AN2)
{
  math::Quaternion q1 =(*m_cross_field)[AN1.getID()];
  math::Quaternion q2 =(*m_cross_field)[AN2.getID()];

  return q1.angle(q2);

}
/*----------------------------------------------------------------------------*/
double AdvancingFrontFieldCommon::
computeStability(Node& ANode, const int AMark)
{
  //ids of the node, which are in the distance ball of ANode
  std::vector<TCellID>& ball = m_ball[ANode.getID()];
  std::vector<double>& distances = m_distance_ball[ANode.getID()];
 // std::vector<double>& angles = m_smooth_ball[ANode.getID()];

  //as we compute the stability of ANode, a quaternion is defined on it
  math::Quaternion quat_ref = (*m_cross_field)[ANode.getID()];
  //==================================================================
  // Among the nodes in the distance ball, we keep the ones where 
  // a quaternion is defined (i.e. marked with AMark)
  //==================================================================
  std::vector<TCellID> defined_node_ids;
  std::vector<math::Quaternion> defined_quat;
  for (unsigned int i = 0; i < ball.size(); i++)
    {
      Node ni = m_mesh->get<Node>(ball[i]);
      if (m_mesh->isMarked(ni, AMark))
	{
	  defined_node_ids.push_back(i);
	}
    }
  //==================================================================
  // Now we order by growing distance the available contributions
  //==================================================================

  std::set<Contribution> ordered_contrib;
  for (unsigned int i = 0; i < defined_node_ids.size(); i++)
    {
      Contribution ci;
      ci.distance = distances[defined_node_ids[i]];
      // a defined node is alive, so a quaternion is associated to it
      math::Quaternion quat_i = (*m_cross_field)[defined_node_ids[i]];
      ci.angle = quat_ref.angle(quat_i);
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
double  AdvancingFrontFieldCommon::
computeStability(std::set<Contribution>& AContributions)
{
  double res = 0.0;
  double lastDist = 0.0;
  double delta = 0.0;

  std::set<Contribution>::const_iterator it = AContributions.begin();
  const std::set<Contribution>::const_iterator it_end = AContributions.end();
  while (it != it_end)
    {
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
/*---------------------------------------------------------------------------*/
bool AdvancingFrontFieldCommon::belongToTheSameFace(
						    gmds::Node& ATo,
						    gmds::Node& AFrom1,
						    gmds::Node& AFrom2,
						    gmds::Face& AFace)
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
/*---------------------------------------------------------------------------*/
void AdvancingFrontFieldCommon::extrapolateQuaternion(
						      gmds::Node& ATo,
						      gmds::Node& AFrom1,
						      gmds::Node& AFrom2,
						      gmds::Face& AFace)
{
  //===================================================================
  //We get quaternions
  //===================================================================
  math::Quaternion q1 = (*m_cross_field)[AFrom1.getID()];
  math::Quaternion q2 = (*m_cross_field)[AFrom2.getID()];
  math::Vector normal = AFace.normal();

  //q1 = q1.alignWith(normal);
  //q2 = q2.alignWith(normal);

  std::vector<math::Quaternion> qs;
  qs.push_back(q1);
  qs.push_back(q2);

  //===================================================================
  //We get coeff (inverse of distance)
  //===================================================================
  std::vector<TCoord> cs;
  math::Point p0 = ATo.getPoint();
  math::Point p1 = AFrom1.getPoint();
  math::Point p2 = AFrom2.getPoint();

  math::Vector v01(p0, p1);
  math::Vector v02(p0, p2);
  cs.push_back(1.0/v01.norm());
  cs.push_back(1.0/v02.norm());

  //===================================================================
  //quaternion extrapolated in the face AFace plane
  //===================================================================
  math::Quaternion q = math::Quaternion::mean(qs, cs);

  //===================================================================
  //We compute the normal to the surface now
  //===================================================================
  std::vector<Face> adj_faces = ATo.get<Face>();
  std::vector<Face> adj_bnd_faces;
  for (unsigned int i = 0; i < adj_faces.size(); i++)
    {
      Face fi = adj_faces[i];
      if (m_mesh->isMarked(fi, m_markFaceOnSurf))
	adj_bnd_faces.push_back(fi);
    }
  math::Vector out_vec;
  for (unsigned int i = 0; i < adj_bnd_faces.size(); i++)
    {
      Face fi = adj_bnd_faces[i];
      out_vec = out_vec + getOutputNormalOfABoundaryFace(fi);
    }
  out_vec = out_vec / adj_bnd_faces.size();

  //===================================================================
  //We assign the quaternion align with out_vec
  //===================================================================
  math::Quaternion res = q.alignWith(out_vec);
  (*m_cross_field)[ATo.getID()] = res;

}
/*----------------------------------------------------------------------------*/
void AdvancingFrontFieldCommon::extrapolateQuaternion(
						      gmds::Node& ATo,
						      std::vector<gmds::Node>& AFrom)
{
  if (m_mesh->isMarked(ATo, m_markNodeOnBnd))
    extrapolateQuaternionOnSurf(ATo, AFrom);
  else
    extrapolateQuaternionInVol(ATo, AFrom);

}
/*----------------------------------------------------------------------------*/
void AdvancingFrontFieldCommon::extrapolateQuaternionOnSurf(
							    gmds::Node& ATo,
							    std::vector<gmds::Node>& AFrom)
{
  math::Point p0 = ATo.getPoint();

  //===================================================================
  //We compute the normal to the surface at ATo
  //===================================================================
  std::vector<Face> adj_faces = ATo.get<Face>();
  std::vector<Face> adj_bnd_faces;
  for (unsigned int i = 0; i < adj_faces.size(); i++)
    {
      Face fi = adj_faces[i];
      if (m_mesh->isMarked(fi, m_markFaceOnSurf))
	adj_bnd_faces.push_back(fi);
    }
  math::Vector out_vec;
  for (unsigned int i = 0; i < adj_bnd_faces.size(); i++)
    {
      Face fi = adj_bnd_faces[i];
      out_vec = out_vec + getOutputNormalOfABoundaryFace(fi);
    }
  out_vec = out_vec / adj_bnd_faces.size();

  //===================================================================
  //We get quaternions and  coeff (inverse of distance)
  //===================================================================
  std::vector<math::Quaternion> qs;
  std::vector<TCoord> cs;
  for (unsigned int i = 0; i < AFrom.size(); i++)
    {
      Node current = AFrom[i];
      //quaternion 
      math::Quaternion q = (*m_cross_field)[current.getID()];
      math::Quaternion q_aligned = q.alignWith(out_vec);
      qs.push_back(q);
      //coeef
      math::Point p1 = current.getPoint();
      math::Vector v01(p0, p1);
      cs.push_back(1.0/v01.norm());
    }

  //===================================================================
  //quaternion extrapolated in the face AFace plane
  //===================================================================
  math::Quaternion q = math::Quaternion::mean(qs, cs);

  //===================================================================
  //We assign the quaternion align with out_vec
  //===================================================================

  math::Quaternion res = q.alignWith(out_vec);
  (*m_cross_field)[ATo.getID()] = res;

}
/*----------------------------------------------------------------------------*/
void AdvancingFrontFieldCommon::extrapolateQuaternionInVol(
							   gmds::Node& ATo,
							   std::vector<gmds::Node>& AFrom)
{
  math::Point p0 = ATo.getPoint();
  //===================================================================
  //We get quaternions and  coeff (inverse of distance)
  //===================================================================
  std::vector<math::Quaternion> qs;
  std::vector<TCoord> cs;
  for (unsigned int i = 0; i < AFrom.size(); i++)
    {
      Node current = AFrom[i];
      //quaternion 
      math::Quaternion q = (*m_cross_field)[current.getID()];
      qs.push_back(q);
      //coeef
      math::Point p1 = current.getPoint();
      math::Vector v01(p0, p1);
      cs.push_back(1.0/v01.norm());
    }

  //===================================================================
  //quaternion extrapolated 
  //===================================================================
  math::Quaternion q = math::Quaternion::mean(qs, cs);
  (*m_cross_field)[ATo.getID()] = q;

}
/*----------------------------------------------------------------------------*/
std::vector<Node> AdvancingFrontFieldCommon::
getAdjacentNodes(Node& ANode, const int AMark)
{
  std::set<Node> adj_nodes;
  std::vector<Face>  adj_faces = ANode.get<Face>();

  for (unsigned int j = 0; j < adj_faces.size(); j++)
    {
      Face f_j = adj_faces[j];

      //only boundary faces have to be taken into account
      if (!m_mesh->isMarked(f_j, m_markFaceOnSurf))
	continue;
      std::vector<Node> nodes_fj = f_j.get<Node>();
      for (unsigned int j = 0; j < nodes_fj.size(); j++)
	{
	  Node n_j = nodes_fj[j];
	  if (m_mesh->isMarked(n_j, AMark) && n_j.getID() != ANode.getID())
	    {
	      adj_nodes.insert(n_j);
	    }
	}
    }
  std::vector<Node> adj;
  adj.insert(adj.end(), adj_nodes.begin(), adj_nodes.end());

  return adj;
}

/*----------------------------------------------------------------------------*/
std::vector<Node> AdvancingFrontFieldCommon::
getAdjacentNodesByEdge(Node& ANode, const int AMark)
{
  std::set<Node> adj_nodes;
  std::vector<Edge>  adj_edges = ANode.get<Edge>();

  for (unsigned int j = 0; j < adj_edges.size(); j++)
    {
      Edge e_j = adj_edges[j];
      std::vector<Node> nodes_ej = e_j.get<Node>();
      for (unsigned int j = 0; j < nodes_ej.size(); j++)
	{
	  Node n_j = nodes_ej[j];
	  if (m_mesh->isMarked(n_j, AMark) && n_j.getID() != ANode.getID())
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
std::vector<Edge>  AdvancingFrontFieldCommon::
getEdgesOnCurve(Node& ANode) const
{
    std::vector<Edge> edges_on_curve;
    std::vector<Edge> adj_edges = ANode.get<Edge>();
    for (unsigned int i = 0; i < adj_edges.size(); i++)
    {
        Edge ei = adj_edges[i];
        if (m_mesh->isMarked(ei, m_markEdgeOnCurv))
            edges_on_curve.push_back(ei);
    }
    return edges_on_curve;
}
/*---------------------------------------------------------------------------*/
std::vector<Face>  AdvancingFrontFieldCommon::
getFacesOnSurface(Edge& AEdge) const
{
  std::vector<Face> faces_on_surf;
  std::vector<Face> adj_faces = AEdge.get<Face>();
  for (unsigned int i = 0; i < adj_faces.size(); i++)
    {
      Face fi = adj_faces[i];
      if (m_mesh->isMarked(fi, m_markFaceOnSurf))
	faces_on_surf.push_back(fi);
    }
  return faces_on_surf;
}
/*---------------------------------------------------------------------------*/
Node  AdvancingFrontFieldCommon::getNeighboorOn(Node& AN, Edge& AE) const
{
  std::vector<Node> nodes = AE.get<Node>();
  if (nodes[0].getID() == AN.getID())
    return nodes[1];

  return nodes[0];
}
/*---------------------------------------------------------------------------*/
void  AdvancingFrontFieldCommon::initQuaternionOnCurves()
{
  //for each node on a geometrical curve, we compute
  //its associated quaternion
  IGMesh::node_iterator it_node = m_mesh->nodes_begin();
  for (; !it_node.isDone(); it_node.next())
    {
      Node current_node = it_node.value();
	
      math::Chart current_chart;
      if (m_mesh->isMarked(current_node, m_markNodeOnCurv) &&
	  !m_mesh->isMarked(current_node, m_markNodeOnPnt))
	{
	  //current_node is on a geometric curve
	  //We get adjacent edges that are classified onto
	  //a geometric curve. We must get one or two edges
	  //at most.
	  std::vector<Edge> ridges = getEdgesOnCurve(current_node);
	  math::Vector v1, newN1, newN2;
	  if (ridges.size() == 1)
	    {
	      Edge current_edge = ridges[0];
	      Node node1 = getNeighboorOn(current_node, current_edge);
	      Node node2 = current_node;

	      //we build the direction vector of the current edge
	      math::Point p1 = node1.getPoint();
	      math::Point p2 = node2.getPoint();
	      v1 = math::Vector(p1, p2);
	      v1.normalize();
	      //We retrieve adjacent boundary faces and we get their normal vectors
	      std::vector<Face> boundary_faces = getFacesOnSurface(current_edge);

	      math::Vector n1 = boundary_faces[0].normal();
	      math::Vector n2 = boundary_faces[1].normal();
	      n1.normalize();
	      n2.normalize();

	      //They are projected onto the tangent vector
	      n1 = n1 - (v1.dot(n1) * v1);
	      n2 = n2 - (v1.dot(n2) * v1);
	      n1.normalize();
	      n2.normalize();

	      newN1 = n1;
	      newN2 = n2;
	    }
	  else if (ridges.size() == 2)
	    {
	      //With 2 adajcent edges on the curve, we compute average values

	      Edge edge1 = ridges[0];
	      Edge edge2 = ridges[1];

	      Node node1 = getNeighboorOn(current_node, edge1);
	      Node node2 = getNeighboorOn(current_node, edge2);
	      //we build the direction vector linking the two adjacent nodes
	      math::Point p1 = node1.getPoint();
	      math::Point p2 = node2.getPoint();
	      v1 = math::Vector(p1, p2);
	      v1.normalize();

	      //Normal are computed in the same way
	      std::vector<Face> faces1 = getFacesOnSurface(edge1);
	      std::vector<Face> faces2 = getFacesOnSurface(edge2);
	      Face SelectedFace[4];
	      SelectedFace[0] = faces1[0];
	      SelectedFace[1] = faces1[1];
	      SelectedFace[2] = faces2[0];
	      SelectedFace[3] = faces2[1];

	      math::Vector n1 = SelectedFace[0].normal();
	      math::Vector n2 = SelectedFace[1].normal();
	      math::Vector n3 = SelectedFace[2].normal();
	      math::Vector n4 = SelectedFace[3].normal();

	      n1.normalize();
	      n2.normalize();
	      n3.normalize();
	      n4.normalize();


	      //Calcul des produits scalaires
	      const TCoord prodScal1 = abs(n1.dot(n3));
	      const TCoord prodScal2 = abs(n1.dot(n4));

	      //Correction
	      if (prodScal1 < prodScal2)
		{
		  Face fTmp = SelectedFace[2];
		  SelectedFace[2] = SelectedFace[3];
		  SelectedFace[3] = fTmp;
		}

	      //Recuperation des vecteurs normaux.

	      math::Vector normal1 = SelectedFace[0].normal();
	      math::Vector normal2 = SelectedFace[1].normal();
	      {
		math::Vector normal3 = SelectedFace[2].normal();
		math::Vector normal4 = SelectedFace[3].normal();
		//Mise dans le mememe sens des normales
		TCoord scal;
		scal = normal1.dot(normal3);
		if (scal < 0.0)
		  normal3 = normal3.opp();
		scal = normal2.dot(normal4);
		if (scal < 0.0)
		  normal4 = normal4.opp();

		//Moyenne des normales
		normal1 = normal1 + normal3;
		normal2 = normal2 + normal4;

		normal1.normalize();
		normal2.normalize();
	      }
	      //Projection des normales sur la tangente
	      normal1 = normal1 - (v1.dot(normal1)* v1);
	      normal2 = normal2 - (v1.dot(normal2)* v1);
	      //Normalisation des vecteurs projetes
	      normal1.normalize();
	      normal2.normalize();

	      newN1 = normal1;
	      newN2 = normal2;
	    }
	  else{
	    std::cout << "Nb ridges for an edge adjacent to node "
		      << current_node.getID() << ": " << ridges.size() << std::endl;
	    throw GMDSException("A ridge node has an illegal number of edges.");
	  }
	  //Creation of N1 prime
	  math::Vector N1Prime = v1.cross(newN1);
	  if (N1Prime.dot(newN2) < 0.0)
	    N1Prime = N1Prime.opp();

	  //Computation of v2 from N1' and N2
	  math::Vector v2 = N1Prime + newN2;
	  v2.normalize();

	  //Computation of v3 from v1 and v2.
	  const math::Vector v3 = v1.cross(v2);

	  current_chart = math::Chart(v1, v2, v3);
	  math::Quaternion q(current_chart);

	  (*m_cross_field)[current_node.getID()] = q;
	}//if (m_mesh->isMarked(current_node, m_markNodeOnCurv))

    }//for (; !it_node.isDone(); it_node.next())




}
/*----------------------------------------------------------------------------*/
math::Vector AdvancingFrontFieldCommon::
getOutputNormal(Face& AFace, Region& ARegion)
{
  std::vector<Node> region_nodes = ARegion.get<Node>();
  std::vector<Node> face_nodes = AFace.get<Node>();

  if (region_nodes.size() != 4)
    throw GMDSException("SingularityGraphBuilder::getOutputNormal can only be used on tetrahedral regions");
  if (face_nodes.size() != 3)
    throw GMDSException("SingularityGraphBuilder::getOutputNormal can only be used on triangular faces");

  //we go through all the nodes of ARegion to find the one that do not belong
  //to AFAce
  for (unsigned int i = 0; i<region_nodes.size(); i++)
    {
      Node  n = region_nodes[i];
      if (n != face_nodes[0] && n != face_nodes[1] && n != face_nodes[2])
	{
	  //n is the node opposite to the face AFace
	  Node n0 = face_nodes[0];
	  Node n1 = face_nodes[1];
	  Node n2 = face_nodes[2];
	  math::Vector normal_to_face = AFace.normal();
	  math::Vector in_vector(n0.getPoint(), n.getPoint());
	  if (normal_to_face.dot(in_vector)>0.0)
	    {
	      return math::Vector(-normal_to_face.get(0),
				  -normal_to_face.get(1),
				  -normal_to_face.get(2));
	    }
	  else
	    {
	      return normal_to_face;
	    }

	} //if (n != face_nodes[0] && n != face_nodes[1] && n != face_nodes[2])
	
    }//for (unsigned int i = 0; i<region_nodes.size(); i++)
  return math::Vector(0, 0, 0);
}
/*----------------------------------------------------------------------------*/
math::Vector AdvancingFrontFieldCommon::
getOutputNormalOfABoundaryFace(Face& AFace)
{
  std::vector<Region> adj_regions = AFace.get<Region>();
  if (adj_regions.size() != 1)
    throw GMDSException("A boundary face must be adjacent to only 1 region!!!");

  return getOutputNormal(AFace, adj_regions[0]);
}

/*---------------------------------------------------------------------------*/
void AdvancingFrontFieldCommon::writeForDebug(const std::string AFileName)
{
  static int nb_file = 0;
  VTKWriter<IGMesh> writer(*m_mesh);
  if (AFileName != "")
    {
      std::cout<<"Write in "<<AFileName<<std::endl;
      writer.write(AFileName, DIM3 | R | F | N);

    }
  else
    {
      std::stringstream file_name;
      file_name <<m_output_directory_name<<"/FFG_Debug_" << nb_file;
      std::cout<<"Write in "<<file_name<<std::endl;
      writer.write(file_name.str(), DIM3 | R | F | N);
    }

  IGMesh::node_iterator it = m_mesh->nodes_begin();
  double x_min = 100000;
  double y_min = 100000;
  double z_min = 100000;
  double x_max = -100000;
  double y_max = -100000;
  double z_max = -100000;
  for (; !it.isDone(); it.next())
    {
      Node n = it.value();
      math::Point p = n.getPoint();
      if (p.X() < x_min)
	x_min = p.X();
      if (p.X() > x_max)
	x_max = p.X();

      if (p.Y() < y_min)
	y_min = p.Y();
      if (p.Y() > y_max)
	y_max = p.Y();

      if (p.Z() < z_min)
	z_min = p.Z();
      if (p.Z() > z_max)
	z_max = p.Z();
    }
  double dist_x = x_max - x_min;
  double dist_y = y_max - y_min;
  double dist_z = z_max - z_min;

  double cube_size = 0;
  if (dist_x <= dist_y && dist_x <= dist_z){
    cube_size = dist_x;
  }
  else if (dist_y <= dist_x && dist_y <= dist_z){
    cube_size = dist_y;
  }
  else
    cube_size = dist_z;

  cube_size /= 20;

  MeshModel model_cube(DIM3 | R | N | R2N);
  IGMesh mesh_cube(model_cube);
  for (it = m_mesh->nodes_begin(); !it.isDone(); it.next())
    {
      Node n = it.value();
      math::Point center = n.getPoint();
      if (m_mesh->isMarked(n, m_mark_alive)){
	math::Quaternion q = (*m_cross_field)[n.getID()];
	math::Chart t(q);
	math::Vector vx = t.X().getNormalize();
	math::Vector vy = t.Y().getNormalize();
	math::Vector vz = t.Z().getNormalize();
	math::Point p1 = center + (vx + vy - vz)*cube_size;
	Node n1 = mesh_cube.newNode(p1);
	math::Point p2 = center + (vx - vy - vz)*cube_size;
	Node n2 = mesh_cube.newNode(p2);
	math::Point p3 = center + (vx + vy + vz).opp()*cube_size;
	Node n3 = mesh_cube.newNode(p3);
	math::Point p4 = center + (vy - vx - vz)*cube_size;
	Node n4 = mesh_cube.newNode(p4);

	math::Point p5 = center + (vx + vy + vz)*cube_size;
	Node n5 = mesh_cube.newNode(p5);
	math::Point p6 = center + (vx - vy + vz)*cube_size;
	Node n6 = mesh_cube.newNode(p6);
	math::Point p7 = center + (vx + vy - vz).opp()*cube_size;
	Node n7 = mesh_cube.newNode(p7);
	math::Point p8 = center + (vy - vx + vz)*cube_size;
	Node n8 = mesh_cube.newNode(p8);

	mesh_cube.newHex(n1, n2, n3, n4, n5, n6, n7, n8);
      }
    }
  VTKWriter<IGMesh> writer_cube(mesh_cube);

  std::stringstream file_name_cube;
  file_name_cube <<m_output_directory_name<<"/FFG_Debug_Cube_" << nb_file;
  writer_cube.write(file_name_cube.str(), DIM3 | R | F | N);

  nb_file++;
}
/*---------------------------------------------------------------------------*/
