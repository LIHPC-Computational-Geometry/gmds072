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
/*---------------------------------------------------------------------------*/
#include <iostream>
/*---------------------------------------------------------------------------*/
#include <GMDS/Algo/BoundaryOperator.h>
#include <GMDS/Algo/DistanceFieldBuilder2D.h>
#include <GMDS/Algo/DistanceFieldBuilder3D.h>
/*---------------------------------------------------------------------------*/
#include <GMDS/IO/VTKWriter.h>
#include <sstream>
#include <map>
/*---------------------------------------------------------------------------*/
#include "AdvancingFrontFieldGen3D.h"
#include "SmoothingHLBFGS.h"
#include "FrameFieldSmoother.h"
#include "FrameFieldLaplacianSmoothing.h"
/*---------------------------------------------------------------------------*/
using namespace gmds;
/*---------------------------------------------------------------------------*/
AdvancingFrontFieldGen3D::AdvancingFrontFieldGen3D(gmds::IGMesh* AMesh)
  : AdvancingFrontFieldCommon(AMesh),
    m_distance_method(EDistanceField_DisconnectedSurfAndVol)
{
  m_distance_field_2D = 0;
  m_distance_field_3D = 0;
  m_max_distance_2D = 0;
  m_max_distance_3D = 0;
  m_max_distance_2D_div2 = 0;
  m_max_distance_3D_div2 = 0;
  m_distanceBuilder2D = DistanceFieldBuilder2D(m_mesh);
  m_distanceBuilder3D = DistanceFieldBuilder3D(m_mesh);
		  
}
/*---------------------------------------------------------------------------*/
void AdvancingFrontFieldGen3D::execute()
{
  std::cout<<"============ Start frame field generation ==========="<<std::endl;
   std::cout << "Mark' initialization " << std::endl; 
   initMarks(); //SHARED BEHAVIOUR
   std::cout << "  -->  DONE" << std::endl;
  //==================================================================
  //STEP 1 - We mark all the nodes, edges, faces classified on
  // boundary entities (points, curves, surfaces)
  //==================================================================
  std::cout << "=============================================================" << std::endl;
  std::cout << "Mark boundary " << std::endl;
  markBoundaryCells();  //SHARED BEHAVIOUR
  std::cout << "  -->  DONE" << std::endl;
  //==================================================================
  //STEP 2 - We store all the nodes classified on curves and surfaces
  // in STL vectors
  //=================================================================
  std::cout << "=============================================================" << std::endl;
  std::cout << " Storage of nodes classified on curves and surfaces" << std::flush;
  IGMesh::node_iterator it_nodes = m_mesh->nodes_begin();
  for (; !it_nodes.isDone(); it_nodes.next())
    {
      Node n = it_nodes.value();
      if (m_mesh->isMarked(n, m_markNodeOnCurv) ||
	  m_mesh->isMarked(n, m_markNodeOnPnt))
	{
	  m_curve_nodes.push_back(n);
	  m_boundary_nodes.push_back(n);
	}
      if (m_mesh->isMarked(n, m_markNodeOnSurf))
	{
	  m_surf_nodes.push_back(n);
	  m_boundary_nodes.push_back(n);
	}
    }
  std::cout << "    DONE" << std::endl;

  //==================================================================
  //STEP 2 - For nodes on curves, we compute quaternion from the 
  // geometric information we have
  //==================================================================
  std::cout << "=============================================================" << std::endl;
  std::cout << "Initialization of quaternions on curves" << std::flush;
  initQuaternionOnCurves(); //SHARED BEHAVIOUR

  writeForDebug();
  std::cout << "    DONE" << std::endl;

  //==================================================================
  //STEP 3 - A quaternion is associated to each node classified on 
  // a geometric point
  //==================================================================
  std::cout << "=============================================================" << std::endl;
  std::cout << "Initialization of quaternions on points" << std::flush;
  initQuaternionOnPoints();  //SHARED BEHAVIOUR
  writeForDebug();
  std::cout << "    DONE" << std::endl;


  //==================================================================
  //STEP 4- Computes the distance field(s)
  //==================================================================
  std::cout << "=============================================================" << std::endl;
  std::cout << "Computation of distance fields" << std::endl;
  computeDistanceFields(); //CHILDREN'S DELEGATED BEHAVIOUR
  writeForDebug();
  std::cout << "    DONE" << std::endl;



  //==================================================================
  //STEP 5 - buid and compute the distance to the approximated medial
  // axis - NOT DONE NOW
  //==================================================================
  //std::vector<gmds::Node>  nodes_of_approximate_MA =
  //	computeApproximateMedialAxis(surf_nodes);
  //m_medial_axis_distance_field =
  //			dist_computer.computeDistance(nodes_of_approximate_MA,
  //											surface_nodes_ids);



  //==================================================================
  //STEP 6 - Marks nodes where a quaternion is defined on
  //==================================================================
  std::cout << "=============================================================" << std::endl;
  std::cout << "Alive nodes initialization" << std::flush;
  initAliveNodes();
  writeForDebug();
  std::cout << "    DONE" << std::endl;

  //==================================================================
  //STEP 6 - We compute the normal along the boundary
  //==================================================================
  std::cout << "=============================================================" << std::endl;
  std::cout << "Computes normals" << std::flush;
  initNormalOnBoundarySurfaceNodes();
  writeForDebug();
  std::cout << "    DONE" << std::endl;

  //==================================================================
  //STEP 8 - We get in an ordered way all the nodes, which are closer 
  // to the boundary than to the MA. In this fist version we only uses
  // a proportion of the nodes added during the distance computation
  //==================================================================
  //std::cout << "=============================================================" << std::endl;
  //std::cout << "Boundary snapping" << std::flush;
  //initQuaternionsBySnapping();
  //writeForDebug();
  //std::cout << "    DONE" << std::endl;


  //==================================================================
  //STEP 9 - Now, we assign a quaternion to every  nodes !!!
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
  std::cout << "=============================================================" << std::endl;
  buildFrameField();

  math::Vector aleat1(1, 1, 1);
  it_nodes = m_mesh->nodes_begin();
  for (; !it_nodes.isDone(); it_nodes.next())
    {
      Node n = it_nodes.value();
      if (!m_mesh->isMarked(n, m_markNodeOnCurv) &&
	  !m_mesh->isMarked(n, m_markNodeOnSurf) &&
	  !m_mesh->isMarked(n, m_markNodeOnPnt))
	{
	  /*	math::Quaternion q = (*m_cross_field)[n.getID()];
		q = q.alignWith(aleat1);
		(*m_cross_field)[n.getID()] = q;*/
	}
    }
  writeForDebug();
  std::cout << "    DONE" << std::endl;
  //==================================================================
  // IMPROVEMENT STEP - Smoothing ???
  //==================================================================
  std::cout << "=============================================================" << std::endl;
  std::cout << " smoothing" << std::flush;
  smoothAll();
  std::cout << "    done" << std::endl;

    
  //==================================================================
  // DEBUGGING STEP - Tet coloring
  //==================================================================
  std::cout << "=============================================================" << std::endl;
  std::cout << " Coloring simplices" << std::endl;
  colorSimplices();
  writeForDebug();
  std::cout << "    DONE" << std::endl;

  //==================================================================
  // CLEANING - Boolean marks are cleaned
  //==================================================================
  cleanMarks();
}

/*---------------------------------------------------------------------------*/
void AdvancingFrontFieldGen3D::
setDistanceFieldMethod(const EDistanceFieldMethod& AMethod)
{
	m_distance_method = AMethod;
}
/*---------------------------------------------------------------------------*/
double AdvancingFrontFieldGen3D::
getDistance(const Node& ANode)
{
	double val = 0;
	if (m_distance_method == EDistanceField_VolFromEdge)
	{
		val = (*m_distance_field_3D)[ANode.getID()];
	}
	else
	{
		//We can use the surface and the volume distance
		if (m_mesh->isMarked(ANode, m_markNodeOnSurf) ||
			m_mesh->isMarked(ANode, m_markNodeOnCurv) ||
			m_mesh->isMarked(ANode, m_markNodeOnPnt))
		{
			val = (*m_distance_field_2D)[ANode.getID()];
		}
		else //ANode is in the volume
		{
			val = (*m_distance_field_3D)[ANode.getID()];
		}
	}
	return val;
}
/*---------------------------------------------------------------------------*/
void AdvancingFrontFieldGen3D::computeDistanceFields()
{
	//==================================================================
	//we build the inner_node collection
	//==================================================================
	m_inner_nodes.clear();
	IGMesh::node_iterator it_nodes = m_mesh->nodes_begin();

	for (; !it_nodes.isDone(); it_nodes.next())
	{
		Node n = it_nodes.value();
		if (!m_mesh->isMarked(n, m_markNodeOnCurv) &&
			!m_mesh->isMarked(n, m_markNodeOnPnt) &&
			!m_mesh->isMarked(n, m_markNodeOnSurf))
		{
			m_inner_nodes.push_back(n);
		}
	}
	//==================================================================
	// We compute the distance fields 
	//==================================================================
	if (m_distance_method == EDistanceField_DisconnectedSurfAndVol)
	{
		// Computes the 2D distance field
		DistanceFieldBuilder2D distanceBuilder2D(m_mesh);
		if (!distanceBuilder2D.isValid())
		{
			std::cout << "Invalid model for distance computation" << std::endl;
			throw GMDSException("Invalid model for distance computation");
		}
		m_distance_field_2D = distanceBuilder2D.computeDistance(
			m_curve_nodes, m_surf_nodes, m_markFaceOnSurf);

		// Computes the 3D distance field
		DistanceFieldBuilder3D distanceBuilder3D(m_mesh);
		if (!distanceBuilder3D.isValid())
		{
			std::cout << "Invalid model for distance computation" << std::endl;
			throw GMDSException("Invalid model for distance computation");
		}
		m_distance_field_3D = distanceBuilder3D.computeDistance(m_boundary_nodes, m_inner_nodes);


		//==================================================================
		//We get the maximum distance 
		//==================================================================
		it_nodes = m_mesh->nodes_begin();
		for (; !it_nodes.isDone(); it_nodes.next())
		{
			Node n = it_nodes.value();
			double dist2D = (*m_distance_field_2D)[n.getID()];
			double dist3D = (*m_distance_field_3D)[n.getID()];
			if (dist2D > m_max_distance_2D)
				m_max_distance_2D = dist2D;
			if (dist3D > m_max_distance_3D)
				m_max_distance_3D = dist3D;
		}


	}
	else
	{
		// Computes the 3D distance field from curves
		DistanceFieldBuilder3D distanceBuilder3D(m_mesh);
		if (!distanceBuilder3D.isValid())
		{
			std::cout << "Invalid model for distance computation" << std::endl;
			throw GMDSException("Invalid model for distance computation");
		}

		std::vector<Node> to_propag;
		to_propag.insert(to_propag.end(), m_surf_nodes.begin(), m_surf_nodes.end());
		to_propag.insert(to_propag.end(), m_inner_nodes.begin(), m_inner_nodes.end());
		m_distance_field_3D = distanceBuilder3D.computeDistance(m_curve_nodes, to_propag);
	}

	std::cout << "Distance max 2D: " << m_max_distance_2D << std::endl;
	std::cout << "Distance max 3D: " << m_max_distance_3D << std::endl;
	m_max_distance_2D_div2 = m_max_distance_2D;// / 2.0;
	m_max_distance_3D_div2 = m_max_distance_3D;// / 2.0;

}
/*---------------------------------------------------------------------------*/
void AdvancingFrontFieldGen3D::initQuaternionsBySnapping()
{
	std::cout << "NOT IMPLEMENTED: AdvancingFrontFieldGen3D::initQuaternionsBySnapping" << std::endl;
	throw GMDSException("NOT IMPLEMENTED: AdvancingFrontFieldGen3D::initQuaternionsBySnapping");
	//std::vector<TCellID> ordered_node_ids =
	//	m_distanceBuilder2D.getInsertionOrder();
	//int nb_direct_snap = 0;// ordered_node_ids.size() / 3;

	//for (int i = 0; i < nb_direct_snap; i++)
	//{
	//	//As these nodes were inserted during the distance field
	//	//insertion process, at least two adjacent nodes have a
	//	//defined quaternion
	//	gmds::Node n = m_mesh->get<Node>(ordered_node_ids[i]);
	//	bool done = computeLocalFrame(n, m_mark_alive);

	//}

}

/*----------------------------------------------------------------------------*/
std::vector<Node> AdvancingFrontFieldGen3D::
getAdjacentNodesByRegion(Node& ANode)
{
	std::set<Node> adj_nodes;
	std::vector<Region>  adj_regions = ANode.get<Region>();

	for (unsigned int j = 0; j < adj_regions.size(); j++)
	{
		Region r_j = adj_regions[j];
		std::vector<Node> nodes_rj = r_j.get<Node>();
		for (unsigned int j = 0; j < nodes_rj.size(); j++)
		{
			Node n_j = nodes_rj[j];
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
/*----------------------------------------------------------------------------*/
void AdvancingFrontFieldGen3D::getNodesInTheSameTetrahedron(
	const std::vector<Node>& AIN,
	std::vector<std::vector<Node> >& AOUT)
{

	for (unsigned int i = 0; i < AIN.size(); i++)
	{
		Node ni = AIN[i];
		std::vector<TCellID> regions_i = ni.getIDs<Region>();
		for (unsigned int j = i + 1; j < AIN.size(); j++)
		{
			Node nj = AIN[j];
			std::vector<TCellID> regions_j = nj.getIDs<Region>();

			//we look for all the regions that would be both in fsi and fsj
			std::vector<TCellID> regions_ij;
			for (unsigned int ki = 0; ki < regions_i.size(); ki++)
			{
				TCellID id1 = regions_i[ki];
				for (unsigned int kj = 0; kj < regions_j.size(); kj++)
				{
					TCellID id2 = regions_j[kj];
					if (id1 == id2)
					{
						regions_ij.push_back(id1);
					}
				}
			}
			if (!regions_ij.empty())
			{
				//for each region, we look if a third node of AIN would be adjacent to
				for (unsigned int ir = 0; ir < regions_ij.size(); ir++)
				{
					Region r = m_mesh->get<Region>(regions_ij[ir]);
					std::vector<TCellID> r_nodes = r.getIDs<Node>();
					for (unsigned int k = j + 1; k < AIN.size(); k++)
					{
						TCellID nk_id = AIN[k].getID();
						bool found = false;
						for (unsigned l = 0; l < r_nodes.size(); l++)
						{
							if (r_nodes[l] == nk_id)
								found = true;
						}
						if (found){
							std::vector<Node> triplet;
							triplet.resize(3);
							triplet[0] = ni;
							triplet[1] = nj;
							triplet[2] = AIN[k];
							AOUT.push_back(triplet);
						}
					}//for (unsigned int k = j + 1; k < AIN.size(); k++)

				}//for (unsigned int ir = 0; ir < regions_ij.size(); ir++)

			}//if (!regions_ij.empty())

		}//for (unsigned int j = i + 1; j < AIN.size(); j++)

	}//	for (unsigned int i = 0; i < AIN.size(); i++)

}



/*----------------------------------------------------------------------------*/
double AdvancingFrontFieldGen3D::
distancePonderation(const double AVal, const double ADist)
{
	if (ADist == 0)
		return 0;
	else
		return AVal;// *ADist;// *log(ADist);
}
/*----------------------------------------------------------------------------*/
void AdvancingFrontFieldGen3D::
initNarrowBandSurf(std::list<StabilityInfo>& ANarrowBand)
{
	IGMesh::node_iterator it_node = m_mesh->nodes_begin();
	for (; !it_node.isDone(); it_node.next())
	{
		Node current_node = it_node.value();

		if (m_mesh->isMarked(current_node, m_mark_alive))
		{
			//current_node is alive, so classified on a curve
			//We get adjacent nodes that are on the surface
			std::vector<Node> adj_nodes =
				getAdjacentNodes(current_node, m_markNodeOnSurf);

			for (unsigned int j = 0; j < adj_nodes.size(); j++)
			{
				Node n_j = adj_nodes[j];
				// we look for nodes that are not yet in the narrow band or alive
				if (m_mesh->isMarked(n_j, m_mark_narrow) ||
					m_mesh->isMarked(n_j, m_mark_alive))
					continue;

				std::vector<Node> adj_nodes_j =
					getAdjacentNodes(n_j, m_markNodeOnSurf);

				std::vector<Node> ext_nodes;
				int nb_alive_adj_nodes = 0;
				for (unsigned int k = 0; k < adj_nodes_j.size(); k++)
				{
					Node n_k = adj_nodes_j[k];
					if (m_mesh->isMarked(n_k, m_mark_alive) &&
						!m_mesh->isMarked(n_k, m_markNodeOnPnt))
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
						m_mesh->mark(n_j, m_mark_narrow);
						// We can now extrapolate a quaternion in this node... 
						extrapolateQuaternion(n_j, ext_nodes[0], ext_nodes[1], commonFace);
						//we keep in mind the order in which quaternions are created
						m_nodes_by_creation_order.push_back(n_j.getID());
						// then compute a stability value in n_j
						double local_stab = computeStability(n_j, m_mark_alive);
						//weighted by the distance to the curves
						double stab = distancePonderation(local_stab,
							(*m_distance_field_2D)[n_j.getID()]);
						m_stability_value[n_j.getID()] = stab;

						ANarrowBand.push_back(StabilityInfo(n_j, stab));
					}

				}//if (nb_alive_adj_nodes == 2)
				else if (nb_alive_adj_nodes > 2)
				{
					m_mesh->mark(n_j, m_mark_narrow);
					// We can now extrapolate a quaternion in this node... 
					extrapolateQuaternion(n_j, ext_nodes);
					//we keep in mind the order in which quaternions are created
					m_nodes_by_creation_order.push_back(n_j.getID());
					// then compute a stability value in n_j
					double local_stab = computeStability(n_j, m_mark_alive);
					//weighted by the distance to the curves
					double stab = distancePonderation(local_stab,
						(*m_distance_field_2D)[n_j.getID()]);
					m_stability_value[n_j.getID()] = stab;
					ANarrowBand.push_back(StabilityInfo(n_j.getID(), stab));
				}

			}//for (unsigned int j = 0; j < adj_nodes.size(); j++)

		}//if (m_mesh->isMarked(current_node, AMarkAlive))

	}//for (; !it_node.isDone(); it_node.next())

	std::cout << "Surface Narrow Band Size: " << ANarrowBand.size() << std::endl;
}

/*---------------------------------------------------------------------------*/
AdvancingFrontFieldGen3D::StabilityInfo AdvancingFrontFieldGen3D::
extractElectFrom(std::list<StabilityInfo>& AL)
{
	//==================================================================
	//we look for the most stable one
	//==================================================================
	std::list<StabilityInfo>::iterator it = AL.begin();
	StabilityInfo elect = *it;
	for (; it != AL.end(); it++)
	{
		if (it->stability < elect.stability)
			elect = *it;
	}
	//==================================================================
	// we remove it from AL
	//==================================================================
	AL.remove(elect);

	//==================================================================
	//we return it
	//==================================================================
	return elect;
}
/*---------------------------------------------------------------------------*/
void AdvancingFrontFieldGen3D::
propagateQuaternionsSurf(std::list<StabilityInfo>& ANarrowBand)
{

	int i = 0;
	int step = m_boundary_nodes.size() / 20;
	while (!ANarrowBand.empty())
	{
		i++;
		if (i == step)
		{
			i = 0;
			writeForDebug();
		}
		//we extract the node with the smallest distance computed in
		StabilityInfo trial_info = extractElectFrom(ANarrowBand);
		Node trial_node = m_mesh->get<Node>(trial_info.nodeID);

		//the trial node becomes alive ...
		m_mesh->mark(trial_node, m_mark_alive);

		//... and is removed from the narrow band
		m_mesh->unmark(trial_node, m_mark_narrow);

		/*	std::cout << "Trial: " << trial_node.getID()
				<< " with stab: " << trial_info.stability << std::endl;*/

		//the quaternion extrapolated in trial_node has been already assigned
		//so, nothing to do...? NO, as trial_node becomes alive, all the nodes
		// in his ball must update their stability info.
		std::vector<TCellID> trial_ball = m_ball[trial_node.getID()];
		for (unsigned int i = 0; i < trial_ball.size(); i++)
		{
			//stability is only computed from nodes in the front
			Node current_node = m_mesh->get<Node>(trial_ball[i]);
			updateStabilityInfo(current_node, m_mark_alive);
		}



		//adjacent nodes that are not alive or in the narrow band are added
		//in the narrow band if they are now adjacent to 2 alive nodes
		std::vector<Node> trial_adj_nodes =
			getAdjacentNodes(trial_node, m_markNodeOnSurf);

		for (unsigned int j = 0; j < trial_adj_nodes.size(); j++)
		{
			Node n_j = trial_adj_nodes[j];

			// we look for nodes that are not yet in the narrow band or alive
			if (m_mesh->isMarked(n_j, m_mark_narrow) ||
				m_mesh->isMarked(n_j, m_mark_alive))
				continue;
			//we check the mark of the nodes adjacent to n_j
			std::vector<Node> adj_nodes_j = getAdjacentNodes(n_j, m_markNodeOnSurf);
			std::vector<Node> ext_nodes;
			int nb_alive_adj_nodes = 0;
			for (unsigned int k = 0; k < adj_nodes_j.size(); k++)
			{
				Node n_k = adj_nodes_j[k];
				if (m_mesh->isMarked(n_k, m_mark_alive) &&
					!m_mesh->isMarked(n_k, m_markNodeOnPnt))
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
					m_mesh->mark(n_j, m_mark_narrow);
					// We can now extrapolate a quaternion in this node... 
					extrapolateQuaternion(n_j, ext_nodes[0], ext_nodes[1], commonFace);


					//we keep in mind the order in which quaternions are created
					m_nodes_by_creation_order.push_back(n_j.getID());
					// then compute a stability value in n_j
					double local_stab = computeStability(n_j, m_mark_alive);
					//weighted by the distance to the curves
					double stab = distancePonderation(local_stab,
						(*m_distance_field_2D)[n_j.getID()]);
					m_stability_value[n_j.getID()] = stab;
					ANarrowBand.push_back(StabilityInfo(n_j.getID(), stab));
				}
			}
			else if (nb_alive_adj_nodes > 2)
			{
				m_mesh->mark(n_j, m_mark_narrow);
				// We can now extrapolate a quaternion in this node... 
				extrapolateQuaternion(n_j, ext_nodes);

				//smoothNodeAndBall(n_j);

				//we keep in mind the order in which quaternions are created
				m_nodes_by_creation_order.push_back(n_j.getID());
				// then compute a stability value in n_j
				double local_stab = computeStability(n_j, m_mark_alive);
				//weighted by the distance to the curves
				double stab = distancePonderation(local_stab,
					(*m_distance_field_2D)[n_j.getID()]);
				m_stability_value[n_j.getID()] = stab;
				ANarrowBand.push_back(StabilityInfo(n_j.getID(), stab));
			}
		}//for (unsigned int k = 0; k < adj_nodes_j.size(); k++)

	}
}
/*----------------------------------------------------------------------------*/
void AdvancingFrontFieldGen3D::buildFrameField()
{
	//==================================================================
	// Step 0 - We initialize the distance ball only for surf elements
	//==================================================================
	std::cout << std::endl << "2D init distance balls" << std::endl;
	initDistanceBallsOnSurf();
	//==================================================================
	// Step 1 - We build the frames on the surfaces
	//==================================================================
	// A node of the narrow band must be adjacent to 2 alive nodes, 
	// while being not alive.
	std::cout << "Init Narrow Band Surf" << std::endl;
	std::list<StabilityInfo>narrow_band;

	initNarrowBandSurf(narrow_band);
	writeForDebug();
	std::cout << "    DONE" << std::endl;


	//======================================================================
	// Step 2 - Advancing front loop on the surface
	//======================================================================
	std::cout << "=============================================================" << std::endl;
	std::cout << " Advancing-front progress of the narrow band on the surface" << std::endl;
	propagateQuaternionsSurf(narrow_band);
	writeForDebug();
	std::cout << "    DONE" << std::endl;

	////======================================================================
	//// Step 2 bis -we smooth the boundary
	////======================================================================
	std::cout << "=============================================================" << std::endl;
	std::cout << "final surface smoothing" << std::endl;
	surfaceSmoothing();
	writeForDebug();
	std::cout << "    DONE" << std::endl;
	//==================================================================
	// Step 0 - We initialize the distance ball for the volume elts
	//==================================================================
	std::cout << std::endl << "3D init distance balls" << std::endl;
	initDistanceBallsInVol();

	//==================================================================
	// Step 3 - We build the frames in the volume
	//==================================================================
	// A node of the narrow band must be adjacent to 2 alive nodes, 
	// while being not alive.
	std::cout << "Init Narrow Band Vol" << std::endl;
	narrow_band.clear();
	initNarrowBand(narrow_band);
	writeForDebug();
	std::cout << "    DONE" << std::endl;


	//======================================================================
	// Step 4 - Advancing front loop in the volume
	//======================================================================
	std::cout << "=============================================================" << std::endl;
	std::cout << " Advancing-front progress of the narrow band in the volume" << std::endl;
	propagateQuaternions(narrow_band);
	writeForDebug();
	std::cout << "    DONE" << std::endl;


	//std::cout << "=============================================================" << std::endl;
	//std::cout << "final volume smoothing" << std::endl;

	//for (int loop = 0; loop < 5; loop++)
	//{

	//	IGMesh::node_iterator it_node = m_mesh->nodes_begin();
	//	for (; !it_node.isDone(); it_node.next())
	//	{
	//		Node current_node = it_node.value();
	//		if (m_mesh->isMarked(current_node, m_mark_alive))
	//			internalSmooth(current_node);

	//	}
	//}
	//writeForDebug();
	//std::cout << "    DONE" << std::endl;
}

/*----------------------------------------------------------------------------*/
void AdvancingFrontFieldGen3D::
surfaceSmoothing()
{
	int mark_smooth = m_mesh->getNewMark<Node>();
	std::vector<Node> smooth_nodes;
	IGMesh::node_iterator it_nodes = m_mesh->nodes_begin();
	for (; !it_nodes.isDone(); it_nodes.next())
	{
		Node n = it_nodes.value();
		if ( m_mesh->isMarked(n, m_markNodeOnSurf) && 
			!m_mesh->isMarked(n, m_markNodeOnCurv) &&
			!m_mesh->isMarked(n, m_markNodeOnPnt))
		{
			smooth_nodes.push_back(n);
			m_mesh->mark(n, mark_smooth);
		}

	}

	smooth(mark_smooth);
	for (unsigned int i = 0; i < smooth_nodes.size(); i++)
	{
		m_mesh->unmark(smooth_nodes[i], mark_smooth);

	}
	m_mesh->freeMark<Node>(mark_smooth);

	//std::map<TCellID, math::Quaternion> quat[2];

	////initialization of quaternions collections
	//for (unsigned int i = 0; i < m_boundary_nodes.size(); i++)
	//{
	//	Node current_node = m_boundary_nodes[i];
	//	quat[0][current_node.getID()] = (*m_cross_field)[current_node.getID()];
	//	quat[1][current_node.getID()] = (*m_cross_field)[current_node.getID()];
	//}

	//int current = 0, next = 1;
	//int nb_loop = 5;
	//for (int loop = 0; loop < nb_loop; loop++)
	//{

	//	for (unsigned int i = 0; i < m_boundary_nodes.size(); i++)
	//	{
	//		Node n = m_boundary_nodes[i];
	//		TCellID n_id = n.getID();
	//		math::Vector normal = (*m_surf_normal)[n_id];

	//		math::Quaternion qref = quat[current][n_id];

	//		if (m_mesh->isMarked(n, m_markNodeOnPnt) ||
	//			m_mesh->isMarked(n, m_markNodeOnCurv))
	//		{
	//			quat[next][n_id] = quat[current][n_id];
	//		}
	//		else
	//		{
	//			//we get the adjacent quaternions associated to nodes that are not
	//			//in the volume
	//			std::vector<Node> adj_nodes;
	//			adj_nodes = getAdjacentNodes(n, m_markNodeOnBnd);
	//			std::vector<math::Quaternion> local_quats;
	//			std::vector<TCoord> local_coefs;

	//			for (int i = 0; i < adj_nodes.size(); i++)
	//			{
	//				Node n = adj_nodes[i];
	//				if (!m_mesh->isMarked(n, m_markNodeOnPnt) &&
	//					m_mesh->isMarked(n, m_mark_alive))
	//				{
	//					math::Quaternion q = quat[current][n.getID()];
	//					q = q.alignWith(normal);
	//					local_quats.push_back(q);
	//					local_coefs.push_back(1.0);
	//				}
	//			}
	//			math::Quaternion res = math::Quaternion::mean(local_quats, local_coefs);
	//			res = res.alignWith(normal);
	//			quat[next][n_id] = res;
	//		}


	//	} //for (unsigned int i = 0; i < m_boundary_nodes.size(); i++)

	//	int tmp_index = current;
	//	current = next;
	//	next = tmp_index;
	//}//	for (int loop = 0; loop < nb_loop; loop++)


	//for (unsigned int i = 0; i < m_boundary_nodes.size(); i++)
	//{
	//	Node current_node = m_boundary_nodes[i];
	//	(*m_cross_field)[current_node.getID()] =
	//		quat[current][current_node.getID()];
	//}
}
/*----------------------------------------------------------------------------*/
void AdvancingFrontFieldGen3D::
initNarrowBand(std::list<StabilityInfo>& ANarrowBand)
{
	ANarrowBand.clear();
	//======================================================================
	// We initialize a quaternion value for each node put in the 
	// narrow band. This value is extrapolated only from the incident alive
	// nodes
	//======================================================================
	IGMesh::node_iterator it_node = m_mesh->nodes_begin();
	for (; !it_node.isDone(); it_node.next())
	{
		Node current_node = it_node.value();
		if (m_mesh->isMarked(current_node, m_mark_alive))
		{
			//We get adjacent nodes that are on the surface
			std::vector<Node>  adj_nodes =
				getAdjacentNodesByRegion(current_node);

			for (unsigned int j = 0; j < adj_nodes.size(); j++)
			{
				Node n_j = adj_nodes[j];
				// we look for nodes that are not yet in the narrow band or alive
				if (m_mesh->isMarked(n_j, m_mark_narrow) ||
					m_mesh->isMarked(n_j, m_mark_alive))
					continue;

				std::vector<Node> adj_nodes_j =
					getAdjacentNodesByRegion(n_j);

				std::vector<Node> ext_nodes;
				int nb_alive_adj_nodes = 0;
				for (unsigned int k = 0; k < adj_nodes_j.size(); k++)
				{
					Node n_k = adj_nodes_j[k];
					if (m_mesh->isMarked(n_k, m_mark_alive) &&
						!m_mesh->isMarked(n_k, m_markNodeOnPnt))
					{
						nb_alive_adj_nodes++;
						ext_nodes.push_back(n_k);
					}
				}
				if (nb_alive_adj_nodes == 3)
				{
					// we have to be sure that n_j, ext_nodes[0], ext_nodes[1],
					// ext_nodes[2] belong to the same face. Iy yes we can put 
					// n_j in the narrow band. 
					Region commonRegion;
					if (belongToTheSameRegion(n_j, ext_nodes[0], ext_nodes[1], ext_nodes[2], commonRegion))
					{
						m_mesh->mark(n_j, m_mark_narrow);
						// We can now extrapolate a quaternion in this node... 
						extrapolateQuaternionInTet(n_j, ext_nodes);
						//we keep in mind the order in which quaternions are created
						m_nodes_by_creation_order.push_back(n_j.getID());
						// then compute a stability value in n_j
						double local_stab = computeStability(n_j, m_mark_alive);
						//weighted by the distance to the boundary
						double stab = distancePonderation(local_stab,
							(*m_distance_field_3D)[n_j.getID()]);
						m_stability_value[n_j.getID()] = stab;
						ANarrowBand.push_back(StabilityInfo(n_j.getID(), stab));
					}

				}//if (nb_alive_adj_nodes == 3)
				else if (nb_alive_adj_nodes > 3)
				{
					std::vector< std::vector<Node> > node_triplets;
					getNodesInTheSameTetrahedron(ext_nodes, node_triplets);
					if (!node_triplets.empty())
					{
						m_mesh->mark(n_j, m_mark_narrow);

						extrapolateQuaternion(n_j, ext_nodes);
						//we keep in mind the order in which quaternions are created
						m_nodes_by_creation_order.push_back(n_j.getID());

						// then compute a stability value in n_j
						double local_stab = computeStability(n_j, m_mark_alive);
						//weighted by the distance to the boundary
						double stab = distancePonderation(local_stab,
							(*m_distance_field_3D)[n_j.getID()]);
						m_stability_value[n_j.getID()] = stab;
						ANarrowBand.push_back(StabilityInfo(n_j.getID(), stab));

					}//if (!node_triplets.empty())

				}//else if (nb_alive_adj_nodes > 3)

			}//for (unsigned int j = 0; j < adj_nodes.size(); j++)

		}//if (m_mesh->isMarked(current_node, AMarkAlive))

	}//for (; !it_node.isDone(); it_node.next())

	std::cout << "Narrow Band Size: " << ANarrowBand.size() << std::endl;
}

/*----------------------------------------------------------------------------*/
void AdvancingFrontFieldGen3D::initDistanceBallsOnSurf()
{

	for (unsigned int i = 0; i < m_boundary_nodes.size(); i++)
	{
		Node ni = m_boundary_nodes[i];
		m_ball[ni.getID()].clear();
		m_ball_location[ni.getID()].clear();
		m_distance_ball[ni.getID()].clear();
		m_smooth_ball[ni.getID()].clear();
	}

	for (unsigned int i = 0; i < m_boundary_nodes.size(); i++)
	{
		Node ni = m_boundary_nodes[i];
		math::Point pi = ni.getPoint();
		double di = (*m_distance_field_2D)[ni.getID()];
	/*	if (di>m_max_distance_2D_div2)
			di = m_max_distance_2D_div2;*/

		double di2 = di*di;
		int ball_size = 0;
		for (unsigned int j = i + 1; j < m_boundary_nodes.size(); j++)
		{
			//i!=j by definition
			Node nj = m_boundary_nodes[j];
			math::Point pj = nj.getPoint();
			double dj = (*m_distance_field_2D)[nj.getID()];
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
		if (ball_size == 0)
		{
			std::vector<Node> adj_nodes = getAdjacentNodes(ni, m_markNodeOnBnd);
			for (unsigned int k = 0; k < adj_nodes.size(); k++)
			{
				putInBalls(ni, adj_nodes[k]);

			}
		}

	}
}
/*----------------------------------------------------------------------------*/
void AdvancingFrontFieldGen3D::putInBalls(Node& AN1, Node& AN2)
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

	// Warning, we initialize the contribution of AN1 
	// to AN2 to -10! It is just a way to field the
	// values at the initialization, but this value
	// is ill-defined. It will have to be updated when
	// a quaternion will be assigned to AN1 and AN2.

	m_smooth_ball[AN1.getID()].push_back(-10);
	m_smooth_ball[AN2.getID()].push_back(-10);

}
/*----------------------------------------------------------------------------*/
void AdvancingFrontFieldGen3D::initDistanceBallsInVol()
{
	std::vector<Node> allNodes;
	allNodes.reserve(m_mesh->getNbNodes());

	IGMesh::node_iterator it = m_mesh->nodes_begin();
	for (; !it.isDone(); it.next())
	{
		allNodes.push_back(it.value());
	}

	for (unsigned int i = 0; i < allNodes.size(); i++)
	{
		Node ni = allNodes[i];
		m_ball[ni.getID()].clear();
		m_ball_location[ni.getID()].clear();
		m_distance_ball[ni.getID()].clear();
		m_smooth_ball[ni.getID()].clear();
	}
	//TODO, Uses ANN structure to fill the ball, it will drastically improve
	// performances of this step
	for (unsigned int i = 0; i < allNodes.size(); i++)
	{
		Node ni = allNodes[i];
		//We only consider nodes that are inside the volume
		if (m_mesh->isMarked(ni, m_markNodeOnBnd))
			continue;

		math::Point pi = ni.getPoint();
		double di = (*m_distance_field_3D)[ni.getID()];

		if (di > m_max_distance_3D_div2)
			di = m_max_distance_3D_div2;

		int ball_size = 0;
		for (unsigned int j = i + 1; j < allNodes.size(); j++)
		{
			Node nj = allNodes[j];
			//i!=j by definition
			math::Point pj = nj.getPoint();
			math::Vector vij(pi, pj);

			if (m_mesh->isMarked(nj, m_markNodeOnBnd))
			{ //we have a boundary node. It only
				if (vij.norm() <= di){
					putInBalls(ni, nj);
					ball_size++;
				}
			}
			else
			{
				double dj = (*m_distance_field_3D)[nj.getID()];
				if (dj > m_max_distance_3D_div2)
					dj = m_max_distance_3D_div2;

				double d = (di > dj) ? di : dj;
				if (vij.norm() <= d)
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
		}//for (unsigned int j = i + 1; j < allNodes.size(); j++)
		if (ball_size == 0){
			std::vector<Node> adj_nodes = getAdjacentNodesByRegion(ni);
			for (unsigned int k = 0; k < adj_nodes.size(); k++)
			{
				putInBalls(ni, adj_nodes[k]);

			}
		}
	}

}

/*---------------------------------------------------------------------------*/
bool AdvancingFrontFieldGen3D::
mustUpdateStabilityInfo(const int ANarrowBandSize)
{
	return true;
}
/*---------------------------------------------------------------------------*/
void AdvancingFrontFieldGen3D::
propagateQuaternions(std::list<StabilityInfo>& ANarrowBand)
{
	int i = 0;
	int step = m_mesh->getNbNodes() / 20;
	while (!ANarrowBand.empty())
	{
		i++;
		if (i == step)
		{
			i = 0;
			writeForDebug();
		}
		//we extract the node with the smallest stability computed in
		StabilityInfo trial_info = extractElectFrom(ANarrowBand);
		Node trial_node = m_mesh->get<Node>(trial_info.nodeID);

		// The quaternion extrapolated in trial_node has been already assigned
		// so, nothing to do...? NO, as trial_node becomes alive, all the nodes
		// in his ball must update their stability info.
		// But this computation is quite expensive, so we not do it too often.
		if (mustUpdateStabilityInfo(ANarrowBand.size()))
		{
			std::vector<TCellID> trial_ball = m_ball[trial_node.getID()];
			for (unsigned int i = 0; i < trial_ball.size(); i++)
			{
				//stability is only computed from nodes in the front
				Node current_node = m_mesh->get<Node>(trial_ball[i]);
				updateStabilityInfo(current_node, m_mark_alive);
			}
		}
		//the trial node becomes alive ...
		m_mesh->mark(trial_node, m_mark_alive);

		//... and is removed from the narrow band
		m_mesh->unmark(trial_node, m_mark_narrow);

		//adjacent nodes that are not alive or in the narrow band are added
		//in the narrow band if they are now adjacent to 2 alive nodes
		std::vector<Node> trial_adj_nodes =
			getAdjacentNodesByRegion(trial_node);

		for (unsigned int j = 0; j < trial_adj_nodes.size(); j++)
		{
			Node n_j = trial_adj_nodes[j];
			// we look for nodes that are not yet in the narrow band or alive
			if (m_mesh->isMarked(n_j, m_mark_narrow) ||
				m_mesh->isMarked(n_j, m_mark_alive))
				continue;
			//we check the mark of the nodes adjacent to n_j
			std::vector<Node> adj_nodes_j = getAdjacentNodesByRegion(n_j);
			std::vector<Node> ext_nodes;
			int nb_alive_adj_nodes = 0;
			for (unsigned int k = 0; k < adj_nodes_j.size(); k++)
			{
				Node n_k = adj_nodes_j[k];
				if (m_mesh->isMarked(n_k, m_mark_alive) &&
					!m_mesh->isMarked(n_k, m_markNodeOnPnt))
				{
					nb_alive_adj_nodes++;
					ext_nodes.push_back(n_k);
				}
			}

			if (nb_alive_adj_nodes == 3)
			{
				// we have to be sure that n_j, ext_nodes[0], ext_nodes[1],
				// ext_nodes[2] belong to the same face. Iy yes we can put 
				// n_j in the narrow band.
				Region commonRegion;
				if (belongToTheSameRegion(n_j, ext_nodes[0], ext_nodes[1], ext_nodes[2], commonRegion))
				{
					m_mesh->mark(n_j, m_mark_narrow);
					// We can now extrapolate a quaternion in this node... 
					extrapolateQuaternionInTet(n_j, ext_nodes);


					//we keep in mind the order in which quaternions are created
					m_nodes_by_creation_order.push_back(n_j.getID());

					// then compute a stability value in n_j
					double local_stab = computeStability(n_j, m_mark_alive);
					//weighted by the distance to the boundary
					double stab = distancePonderation(local_stab,
						(*m_distance_field_3D)[n_j.getID()]);
					m_stability_value[n_j.getID()] = stab;

					ANarrowBand.push_back(StabilityInfo(n_j.getID(), stab));
				}

			}//if (nb_alive_adj_nodes == 3)
			else if (nb_alive_adj_nodes > 3)
			{

				//we look for 2 nodes of ext_nodes that belongs to the same face
				//we extrapolate the distance in n_j
				std::vector< std::vector<Node> > node_triplets;

				getNodesInTheSameTetrahedron(ext_nodes, node_triplets);
				if (!node_triplets.empty())
				{
					m_mesh->mark(n_j, m_mark_narrow);
					extrapolateQuaternion(n_j, ext_nodes);

					//smoothNodeAndBall(n_j);

					//we keep in mind the order in which quaternions are created
					m_nodes_by_creation_order.push_back(n_j.getID());

					// then compute a stability value in n_j
					double local_stab = computeStability(n_j, m_mark_alive);
					//weighted by the distance to the boundary
					double stab = distancePonderation(local_stab,
						(*m_distance_field_3D)[n_j.getID()]);
					m_stability_value[n_j.getID()] = stab;
					ANarrowBand.push_back(StabilityInfo(n_j.getID(), stab));

				}//if (!node_triplets.empty())

			}//else if (nb_alive_adj_nodes > 3)

		}//for (unsigned int k = 0; k < adj_nodes_j.size(); k++)

	}
}
/*----------------------------------------------------------------------------*/
void AdvancingFrontFieldGen3D::smoothAll()
{
	int mark_smooth = m_mesh->getNewMark<Node>();
	IGMesh::node_iterator it_nodes = m_mesh->nodes_begin();
	for (; !it_nodes.isDone(); it_nodes.next())
	{
		Node n = it_nodes.value();
		if (!m_mesh->isMarked(n, m_markNodeOnPnt) &&
			!m_mesh->isMarked(n, m_markNodeOnCurv))
			m_mesh->mark(n, mark_smooth);

	}
	smooth(mark_smooth);

	  m_mesh->unmarkAll<Node>(mark_smooth);
	    m_mesh->freeMark<Node>(mark_smooth);
	      
	      
	      
	      }
/*----------------------------------------------------------------------------*/
	      void AdvancingFrontFieldGen3D::colorSimplices()
	      {
		Variable<int>* var_sing = 0;
		  try{
		    var_sing = m_mesh->getVariable<int>(GMDS_REGION, "sing_tet");
		      }
		  catch (GMDSException& e){
		    var_sing = m_mesh->newVariable<int>(GMDS_REGION, "sing_tet");
		      }
		    
		    IGMesh::region_iterator it = m_mesh->regions_begin();
		      
	//	      int nbOfCluster = 0;
			int nbColoredTet = 0;
			  //=========================================================================
			  // INTERN SKELETON CREATION
			  //=========================================================================
			  // FIRST LOOP ON REGIONS TO GET ALL THE 3-SING. TETS
			  //=========================================================================
			  for (; !it.isDone(); it.next()){
			    Region current_region = it.value();
			      std::vector<Node> nodes = current_region.get<Node>();
				bool onPnt = false;
				  for (unsigned int i_node = 0; i_node < nodes.size(); i_node++)
				    {
				      Node ni = nodes[i_node];
					if (m_mesh->isMarked(ni, m_markNodeOnPnt))
					  onPnt = true;

					  }
				  if (onPnt)
				    {
				      (*var_sing)[current_region.getID()] = 0;
					}
				  else
				    {
				      std::vector<TCellID> nodeIDs = current_region.getIDs<Node>();
					int ID1 = nodeIDs[0];
					  int ID2 = nodeIDs[1];
					    int ID3 = nodeIDs[2];
					      int ID4 = nodeIDs[3];
						
						int singTypeTmp = math::Quaternion::testSingularity((*m_cross_field)[ID1],
													 (*m_cross_field)[ID2],
													 (*m_cross_field)[ID3],
													 (*m_cross_field)[ID4]);
						  if (singTypeTmp != 0)
						    nbColoredTet++;
						    
						    (*var_sing)[current_region.getID()] = singTypeTmp;
						      }
				    }
			    std::cout << "Nb. colored tetrahedra: " << nbColoredTet << std::endl;
			      }

/*---------------------------------------------------------------------------*/
std::vector<gmds::Node> AdvancingFrontFieldGen3D::
computeApproximateMedialAxis(const std::vector<gmds::Node>& ANodes)
{
	std::vector<gmds::Node> medial_axis_nodes;

	std::vector<gmds::Node>::const_iterator it = ANodes.begin();
	for (; it != ANodes.end(); it++)
	{

	}

	return medial_axis_nodes;
}



/*---------------------------------------------------------------------------*/
bool AdvancingFrontFieldGen3D::
computeLocalFrame(gmds::Node& ANode, const int& AMark)
{

	std::vector<Node> adj_nodes;
	std::vector<Face>  adj_faces = ANode.get<Face>();
	bool found_face = false;
	Face current_face;
	for (unsigned int j = 0; !found_face && j < adj_faces.size(); j++)
	{
		adj_nodes.clear();
		Face f_j = adj_faces[j];
		std::vector<Node> nodes_fj = f_j.get<Node>();
		for (unsigned int j = 0; j < nodes_fj.size(); j++)
		{
			Node n_j = nodes_fj[j];
			if (m_mesh->isMarked(n_j, AMark) && n_j.getID() != ANode.getID())
			{
				adj_nodes.push_back(n_j);
			}
		}
		if (adj_nodes.size() == 2)
		{
			found_face = true;
			current_face = f_j;
		}
	}
	if (!found_face)
		return false;
	//		throw GMDSException("ERROR during the advancing front algoritm (on SURFACES)");

	//we mark the node ANode as treated
	m_mesh->mark(ANode, AMark);
	//we compute the quaternion in ANode
	//TODO	extrapolateQuaternion(ANode, adj_nodes[0], adj_nodes[1], current_face);
	return true;
}
/*---------------------------------------------------------------------------*/
bool AdvancingFrontFieldGen3D::belongToTheSameRegion(
	gmds::Node& ATo,
	gmds::Node& AFrom1,
	gmds::Node& AFrom2,
	gmds::Node& AFrom3,
	gmds::Region& ARegion)
{
	std::vector<TCellID> common_regions_12;
	std::vector<TCellID> common_regions_123;
	std::vector<TCellID> common_regions_1234;
	std::vector<TCellID> rs1 = ATo.getIDs<Region>();
	std::vector<TCellID> rs2 = AFrom1.getIDs<Region>();
	for (unsigned int i = 0; i < rs1.size(); i++)
	{
		TCellID id1 = rs1[i];
		bool found = false;
		for (unsigned int j = 0; !found && j < rs2.size(); j++)
		{
			TCellID id2 = rs2[j];
			if (id1 == id2)
			{
				found = true;
				common_regions_12.push_back(id1);
			}
		}
	}
	rs2 = AFrom2.getIDs<Region>();
	for (unsigned int i = 0; i < common_regions_12.size(); i++)
	{
		TCellID id1 = common_regions_12[i];
		bool found = false;
		for (unsigned int j = 0; !found && j < rs2.size(); j++)
		{
			TCellID id2 = rs2[j];
			if (id1 == id2)
			{
				found = true;
				common_regions_123.push_back(id1);
			}
		}
	}
	rs2 = AFrom3.getIDs<Region>();
	for (unsigned int i = 0; i < common_regions_123.size(); i++)
	{
		TCellID id1 = common_regions_123[i];
		bool found = false;
		for (unsigned int j = 0; !found && j < rs2.size(); j++)
		{
			TCellID id2 = rs2[j];
			if (id1 == id2)
			{
				found = true;
				common_regions_1234.push_back(id1);
			}
		}
	}
	if (common_regions_1234.size() == 1){
		ARegion = m_mesh->get<Region>(common_regions_1234[0]);
		return true;
	}
	return false;
}
/*----------------------------------------------------------------------------*/
void AdvancingFrontFieldGen3D::extrapolateQuaternionInTet(
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
		//coeff
		math::Point p1 = current.getPoint();
		math::Vector v01(p0, p1);
		cs.push_back(1.0/v01.norm());
	}

	//===================================================================
	//quaternion extrapolated
	//===================================================================
	math::Quaternion q = math::Quaternion::mean(qs, cs);

	//===================================================================
	//We assign the quaternion align with out_vec
	//===================================================================
	(*m_cross_field)[ATo.getID()] = q;

}
/*----------------------------------------------------------------------------*/
void AdvancingFrontFieldGen3D::smooth(const int AMark)
{

	FrameFieldSmoother smoother(m_mesh, m_cross_field, m_surf_normal);
	smoother.initBoundaryMarks(m_markNodeOnPnt, m_markNodeOnCurv, m_markNodeOnSurf);
	smoother.selectNodes(AMark);
	smoother.execute();

}

/*----------------------------------------------------------------------------*/
void AdvancingFrontFieldGen3D::smoothNodeAndBall(Node& ANode)
{
	int mark_smooth = m_mesh->getNewMark<Node>();
	std::vector<TCellID>& ball = m_ball[ANode.getID()];
	std::vector<Node> smooth_nodes;

	for (unsigned int i = 0; i < ball.size(); i++)
	{
		Node current_node = m_mesh->get<Node>(ball[i]);
		if (m_mesh->isMarked(current_node, m_mark_narrow))
		{
			smooth_nodes.push_back(current_node);
			m_mesh->mark(current_node, mark_smooth);
		}
	}

	smooth(mark_smooth);
	for (unsigned int i = 0; i < smooth_nodes.size(); i++)
	{
		m_mesh->unmark(smooth_nodes[i], mark_smooth);

	}
	m_mesh->freeMark<Node>(mark_smooth);
}
/*----------------------------------------------------------------------------*/
void AdvancingFrontFieldGen3D::internalSmooth(Node& ANode)
{
	bool onSurface = false;
	math::Vector normal;
	if (m_mesh->isMarked(ANode, m_markNodeOnPnt) ||
		m_mesh->isMarked(ANode, m_markNodeOnCurv))
		return;

	if (m_mesh->isMarked(ANode, m_markNodeOnSurf))
	{
		onSurface = true;
		normal = (*m_surf_normal)[ANode.getID()];
	}
	//we get the adjacent quaternions associated to nodes that are not
	//in the volume
	std::vector<Node> adj_nodes;
	if (onSurface)
		adj_nodes = getAdjacentNodes(ANode, m_markNodeOnBnd);
	else
		adj_nodes = getAdjacentNodesByRegion(ANode);

	std::vector<math::Quaternion> quats;
	std::vector<TCoord> coefs;

	for (int i = 0; i < adj_nodes.size(); i++)
	{
		Node n = adj_nodes[i];
		if (!m_mesh->isMarked(n, m_markNodeOnPnt) &&
			m_mesh->isMarked(n, m_mark_alive))
		{
			math::Quaternion q = (*m_cross_field)[n.getID()];
			if (onSurface)
				q = q.alignWith(normal);
			quats.push_back(q);
			coefs.push_back(1.0);
		}
	}
	math::Quaternion res = math::Quaternion::mean(quats, coefs);
	if (onSurface)
		res = res.alignWith(normal);
	(*m_cross_field)[ANode.getID()] = res;

}
/*---------------------------------------------------------------------------*/
