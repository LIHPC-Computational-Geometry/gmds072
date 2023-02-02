/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * The GMDS library is a computer program whose purpose is to provide a set of
 * functionnalities to represent and handle any type of meshes (2D, 3D,
 * triangles, tetrahedra, quad, hexa, polygons, polyhedra, etc.) and write
 * meshing algorithms. So it gathers many mathematical objects like points,
 * segment, quaternions, etc. and basic algorithms useful to build more evolved
 * ones.
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
/** \file    DistanceFieldBuilder3D.cpp
 *  \author  F. LEDOUX
 *  \date    08/31/2010
 */
/*----------------------------------------------------------------------------*/
#include <map>
#include <set>
/*----------------------------------------------------------------------------*/
#include <GMDS/Algo/DistanceFieldBuilder3D.h>
#include <GMDS/Math/Matrix.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
DistanceFieldBuilder3D::DistanceFieldBuilder3D(IGMesh* AMesh)
:m_mesh(AMesh)
{
}
/*----------------------------------------------------------------------------*/
DistanceFieldBuilder3D::~DistanceFieldBuilder3D()
{}
/*----------------------------------------------------------------------------*/
bool DistanceFieldBuilder3D::isValid()
{
	//R2N and N2R connectivites are required
	MeshModel model = m_mesh->getModel();
	if (!model.has(N2R))
	{
		std::cout << "Error in the32D distance field computation, N2R connectivity must be available" << std::endl;
		return false;
	}

	//The mesh regions must be tetrahedra only
	IGMesh::region_iterator it_r = m_mesh->regions_begin();
	for (; !it_r.isDone(); it_r.next()){
		Region r = it_r.value();
		if (r.getType() != GMDS_TETRA)
		{
			std::cout << "Error in the 3D distance field computation, ";
			std::cout << "only tetrahedral regions are supported" << std::endl;
			return false;
		}
	}
	return true;
}
/*----------------------------------------------------------------------------*/
std::vector<Node> DistanceFieldBuilder3D::
getAdjacentNodes(Node& ANode, const int AMark)
{
	std::set<Node>      adj_nodes;
	std::vector<Region> adj_regions = ANode.get<Region>();

	for (unsigned int j = 0; j < adj_regions.size(); j++)
	{
		Region r_j = adj_regions[j];
		std::vector<Node> nodes_rj = r_j.get<Node>();
		for (unsigned int j = 0; j < nodes_rj.size(); j++)
		{
			Node n_j = nodes_rj[j];
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
std::vector<TCellID> DistanceFieldBuilder3D::getInsertionOrder() const
{
	return m_insertion_order;
}
/*----------------------------------------------------------------------------*/
double DistanceFieldBuilder3D::
extrapolateDistance(Node& AN, Node& AN1, Node& AN2, Node& AN3,
Variable<double>*& ADist)
{
	//std::cout << "# Compute 3D distance in " << AN.getID() << " from "
	//	<< AN1.getID() << ", "
	//	<< AN2.getID() << " and "
	//	<< AN3.getID();
	double d1 = (*ADist)[AN1.getID()];
	double d2 = (*ADist)[AN2.getID()];
	double d3 = (*ADist)[AN3.getID()];
//	std::cout << " with " << d1 << ", " << d2 << ", " << d3 << std::endl;
	double J_tab[3][3];

	math::Point p4 = AN.getPoint();
	math::Point p1 = AN1.getPoint();
	math::Point p2 = AN2.getPoint();
	math::Point p3 = AN3.getPoint();
	//first line
	J_tab[0][0] = p2.X() - p1.X();
	J_tab[0][1] = p3.X() - p1.X();
	J_tab[0][2] = p4.X() - p1.X();
	//second line
	J_tab[1][0] = p2.Y() - p1.Y();
	J_tab[1][1] = p3.Y() - p1.Y();
	J_tab[1][2] = p4.Y() - p1.Y();
	//third line
	J_tab[2][0] = p2.Z() - p1.Z();
	J_tab[2][1] = p3.Z() - p1.Z();
	J_tab[2][2] = p4.Z() - p1.Z();

	// matrix inversion
	math::Matrix<3, 3,double> J(J_tab);
	math::Matrix<3, 3,double> J_inv = J.inverse();
	//first column
	double a2 = J_inv.get(0, 0);
	double a3 = J_inv.get(1, 0);
	double a4 = J_inv.get(2, 0);
	//second column
	double b2 = J_inv.get(0, 1);
	double b3 = J_inv.get(1, 1);
	double b4 = J_inv.get(2, 1);
	//third column
	double c2 = J_inv.get(0, 2);
	double c3 = J_inv.get(1, 2);
	double c4 = J_inv.get(2, 2);
	double d2_bar = d2 - d1;
	double d3_bar = d3 - d1;

	double gamma_x = a2*d2_bar + a3*d3_bar;
	double gamma_y = b2*d2_bar + b3*d3_bar;
	double gamma_z = c2*d2_bar + c3*d3_bar;

	double A = a4*a4 + b4*b4 + c4*c4;
	double B = 2 * (gamma_x*a4 + gamma_y*b4 + gamma_z*c4);
	double C = gamma_x*gamma_x + gamma_y*gamma_y + gamma_z*gamma_z - 1;

	double A2 = 2 * A;
	double DELTA_SQ = B*B - 4 * A*C;
	if (DELTA_SQ < 1e-5)
		DELTA_SQ = 0;

	double d4 = 0;
	if (DELTA_SQ<0)
	{
		d4 = -1;
	}
	else
	{
		double DELTA = sqrt(DELTA_SQ);
		double root1 = (-B + DELTA) / A2;
		double root2 = (-B - DELTA) / A2;
		double d4_bar = 0;
		if (root1 > root2)
		{
			d4_bar = root1;
		}
		else
		{
			d4_bar = root2;
		}
		//and now the absolute distance is
		d4 = d1 + d4_bar;
	}
//	std::cout << "      ==> " << d4 << std::endl;
	return d4;
}
/*---------------------------------------------------------------------------*/
bool DistanceFieldBuilder3D::belongToTheSameTetrahedron(
	Node& AN1,
	Node& AN2,
	Node& AN3,
	Node& AN4,
	Region& AR)
{
	std::vector<TCellID> common_12;
	std::vector<TCellID> common_123;
	std::vector<TCellID> common_1234;
	std::vector<TCellID> rs1 = AN1.getIDs<Region>();
	std::vector<TCellID> rs2 = AN2.getIDs<Region>();
	std::vector<TCellID> rs3 = AN3.getIDs<Region>();
	std::vector<TCellID> rs4 = AN4.getIDs<Region>();

	//=============================================================
	// STEP 1 - we get the regions common to AN1 and AN2
	//=============================================================
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
				common_12.push_back(id1);
			}
		}
	}
	if (common_12.empty())
		return false;

	//=============================================================
	// STEP 2 - we get the regions common to AN1, AN2 and AN3
	//=============================================================
	for (unsigned int i = 0; i < common_12.size(); i++)
	{
		TCellID id1 = common_12[i];
		bool found = false;
		for (unsigned int j = 0; !found && j < rs3.size(); j++)
		{
			TCellID id2 = rs3[j];
			if (id1 == id2)
			{
				found = true;
				common_123.push_back(id1);
			}
		}
	}
	if (common_123.empty())
		return false;


	//=============================================================
	// STEP 3 - we get the regions common to AN1, AN2, AN3 and AN4
	//=============================================================
	for (unsigned int i = 0; i < common_123.size(); i++)
	{
		TCellID id1 = common_123[i];
		bool found = false;
		for (unsigned int j = 0; !found && j < rs4.size(); j++)
		{
			TCellID id2 = rs4[j];
			if (id1 == id2)
			{
				found = true;
				common_1234.push_back(id1);
			}
		}
	}
	if (common_1234.empty())
		return false;

	if (common_1234.size() == 1){
		AR = m_mesh->get<Region>(common_1234[0]);
		return true;
	}
	return false;
}
/*----------------------------------------------------------------------------*/
void DistanceFieldBuilder3D::getNodesInTheSameTetrahedron(
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
void DistanceFieldBuilder3D::initNarrowBand(
	std::vector<Node>& AFrom,
	std::multimap<double,Node>& ANarrowBand,
	std::map<TCellID, double>& AExtrapolateDistance,
	const int AMarkNodeToWorkOn,
	const int AMarkNodeAlive,
	const int AMarkNodeNarrow)
{


	for (unsigned int i = 0; i < AFrom.size(); i++)
	{
		Node n_i = AFrom[i];
		std::vector<Node> adj_nodes = getAdjacentNodes(n_i, AMarkNodeToWorkOn);
		for (unsigned int j = 0; j < adj_nodes.size(); j++)
		{
			Node n_j = adj_nodes[j];
			// we look for nodes that are not yet in the narrow band or alive
			if (m_mesh->isMarked(n_j, AMarkNodeNarrow) || m_mesh->isMarked(n_j, AMarkNodeAlive))
				continue;

			std::vector<Node> adj_nodes_j = getAdjacentNodes(n_j, AMarkNodeToWorkOn);
			std::vector<Node> ext_nodes;
			int nb_alive_adj_nodes = 0;
			for (unsigned int k = 0; k < adj_nodes_j.size(); k++)
			{
				Node n_k = adj_nodes_j[k];
				if (m_mesh->isMarked(n_k, AMarkNodeAlive))
				{
					nb_alive_adj_nodes++;
					ext_nodes.push_back(n_k);
				}
			}
			if (nb_alive_adj_nodes == 3)
			{
				Region commonTet;
				if (belongToTheSameTetrahedron(n_j, ext_nodes[0], ext_nodes[1],
					ext_nodes[2], commonTet))
				{
					m_mesh->mark(n_j, AMarkNodeNarrow);
					//we extrapolate the distance in n_j
					double d_j = extrapolateDistance(n_j, ext_nodes[0], ext_nodes[1],
						ext_nodes[2], m_distance);
					if (d_j < 0)
					{
						d_j = 0;
						std::cout << "Warning, we get a negative distance"
							<< std::endl;
						AExtrapolateDistance[n_j.getID()] = 10000;
						ANarrowBand.insert(std::pair<double, Node>(d_j, n_j));
					}
					else
					{
						AExtrapolateDistance[n_j.getID()] = d_j;
						ANarrowBand.insert(std::pair<double, Node>(d_j, n_j));
					}
				}
			} //if (nb_alive_adj_nodes == 3)
			else if (nb_alive_adj_nodes > 3)
			{

				//we look for 2 nodes of ext_nodes that belongs to the same face
				//we extrapolate the distance in n_j
				std::vector< std::vector<Node> > node_triplets;

				getNodesInTheSameTetrahedron(ext_nodes, node_triplets);
				if (!node_triplets.empty())
				{
					m_mesh->mark(n_j, AMarkNodeNarrow);

					Node from1 = node_triplets[0][0];
					Node from2 = node_triplets[0][1];
					Node from3 = node_triplets[0][2];

					double d_j = extrapolateDistance(n_j, from1, from2, from3, m_distance);
					if (d_j < 0)
					{
						d_j = 0;
						std::cout << "Warning, we get a negative distance"
							<< std::endl;
						AExtrapolateDistance[n_j.getID()] = 10000;
						ANarrowBand.insert(std::pair<double, Node>(d_j, n_j));
					}
					else
					{
						AExtrapolateDistance[n_j.getID()] = d_j;
						ANarrowBand.insert(std::pair<double, Node>(d_j, n_j));
					}
				}//if (!node_triplets.empty())

			}//else if (nb_alive_adj_nodes > 3)
		}
	}
	//	std::cout << "Narrow Band Size: " << ANarrowBand.size() << std::endl;
}

/*----------------------------------------------------------------------------*/
void DistanceFieldBuilder3D::advanceDistanceFront(
	std::multimap<double,Node>& ANarrowBand,
	std::map<TCellID, double>& AExtrapolateDistance,
	const int AMarkNodeToWorkOn,
	const int AMarkNodeAlive,
	const int AMarkNodeNarrow)
{
	while (!ANarrowBand.empty())
	{
		//we take the node with the smallest distance computed in
		std::multimap<double, Node>::iterator it = ANarrowBand.begin();
		Node trial_node = it->second;
		double trial_dist = it->first;
		(*m_distance)[trial_node.getID()] = trial_dist;
		m_insertion_order.push_back(trial_node.getID());
		ANarrowBand.erase(it);

		//the trial node becomes alive ...
		m_mesh->mark(trial_node, AMarkNodeAlive);
		//... and is removed from the narrow band
		m_mesh->unmark(trial_node, AMarkNodeNarrow);

		//adjacent nodes that are not alive or in tha narrow band are added
		//in the narrow band if they are now adjacent to 2 alive nodes
		std::vector<Node> trial_adj_nodes = getAdjacentNodes(trial_node, AMarkNodeToWorkOn);

		for (unsigned int j = 0; j < trial_adj_nodes.size(); j++)
		{
			Node n_j = trial_adj_nodes[j];
			// we look for nodes that are not yet in the narrow band or alive
			if (m_mesh->isMarked(n_j, AMarkNodeNarrow) || m_mesh->isMarked(n_j, AMarkNodeAlive))
				continue;

			std::vector<Node> adj_nodes_j = getAdjacentNodes(n_j, AMarkNodeToWorkOn);
			std::vector<Node> ext_nodes;
			int nb_alive_adj_nodes = 0;
			for (unsigned int k = 0; k < adj_nodes_j.size(); k++)
			{
				Node n_k = adj_nodes_j[k];
				if (m_mesh->isMarked(n_k, AMarkNodeAlive))
				{
					nb_alive_adj_nodes++;
					ext_nodes.push_back(n_k);
				}
			}

			if (nb_alive_adj_nodes == 3)
			{
				// we have to be sure that n_j, ext_nodes[0], ext_nodes[1]
				// belong to the same face. Iy yes we can put n_j in the 
				// narrow band.
				Region commonTet;
				if (belongToTheSameTetrahedron(n_j, ext_nodes[0], ext_nodes[1],
					ext_nodes[2], commonTet))
				{
					m_mesh->mark(n_j, AMarkNodeNarrow);

					//we extrapolate the distance in n_j
					//std::cout << "In tet: " << commonTet.getID() << std::endl;
					double d_j = extrapolateDistance(n_j, ext_nodes[0], ext_nodes[1],
						ext_nodes[2], m_distance);
					if (d_j < 0)
					{
						d_j = 0;
						std::cout << "Warning, we get a negative distance"
							<< std::endl;
						AExtrapolateDistance[n_j.getID()] = 10000;
						ANarrowBand.insert(std::pair<double, Node>(d_j, n_j));
					}
					else
					{
						AExtrapolateDistance[n_j.getID()] = d_j;
						ANarrowBand.insert(std::pair<double, Node>(d_j, n_j));
					}

				}
			} //if (nb_alive_adj_nodes == 3)
			else if (nb_alive_adj_nodes > 3)
			{

				//we look for 2 nodes of ext_nodes that belongs to the same face
				//we extrapolate the distance in n_j
				std::vector< std::vector<Node> > node_triplets;

				getNodesInTheSameTetrahedron(ext_nodes, node_triplets);
				if (!node_triplets.empty())
				{
					m_mesh->mark(n_j, AMarkNodeNarrow);

					Node from1 = node_triplets[0][0];
					Node from2 = node_triplets[0][1];
					Node from3 = node_triplets[0][2];

					double d_j = extrapolateDistance(n_j, from1, from2, from3, m_distance);
					if (d_j < 0)
					{
						d_j = 0;
						std::cout << "Warning, we get a negative distance"
							<< std::endl;
						AExtrapolateDistance[n_j.getID()] = 10000;
						ANarrowBand.insert(std::pair<double, Node>(d_j, n_j));
					}
					else
					{
						AExtrapolateDistance[n_j.getID()] = d_j;
						ANarrowBand.insert(std::pair<double, Node>(d_j, n_j));
					}

				}//if (!node_triplets.empty())

			}//else if (nb_alive_adj_nodes > 3)
		}
	}
}
/*----------------------------------------------------------------------------*/
Variable<double>* DistanceFieldBuilder3D::
computeDistance(std::vector<Node>& AFrom,
std::vector<Node>& AToCompute)
{
	m_insertion_order.clear();
	m_distance = m_mesh->newVariable<double>(GMDS_NODE, "distance3D");

	//======================================================================
	// Boolean marks definition
	//======================================================================
	// a mark is given to the node to be computed
	int node_to_use = m_mesh->getNewMark<Node>();
	//all the nodes a value is computed are marked node_alive
	int node_alive = m_mesh->getNewMark<Node>();
	//all the nodes that are in the narrow band
	int node_narrow = m_mesh->getNewMark<Node>();

	//======================================================================
	// Initialization - ALIVE NODES
	//======================================================================
	//we mark all the nodes that we want to compute a distance on
	for (unsigned int i = 0; i < AToCompute.size(); i++)
	{
		Node n_i = AToCompute[i];
		m_mesh->mark(n_i, node_to_use);
	}

	for (unsigned int i = 0; i < AFrom.size(); i++)
	{
		Node n_i = AFrom[i];
		//initial nodes are evaluated to zero
		(*m_distance)[n_i.getID()] = 0;
		// the node is also marked as done
		m_mesh->mark(n_i, node_alive);
		m_mesh->mark(n_i, node_to_use);
	}
	// now, all the nodes where a distance is assigned are marked alive,
	// and all the others are not marked alive.

	//======================================================================
	// Initialization - NARROW BAND NODES
	//======================================================================
	// We build the narrow band list, which contains all the nodes we are going
	// to work on. They are adjacent to 3 alive nodes, while being
	// not alive.
	std::multimap<double,Node> narrow_band;
	// we keep in mind the extrapolate distance for any node in the narrow
	// band
	std::map<TCellID, double> extrapolate_distance;

	initNarrowBand(AFrom, narrow_band, extrapolate_distance,
		node_to_use, node_alive, node_narrow);

	//======================================================================
	// Advancing front loop
	//======================================================================
	advanceDistanceFront(narrow_band, extrapolate_distance,
		node_to_use, node_alive, node_narrow);


	//======================================================================
	// Cleaning of the Boolean marks
	//======================================================================
	m_mesh->unmarkAll<Node>(node_alive);
	m_mesh->unmarkAll<Node>(node_narrow);
	m_mesh->unmarkAll<Node>(node_to_use);
	m_mesh->freeMark<Node>(node_alive);
	m_mesh->freeMark<Node>(node_narrow);
	m_mesh->freeMark<Node>(node_to_use);

	return m_distance;
}
/*----------------------------------------------------------------------------*/
