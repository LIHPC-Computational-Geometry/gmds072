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
// STL Header files
#include <iostream>
#include <map>
#include <set>
/*---------------------------------------------------------------------------*/
#include "SmoothingHLBFGS.h"
/*---------------------------------------------------------------------------*/
// GMDS Header files
#include "GMDS/Math/Matrix.h"/*---------------------------------------------------------------------------*/
// HLBFGS Header files
#include "HLBFGS.h"
/*---------------------------------------------------------------------------*/
using namespace gmds;

/*----------------------------------------------------------------------------*/
/*
* IN N -> nb de variables a traiter (3*nb sommets pour nous)
* IN x -> taille N contient les angles d'Euler pour chaque sommet à traiter
*  		initialisé avant
* IN prev_x -> NON UTILISE valeur précédente de x sinon
* OUT func-> pointeur sur la fonction calculée
* OUT grad -> grad(func)
*/
void SmoothingHLBFGS::
evalF(int N, double* x, double *prev_x, double* func, double* grad)
{
	size_t nvars = N;

	//grad est deja de taille N car pointe sur un vecteur de double initialisé
	//avec la bonne taille dans HLBFGS

	//INIT DES FONCTIONS TRIGO POUR CALCUL DIFFERENTIELS A VENIR
	std::vector<TCoord> trig_cos_vectors, trig_sin_vectors;
	trig_cos_vectors.resize(N);
	trig_sin_vectors.resize(N);

	for (unsigned int i = 0; i < N; i++){
		trig_cos_vectors[i] = cos(x[i]);
		trig_sin_vectors[i] = sin(x[i]);
	}

	*func = 0.0;

	// DOIT ON INITIALISER g?? RIEN DE DIT DANS HLBFGS TODO
		for (unsigned int i = 0;i < N;i++)
			grad[i] = 0.0;


	// PARCOURS DES ARETES POUR INITIALISER LES MATRICES
	for (unsigned int i = 0; i < m_edges.size(); i++)
	{
		std::vector<Node> edge_nodes = m_edges[i].get<Node>();
		Node n0 = edge_nodes[0];
		Node n1 = edge_nodes[1];

		int index0 = m_candidates_index[n0.getID()];
		int index1 = m_candidates_index[n1.getID()];

		//		std::cout<<"premier vertex en pos "<<vi->getX()<<" "<<vi->getY()<<" "<<vi->getZ()<<std::endl;
		//		std::cout<<"second vertex en pos "<<vj->getX()<<" "<<vj->getY()<<" "<<vj->getZ()<<std::endl;



		gmds::math::Matrix<3, 3, double> Mix, Miy, Miz;
		gmds::math::Matrix<3, 3, double> Mjx, Mjy, Mjz;
		gmds::math::Matrix<3, 3, double> DMix, DMiy, DMiz;
		gmds::math::Matrix<3, 3, double> DMjx, DMjy, DMjz;
		gmds::math::Matrix<3, 3, double> M, Mi, Mj, Mi1, Mi2, Mi3;
		// toutes les matrices 3x3 precedentes ne contiennent que des 0
		TCoord diff1[3][3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
		TCoord diff2[3][3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
		TCoord diffi[3] = { 0, 0, 0 };
		TCoord diffj[3] = { 0, 0, 0 };
		TCoord e[18];

		// Matrix computation for node n0
		initMatrix(n0, Mix, Miy, Miz, DMix, DMiy, DMiz, trig_cos_vectors, trig_sin_vectors);

		Mi = Mix * Miy * Miz;
		Mi1 = DMix * Miy * Miz;
		Mi2 = Mix * DMiy * Miz;
		Mi3 = Mix * Miy * DMiz;

		// Matrix computation for node n1
		initMatrix(n1, Mjx, Mjy, Mjz, DMjx, DMjy, DMjz, trig_cos_vectors, trig_sin_vectors);

		Mj = Mjx * Mjy * Mjz;

		//ATTENTION, il est possible que ce soit Mi a transposer dans les lignes
		//suivantes. Idem Mi1, Mi2, Mi3
		//		std::cout<<"Mi "<<Mi<<std::endl;
		//		std::cout<<"Mj "<<Mj<<std::endl;
		//M = Mi * Mj.transpose();
		M = Mi.transpose() * Mj;

		//		M = Mj.transpose() * Mi;
		//		std::cout<<"M "<<M<<std::endl;

		//		gmds::math::Matrix<reel,3,3> tmpMat = (Mi * (DMjx * Mjy * Mjz).transpose());
		gmds::math::Matrix<3, 3, double> tmpMat = (Mi.transpose() * (DMjx * Mjy * Mjz));
		//		gmds::math::Matrix<reel,3,3> tmpMat = ((DMjx * Mjy * Mjz).transpose() * Mi);
		tmpMat.getTab(diff2[0]);

		//		tmpMat = ((Mjx * DMjy * Mjz).transpose() * Mi);
		tmpMat = (Mi.transpose() * (Mjx * DMjy * Mjz));
		//		tmpMat = (Mi * (Mjx * DMjy * Mjz).transpose());
		tmpMat.getTab(diff2[1]);

		//		tmpMat = (Mi * (Mjx * Mjy * DMjz).transpose());
		tmpMat = (Mi.transpose() * (Mjx * Mjy * DMjz));
		//		tmpMat = ((Mjx * Mjy * DMjz).transpose() * Mj);
		tmpMat.getTab(diff2[2]);

		//		(Mj.transpose() * Mi1).getTab(diff1[0]);
		//		(Mj.transpose() * Mi2).getTab(diff1[1]);
		//		(Mj.transpose() * Mi3).getTab(diff1[2]);
		(Mi1.transpose() * Mj).getTab(diff1[0]);
		(Mi2.transpose() * Mj).getTab(diff1[1]);
		(Mi3.transpose() * Mj).getTab(diff1[2]);
		//		(Mi1 * Mj.transpose()).getTab(diff1[0]);
		//		(Mi2 * Mj.transpose()).getTab(diff1[1]);
		//		(Mi3 * Mj.transpose()).getTab(diff1[2]);

		//		diff2[0] = Mi * (DMjx * Mjy * Mjz).transpose();
		//		diff2[1] = Mi * (Mjx * DMjy * Mjz).transpose();
		//		diff2[2] = Mi * (Mjx * Mjy * DMjz).transpose();
		//		diff1[0] = Mi1 * Mj.transpose();
		//		diff1[1] = Mi2 * Mj.transpose();
		//		diff1[2] = Mi3 * Mj.transpose();

		e[0] = M.get(0, 0) * M.get(0, 1); e[1] = M.get(0, 0) * M.get(0, 2); e[2] = M.get(0, 1) * M.get(0, 2);
		e[3] = M.get(1, 0) * M.get(1, 1); e[4] = M.get(1, 0) * M.get(1, 2); e[5] = M.get(1, 1) * M.get(1, 2);
		e[6] = M.get(2, 0) * M.get(2, 1); e[7] = M.get(2, 0) * M.get(2, 2); e[8] = M.get(2, 1) * M.get(2, 2);
		e[9] = M.get(0, 0) * M.get(1, 0); e[10] = M.get(0, 0) * M.get(2, 0); e[11] = M.get(1, 0) * M.get(2, 0);
		e[12] = M.get(0, 1) * M.get(1, 1); e[13] = M.get(0, 1) * M.get(2, 1); e[14] = M.get(1, 1) * M.get(2, 1);
		e[15] = M.get(0, 2) * M.get(1, 2); e[16] = M.get(0, 2) * M.get(2, 2); e[17] = M.get(1, 2) * M.get(2, 2);

		double m_sigma = 1.0;

		double local_f = 0.0;


		static const int eflag[9][2][2] = {
			{ { 0, 0 }, { 0, 1 } }, { { 0, 0 }, { 0, 2 } }, { { 0, 1 }, { 0, 2 } },
			{ { 1, 0 }, { 1, 1 } }, { { 1, 0 }, { 1, 2 } }, { { 1, 1 }, { 1, 2 } },
			{ { 2, 0 }, { 2, 1 } }, { { 2, 0 }, { 2, 2 } }, { { 2, 1 }, { 2, 2 } } };


		for (int k = 0; k < 9; k++)
		{
			double v2 = exp(e[k + 9] * e[k + 9] / m_sigma);

			local_f += 2 * (v2 - 1);
			for (int l = 0; l < 3; l++)
			{

				diffi[l] += 2 * (e[k + 9] * (diff1[l][eflag[k][0][1]][eflag[k][0][0]] * M.get(eflag[k][1][1], eflag[k][1][0]) + diff1[l][eflag[k][1][1]][eflag[k][1][0]] * M.get(eflag[k][0][1], eflag[k][0][0])) * v2 / m_sigma);
				diffj[l] += 2 * (e[k + 9] * (diff2[l][eflag[k][0][1]][eflag[k][0][0]] * M.get(eflag[k][1][1], eflag[k][1][0]) + diff2[l][eflag[k][1][1]][eflag[k][1][0]] * M.get(eflag[k][0][1], eflag[k][0][0])) * v2 / m_sigma);
			}
		}


		*func += (local_f / 2.0);
		for (int k = 0; k < 3; k++)
		{
			grad[3 * index0 + k] += diffi[k];
			grad[3 * index1 + k] += diffj[k];
		}

	}//FIN DE LA BOUCLE SUR LES ARETES
	//return func/2;
}
/*---------------------------------------------------------------------------*/
void SmoothingHLBFGS::initMatrix(Node& ANode,
	gmds::math::Matrix<3, 3, double>& Mix,
	gmds::math::Matrix<3, 3, double>& Miy,
	gmds::math::Matrix<3, 3, double>& Miz,
	gmds::math::Matrix<3, 3, double>& DMix,
	gmds::math::Matrix<3, 3, double>& DMiy, 
	gmds::math::Matrix<3, 3, double>& DMiz,
	std::vector<TCoord>& trig_cos_vectors,
	std::vector<TCoord>& trig_sin_vectors)
{
	Node n = ANode;
	int index = m_candidates_index[n.getID()];
	// contribution du sommet 0
	if (m_mesh->isMarked(n, m_mark_curve))
	{
		//point sur une arete vive
		Mix.set(0, 0, 1);
		Mix.set(1, 1, trig_cos_vectors[3 * index]);
		Mix.set(2, 2, trig_cos_vectors[3 * index]);
		Mix.set(1, 2, -trig_sin_vectors[3 * index]);
		Mix.set(2, 1, trig_sin_vectors[3 * index]);

		Miy.set(1, 1, 1);
		Miy.set(0, 0, trig_cos_vectors[3 * index + 1]);
		Miy.set(2, 2, trig_cos_vectors[3 * index + 1]);
		Miy.set(0, 2, trig_sin_vectors[3 * index + 1]);
		Miy.set(2, 0, -trig_sin_vectors[3 * index + 1]);

		Miz.set(2, 2, 1);
		Miz.set(0, 0, trig_cos_vectors[3 * index + 2]);
		Miz.set(1, 1, trig_cos_vectors[3 * index + 2]);
		Miz.set(0, 1, -trig_sin_vectors[3 * index + 2]);
		Miz.set(1, 0, trig_sin_vectors[3 * index + 2]);

		// valeurs stationnaires pour le processus donc derivees nulles
	}
	else if (m_mesh->isMarked(n, m_mark_surface))
	{
		Mix.set(0, 0, 1.0);
		Mix.set(1, 1, 1.0);
		Mix.set(2, 2, 1.0);

		Miy.set(0, 0, 1.0);
		Miy.set(1, 1, 1.0);
		Miy.set(2, 2, 1.0);

		double VX[3];
		double VY[3];
		math::Chart t = m_boundary_triad[n.getID()];
		math::Vector tX = t.X();
		math::Vector tY = t.Y();
		math::Vector tZ = t.Z();

		VX[0] = trig_cos_vectors[3 * index + 2] * tY[0] + trig_sin_vectors[3 * index + 2] * tZ[0];
		VX[1] = trig_cos_vectors[3 * index + 2] * tY[1] + trig_sin_vectors[3 * index + 2] * tZ[1];
		VX[2] = trig_cos_vectors[3 * index + 2] * tY[2] + trig_sin_vectors[3 * index + 2] * tZ[2];

		VY[0] = -trig_sin_vectors[3 * index + 2] * tY[0] + trig_cos_vectors[3 * index + 2] * tZ[0];
		VY[1] = -trig_sin_vectors[3 * index + 2] * tY[1] + trig_cos_vectors[3 * index + 2] * tZ[1];
		VY[2] = -trig_sin_vectors[3 * index + 2] * tY[2] + trig_cos_vectors[3 * index + 2] * tZ[2];

		Miz.set(0, 0, VX[0]);
		Miz.set(1, 0, VX[1]);
		Miz.set(2, 0, VX[2]);

		Miz.set(0, 1, VY[0]);
		Miz.set(1, 1, VY[1]);
		Miz.set(2, 1, VY[2]);

		Miz.set(0, 2, tX[0]);
		Miz.set(1, 2, tX[1]);
		Miz.set(2, 2, tX[2]);

		DMiz.set(0, 0, -trig_sin_vectors[3 * index + 2] * tY[0] + trig_cos_vectors[3 * index + 2] * tZ[0]);
		DMiz.set(1, 0, -trig_sin_vectors[3 * index + 2] * tY[1] + trig_cos_vectors[3 * index + 2] * tZ[1]);
		DMiz.set(2, 0, -trig_sin_vectors[3 * index + 2] * tY[2] + trig_cos_vectors[3 * index + 2] * tZ[2]);

		DMiz.set(0, 1, -trig_cos_vectors[3 * index + 2] * tY[0] - trig_sin_vectors[3 * index + 2] * tZ[0]);
		DMiz.set(1, 1, -trig_cos_vectors[3 * index + 2] * tY[1] - trig_sin_vectors[3 * index + 2] * tZ[1]);
		DMiz.set(2, 1, -trig_cos_vectors[3 * index + 2] * tY[2] - trig_sin_vectors[3 * index + 2] * tZ[2]);

		DMiz.set(0, 2, 0.0);
		DMiz.set(1, 2, 0.0);
		DMiz.set(2, 2, 0.0);

	}
	else{ //CAS GENERAL
		//point pas sur le bord, remplissage classique des matrices
		//			std::cout<<"vi inside"<<std::endl;
		Mix.set(0, 0, 1);
		Mix.set(1, 1, trig_cos_vectors[3 * index]);
		Mix.set(2, 2, trig_cos_vectors[3 * index]);
		Mix.set(1, 2, -trig_sin_vectors[3 * index]);
		Mix.set(2, 1, trig_sin_vectors[3 * index]);

		DMix.set(1, 1, -trig_sin_vectors[3 * index]);
		DMix.set(2, 2, -trig_sin_vectors[3 * index]);
		DMix.set(1, 2, -trig_cos_vectors[3 * index]);
		DMix.set(2, 1, trig_cos_vectors[3 * index]);

		Miy.set(1, 1, 1);
		Miy.set(0, 0, trig_cos_vectors[3 * index + 1]);
		Miy.set(2, 2, trig_cos_vectors[3 * index + 1]);
		Miy.set(0, 2, trig_sin_vectors[3 * index + 1]);
		Miy.set(2, 0, -trig_sin_vectors[3 * index + 1]);

		DMiy.set(0, 0, -trig_sin_vectors[3 * index + 1]);
		DMiy.set(2, 2, -trig_sin_vectors[3 * index + 1]);
		DMiy.set(0, 2, trig_cos_vectors[3 * index + 1]);
		DMiy.set(2, 0, -trig_cos_vectors[3 * index + 1]);

		Miz.set(2, 2, 1);
		Miz.set(0, 0, trig_cos_vectors[3 * index + 2]);
		Miz.set(1, 1, trig_cos_vectors[3 * index + 2]);
		Miz.set(0, 1, -trig_sin_vectors[3 * index + 2]);
		Miz.set(1, 0, trig_sin_vectors[3 * index + 2]);

		DMiz.set(0, 0, -trig_sin_vectors[3 * index + 2]);
		DMiz.set(1, 1, -trig_sin_vectors[3 * index + 2]);
		DMiz.set(0, 1, -trig_cos_vectors[3 * index + 2]);
		DMiz.set(1, 0, trig_cos_vectors[3 * index + 2]);
	}
}
/*----------------------------------------------------------------------------*/
void SmoothingHLBFGS::rebuildQuaternions(double*& eulerAngles)
{

	for (unsigned int i = 0; i < m_candidates.size(); i++)
	{


		Node n = m_candidates[i];
		TCoord angleX = eulerAngles[3 * i];
		TCoord angleY = eulerAngles[3 * i + 1];
		TCoord angleZ = eulerAngles[3 * i + 2];

		math::Quaternion q;
		if (m_mesh->isMarked(n,m_mark_curve))
		{
			//sur une courbe
			q.setFromEulerAngle(angleX, angleY, angleZ);
		}
		else if (m_mesh->isMarked(n, m_mark_surface))
		{
			//sur une surface
			math::Chart t = m_boundary_triad[n.getID()];

			math::Vector newVY = cos(angleZ)*t.Y();
			math::Vector newVZ = sin(angleZ)*t.Z();
			math::Vector tmp1 = newVY+newVZ;

			math::Vector tmp2 = t.X().cross(tmp1);

			math::Quaternion q = math::Quaternion(math::Chart(t.X(), tmp1, tmp2));
			//TODO : change si on repasse a Z = normale
		}
		else {
			//dans le volume
			q.setFromEulerAngle(angleX, angleY, angleZ);
		}

		//mise a jour du Quaternion
		if (m_mesh->isMarked(n,m_mark_surface))
		{
			math::Vector normal = (*m_surf_normal)[n.getID()];
			q=q.alignWith(normal);
		}
	//	std::cout << n.getID() << ": "<<(*m_frame_field)[n.getID()] << " -> " << q << std::endl;
		(*m_frame_field)[n.getID()] = q;

	}
}

/*---------------------------------------------------------------------------*/
gmds::IGMesh* SmoothingHLBFGS::m_mesh = 0;
int SmoothingHLBFGS::m_mark_point = 0;
int SmoothingHLBFGS::m_mark_curve = 0;
int SmoothingHLBFGS::m_mark_surface = 0;
std::vector<gmds::Node> SmoothingHLBFGS::m_candidates;
std::vector<gmds::Edge> SmoothingHLBFGS::m_edges;
std::map<gmds::TCellID, int> SmoothingHLBFGS::m_candidates_index;
std::map<gmds::TCellID, gmds::math::Chart> SmoothingHLBFGS::m_boundary_triad;
/*---------------------------------------------------------------------------*/
SmoothingHLBFGS::SmoothingHLBFGS(
	gmds::IGMesh* AMesh,
	gmds::Variable<gmds::math::Quaternion>*& AField,
	gmds::Variable<gmds::math::Vector>*& ANormal)
	: m_frame_field(AField), m_surf_normal(ANormal),
	m_mark_candidates(0), m_candidate_mark_initialized(false),
	m_boundary_marks_initialized(false)
{
	m_mesh = AMesh;
}
/*---------------------------------------------------------------------------*/
void SmoothingHLBFGS::selectNodes(const int AMark)
{
	m_mark_candidates = AMark;
	m_candidate_mark_initialized = true;
}
/*---------------------------------------------------------------------------*/
void SmoothingHLBFGS::initBoundaryMarks(
	const int AMarkPnt, const int AMarkCurve,
	const int AMarkSurf)
{
	m_mark_point = AMarkPnt;
	m_mark_curve = AMarkCurve;
	m_mark_surface = AMarkSurf;
	m_boundary_marks_initialized = true;
}
/*---------------------------------------------------------------------------*/
void SmoothingHLBFGS::execute()
{
	if (!m_candidate_mark_initialized)
		throw GMDSException("Candidate mark is not initialized");
	if (!m_boundary_marks_initialized)
		throw GMDSException("Boundary marks are not initialized");

	finalSmoothing();
}

/*---------------------------------------------------------------------------*/
void SmoothingHLBFGS::initCandidates()
{
	m_candidates.clear();
	IGMesh::node_iterator it_nodes = m_mesh->nodes_begin();
	for (; !it_nodes.isDone(); it_nodes.next())
	{
		Node n = it_nodes.value();
		if (m_mesh->isMarked(n,m_mark_candidates))
			m_candidates.push_back(it_nodes.value());
	}
	std::cout << std::endl << "Nb candidates (" << m_candidates.size()
		<< " / " << m_mesh->getNbNodes() <<")"<< std::endl;;
}
/*---------------------------------------------------------------------------*/
int SmoothingHLBFGS::initOptimizationData(double*& eulerAngles)
{
	int nb_nodes_on_curves = 0;
	int nb_nodes_on_surfaces = 0;
	int nb_nodes_in_volumes = 0;

	// on calcule des angles d'Euler en tout sommet défini (ridge y compris)
	for (unsigned int i = 0; i < m_candidates.size(); i++)
	{
		Node n = m_candidates[i];
		TCellID n_id = n.getID();
		math::Quaternion q = (*m_frame_field)[n_id];
		//we store the candidate index of n_i
		m_candidates_index[n.getID()] = i;

		if (m_mesh->isMarked(n, m_mark_curve))
		{ //node on an edge??
			q.toEulerAngle(eulerAngles[3 * i],
				eulerAngles[3 * i + 1],
				eulerAngles[3 * i + 2]);
			nb_nodes_on_curves++;
		}
		else if (m_mesh->isMarked(n, m_mark_surface))
		{

			nb_nodes_on_surfaces++;
			math::Chart t(q);
			math::Vector X[6];
			X[0] = t.X(); X[1] = t.Y(); X[2] = t.Z();
			X[3] = X[0].opp();
			X[4] = X[1].opp();
			X[5] = X[2].opp();

			math::Vector normal = (*m_surf_normal)[n_id];

			// We get among X, the vector that is the most aligned with
			// the vertex normal
			int indexNormal = 0;
			TCoord maxVal = normal.dot(X[0]);
			for (unsigned int j = 1; j<6; j++){
			  TCoord currentVal = normal.dot(X[j]);
				if (currentVal>maxVal){
					indexNormal = j;
					maxVal = currentVal;
				}
			}

			// On construit la Chart externe
			math::Vector ChartRes[3];
			ChartRes[0] = X[indexNormal];
			if (indexNormal == 0){
				ChartRes[1] = X[1]; ChartRes[2] = X[2];
			}
			else if (indexNormal == 1){
				ChartRes[1] = X[2]; ChartRes[2] = X[0];
			}
			else if (indexNormal == 2){
				ChartRes[1] = X[0]; ChartRes[2] = X[1];
			}
			else if (indexNormal == 3){
				ChartRes[1] = X[2]; ChartRes[2] = X[1];
			}
			else if (indexNormal == 4){
				ChartRes[1] = X[0]; ChartRes[2] = X[2];
			}
			else if (indexNormal == 5){
				ChartRes[1] = X[0]; ChartRes[2] = X[4];
			}

			m_boundary_triad[n_id] =
				math::Chart(ChartRes[0], ChartRes[1], ChartRes[2]);
			//TODO : change si on repasse a Z = normale
			//			std::cout<<"au point "<<v->getCoord()[0]<<" "<<v->getCoord()[1]<<" "<<v->getCoord()[2]<<" on a une Chart qui vaut "<<std::endl;
			//			std::cout<<"          "<<ChartRes[0]<<std::endl;
			//			std::cout<<"          "<<ChartRes[1]<<std::endl;
			//			std::cout<<"          "<<ChartRes[2]<<std::endl;
			//			boundary_triad[v] = Chart(ChartRes[1],ChartRes[2],ChartRes[0]);
			
			math::Quaternion qTmp(m_boundary_triad[n_id]);

			qTmp.toEulerAngle(eulerAngles[3 * i],
				eulerAngles[3 * i + 1],
				eulerAngles[3 * i + 2]);
			// TOUT A ZERO SI VOLUMIQUE????
	/*		eulerAngles[3 * i] = 0;
			eulerAngles[3 * i + 1] = 0;
			eulerAngles[3 * i + 2] = 0;
*/
		}
		else {//cas general dans le volume
			q.toEulerAngle(eulerAngles[3 * i],
				eulerAngles[3 * i + 1],
				eulerAngles[3 * i + 2]);

			nb_nodes_in_volumes++;
		}
	}
	 
	//Definition of the support Edges
	m_edges.clear();
	for (unsigned int i = 0; i < m_candidates.size(); i++){
		Node n = m_candidates[i];
		std::vector<Edge> adj_edges = n.get<Edge>();

		for (unsigned int j = 0; j < adj_edges.size(); j++)
		{
			Edge e_j = adj_edges[j];

			std::vector<Node> e_j_nodes = e_j.get<Node>();
			Node other_node =
				(e_j_nodes[0].getID() == n.getID()) ? e_j_nodes[1] : e_j_nodes[0];

			if (m_mesh->isMarked(other_node, m_mark_candidates))
			{
				if (other_node.getID() < n.getID()) //to store the edge only once
					m_edges.push_back(e_j);
			}
		}
	}

	return nb_nodes_in_volumes;
}
/*----------------------------------------------------------------------------*/
void SmoothingHLBFGS::
newiteration(int iter, int call_iter, double *x, double* f, double *g, double* gnorm)
{
	std::cout << iter << ": " << call_iter << " " << *f << " " << *gnorm << std::endl;
}
/*---------------------------------------------------------------------------*/
void SmoothingHLBFGS::finalSmoothing()
{
	//We collect the nodes we want to smooth the Quaternion on
	initCandidates();


	double* euler_angles = new double[3 * m_candidates.size()];

	int nb_candidates = initOptimizationData(euler_angles); //new de euler_angle dans la fonction init


	std::cout << "Global Euler Smoothing with HLBFGS"
		<< " (nb candidates= " << m_candidates.size() <<")"<< std::endl;

	//Init. de HLBFGS
	double parameter[20];
	int info[20];
	int N = 3 * m_candidates.size(); // NOMBRE DE VARIABLES: 3 ANGLES PAR SOMMET
	int M = 7;//CHOIX DE L'ALGO D'OPTIMISATION???
	int T = 0;
	int num_iter = 1000;
	bool with_hessian = false;
	INIT_HLBFGS(parameter, info);
	// CHOIX DES PARAMETRES A VERIFIER ??? Ici valeur des exemples du tutoriel de HLBGS
	info[0] = 20;
	info[4] = num_iter;
	info[6] = T;
	info[7] = with_hessian ? 1 : 0;
	info[10] = 0;
	info[11] = 1;

	//std::cout << "AV" << std::endl;
	//for (unsigned int ii = 0; ii<N; ii++){
	//	if ((euler_angles[ii] / 3.14159265358979323846<10e-10) && (euler_angles[ii] / 3.14159265358979323846> (-10e-10)))
	//		std::cout << 0 << " ";
	//	else
	//		std::cout << euler_angles[ii] / 3.14159265358979323846 << " ";
	//}
	//std::cout << std::endl;
	std::cout << "=== Avant le lissage" << std::endl;
	HLBFGS(N, M, &euler_angles[0], (this->evalF), 
		0, HLBFGS_UPDATE_Hessian, 
		this->newiteration, parameter, info);
	std::cout << "=== Apres le lissage, nb iter : " << info[2] 
		<< " pour un nb deval de " << info[1] << std::endl;
	/*

	std::cout << "AP" << std::endl;
	for (unsigned int ii = 0; ii<N; ii++){
		if ((euler_angles[ii] / 3.14159265358979323846<10e-10) && (euler_angles[ii] / 3.14159265358979323846> (-10e-10)))
			std::cout << 0 << " ";
		else
			std::cout << euler_angles[ii] / 3.14159265358979323846 << " ";
	}*/
	std::cout << "Quaternions rebuild" << std::endl;
	rebuildQuaternions(euler_angles);
	std::cout << "Quaternions rebuild done" << std::endl;
	

	delete[] euler_angles;
    
 
}

