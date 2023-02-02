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
/** \file    MatrixMarketWriter.h
 *  \author  F. LEDOUX
 *  \date    16 may 2014
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MATRIXMARKETWRITER_H_
#define GMDS_MATRIXMARKETWRITER_H_
/*----------------------------------------------------------------------------*/
// headers of STL files
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
#include <GMDS/IO/MatrixMarketWriter_def.h>
/*----------------------------------------------------------------------------*/
template<typename TMesh>
MatrixMarketWriter<TMesh>::MatrixMarketWriter(TMesh& AMesh)
:m_mesh(AMesh)
{
	m_mesh_dimension = m_mesh.getDim();
	MeshModel model = m_mesh.getModel();
	if(model.has(R))
		m_max_cell_dim = 3;
	else if(model.has(F))
		m_max_cell_dim = 2;
	else
		throw GMDSException("Only meshes with face or regions can be exported");

	if(m_max_cell_dim==3 && !model.has(R2R))
		throw GMDSException("Only volume meshes with R2R can be exported");

	if(m_max_cell_dim==2 && !model.has(F2F))
		throw GMDSException("Only surface meshes with F2F can be exported");
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
MatrixMarketWriter<TMesh>::~MatrixMarketWriter()
{}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void MatrixMarketWriter<TMesh>::write(
		const std::string& AFileName,
		const bool& AWithCoords)
{
	std::stringstream fileName;
	fileName<<AFileName.c_str()<<".mtx";
	std::ofstream output(fileName.str().c_str(),std::ios::out);
	if(!output)
		throw GMDSException("Impossible to create a Matrix Market File");

	output<<"%%MatrixMarket matrix coordinate real general\n";
	output<<"%-------------------------------------------------------------------------------\n";
	output<<"% Generated from the GMDS library\n";
	output<<"%-------------------------------------------------------------------------------\n";

	/// The next line contains the number of lines and columns and the number of non-zero values.
	TInt nb_lines_columns = 0;
	TInt nb_entries = 0;
	if(m_max_cell_dim==3){
		nb_lines_columns = m_mesh.getNbRegions();
		nb_entries = build3DGraph();
		output<<nb_lines_columns<<" "<<nb_lines_columns<<" "<<nb_entries<<"\n";
		drop3DGraph(output);
	}
	else {
		nb_lines_columns = m_mesh.getNbFaces();
		nb_entries = build2DGraph();
		output<<nb_lines_columns<<" "<<nb_lines_columns<<" "<<nb_entries<<"\n";
		drop2DGraph(output);
	}

	if(AWithCoords){
		std::stringstream fileNameCoord;
		fileNameCoord<<AFileName.c_str()<<".mtx_coord";
		std::ofstream output_coord(fileNameCoord.str().c_str(),std::ios::out);
		output_coord<<"%%MatrixMarket matrix coordinate real general\n";
		output_coord<<"%-------------------------------------------------------------------------------\n";
		output_coord<<"% Generated from the GMDS library\n";
		output_coord<<"% Storage of graph coordinates\n";
		output_coord<<"%-------------------------------------------------------------------------------\n";
		output_coord<<nb_lines_columns<<" 3 "<<nb_lines_columns<<"\n";
		if(m_max_cell_dim==3)
			write3DCoord(output_coord);
		else
			write2DCoord(output_coord);

		output_coord.close();
	}

	output.close();
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
TInt MatrixMarketWriter<TMesh>::build2DGraph()
{
	TInt nb_entries =0;
	typename TMesh::face_iterator it = m_mesh.faces_begin();
	while(!it.isDone()){

		Face f = it.value();
		TCellID f_id = f.getID();
		std::vector<TCellID> adj_faces_ids = f.getIDs<Face>();
		for(unsigned int i=0;i<adj_faces_ids.size();i++){
			m_graph[f_id].push_back(adj_faces_ids[i]);
			nb_entries++;
		}
		it.next();
	}
	return nb_entries;
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
TInt MatrixMarketWriter<TMesh>::build3DGraph()
{
	TInt nb_entries =0;
	typename TMesh::region_iterator it = m_mesh.regions_begin();
	while(!it.isDone()){
		Region r = it.value();
		TCellID r_id = r.getID();
		std::vector<TCellID> adj_regions_ids = r.getIDs<Region>();
		for(unsigned int i=0;i<adj_regions_ids.size();i++){
			m_graph[r_id].push_back(adj_regions_ids[i]);
			nb_entries++;
		}
		it.next();
	}
	return nb_entries;
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void MatrixMarketWriter<TMesh>::drop2DGraph(std::ofstream& AOut)
{
	typename TMesh::face_iterator it = m_mesh.faces_begin();
	while(!it.isDone()){
		Face f = it.value();
		TCellID f_id = f.getID();
		std::vector<TCellID> adj_ids = m_graph[f_id];
		for(unsigned int i=0;i<adj_ids.size();i++)
			AOut<<f_id<<" "<<adj_ids[i]<<" 1.0\n";
		it.next();
	}
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void MatrixMarketWriter<TMesh>::drop3DGraph(std::ofstream& AOut)
{
	typename TMesh::region_iterator it = m_mesh.regions_begin();
	while(!it.isDone()){
		Region r = it.value();
		TCellID r_id = r.getID();
		std::vector<TCellID> adj_ids = m_graph[r_id];
		for(unsigned int i=0;i<adj_ids.size();i++)
			AOut<<r_id<<" "<<adj_ids[i]<<" 1.0\n";
		it.next();
	}
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void MatrixMarketWriter<TMesh>::write2DCoord(std::ofstream& AOut)
{
	typename TMesh::face_iterator it = m_mesh.faces_begin();
	while(!it.isDone()){
		Face f = it.value();
		std::vector<Node> nodes = f.get<Node>();
		TCoord x=0, y=0, z=0;
		for(unsigned int i=0;i<nodes.size();i++){
			Node current = nodes[i];
			x+=current.X();
			y+=current.Y();
			z+=current.Z();
		}
		x=x/nodes.size();
		y=y/nodes.size();
		z=z/nodes.size();
		AOut<<x<<" "<<y<<" "<<z<<"\n";
		it.next();
	}
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void MatrixMarketWriter<TMesh>::write3DCoord(std::ofstream& AOut)
{
	typename TMesh::region_iterator it = m_mesh.regions_begin();
	while(!it.isDone()){
		Region r = it.value();
		std::vector<Node> nodes = r.get<Node>();
		TCoord x=0, y=0, z=0;
		for(unsigned int i=0;i<nodes.size();i++){
			Node current = nodes[i];
			x+=current.X();
			y+=current.Y();
			z+=current.Z();
		}
		x=x/nodes.size();
		y=y/nodes.size();
		z=z/nodes.size();
		AOut<<x<<" "<<y<<" "<<z<<"\n";
		it.next();
	}
}
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MATRIXMARKETWRITER_H_ */


