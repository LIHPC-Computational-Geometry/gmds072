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
/** \file    MeditReader.t.h
 *  \author  F. LEDOUX
 *  \date    09/11/2008
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_MEDIT_READER_H_
#define GMDS_MEDIT_READER_H_
/*----------------------------------------------------------------------------*/
// headers of GMDS files
#include <GMDS/IG/IG.h>
/*----------------------------------------------------------------------------*/
// headers of STL files
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
#include <GMDS/IO/MeditReader_def.h>
/*----------------------------------------------------------------------------*/
template<typename TMesh>
MeditReader<TMesh>::MeditReader(TMesh& AMesh)
:m_mesh(AMesh)
{}
/*----------------------------------------------------------------------------*/
template<typename TMesh>MeditReader<TMesh>::~MeditReader()
{}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void MeditReader<TMesh>::read(const std::string& AFileName, int AMask,
							  bool useLabels)
{
	std::ifstream input(AFileName.c_str(),std::ios::in);
    if(!input){
        std::string mess ="Impossible to read this Medit file: "+AFileName;
        throw GMDSException(mess);
    }
	const int max_length = 60;
	char ch[max_length];
	int lg;
	do{
		input.getline(ch,max_length);
		lg = input.gcount();
	} while(std::string(ch).find("Dimension")>std::string(ch).size()
		&& lg>0);

	if(lg==0)
		throw GMDSException("Impossible to read the dimension in the Medit file");


	input>>m_mesh_dimension;
	do{
		input.getline(ch,max_length);
		lg = input.gcount();
	} while(std::string(ch).find("Vertices")>std::string(ch).size() && lg>0);

	readNodes(input,useLabels);


	do{

		do{
			input.getline(ch,max_length);
			lg = input.gcount();
		} while((std::string(ch).find("Triangles")>std::string(ch).size() ) &&
				(std::string(ch).find("Edges")>std::string(ch).size() ) &&
				(std::string(ch).find("Quadrilaterals")>std::string(ch).size() ) &&
				(std::string(ch).find("Tetrahedra")>std::string(ch).size() ) &&
				(std::string(ch).find("Hexahedra")>std::string(ch).size() ) &&
				 lg>0);

		if ((std::string(ch).find("Edges")!=std::string::npos &&
			 std::string(ch).find("Edges")<=std::string(ch).size() ) &&
			 m_mesh.getModel().has(E) && (AMask|E)==AMask)
				readEdges(input,useLabels);
		else if ((std::string(ch).find("Triangles")!=std::string::npos &&
				 std::string(ch).find("Triangles")<=std::string(ch).size() ) &&
				 (AMask|F)==AMask &&  m_mesh.getModel().has(F) )
				readTriangles(input,useLabels);
		else if ((std::string(ch).find("Quadrilaterals")!=std::string::npos &&
				 std::string(ch).find("Quadrilaterals")<=std::string(ch).size() )
				 && (AMask|F)==AMask &&  m_mesh.getModel().has(F) )
				readQuadrilaterals(input,useLabels);
		else if ((std::string(ch).find("Tetrahedra")!=std::string::npos &&
				 std::string(ch).find("Tetrahedra")<=std::string(ch).size() )
				 && (AMask|R)==AMask &&  m_mesh.getModel().has(R) )
				readTetrahedra(input,useLabels);
		else if ((std::string(ch).find("Hexahedra")!=std::string::npos &&
				 std::string(ch).find("Hexahedra")<=std::string(ch).size() ) &&
				 (AMask|R)==AMask &&  m_mesh.getModel().has(R) )
				readHexahedra(input,useLabels);
	}
	while(lg>0);

	input.close();

}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void MeditReader<TMesh>::readNodes(std::ifstream& AIn, bool useLabels)
{
	TInt nb_nodes;
	AIn>>nb_nodes;
	std::map<int,std::string> labels;

	for(unsigned int i = 0; i < nb_nodes; i++)
	{
		  TCoord x, y, z;
		  int ref;
		  Node n;
		  if(m_mesh_dimension==2)
		  {
			  AIn>>x>>y>>ref;
			  n = m_mesh.newNode(x,y);
		  }
		  else
		  {
			  AIn>>x>>y>>z>>ref;
			  n = m_mesh.newNode(x,y,z);
		  }

		  if(useLabels)
		  {
//			  if(labels.find(ref)!=labels.end()){
//				  m_mesh.getCloud(labels[ref]).add(n);
//			  }
//			  else{
//				  std::ostringstream stream;
//				  stream << ref;
//				  typename TMesh::cloud& cl = m_mesh.newCloud(stream.str());
//				  cl.add(n);
//				  labels.insert(std::pair<int,std::string>(ref,stream.str()));
//			  }
		  }
	}
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void MeditReader<TMesh>::readEdges(std::ifstream& AIn, bool useLabels)
{
	TInt nb_edges;
	AIn>>nb_edges;
	std::map<int,std::string> labels;
	for(unsigned int i = 0; i < nb_edges; i++)
	{
		  TCellID x1, x2;
		  int ref;

		  AIn>>x1>>x2>>ref;
		  Node n1 = m_mesh.template get<Node>(x1-1);
		  Node n2 = m_mesh.template get<Node>(x2-1);
		  Edge e = m_mesh.newEdge(n1,n2);

		  if(useLabels)
		  {
//			  if(labels.find(ref)!=labels.end()){
//				  m_mesh.getline(labels[ref]).add(e);
//			  }
//			  else{
//				  std::ostringstream stream;
//				  stream << ref;
//				  typename Mesh<TMesh>::line& li = m_mesh.newLine(stream.str());
//				  li.add(e);
//				  labels.insert(std::pair<int,std::string>(ref,stream.str()));
//			  }
		  }

	}
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void MeditReader<TMesh>::storeFaceLabel(Face &f, const int& lbl)
{
	std::string s;
	if(m_2D_labels.find(lbl)!=m_2D_labels.end())
		s = m_2D_labels[lbl];
	else{
		std::ostringstream stream;
		stream << lbl;
		s =stream.str();
		m_mesh.newSurface(s);
		m_2D_labels.insert(std::pair<int,std::string>(lbl,s));
	}
	m_mesh.getSurface(s).add(f);
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void MeditReader<TMesh>::storeRegionLabel(Region &r, const int& lbl)
{
	std::string v;
	if(m_3D_labels.find(lbl)!=m_3D_labels.end())
		v = m_3D_labels[lbl];
	else{
		std::ostringstream stream;
		stream << lbl;
		v = stream.str();
		m_mesh.newVolume(v);
		m_3D_labels.insert(std::pair<int,std::string>(lbl,v));
	}
	m_mesh.getVolume(v).add(r);

}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void MeditReader<TMesh>::readTriangles(std::ifstream& AIn, bool useLabels)
{
	TInt nb_triangles;
	AIn>>nb_triangles;
	for(unsigned int i = 0; i < nb_triangles; i++)
	{
		  TCellID x1, x2, x3;
		  int ref;

		  AIn>>x1>>x2>>x3>>ref;
		  Node n1 = m_mesh.template get<Node>(x1-1);
		  Node n2 = m_mesh.template get<Node>(x2-1);
		  Node n3 = m_mesh.template get<Node>(x3-1);
		  Face f = m_mesh.newTriangle(n1,n2,n3);
		  if(useLabels)
			  storeFaceLabel(f,ref);
	}
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void MeditReader<TMesh>::readQuadrilaterals(std::ifstream& AIn, bool useLabels)
{
	TInt nb_quads;
	AIn>>nb_quads;
	for(unsigned int i = 0; i < nb_quads; i++)
	{
		  TCellID x1, x2, x3, x4;
		  int ref;

		  AIn>>x1>>x2>>x3>>x4>>ref;
		  Node n1 = m_mesh.template get<Node>(x1-1);
		  Node n2 = m_mesh.template get<Node>(x2-1);
		  Node n3 = m_mesh.template get<Node>(x3-1);
		  Node n4 = m_mesh.template get<Node>(x4-1);
		  Face f  = m_mesh.newQuad(n1,n2,n3,n4);
		  if(useLabels)
			  storeFaceLabel(f,ref);
	}
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void MeditReader<TMesh>::readTetrahedra(std::ifstream& AIn, bool useLabels)
{
	TInt nb_tets;
	AIn>>nb_tets;
	for(unsigned int i = 0; i < nb_tets; i++)
	{
		  TCellID x1, x2, x3, x4;
		  int ref;

		  AIn>>x1>>x2>>x3>>x4>>ref;
		  Node n1 = m_mesh.template get<Node>(x1-1);
		  Node n2 = m_mesh.template get<Node>(x2-1);
		  Node n3 = m_mesh.template get<Node>(x3-1);
		  Node n4 = m_mesh.template get<Node>(x4-1);
		  Region r=m_mesh.newTet(n1,n2,n3,n4);
		  if(useLabels)
			  storeRegionLabel(r,ref);
	}
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void MeditReader<TMesh>::readHexahedra(std::ifstream& AIn, bool useLabels)
{
	TInt nb_hex;
	AIn>>nb_hex;
	for(unsigned int i = 0; i < nb_hex; i++)
	{
		  TCellID x1, x2, x3, x4, x5, x6, x7, x8;
		  int ref;

		  AIn>>x1>>x2>>x3>>x4>>x5>>x6>>x7>>x8>>ref;
		  Node n1 = m_mesh.template get<Node>(x1-1);
		  Node n2 = m_mesh.template get<Node>(x2-1);
		  Node n3 = m_mesh.template get<Node>(x3-1);
		  Node n4 = m_mesh.template get<Node>(x4-1);
		  Node n5 = m_mesh.template get<Node>(x5-1);
		  Node n6 = m_mesh.template get<Node>(x6-1);
		  Node n7 = m_mesh.template get<Node>(x7-1);
		  Node n8 = m_mesh.template get<Node>(x8-1);
		  Region r = m_mesh.newHex(n1,n2,n3,n4,n5,n6,n7,n8);
		  if(useLabels)
			  storeRegionLabel(r,ref);
	}
}
/*----------------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------------*/
#endif /* GMDS_MEDIT_READER_H_ */
/*----------------------------------------------------------------------------*/
