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
/** \file    CubitFacReader.t.h
 *  \author  legoff
 *  \date    16/10/2014
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_CUBITFACREADER_H_
#define GMDS_CUBITFACREADER_H_
/*----------------------------------------------------------------------------*/
// headers of GMDS files
#include <GMDS/Utils/CommonTypes.h>
#include <GMDS/IG/Node.h>
#include <GMDS/IG/Face.h>
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
#include <GMDS/IO/CubitFacReader_def.h>
/*----------------------------------------------------------------------------*/
template<typename TMesh>
CubitFacReader<TMesh>::CubitFacReader(TMesh& AMesh)
:m_mesh(AMesh)
{}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
CubitFacReader<TMesh>::~CubitFacReader()
{}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void CubitFacReader<TMesh>::read(const std::string& AFileName)
{
	std::ifstream input(AFileName.c_str(),std::ios::in);
	if(!input)
		throw GMDSException("Impossible to read this FAC file");

	readNodes(input);
	readTriangles(input);

	input.close();

}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void CubitFacReader<TMesh>::readNodes(std::ifstream& AIn)
{
	TInt nb_nodes;
	AIn>>nb_nodes;
	for(unsigned int iNode = 0; iNode < nb_nodes; iNode++)
	{
		  TCoord x, y, z;
		  int ref;
		  AIn>>ref>>x>>y>>z;
		  Node n = this->m_mesh.newNode(x,y,z);
		  m_fac2GMDSNodes[ref] = n;
	}
}
/*----------------------------------------------------------------------------*/
template<typename TMesh>
void CubitFacReader<TMesh>::readTriangles(std::ifstream& AIn)
{
	TInt nb_triangles;
	AIn>>nb_triangles;
	for(unsigned int i = 0; i < nb_triangles; i++)
	{
		  TCellID x1, x2, x3;
		  int ref;

		  AIn>>ref>>x1>>x2>>x3;
		  Face f = m_mesh.newTriangle(m_fac2GMDSNodes[x1],
					      m_fac2GMDSNodes[x2],
					      m_fac2GMDSNodes[x3]);
	}
}
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
#endif /* GMDS_CUBITFACREADER_H_ */
/*----------------------------------------------------------------------------*/
