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
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IG.h>
#include <GMDS/IO/MatrixMarketWriter.h>
#include <GMDS/Math/Vector.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
class WriterTest: public ::testing::Test {

  protected:
	WriterTest(){;}
    virtual ~WriterTest(){;}
};
/*----------------------------------------------------------------------------*/
TEST_F(WriterTest ,matrixMarketWrite2Triangles3D)
{
	MeshModel mod = DIM3|N|F|F2N|F2F;
	IGMesh mesh(mod);

	Node t[4];
	Face f[2];

	t[0] = mesh.newNode(0.,0.,0.);
	t[1] = mesh.newNode(0.,1.,0.);
	t[2] = mesh.newNode(1.,1.,0.);
	t[3] = mesh.newNode(1.,0.,0.);

	mesh.newTriangle(t[0],t[1],t[2]);
	mesh.newTriangle(t[0],t[2],t[3]);

	IGMeshDoctor doc(&mesh);
	doc.buildF2F(mod);
	MatrixMarketWriter<IGMesh> writer(mesh);
	writer.write("matrixMarked1",true);
}

/*----------------------------------------------------------------------------*/
