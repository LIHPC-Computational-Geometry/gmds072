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
#include <gtest/gtest.h>
#include <GMDS/IG/IG.h>

#include <GMDS/IO/MeditReader.h>
class ModelChangeTest: public ::testing::Test {

  protected:
	ModelChangeTest(){;}
    virtual ~ModelChangeTest(){;}
};
/*----------------------------------------------------------------------------*/
TEST_F(ModelChangeTest,intersectionModel)
{
	gmds::MeshModel m1 = gmds::F|gmds::N|gmds::E;
	gmds::MeshModel m2 = gmds::N|gmds::F|gmds::R;

	gmds::MeshModel m = gmds::MeshModel::intersection(m1,m2);
	EXPECT_TRUE (m.has(gmds::N));
	EXPECT_TRUE (m.has(gmds::F));
	EXPECT_FALSE(m.has(gmds::E));
	EXPECT_FALSE(m.has(gmds::R));
}
/*----------------------------------------------------------------------------*/
TEST_F(ModelChangeTest,incrementModel)
{
	gmds::MeshModel m = gmds::F|gmds::N|gmds::F2N;

	m.add(gmds::E|gmds::E2N);

	EXPECT_TRUE(m.has(gmds::N));
	EXPECT_TRUE(m.has(gmds::F));
	EXPECT_TRUE(m.has(gmds::E));
	EXPECT_TRUE(m.has(gmds::E2N));
}
/*----------------------------------------------------------------------------*/
TEST_F(ModelChangeTest,meshChange2D_1)
{
	gmds::MeshModel mod = gmds::F|gmds::N|gmds::F2N;
	gmds::IGMesh m(mod);

	gmds::Node n0 = m.newNode(0,0,0);
	gmds::Node n1 = m.newNode(1,0,0);
	gmds::Node n2 = m.newNode(1,1,0);
	gmds::Node n3 = m.newNode(0,1,0);
	gmds::Node n4 = m.newNode(2,0,0);
	gmds::Node n5 = m.newNode(2,1,0);

	gmds::Face f1 = m.newQuad(n0,n1,n2,n3);
	gmds::Face f2 = m.newQuad(n1,n4,n5,n2);

	EXPECT_EQ(4, f1.getNbNodes());
	EXPECT_EQ(0, f1.getNbFaces());

	gmds::MeshModel mod2 = gmds::F|gmds::N|gmds::F2N|gmds::F2F;
	m.changeModel(mod2);

	gmds::IGMesh::face_iterator it_f = m.faces_begin();
	for(;!it_f.isDone();it_f.next()){
		gmds::Face f = it_f.value();
		EXPECT_EQ(1, f.get<gmds::Face>().size());
	}
	EXPECT_EQ(1, f1.get<gmds::Face>().size());

}
/*----------------------------------------------------------------------------*/
TEST_F(ModelChangeTest,meshChange2D_2)
{
	gmds::MeshModel mod = gmds::F|gmds::N|gmds::F2N;
	gmds::IGMesh m(mod);

	gmds::Node n0 = m.newNode(0,0,0);
	gmds::Node n1 = m.newNode(1,0,0);
	gmds::Node n2 = m.newNode(1,1,0);
	gmds::Node n3 = m.newNode(0,1,0);


	m.newTriangle(n0,n1,n2);
	m.newTriangle(n3,n1,n2);

	EXPECT_EQ(0, m.getNbEdges());

	gmds::MeshModel mod2 = gmds::F|gmds::N|gmds::F2N|gmds::E|gmds::E2N;
	m.changeModel(mod2);

	EXPECT_EQ(5, m.getNbEdges());

}
/*----------------------------------------------------------------------------*/
TEST_F(ModelChangeTest,meshChange3D_1)
{
	gmds::MeshModel mod = gmds::R|gmds::N|gmds::R2N;
	gmds::IGMesh m(mod);

	gmds::Node n0 = m.newNode(0,0,0);
	gmds::Node n1 = m.newNode(1,0,0);
	gmds::Node n2 = m.newNode(1,1,0);
	gmds::Node n3 = m.newNode(0,1,0);
	gmds::Node n4 = m.newNode(0,1,1);

	gmds::Region r = m.newTet(n0,n1,n2,n3);
	gmds::Region r3 = m.newTet(n4,n1,n2,n3);

	EXPECT_EQ(0, m.getNbFaces());

	gmds::MeshModel mod2 = gmds::R|gmds::N|gmds::R2N|gmds::R2F|gmds::F|gmds::F2N;
	m.changeModel(mod2);

	EXPECT_EQ(7, m.getNbFaces());

}

/*----------------------------------------------------------------------------*/
TEST_F(ModelChangeTest,meshChange2D_3)
{
	gmds::MeshModel mod = DIM3|F|N|F2N;
	gmds::IGMesh mesh(mod);


	MeditReader<IGMesh> reader(mesh);
	reader.read("Samples/carre0.mesh",N|F);


	EXPECT_NE(0, mesh.getNbTriangles());
	EXPECT_EQ(0, mesh.getNbQuadrilaterals());

	mod.add(E|E2N);
	mesh.changeModel(mod);

	EXPECT_NE(0, mesh.getNbEdges());

}
/*----------------------------------------------------------------------------*/
TEST_F(ModelChangeTest,meshChange2D_4)
{
	gmds::MeshModel mod = DIM3|F|N|F2N;
	gmds::IGMesh mesh(mod);


	MeditReader<IGMesh> reader(mesh);
	reader.read("Samples/carre0.mesh",N|F);


	EXPECT_NE(0, mesh.getNbTriangles());
	EXPECT_EQ(0, mesh.getNbQuadrilaterals());

	mod.add(E|E2N|F2E);
	mesh.changeModel(mod);

	EXPECT_NE(0, mesh.getNbEdges());
}
/*----------------------------------------------------------------------------*/
TEST_F(ModelChangeTest,meshChange2D_5)
{
	gmds::MeshModel mod = DIM3|F|N|F2N;
	gmds::IGMesh mesh(mod);


	MeditReader<IGMesh> reader(mesh);
	reader.read("Samples/carre0.mesh",N|F);


	EXPECT_NE(0, mesh.getNbTriangles());
	EXPECT_EQ(0, mesh.getNbQuadrilaterals());

	mod.add(E|E2N|F2E|E2F);
	mesh.changeModel(mod);

	EXPECT_NE(0, mesh.getNbEdges());
	gmds::IGMesh::edge_iterator it = mesh.edges_begin();
	for(;!it.isDone();it.next())
	{
		gmds::Edge e = it.value();
		std::vector<gmds::Face> faces = e.get<gmds::Face>();
		EXPECT_NE(0, faces.size());
	}
}

/*----------------------------------------------------------------------------*/
