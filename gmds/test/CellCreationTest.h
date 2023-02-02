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
#include <GMDS/IG/IG.h>
/*----------------------------------------------------------------------------*/
class CellCreationTest: public ::testing::Test {
protected:
	CellCreationTest(){;}
	virtual ~CellCreationTest(){;}
};

/*----------------------------------------------------------------------------*/
TEST_F(CellCreationTest,meshNewNode) {
	gmds::MeshModel mod = gmds::F|gmds::N;
	gmds::IGMesh m(mod);
	gmds::Node n = m.newNode(0,0,0);

	EXPECT_EQ(0.0, n.X());
	EXPECT_EQ(0.0, n.Y());
	EXPECT_EQ(0.0, n.Z());
	EXPECT_EQ(0, n.getID());
}
/*----------------------------------------------------------------------------*/
TEST_F(CellCreationTest,meshNewEdge) {
	gmds::MeshModel mod = gmds::F|gmds::N|gmds::E2N;
	gmds::IGMesh m(mod);
	gmds::Node n1 = m.newNode(0,0,0);
	gmds::Node n2 = m.newNode(1,0,0);
	gmds::Edge e = m.newEdge(n1,n2);

	EXPECT_EQ(0, e.getID());

	std::vector<gmds::Node> nodes = e.get<gmds::Node>();
	EXPECT_EQ(n1.getID(),nodes[0].getID());
	EXPECT_EQ(n2.getID(),nodes[1].getID());
}
/*----------------------------------------------------------------------------*/
TEST_F(CellCreationTest,hascell) {
	gmds::MeshModel mod = gmds::F|gmds::N|gmds::F2N;
	gmds::IGMesh m(mod);
	gmds::Node n1 = m.newNode(0,0,0);
	gmds::Node n2 = m.newNode(1,0,0);
	gmds::Node n3 = m.newNode(1,0,0);
	gmds::Face f  = m.newTriangle(n1,n2,n3);

	EXPECT_TRUE(m.has<Node>(n1.getID()));
	EXPECT_TRUE(m.has<Node>(n2.getID()));
	EXPECT_TRUE(m.has<Node>(n3.getID()));
	EXPECT_TRUE(m.has<Face>(f.getID()));

	m.deleteNode(n2);
	m.deleteFace(f);


	EXPECT_TRUE(m.has<Node>(n1.getID()));
	EXPECT_FALSE(m.has<Node>(n2.getID()));
	EXPECT_TRUE(m.has<Node>(n3.getID()));
	EXPECT_FALSE(m.has<Face>(f.getID()));
}
/*----------------------------------------------------------------------------*/
TEST_F(CellCreationTest,meshNewDelNode) {
	gmds::MeshModel mod = gmds::F|gmds::N;
	gmds::IGMesh m(mod);
	gmds::Node n = m.newNode(0,0,0);
	m.deleteNode(n);
	EXPECT_EQ(0, m.getNbNodes());
}
/*----------------------------------------------------------------------------*/
TEST_F(CellCreationTest,meshNewDelNode2) {
	gmds::MeshModel mod = gmds::F|gmds::N;
	gmds::IGMesh m(mod);
	gmds::Node n = m.newNode(0,0,0);
	gmds::Node n2 = m.newNode(0,0,0);
	m.deleteNode(n);
	EXPECT_EQ(1, m.getNbNodes());
	gmds::Node n3 = m.newNode(0,0,0);
	EXPECT_EQ(2, m.getNbNodes());
}
/*----------------------------------------------------------------------------*/
TEST_F(CellCreationTest,meshNewDelManyNodes) {
	gmds::MeshModel mod = gmds::F|gmds::N;
	gmds::IGMesh m(mod);
	const int nb = 10000;
	gmds::Node n[nb];

	for(unsigned int i=0;i<nb;i++){
		n[i] = m.newNode(0,0,0);
	}

	for(unsigned int i=0;i<nb;i++){
		m.deleteNode(n[i]);
		EXPECT_EQ(nb-i-1, m.getNbNodes());
	}

}
/*----------------------------------------------------------------------------*/
TEST_F(CellCreationTest,meshNewTriangle) {
	gmds::MeshModel mod = gmds::F|gmds::N| gmds::F2N| gmds::F2F;
	gmds::IGMesh m(mod);
	gmds::Node init_nodes[3];
	init_nodes[0] = m.newNode(0,0,0);
	init_nodes[1] = m.newNode(1,0,0);
	init_nodes[2] = m.newNode(0,1,0);
	gmds::Face f = m.newTriangle(init_nodes[0],init_nodes[1],init_nodes[2]);

	EXPECT_EQ(0, f.getID());
	EXPECT_EQ(3, f.getNbNodes());
	EXPECT_EQ(3, f.getNbFaces());

	std::vector<gmds::Node> nodes;
	f.get<gmds::Node>(nodes);
	for(unsigned int in=0;in<nodes.size();in++){
		gmds::Node ni = nodes[in];
		EXPECT_EQ(init_nodes[in].getID(),ni.getID());
	}

}
/*----------------------------------------------------------------------------*/
TEST_F(CellCreationTest,meshNewTet) {
	gmds::MeshModel mod = gmds::R|gmds::N| gmds::R2N;
	gmds::IGMesh m(mod);
	gmds::Node init_nodes[4];
	init_nodes[0] = m.newNode(0,0,0);
	init_nodes[1] = m.newNode(1,0,0);
	init_nodes[2] = m.newNode(0,1,0);
	init_nodes[3] = m.newNode(0,0,1);
	gmds::Region r = m.newTet(init_nodes[0],init_nodes[1],init_nodes[2],init_nodes[3]);

	EXPECT_EQ(0, r.getID());
	EXPECT_EQ(4, r.getNbNodes());

	std::vector<gmds::Node> nodes;
	r.get<gmds::Node>(nodes);
	for(unsigned int in=0;in<nodes.size();in++){
		gmds::Node ni = nodes[in];
		EXPECT_EQ(init_nodes[in].getID(),ni.getID());
	}
}
/*----------------------------------------------------------------------------*/
TEST_F(CellCreationTest,meshNewHex) {
	gmds::MeshModel mod = gmds::R|gmds::N| gmds::R2N;
	gmds::IGMesh m(mod);
	gmds::Node init_nodes[8];
	init_nodes[0] = m.newNode(0,0,0);
	init_nodes[1] = m.newNode(1,0,0);
	init_nodes[2] = m.newNode(0,1,0);
	init_nodes[3] = m.newNode(1,1,0);
	init_nodes[4] = m.newNode(0,0,1);
	init_nodes[5] = m.newNode(1,0,1);
	init_nodes[6] = m.newNode(0,1,1);
	init_nodes[7] = m.newNode(1,1,1);
	gmds::Region r = m.newHex(init_nodes[0],init_nodes[1],init_nodes[2],init_nodes[3],
			init_nodes[4],init_nodes[5],init_nodes[6],init_nodes[7]);

	EXPECT_EQ(0, r.getID());
	EXPECT_EQ(8, r.getNbNodes());

	std::vector<gmds::Node> nodes;
	r.get<gmds::Node>(nodes);
	for(unsigned int in=0;in<nodes.size();in++){
		gmds::Node ni = nodes[in];
		EXPECT_EQ(init_nodes[in].getID(),ni.getID());
	}
}
/*----------------------------------------------------------------------------*/
//TEST_F(CellCreationTest,meshNewPolygon) {
//	gmds::MeshModel mod = gmds::F|gmds::N| gmds::F2N| gmds::F2F;
//	gmds::IGMesh m(mod);
//	std::vector<gmds::Node> init_nodes;
//	init_nodes.resize(3);
//	init_nodes[0] = m.newNode(0,0,0);
//	init_nodes[1] = m.newNode(1,0,0);
//	init_nodes[2] = m.newNode(0,1,0);
//
//	gmds::Face f = m.newPolygon(init_nodes);
//	EXPECT_EQ(0, f.getID());
//
//	EXPECT_EQ(3, f.getNbNodes());
//
//
//	std::vector<gmds::Node> nodes;
//	f.get<gmds::Node>(nodes);
//	for(unsigned int in=0;in<nodes.size();in++){
//		gmds::Node ni = nodes[in];
//		EXPECT_EQ(init_nodes[in].getID(),ni.getID());
//	}
//}
/*----------------------------------------------------------------------------*/
TEST_F(CellCreationTest,mesh2DGrid) {
	gmds::MeshModel mod = gmds::F|gmds::N| gmds::F2N| gmds::N2F;
	gmds::IGMesh mesh(mod);

	const int si = 2;
	const int sj = 2;
	gmds::Node t[si][sj];

	/* nodes creation */
	for(int i=0;i<si;i++)
		for(int j=0;j<sj;j++){
			t[i][j] = mesh.newNode(i,j);
		}

	for(int i=1;i<si;i++)
		for(int j=1;j<sj;j++){
			gmds::Face f = mesh.newQuad(t[i-1][j-1],t[i-1][j],t[i][j],t[i][j-1]);
			t[i-1][j-1].add<gmds::Face>(f);
			t[i-1][j].add<gmds::Face>(f);
			t[i][j-1].add<gmds::Face>(f);
			t[i][j].add<gmds::Face>(f);
		}
	int nbLocalNodes=0;
	gmds::IGMesh::face_iterator it  = mesh.faces_begin();
	for(;!it.isDone();it.next())
	{
		gmds::Face f=it.value();
		std::vector<gmds::Node> local_nodes = f.get<gmds::Node>();
		nbLocalNodes +=local_nodes.size();
	}
	EXPECT_EQ(si*sj, mesh.getNbNodes());
	EXPECT_EQ((si-1)*(sj-1), mesh.getNbFaces());
	EXPECT_EQ(4*(si-1)*(sj-1),nbLocalNodes);

}
/*----------------------------------------------------------------------------*/
TEST_F(CellCreationTest,meshNewGrid3D) {
	gmds::MeshModel mod = gmds::DIM3|gmds::F|gmds::N| gmds::F2N;
	gmds::IGMesh mesh(mod);

	const int s = 10;
	gmds::Node t[s][s][s];

	/* nodes creation */
	for(int i=0;i<s;i++)
		for(int j=0;j<s;j++)
			for(int k=0;k<s;k++)
				t[i][j][k] = mesh.newNode(i,j,k);

	/* faces creation */
	for(int i=1;i<s;i++)
		for(int j=1;j<s;j++)
			for(int k=1;k<s;k++)
			{
				mesh.newQuad(t[i-1][j-1][k],t[i][j-1][k],t[i][j][k],t[i-1][j][k]);
				mesh.newTriangle(t[i-1][j-1][k],t[i-1][j][k],t[i][j-1][k]);
				mesh.newTriangle(t[i][j-1][k],t[i-1][j][k],t[i][j][k]);
				std::vector<gmds::Node> nodes_for_poly;
				nodes_for_poly.push_back(t[i-1][j-1][k]);
				nodes_for_poly.push_back(t[i-1][j][k]);
				nodes_for_poly.push_back(t[i][j][k]);
				nodes_for_poly.push_back(t[i][j-1][k]);
				mesh.newPolygon(nodes_for_poly);
			}

	int nb_faces = 0;
	gmds::TCellID current_id;
	gmds::IGMesh::face_iterator itf = mesh.faces_begin();
	for(; !itf.isDone();itf.next())
	{
		gmds::Face f = itf.value();
		current_id = f.getID();
		EXPECT_TRUE(current_id!=NullID);
		nb_faces++;
	}


	EXPECT_TRUE(nb_faces!=0);

	mesh.deleteFace(1);
	mesh.deleteFace(10);
	mesh.deleteFace(5);
	mesh.deleteFace(601);
	mesh.deleteFace(602);
	mesh.deleteFace(603);
	mesh.deleteFace(16);
	mesh.deleteFace(154);
	mesh.deleteFace(768);
	mesh.deleteFace(521);
	mesh.deleteFace(996);

	gmds::TInt nb_faces_before = nb_faces;
	nb_faces=0;
	for(itf = mesh.faces_begin(); !itf.isDone();itf.next())
	{
		gmds::Face f = itf.value();
		current_id = f.getID();
		EXPECT_TRUE(current_id!=NullID);
		nb_faces++;
	}
	EXPECT_TRUE(nb_faces!=0);
	EXPECT_TRUE(nb_faces==nb_faces_before-11);
}
/*----------------------------------------------------------------------------*/
