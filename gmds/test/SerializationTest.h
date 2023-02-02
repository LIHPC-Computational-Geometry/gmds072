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
/*
 * SerializationTest.cpp
 *
 *  Created on: 3 août 2010
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IG.h>
/*----------------------------------------------------------------------------*/
#include <fstream>

#include <strstream>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
class TestVarSerialization{
public:
	TestVarSerialization(const int x=1, const double y=0.25):m_x(x),m_y(y){;}

	int m_x;
	double m_y;
};
/*----------------------------------------------------------------------------*/
class SerializationTest: public ::testing::Test {

  protected:
	SerializationTest(){;}
    virtual ~SerializationTest(){;}
};
/*----------------------------------------------------------------------------*/
TEST_F(SerializationTest ,testSerializeIntVariable)
{
	std::fstream file;
	file.open("test.dat", std::ios::out|std::ios::binary);
	gmds::Variable<int> var("Pression");
	var.setDomain(10);
	for(int i=0;i<10;i++)
		var[i]=10*i;

	var.serialize(file);
	file.close();
	file.open("test.dat", std::ios::in|std::ios::binary);
	gmds::Variable<int> var2;
	var2.unserialize(file);


	file.close();
	EXPECT_TRUE(var2.getName()=="Pression");
	for(int i=0;i<10;i++){
		EXPECT_EQ(i*10,var2[i]);
	}

}
/*----------------------------------------------------------------------------*/
TEST_F(SerializationTest ,testSerializeUserClassVariable)
{
	std::fstream file;
	file.open("test.dat", std::ios::out|std::ios::binary);
	gmds::Variable<TestVarSerialization> var("User class");
	var.setDomain(10);
	for(int i=0;i<10;i++){
		var[i].m_x=i;
		var[i].m_y=i+0.5;
	}

	var.serialize(file);
	file.close();
	file.open("test.dat", std::ios::in|std::ios::binary);
	gmds::Variable<TestVarSerialization> var2;
	var2.unserialize(file);

	file.close();
	EXPECT_TRUE(var2.getName()=="User class");
	for(int i=0;i<10;i++){
		EXPECT_EQ(i	   , var2[i].m_x);
		EXPECT_EQ(i+0.5, var2[i].m_y);
	}
}
/*----------------------------------------------------------------------------*/
TEST_F(SerializationTest ,mesh2D)
{
	const int mask = DIM2|N|F|F2N|N2F;
	IGMesh mesh(mask);
	Node t[3][3];
	Face f[2][2];

	/* nodes creation */
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			t[i][j] = mesh.newNode(i,j);

	/* faces creation */
	for(int i=0;i<2;i++)
		for(int j=0;j<2;j++)	{
			f[i][j] = mesh.newQuad(t[i][j],t[i+1][j],t[i+1][j+1],t[i][j+1]);
			t[i  ][j  ].add<Face>(f[i][j]);
			t[i+1][j  ].add<Face>(f[i][j]);
			t[i+1][j+1].add<Face>(f[i][j]);
			t[i  ][j+1].add<Face>(f[i][j]);
		}
	std::vector<TCellID> faces_of_11_bef = t[1][1].getIDs<Face>();

	int mark = mesh.getNewMark<Node>();

	mesh.mark(t[0][0],mark);
	mesh.mark(t[0][1],mark);

	int id_n00 = t[0][0].getID();
	int id_n01 = t[0][1].getID();


	std::vector<TCellID> node_id_bef;
	IGMesh::node_iterator itn = mesh.nodes_begin();


	for(;!itn.isDone();itn.next())
		node_id_bef.push_back(itn.value().getID());

	t[1][1].remove<Face>(f[0][1]);

	mesh.deleteFace(f[0][1]);


	std::vector<TCellID> face_id_bef;
	IGMesh::face_iterator itf = mesh.faces_begin();

	for(;!itf.isDone();itf.next())
		face_id_bef.push_back(itf.value().getID());

	std::fstream file;
	file.open("test_2D.dat", std::ios::out|std::ios::binary);
	mesh.serialize(file);
	file.close();

	// we make the mesh empty
	mesh.clear();

	EXPECT_EQ(0,mesh.getNbNodes());
	EXPECT_EQ(0,mesh.getNbFaces());

	file.open("test_2D.dat", std::ios::in|std::ios::binary);

	mesh.unserialize(file);

	EXPECT_TRUE(mesh.isMarked(mesh.get<Node>(id_n00),mark));
	EXPECT_TRUE(mesh.isMarked(mesh.get<Node>(id_n01),mark));

	EXPECT_EQ(9,mesh.getNbNodes());
	EXPECT_EQ(3,mesh.getNbFaces());

	itn = mesh.nodes_begin();

	int nbN=0;
	for(;!itn.isDone();itn.next()){
		Node n = itn.value();
		EXPECT_EQ(node_id_bef[nbN], n.getID());
		nbN++;
	}

	EXPECT_EQ(9,nbN);


	itf = mesh.faces_begin();

	int nbF=0;
	for(;!itf.isDone();itf.next()){
		Face f = itf.value();
		EXPECT_EQ(face_id_bef[nbF], f.getID());
		nbF++;
	}
	EXPECT_EQ(3,nbF);

}

/*----------------------------------------------------------------------------*/
TEST_F(SerializationTest ,mesh2D_2)
{
	const int mask = DIM2|N|F|F2N|N2F;
	IGMesh mesh(mask);
	Node t[3][3];
	Face f[2][2];

	/* nodes creation */
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			t[i][j] = mesh.newNode(i,j);

	/* faces creation */
	for(int i=0;i<2;i++)
		for(int j=0;j<2;j++)	{
			f[i][j] = mesh.newQuad(t[i][j],t[i+1][j],t[i+1][j+1],t[i][j+1]);
			t[i  ][j  ].add<Face>(f[i][j]);
			t[i+1][j  ].add<Face>(f[i][j]);
			t[i+1][j+1].add<Face>(f[i][j]);
			t[i  ][j+1].add<Face>(f[i][j]);
		}
	std::vector<TCellID> faces_of_11_bef = t[1][1].getIDs<Face>();

	int mark = mesh.getNewMark<Node>();

	mesh.mark(t[0][0],mark);
	mesh.mark(t[0][1],mark);

	int id_n00 = t[0][0].getID();
	int id_n01 = t[0][1].getID();


	std::vector<TCellID> node_id_bef;
	IGMesh::node_iterator itn = mesh.nodes_begin();


	for(;!itn.isDone();itn.next())
		node_id_bef.push_back(itn.value().getID());

	t[1][1].remove<Face>(f[0][1]);

	mesh.deleteFace(f[0][1]);


	std::vector<TCellID> face_id_bef;
	IGMesh::face_iterator itf = mesh.faces_begin();

	for(;!itf.isDone();itf.next())
		face_id_bef.push_back(itf.value().getID());

	std::ostringstream sub_stream;
	mesh.serialize(sub_stream);

	std::string str = sub_stream.str();
	int sub_size   = str.size()*sizeof(char);
	const char* sub_char = str.c_str();
	std::istrstream sub_stream2((char*)sub_char,sub_size);

	IGMesh mesh2(mask);
	mesh2.unserialize(sub_stream2);

	EXPECT_TRUE(mesh2.isMarked(mesh2.get<Node>(id_n00),mark));
	EXPECT_TRUE(mesh2.isMarked(mesh2.get<Node>(id_n01),mark));

	EXPECT_EQ(9,mesh2.getNbNodes());
	EXPECT_EQ(3,mesh2.getNbFaces());

	itn = mesh2.nodes_begin();

	int nbN=0;
	for(;!itn.isDone();itn.next()){
		Node n = itn.value();
		EXPECT_EQ(node_id_bef[nbN], n.getID());
		nbN++;
	}

	EXPECT_EQ(9,nbN);


	itf = mesh2.faces_begin();

	int nbF=0;
	for(;!itf.isDone();itf.next()){
		Face f = itf.value();
		EXPECT_EQ(face_id_bef[nbF], f.getID());
		nbF++;
	}
	EXPECT_EQ(3,nbF);


}

///*----------------------------------------------------------------------------*/
//void SerializationTestSuite::testSerializeMesh3D(){
//	const int mask = DIM3|N|F|R|R2N|N2R|F2N|N2F|R2F|F2R;
//	Mesh<mask> mesh;
//
//	int si = 10;
//	int sj = 10;
//	int sk = 10;
//	Node* t[si][sj][sk];
//
//	/* nodes creation */
//	for(int i=0;i<si;i++)
//		for(int j=0;j<sj;j++)
//			for(int k=0;k<sk;k++)
//			t[i][j][k] = mesh.newNode((TCoord)i,(TCoord)j,(TCoord)k);
//
//	Region* r;
//	for(int i=0;i<si-1;i++)
//		for(int j=0;j<sj-1;j++)
//			for(int k=0;k<sk-1;k++)
//				r = mesh.newHex(t[i][j][k],t[i+1][j][k],t[i+1][j+1][k],t[i][j+1][k],
//							t[i][j][k+1],t[i+1][j][k+1],t[i+1][j+1][k+1],t[i][j+1][k+1]);
//
//
//	MeshDoctor<mask> doctor(mesh);
//	doctor.buildFacesAndR2F();
//	doctor.updateUpwardConnectivity();
//
//	int nb_init_regions = mesh.getNbRegions();
//	int nb_init_faces 	= mesh.getNbFaces();
//	int nb_init_nodes 	= mesh.getNbNodes();
//
//	mesh.deleteFace(mesh.getLFace(2));
//	mesh.deleteFace(mesh.getLFace(11));
//	mesh.deleteFace(mesh.getLFace(45));
//	mesh.deleteFace(mesh.getLFace(7));
//	mesh.deleteFace(mesh.getLFace(8));
//	mesh.deleteFace(mesh.getLFace(9));
//
//	mesh.deleteNode(mesh.getLNode(0));
//	mesh.deleteNode(mesh.getLNode(1));
//	mesh.deleteNode(mesh.getLNode(2));
//	mesh.deleteNode(mesh.getLNode(3));
//	mesh.deleteNode(mesh.getLNode(4));
//	mesh.deleteNode(mesh.getLNode(10));
//	mesh.deleteNode(mesh.getLNode(20));
//	mesh.deleteNode(mesh.getLNode(30));
//	mesh.deleteNode(mesh.getLNode(40));
//	mesh.deleteNode(mesh.getLNode(50));
//
//	std::vector<id> node_id_bef;
//	Mesh<mask>::nodes_iterator itn = mesh.nodes_begin();
//
//
//	for(;!itn->isDone();itn->next())
//		node_id_bef.push_back((itn->currentItem())->getID());
//
//	std::vector<id> face_id_bef;
//	Mesh<mask>::faces_iterator itf = mesh.faces_begin();
//
//
//	for(;!itf->isDone();itf->next())
//		face_id_bef.push_back((itf->currentItem())->getID());
//
//	std::vector<id> region_id_bef;
//	Mesh<mask>::regions_iterator itr = mesh.regions_begin();
//
//
//	for(;!itr->isDone();itr->next())
//		region_id_bef.push_back((itr->currentItem())->getID());
//
//	std::fstream file;
//	file.open("test.dat", std::ios::out|std::ios::binary);
//
//	mesh.serialize(file);
//
//	file.close();
//	mesh.deleteRegion(mesh.getLRegion(0));
//	mesh.deleteRegion(mesh.getLRegion(10));
//	mesh.deleteRegion(mesh.getLRegion(20));
//	mesh.deleteRegion(mesh.getLRegion(21));
//	mesh.deleteRegion(mesh.getLRegion(30));
//	mesh.deleteRegion(mesh.getLRegion(35));
//	mesh.deleteRegion(mesh.getLRegion(36));
//
//	std::fstream file2;
//	file2.open("test.dat", std::ios::in|std::ios::binary);
//	mesh.unserialize(file2);
//
//	CPPUNIT_ASSERT(mesh.getNbNodes()==nb_init_nodes-10);
//	CPPUNIT_ASSERT(mesh.getNbFaces()==nb_init_faces-6);
//	CPPUNIT_ASSERT(mesh.getNbRegions()==nb_init_regions);
//
//	CPPUNIT_ASSERT(mesh.getLNode(0)==0);
//	CPPUNIT_ASSERT(mesh.getLNode(1)==0);
//	CPPUNIT_ASSERT(mesh.getLNode(2)==0);
//	CPPUNIT_ASSERT(mesh.getLNode(3)==0);
//	CPPUNIT_ASSERT(mesh.getLNode(4)==0);
//	CPPUNIT_ASSERT(mesh.getLNode(10)==0);
//	CPPUNIT_ASSERT(mesh.getLNode(20)==0);
//	CPPUNIT_ASSERT(mesh.getLNode(30)==0);
//	CPPUNIT_ASSERT(mesh.getLNode(40)==0);
//	CPPUNIT_ASSERT(mesh.getLNode(50)==0);
//
//	CPPUNIT_ASSERT(mesh.getLFace(2 )==0);
//	CPPUNIT_ASSERT(mesh.getLFace(11)==0);
//	CPPUNIT_ASSERT(mesh.getLFace(45)==0);
//	CPPUNIT_ASSERT(mesh.getLFace(7 )==0);
//	CPPUNIT_ASSERT(mesh.getLFace(8 )==0);
//	CPPUNIT_ASSERT(mesh.getLFace(9 )==0);
//
//	itn = mesh.nodes_begin();
//
//	int nbN=0;
//	for(;!itn->isDone();itn->next()){
//		Node* n = itn->currentItem();
//		CPPUNIT_ASSERT(n->getID()==node_id_bef[nbN]);
//		nbN++;
//	}
//	CPPUNIT_ASSERT(nbN==nb_init_nodes-10);
//
//	itf = mesh.faces_begin();
//
//	int nbF=0;
//	for(;!itf->isDone();itf->next()){
//		Face* f = itf->currentItem();
//		CPPUNIT_ASSERT(f->getID()==face_id_bef[nbF]);
//		nbF++;
//	}
//	CPPUNIT_ASSERT(nbF==nb_init_faces-6);
//
//	itr = mesh.regions_begin();
//
//	int nbR=0;
//	for(;!itr->isDone();itr->next()){
//		Region* r = itr->currentItem();
//		CPPUNIT_ASSERT(r->getID()==region_id_bef[nbR]);
//		nbR++;
//	}
//	CPPUNIT_ASSERT(nbR==nb_init_regions);
//
//}
