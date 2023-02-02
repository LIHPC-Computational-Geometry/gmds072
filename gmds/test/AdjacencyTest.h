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

using namespace gmds;

class AdjacencyTest: public ::testing::Test {

  protected:
    AdjacencyTest(){;}
    virtual ~AdjacencyTest(){;}
};
/*----------------------------------------------------------------------------*/
TEST_F(AdjacencyTest,meshInitModel) {
  MeshModel mod = F|N;
  IGMesh m(mod);
  EXPECT_EQ(true, mod==m.getModel());

}
/*----------------------------------------------------------------------------*/
TEST_F(AdjacencyTest,meshInitSize) {
  MeshModel mod = F|N;
  IGMesh m(mod);
  EXPECT_EQ(0,m.getNbNodes());
  EXPECT_EQ(0,m.getNbFaces());

}
/*----------------------------------------------------------------------------*/

TEST_F(AdjacencyTest,meshF2N_Get) {
	MeshModel mod = F|N| F2N| F2F;
	IGMesh m(mod);
	Node init_nodes[3];
	init_nodes[0] = m.newNode(0,0,0);
	init_nodes[1] = m.newNode(1,0,0);
	init_nodes[2] = m.newNode(0,1,0);
	Face f = m.newTriangle(init_nodes[0],init_nodes[1],init_nodes[2]);

	std::vector<Node> nodes = f.get<Node>();
	EXPECT_EQ(init_nodes[0].getID(),nodes[0].getID());
	EXPECT_EQ(init_nodes[1].getID(),nodes[1].getID());
	EXPECT_EQ(init_nodes[2].getID(),nodes[2].getID());
}
/*----------------------------------------------------------------------------*/
TEST_F(AdjacencyTest,meshF2E_Get_Fail)
{
	MeshModel mod = F| N | F2N;
	IGMesh m(mod);
	Node init_nodes[3];
	init_nodes[0] = m.newNode(0,0,0);
	init_nodes[1] = m.newNode(1,0,0);
	init_nodes[2] = m.newNode(0,1,0);
	Face f = m.newTriangle(init_nodes[0],init_nodes[1],init_nodes[2]);

	EXPECT_THROW({
	std::vector<Edge> edges = f.get<Edge>();
	}, GMDSException);
}
/*----------------------------------------------------------------------------*/

TEST_F(AdjacencyTest,meshF2N_Set) {
	MeshModel mod = F|N| F2N| F2F;
	IGMesh m(mod);
	Node init_nodes[3];
	init_nodes[0] = m.newNode(0,0,0);
	init_nodes[1] = m.newNode(1,0,0);
	init_nodes[2] = m.newNode(0,1,0);
	Face f = m.newTriangle(init_nodes[0],init_nodes[1],init_nodes[2]);

	std::vector<Node> nodes = f.get<Node>();
	EXPECT_EQ(init_nodes[0].getID(),nodes[0].getID());
	EXPECT_EQ(init_nodes[1].getID(),nodes[1].getID());
	EXPECT_EQ(init_nodes[2].getID(),nodes[2].getID());

	Node new_nodes[3];
	new_nodes[0] = m.newNode(0,0,10);
	new_nodes[1] = m.newNode(1,0,10);
	new_nodes[2] = m.newNode(0,1,10);
	std::vector<Node> new_vect_nodes;
	new_vect_nodes.push_back(new_nodes[0]);
	new_vect_nodes.push_back(new_nodes[1]);
	new_vect_nodes.push_back(new_nodes[2]);

	f.set<Node>(new_vect_nodes);
	nodes = f.get<Node>();
	EXPECT_EQ(new_nodes[0].getID(),nodes[0].getID());
	EXPECT_EQ(new_nodes[1].getID(),nodes[1].getID());
	EXPECT_EQ(new_nodes[2].getID(),nodes[2].getID());
}
/*----------------------------------------------------------------------------*/
TEST_F(AdjacencyTest,meshUpdatePolygon)
{
  MeshModel mod = F|N| F2N| F2F;
  IGMesh m(mod);
  std::vector<Node> init_nodes;
  init_nodes.resize(3);
  init_nodes[0] = m.newNode(0,0,0);
  init_nodes[1] = m.newNode(1,0,0);
  init_nodes[2] = m.newNode(0,1,0);

  Face f = m.newPolygon(init_nodes);
  EXPECT_EQ(0, f.getID());

  EXPECT_EQ(3, f.getNbNodes());


  std::vector<Node> nodes;
  f.get<Node>(nodes);
  for(unsigned int in=0;in<nodes.size();in++){
    Node ni = nodes[in];
    EXPECT_EQ(init_nodes[in].getID(),ni.getID());
  }

  init_nodes.resize(4);
  init_nodes[0] = m.newNode(0,0,10);
  init_nodes[1] = m.newNode(1,0,10);
  init_nodes[2] = m.newNode(0,1,10);
  init_nodes[3] = m.newNode(1,1,10);

  f.set<Node>(init_nodes);
  EXPECT_EQ(0, f.getID());

  EXPECT_EQ(4, f.getNbNodes());


  f.get<Node>(nodes);
  for(unsigned int in=0;in<nodes.size();in++){
	  Node ni = nodes[in];
	  EXPECT_EQ(init_nodes[in].getID(),ni.getID());
	  EXPECT_EQ(init_nodes[in].X(),ni.X());
	  EXPECT_EQ(init_nodes[in].Y(),ni.Y());
	  EXPECT_EQ(init_nodes[in].Z(),ni.Z());
  }
}
/*----------------------------------------------------------------------------*/
TEST_F(AdjacencyTest,adjacency3D_1)
{
	MeshModel mod = DIM3|N|E|F|R|R2F|F2E|E2N|N2R;
	IGMesh mesh(mod);
	Node n;

	Node t[4];
	Edge e[6];
	Face f[4];
	Region r;

	/* nodes creation */
	t[0] = mesh.newNode(0,0,0);
	t[1] = mesh.newNode(1,0,0);
	t[2] = mesh.newNode(0,1,0);
	t[3] = mesh.newNode(0,0,1);


	/* edges creation */
	e[0 ] = mesh.newEdge(t[0],t[1]);
	e[1 ] = mesh.newEdge(t[1],t[2]);
	e[2 ] = mesh.newEdge(t[0],t[2]);
	e[3 ] = mesh.newEdge(t[0],t[3]);
	e[4 ] = mesh.newEdge(t[1],t[3]);
	e[5 ] = mesh.newEdge(t[2],t[3]);

	f[0] = mesh.newTriangle();
	std::vector<Edge> edges;
	edges.resize(3);
	edges[0]=e[0];
	edges[1]=e[4];
	edges[2]=e[3];
	f[0].set<Edge>(edges);

	f[1] = mesh.newTriangle();
	edges[0]=e[1];
	edges[1]=e[5];
	edges[2]=e[4];
	f[1].set<Edge>(edges);

	f[2] = mesh.newTriangle();
	edges[0]=e[2];
	edges[1]=e[5];
	edges[2]=e[3];
	f[2].set<Edge>(edges);

	f[3] = mesh.newTriangle();
	edges[0]=e[0];
	edges[1]=e[1];
	edges[2]=e[2];
	f[3].set<Edge>(edges);

	r = mesh.newTet();
	std::vector<Face> faces;
	faces.resize(4);
	faces[0] = f[0];
	faces[1] = f[1];
	faces[2] = f[2];
	faces[3] = f[3];
	r.set<Face>(faces);

	/* connection N->R */
	t[0].add<Region>(r);
	t[1].add<Region>(r);
	t[2].add<Region>(r);
	t[3].add<Region>(r);


	EXPECT_THROW({
		r.get<Edge>().size();
	}, GMDSException);

	EXPECT_THROW({
		r.get<Node>().size();
	}, GMDSException);


	for(int i=0;i<4;i++){
		EXPECT_EQ(1, t[i].get<Region>().size());

		EXPECT_THROW({ t[i].get<Face>().size();}, GMDSException);
		EXPECT_THROW({ t[i].get<Edge>().size();}, GMDSException);
		EXPECT_THROW({ t[i].get<Node>().size();}, GMDSException);
	}
}
/*----------------------------------------------------------------------------*/
TEST_F(AdjacencyTest,adjacency3D_2)
{
	MeshModel mod = DIM3|N|F|R|R2F|F2R|F2N|N2F;
	IGMesh mesh(mod);

	Node t[4];
	Face f[4];
	Region r;

	/* nodes creation */
	t[0] = mesh.newNode(0,0,0);
	t[1] = mesh.newNode(1,0,0);
	t[2] = mesh.newNode(0,1,0);
	t[3] = mesh.newNode(0,0,1);


	f[0] = mesh.newTriangle(t[0],t[1],t[2]);
	f[1] = mesh.newTriangle(t[0],t[1],t[3]);
	f[2] = mesh.newTriangle(t[0],t[2],t[3]);
	f[3] = mesh.newTriangle(t[1],t[2],t[3]);

	r = mesh.newTet();
	std::vector<Face> faces;
	faces.resize(4);
	faces[0] = f[0];
	faces[1] = f[1];
	faces[2] = f[2];
	faces[3] = f[3];
	r.set<Face>(faces);

	/* connection N->F */
	t[0].add<Face>(f[0]);
	t[0].add<Face>(f[1]);
	t[0].add<Face>(f[2]);

	t[1].add<Face>(f[0]);
	t[1].add<Face>(f[1]);
	t[1].add<Face>(f[3]);

	t[2].add<Face>(f[0]);
	t[2].add<Face>(f[2]);
	t[2].add<Face>(f[3]);

	t[3].add<Face>(f[1]);
	t[3].add<Face>(f[2]);
	t[3].add<Face>(f[3]);

	/* connection F->R */
	f[0].add<Region>(r);
	f[1].add<Region>(r);
	f[2].add<Region>(r);
	f[3].add<Region>(r);

	EXPECT_THROW({r.get<Node>().size();}, GMDSException);
	EXPECT_THROW({r.get<Region>().size();}, GMDSException);


	for(int i=0;i<4;i++){
		EXPECT_THROW({f[i].get<Face>().size();}, GMDSException);
		EXPECT_EQ(3, f[i].get<Node>().size());
		EXPECT_EQ(1, f[i].get<Region>().size());
	}

	for(int i=0;i<4;i++){
		EXPECT_THROW({t[i].get<Region>().size();}, GMDSException);
		EXPECT_EQ(3, t[i].get<Face>().size());
		EXPECT_THROW({t[i].get<Node>().size();}, GMDSException);
	}
}
/*----------------------------------------------------------------------------*/
TEST_F(AdjacencyTest,test2D_N_F_F2N_N2F)
{
	MeshModel mod = DIM2|N|F|F2N|N2F;
	IGMesh mesh(mod);

	const int si = 10;
	const int sj = 10;
	Node t[si][sj];
	Face f[si-1][sj-1];

	/* nodes creation */
	for(int i=0;i<si;i++)
		for(int j=0;j<sj;j++){
			t[i][j] = mesh.newNode(i,j);
		}


	for(int i=1;i<si;i++)
		for(int j=1;j<sj;j++){
			f[i-1][j-1] = mesh.newQuad(t[i-1][j-1],t[i-1][j],t[i][j],t[i][j-1]);
			t[i-1][j-1].add<Face>(f[i-1][j-1]);
			t[i][j].add<Face>(f[i-1][j-1]);
			t[i-1][j].add<Face>(f[i-1][j-1]);
			t[i][j-1].add<Face>(f[i-1][j-1]);
		}

	std::vector<TCellID> nodes00 =f[0][0].getIDs<Node>();
	EXPECT_EQ(4, nodes00.size());

	std::vector<TCellID> faces00 =t[0][0].getIDs<Face>();
	EXPECT_EQ(1, faces00.size());

	std::vector<TCellID> faces11=t[1][1].getIDs<Face>();
	EXPECT_EQ(4, faces11.size());


}
/*----------------------------------------------------------------------------*/
TEST_F(AdjacencyTest,nodePoint3D)
{
	MeshModel mod = R|N| R2N| DIM3;
	IGMesh mesh(mod);
	Node n;

	/* nodes creation */
	n = mesh.newNode(4.5,2.0,3.0);
	math::Point p = n.getPoint();

	EXPECT_EQ(4.5, p.X());
	EXPECT_EQ(2.0, p.Y());
	EXPECT_EQ(3.0, p.Z());

	p.setXYZ(1.0,1.0,1.0);

	n.setPoint(p);
	EXPECT_EQ(1.0, p.X());
	EXPECT_EQ(1.0, p.Y());
	EXPECT_EQ(1.0, p.Z());
}
/*----------------------------------------------------------------------------*/
TEST_F(AdjacencyTest,RetNodeModifPoint)
{
	MeshModel mod = R|N| R2N| DIM3;
	IGMesh mesh(mod);
	Node n0 = mesh.newNode(0,0,0);
	Node n1 = mesh.newNode(1,0,0);
	Node n2 = mesh.newNode(0,1,0);
	Node n3 = mesh.newNode(0,0,1);

	Region r = mesh.newTet(n0,n1,n2,n3);

	std::vector<Node> r_nodes = r.get<Node>();

	n0.setZ(1.0);
	n1.setZ(1.0);
	n2.setZ(1.0);

	for(unsigned int i=0;i<r_nodes.size();i++){
		Node current = r_nodes[i];
		EXPECT_EQ(1.0, current.Z());
	}
}/*----------------------------------------------------------------------------*/

