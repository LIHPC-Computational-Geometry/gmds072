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
#include<string>
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IG.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
class MeshVariableTest: public ::testing::Test {

  protected:
	MeshVariableTest(){;}
    virtual ~MeshVariableTest(){;}
};
/*----------------------------------------------------------------------------*/
class Color{
public:
	Color(double r=0, double g=0, double b=0):r_(r),g_(g),b_(b){}

	Color(const Color& c):r_(c.r_),g_(c.g_),b_(c.b_){}

	double r() const {return r_;}
	double g() const {return g_;}
	double b() const {return b_;}

	void setR(const double& r){r_=r;}
	void setG(const double& g){g_=g;}
	void setB(const double& b){b_=b;}
private:
	double r_;
	double g_;
	double b_;
};
/*----------------------------------------------------------------------------*/
TEST_F(MeshVariableTest,initOfIntVariable) {
	MeshModel mod = DIM2|N|F|F2N;
	IGMesh mesh(mod);
	mesh.newNode(0,0);
	mesh.newNode(0,1);
	Variable<int> *v = mesh.newVariable<int>(GMDS_NODE,"var1");

	mesh.newNode(1,0);
	mesh.newNode(1,1);

	IGMesh::node_iterator it = mesh.nodes_begin();

	for(;!it.isDone();it.next()){
		TCellID i = (it.value()).getID();
		int n = (*v)[i];
		EXPECT_EQ(0,n);
		(*v)[i]=1;
	}
	for(it=mesh.nodes_begin();!it.isDone();it.next()){
		TCellID i = (it.value()).getID();
		int n = (*v)[i];
		EXPECT_EQ(1,n);
	}
}
/*----------------------------------------------------------------------------*/
TEST_F(MeshVariableTest,initOfColorVariable) {
	MeshModel mod = DIM2|N|F|F2N;
	IGMesh mesh(mod);

	Variable<Color> *v =mesh.newVariable<Color>(GMDS_FACE,"color");

	Node n0 = mesh.newNode(0,0,0);
	Node n1 = mesh.newNode(1,0,0);
	Node n2 = mesh.newNode(0,1,0);
	Node n3 = mesh.newNode(1,1,0);

	mesh.newQuad(n0,n1,n2,n3);
	mesh.newTriangle(n0,n1,n2);
	mesh.newTriangle(n1,n2,n3);


	IGMesh::face_iterator it = mesh.faces_begin();

	int i=1;
	for(;!it.isDone();it.next()){
		TCellID k = (it.value()).getID();
		Color c = (*v)[k];
		EXPECT_EQ(0,c.r());
		EXPECT_EQ(0,(*v)[k].g());
		EXPECT_EQ(0,(*v)[k].b());
		(*v)[k].setR(i++);
		(*v)[k].setG(i++);
		(*v)[k].setB(i++);
	}
	i=1;
	for(it=mesh.faces_begin();!it.isDone();it.next()){
		TCellID k = (it.value()).getID();
		Color c = (*v)[k];
		EXPECT_EQ(i++,c.r());
		EXPECT_EQ(i++,c.g());
		EXPECT_EQ(i++,c.b());
	}
}
/*----------------------------------------------------------------------------*/
TEST_F(MeshVariableTest,initOfColorVariableAfter) {
	MeshModel mod = DIM2|N|F|F2N;
	IGMesh mesh(mod);


	Node n0 = mesh.newNode(0,0,0);
	Node n1 = mesh.newNode(1,0,0);
	Node n2 = mesh.newNode(0,1,0);
	Node n3 = mesh.newNode(1,1,0);

	mesh.newQuad(n0,n1,n2,n3);
	mesh.newTriangle(n0,n1,n2);
	mesh.newTriangle(n1,n2,n3);
	Variable<Color> *v = mesh.newVariable<Color>(GMDS_FACE,"color");

	EXPECT_EQ(3, v->getNbValues());
	IGMesh::face_iterator it = mesh.faces_begin();

	int i=1;
	for(;!it.isDone();it.next()){
		TCellID k = (it.value()).getID();
		Color c = (*v)[k];
		EXPECT_EQ(0,c.r());
		EXPECT_EQ(0,(*v)[k].g());
		EXPECT_EQ(0,(*v)[k].b());
		(*v)[k].setR(i++);
		(*v)[k].setG(i++);
		(*v)[k].setB(i++);
	}
	i=1;
	for(it=mesh.faces_begin();!it.isDone();it.next()){
		TCellID k = (it.value()).getID();
		Color c = (*v)[k];
		EXPECT_EQ(i++,c.r());
		EXPECT_EQ(i++,c.g());
		EXPECT_EQ(i++,c.b());
	}
}

/*----------------------------------------------------------------------------*/
TEST_F(MeshVariableTest,deleteInt) {
	MeshModel mod = DIM2|N|F|F2N;
	IGMesh mesh(mod);

	mesh.newNode(0,0);
	mesh.newNode(0,1);
	mesh.newNode(1,0);
	mesh.newNode(1,1);

	Variable<int> *v =mesh.newVariable<int>(GMDS_NODE,"int");

	IGMesh::node_iterator it = mesh.nodes_begin();

	for(;!it.isDone();it.next()){
		Node n = it.value();
		int i = (*v)[n.getID()];
		EXPECT_EQ(0,i);
		(*v)[n.getID()]=10;
	}

	mesh.deleteNode(1);
	mesh.deleteNode(3);
	mesh.newNode(2,2);
	mesh.newNode(3,3);

	for(it=mesh.nodes_begin();!it.isDone();it.next()){
		Node n = it.value();
		int i = (*v)[n.getID()];
		if(n.getID()%2==0)
			EXPECT_EQ(10,i);
		else
			EXPECT_EQ(0,i);
	}
}

/*----------------------------------------------------------------------------*/
TEST_F(MeshVariableTest,getInt) {
	MeshModel mod = DIM2|N|F|F2N;
	IGMesh mesh(mod);
	mesh.newNode(0,0);
	mesh.newNode(0,1);
	mesh.newNode(1,0);
	mesh.newNode(1,1);

	Variable<int> *v =mesh.newVariable<int>(GMDS_NODE,"int");

	IGMesh::node_iterator it = mesh.nodes_begin();

	for(;!it.isDone();it.next()){
		Node n = it.value();
		int i = (*v)[n.getID()];
		EXPECT_EQ(0,i);
		(*v)[n.getID()]=10;
	}

	Variable<int> *v2 = mesh.getVariable<int>(GMDS_NODE,"int");
	EXPECT_TRUE(v2==v);


	Variable<double> *v3 = mesh.getVariable<double>(GMDS_NODE,"int");
	EXPECT_TRUE(v3==0);
}
/*----------------------------------------------------------------------------*/
TEST_F(MeshVariableTest,doesExistYES) {
	MeshModel mod = DIM2|N|F|F2N;
	IGMesh mesh(mod);
	mesh.newNode(0,0);
	mesh.newNode(0,1);
	mesh.newNode(1,0);
	mesh.newNode(1,1);

	Variable<int> *v =mesh.newVariable<int>(GMDS_NODE,"int");

	EXPECT_TRUE(mesh.doesVariableExist(GMDS_NODE,"int"));
}
/*----------------------------------------------------------------------------*/
TEST_F(MeshVariableTest,doesExistNO) {
	MeshModel mod = DIM2|N|F|F2N;
	IGMesh mesh(mod);

	mesh.newNode(0,0);
	mesh.newNode(0,1);
	mesh.newNode(1,0);
	mesh.newNode(1,1);

	Variable<int> *v =mesh.newVariable<int>(GMDS_NODE,"int");

	EXPECT_FALSE(mesh.doesVariableExist(GMDS_FACE,"int"));
	EXPECT_FALSE(mesh.doesVariableExist(GMDS_NODE,"int2"));
}
