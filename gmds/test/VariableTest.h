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
#include <GMDS/Utils/Variable.h>
#include <GMDS/Utils/VariableManager.h>
/*----------------------------------------------------------------------------*/
#include <string>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
class VariableTest: public ::testing::Test {

  protected:
	VariableTest(){;}
    virtual ~VariableTest(){;}
};
/*----------------------------------------------------------------------------*/
class TestVar{
public:
	TestVar(const int x=1, const int y=1):m_x(x),m_y(y){;}

	int m_x;
	int m_y;
};
/*----------------------------------------------------------------------------*/
TEST_F(VariableTest,varInt){
	Variable<int> v;
	EXPECT_EQ(0,v.getNbValues());
	v.setDomain(10);
	v[1]=20;
	EXPECT_EQ(20,v[1]);
}
/*----------------------------------------------------------------------------*/
TEST_F(VariableTest,varClassInit){
	Variable<TestVar> v;
	EXPECT_EQ(0,v.getNbValues());
}
/*----------------------------------------------------------------------------*/
TEST_F(VariableTest,varPtrInit){
	Variable<TestVar*> v;
	EXPECT_EQ(0,v.getNbValues());
}
/*----------------------------------------------------------------------------*/
TEST_F(VariableTest,varStringInit){

	Variable<std::string> v;
	v.setDomain(4);
	v[0]="I";
	v[1]="am";
	v[2]="the";
	v[3]="test";
	EXPECT_EQ(4,v.getNbValues());
	v.removeEntry(2);
	EXPECT_EQ(3,v.getNbValues());
	v[2] = "just a";
	EXPECT_EQ(4,v.getNbValues());
}
/*----------------------------------------------------------------------------*/
TEST_F(VariableTest,varManagerTest){

	VariableManager manager;
	Variable<int>* 			v1 = manager.newVariable<int>("var1");
	Variable<TestVar*>* 	v2 = manager.newVariable<TestVar*>("var2");
	Variable<TestVar>* 		v3 = manager.newVariable<TestVar>("var3");
	Variable<std::string>* 	v4 = manager.newVariable<std::string>("var4");

	v1->setDomain(4);
	v1->set(1,13);
	v1->set(2,14);
	v1->set(3,15);
	EXPECT_EQ(3,v1->getNbValues());
	TestVar tv;
	v3->setDomain(10);
	v3->set(1,tv);
	EXPECT_EQ(1,v3->getNbValues());
	v4->setDomain(100);
	v4->set(1,"toto");
	v4->set(2,"titi");
	EXPECT_EQ(2,v4->getNbValues());
	EXPECT_EQ(4,manager.getNbVariables());
	manager.deleteVariable("var1");
	manager.deleteVariable("var2");
	manager.deleteVariable("var3");
	manager.deleteVariable("var4");

}
/*----------------------------------------------------------------------------*/
TEST_F(VariableTest,varManagerDelete){
	VariableManager manager;
	Variable<int>* 			v1 = manager.newVariable<int>("var1");
	Variable<std::string>* 	v2 = manager.newVariable<std::string>("var2");
	v1->setDomain(4);
	v1->set(1,13);
	v1->set(2,14);
	v1->set(3,15);
	v2->setDomain(100);
	v2->set(1,"toto");
	v2->set(2,"titi");
	manager.deleteVariable("var1");
	EXPECT_EQ(1,manager.getNbVariables());
	EXPECT_TRUE("titi"==(*v2)[2]);
	manager.deleteVariable("var2");
}
/*----------------------------------------------------------------------------*/
TEST_F(VariableTest,varManagerDelete2){
	VariableManager manager;
	Variable<int>* 			v1 = manager.newVariable<int>("var1");
	Variable<std::string>* 	v2 = manager.newVariable<std::string>("var2");

	v1->setDomain(4);
	v1->set(1,13);
	v1->set(2,14);
	v1->set(3,15);
	v2->setDomain(100);
	v2->set(1,"toto");
	v2->set(2,"titi");
	manager.deleteVariable("var2");
	EXPECT_EQ(1,manager.getNbVariables());
	manager.deleteVariable("var1");
}
/*----------------------------------------------------------------------------*/
