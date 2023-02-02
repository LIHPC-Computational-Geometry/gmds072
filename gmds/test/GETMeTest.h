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
/*
 * GETMeTest.h
 *
 *  Created on: 20 oct. 2015
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <GMDS/Algo/GETMe.h>
#include <GMDS/IO/VTKWriter.h>
#include <GMDS/Utils/RandomGenerator.h>
/*----------------------------------------------------------------------------*/
class GETMeTest: public ::testing::Test {
protected:
	GETMeTest(){;}
	virtual ~GETMeTest(){;}
};

using namespace gmds;

/*----------------------------------------------------------------------------*/
TEST_F(GETMeTest,quadScaledJacobian) {

        MeshModel mod = DIM2|N|F|F2N;
        IGMesh mesh(mod);

        Node n1 = mesh.newNode(0,0,0);
        Node n2 = mesh.newNode(1,0,0);
        Node n3 = mesh.newNode(1,1,0);
        Node n4 = mesh.newNode(0,1,0);

        Face f = mesh.newQuad(n1,n2,n3,n4);

        double scaledJacobian = f.computeScaledJacobian2D();

        EXPECT_EQ(1.,scaledJacobian);
}
/*----------------------------------------------------------------------------*/
TEST_F(GETMeTest,callerTestQuad) {
        
        MeshModel mod = DIM2|N|F|F2N;
        IGMesh mesh(mod);

	const int nx = 3;
        const int ny = 3;

        gmds::Node nodes[nx+1][ny+1];

        const double xmin = 0.;
        const double ymin = 0.;

        const double dx = 1.;
        const double dy = 1.;

	int markFixedNodes = mesh.getNewMark<gmds::Node>();
	
	gmds::RandomGenerator randGen;
	randGen.init();

        for(int i=0; i<=nx; i++) {
                for(int j=0; j<=ny; j++) {

				double xpos = xmin+i*dx;
				double ypos = ymin+j*dy;

				if(i!= 0 && j!=0 && i!=nx && j!=ny) {
			
					double xrand = (2.*randGen.value()) -1.;
					double yrand = (2.*randGen.value()) -1.;

					xpos += 3*xrand;
					ypos += 3*yrand;
				}
				
                                nodes[i][j] = mesh.newNode(xpos,ypos);

				if(i== 0 || j==0 || i==nx || j==ny) {
					mesh.mark(nodes[i][j],markFixedNodes);
				}
                        }
	}
        

        gmds::Face faces[nx][ny];

        for(int i=0; i<nx; i++) {
                for(int j=0; j<ny; j++) {
		        gmds::Face f = mesh.newQuad(
				                    nodes[i][j]    ,
                                                nodes[i+1][j]     ,
                                                nodes[i+1][j+1]     ,
                                                nodes[i][j+1]
				       );
			}
		}

	gmds::VTKWriter<gmds::IGMesh> writer(mesh);
        writer.write("GETMecallerTestQuad_before",gmds::N|gmds::F);

	gmds::GETMe getMe(mesh);

	double minScaledJacobian = getMe.computeMinScaledJacobian2D();
	std::cout<<"minScaledJacobian "<<minScaledJacobian<<std::endl;

//getMe.exec(1000,0.9,markFixedNodes);
getMe.execSimult2D(     
                        100,
                        0.8,
                        0.5,
                        0.,
                        1.,
                        1.,
                        false,
                        true,
			false,
                        markFixedNodes
        );

	writer.write("GETMecallerTestQuad_after",gmds::N|gmds::F);

	minScaledJacobian = getMe.computeMinScaledJacobian2D();
	std::cout<<"minScaledJacobian "<<minScaledJacobian<<std::endl;

	mesh.unmarkAll<gmds::Node>(markFixedNodes);
	mesh.freeMark<gmds::Node>(markFixedNodes);		

	EXPECT_LT(0.8,minScaledJacobian);
}
/*----------------------------------------------------------------------------*/
TEST_F(GETMeTest,callerTestHex) {
        
        MeshModel mod = DIM3|N|R|R2N;
        IGMesh mesh(mod);

	const int nx = 3;
        const int ny = 3;
        const int nz = 3;

        gmds::Node nodes[nx+1][ny+1][nz+1];

        const double xmin = 0.;
        const double ymin = 0.;
        const double zmin = 0.;

        const double dx = 1.;
        const double dy = 1.;
        const double dz = 1.;

	int markFixedNodes = mesh.getNewMark<gmds::Node>();
	
	gmds::RandomGenerator randGen;
	randGen.init();

        for(int i=0; i<=nx; i++) {
                for(int j=0; j<=ny; j++) {
                        for(int k=0; k<=nz; k++) {

				double xpos = xmin+i*dx;
				double ypos = ymin+j*dy;
				double zpos = zmin+k*dz;

				if(i!= 0 && j!=0 && k!=0 && i!=nx && j!=ny && k!=nz) {
			
					double xrand = (2.*randGen.value()) -1.;
					double yrand = (2.*randGen.value()) -1.;
					double zrand = (2.*randGen.value()) -1.;

					xpos += 3*xrand;
					ypos += 3*yrand;
					zpos += 3*zrand; 
				}
				
                                nodes[i][j][k] = mesh.newNode(xpos,ypos,zpos);

				if(i== 0 || j==0 || k==0 || i==nx || j==ny || k==nz) {
					mesh.mark(nodes[i][j][k],markFixedNodes);
				}
                        }
                }
        }

        gmds::Region regions[nx][ny][nz];

        for(int i=0; i<nx; i++) {
                for(int j=0; j<ny; j++) {
                        for(int k=0; k<nz; k++) {
                                gmds::Region r = mesh.newHex(
                                                nodes[i][j][k+1],
                                                nodes[i][j+1][k+1],
                                                nodes[i+1][j+1][k+1],
                                                nodes[i+1][j][k+1],
                                                nodes[i][j][k],
                                                nodes[i][j+1][k],
                                                nodes[i+1][j+1][k],
                                                nodes[i+1][j][k]
                                );
			}
		}
	}

	gmds::VTKWriter<gmds::IGMesh> writer(mesh);
        writer.write("GETMecallerTestHex_before",gmds::N|gmds::R);

	gmds::GETMe getMe(mesh);

	double minScaledJacobian = getMe.computeMinNormalizedScaledJacobian();
	std::cout<<"minNormalizedScaledJacobian "<<minScaledJacobian<<std::endl;

	getMe.exec(100,10000,0.8,false,markFixedNodes);

	writer.write("GETMecallerTestHex_after",gmds::N|gmds::R);

	minScaledJacobian = getMe.computeMinNormalizedScaledJacobian();
	std::cout<<"minNormalizedScaledJacobian "<<minScaledJacobian<<std::endl;

	mesh.unmarkAll<gmds::Node>(markFixedNodes);
	mesh.freeMark<gmds::Node>(markFixedNodes);		

	EXPECT_LT(0.8,minScaledJacobian);
}
/*----------------------------------------------------------------------------*/
TEST_F(GETMeTest,getme_seq_oneHex) {
        
        MeshModel mod = DIM3|N|R|R2N;
        IGMesh mesh(mod);

	gmds::Node n0 = mesh.newNode( 0., 0., 0.);
	gmds::Node n1 = mesh.newNode( 1., 0., 0.);
	gmds::Node n2 = mesh.newNode( 1., 1., 0.);
	gmds::Node n3 = mesh.newNode( 0., 1., 0.);
	gmds::Node n4 = mesh.newNode( 0., 0., 1.);
	gmds::Node n5 = mesh.newNode( 1., 0., 1.);
	gmds::Node n6 = mesh.newNode( 1., 1., 1.);
	gmds::Node n7 = mesh.newNode( 0., 1., 0.5);
	
	int markFixedNodes = mesh.getNewMark<gmds::Node>();

	mesh.newHex(n0,n1,n2,n3,n4,n5,n6,n7);

	gmds::VTKWriter<gmds::IGMesh> writer(mesh);
        writer.write("GETMe_getme_seq_oneHex_before",gmds::N|gmds::R);

	gmds::GETMe getMe(mesh);

	double minScaledJacobian = getMe.computeMinNormalizedScaledJacobian();
	std::cout<<"minNormalizedScaledJacobian "<<minScaledJacobian<<std::endl;

	//getMe.exec(100,0.8,false,markFixedNodes);
	getMe.execSeq(
			 10000,
			 0.9,
			 0.01,
			 1.0,
			 0.0002,
			 0.0004,
			 0.0002,
			 false,
			 false,
			 false,
			 markFixedNodes
			 );

	writer.write("GETMe_getme_seq_oneHex_after",gmds::N|gmds::R);

	minScaledJacobian = getMe.computeMinNormalizedScaledJacobian();
	std::cout<<"minNormalizedScaledJacobian "<<minScaledJacobian<<std::endl;

	mesh.unmarkAll<gmds::Node>(markFixedNodes);
	mesh.freeMark<gmds::Node>(markFixedNodes);		

	EXPECT_LT(0.9,minScaledJacobian);
}
/*----------------------------------------------------------------------------*/
TEST_F(GETMeTest,getme_seq_hex) {
        
        MeshModel mod = DIM3|N|R|R2N;
        IGMesh mesh(mod);

	const int nx = 3;
        const int ny = 3;
        const int nz = 3;

        gmds::Node nodes[nx+1][ny+1][nz+1];

        const double xmin = 0.;
        const double ymin = 0.;
        const double zmin = 0.;

        const double dx = 1.;
        const double dy = 1.;
        const double dz = 1.;

	int markFixedNodes = mesh.getNewMark<gmds::Node>();
	
	gmds::RandomGenerator randGen;
	randGen.init();

        for(int i=0; i<=nx; i++) {
                for(int j=0; j<=ny; j++) {
                        for(int k=0; k<=nz; k++) {

				double xpos = xmin+i*dx;
				double ypos = ymin+j*dy;
				double zpos = zmin+k*dz;

				if(i!= 0 && j!=0 && k!=0 && i!=nx && j!=ny && k!=nz) {
			
					double xrand = (2.*randGen.value()) -1.;
					double yrand = (2.*randGen.value()) -1.;
					double zrand = (2.*randGen.value()) -1.;

					//xpos += 3*xrand;
					//ypos += 3*yrand;
					//zpos += 3*zrand; 
					xpos += 3*xrand;
					ypos += 3*yrand;
					zpos += 3*zrand;
				}
				
                                nodes[i][j][k] = mesh.newNode(xpos,ypos,zpos);

				if(i== 0 || j==0 || k==0 || i==nx || j==ny || k==nz) {
					mesh.mark(nodes[i][j][k],markFixedNodes);
				}
                        }
                }
        }

        gmds::Region regions[nx][ny][nz];

        for(int i=0; i<nx; i++) {
                for(int j=0; j<ny; j++) {
                        for(int k=0; k<nz; k++) {
                                gmds::Region r = mesh.newHex(
                                                nodes[i][j][k+1],
                                                nodes[i][j+1][k+1],
                                                nodes[i+1][j+1][k+1],
                                                nodes[i+1][j][k+1],
                                                nodes[i][j][k],
                                                nodes[i][j+1][k],
                                                nodes[i+1][j+1][k],
                                                nodes[i+1][j][k]
                                );
			}
		}
	}

	gmds::VTKWriter<gmds::IGMesh> writer(mesh);
        writer.write("GETMe_getme_seq_hex_before",gmds::N|gmds::R);

	gmds::GETMe getMe(mesh);

	double minScaledJacobian = getMe.computeMinNormalizedScaledJacobian();
	std::cout<<"minNormalizedScaledJacobian "<<minScaledJacobian<<std::endl;

	//getMe.exec(100,0.8,false,markFixedNodes);
	getMe.execSeq(
			 10000,
			 0.9,
			 0.01,
			 1.0,
			 0.0002,
			 0.0004,
			 0.0002,
			 false,
			 false,
			 false,
			 markFixedNodes
			 );

	writer.write("GETMe_getme_seq_hex_after",gmds::N|gmds::R);

	minScaledJacobian = getMe.computeMinNormalizedScaledJacobian();
	std::cout<<"minNormalizedScaledJacobian "<<minScaledJacobian<<std::endl;

	mesh.unmarkAll<gmds::Node>(markFixedNodes);
	mesh.freeMark<gmds::Node>(markFixedNodes);		

	EXPECT_LT(0.2,minScaledJacobian);
}
/*----------------------------------------------------------------------------*/
TEST_F(GETMeTest,getme_seq_oneTet) {

        MeshModel mod = DIM3|N|R|R2N;
        IGMesh mesh(mod);

        gmds::Node n0 = mesh.newNode( 0., 0., 0.);
        gmds::Node n1 = mesh.newNode( 1., 0., 0.);
        gmds::Node n2 = mesh.newNode( 1., 1., 0.);
        gmds::Node n3 = mesh.newNode( 1., 0., 0.5);

        int markFixedNodes = mesh.getNewMark<gmds::Node>();

        mesh.newTet(n0,n1,n2,n3);

        gmds::VTKWriter<gmds::IGMesh> writer(mesh);
        writer.write("GETMe_getme_seq_oneTet_before",gmds::N|gmds::R);

        gmds::GETMe getMe(mesh);

        double minScaledJacobian = getMe.computeMinNormalizedScaledJacobian();
        std::cout<<"minNormalizedScaledJacobian "<<minScaledJacobian<<std::endl;

        //getMe.exec(100,0.8,false,markFixedNodes);
        getMe.execSeq(
                         10000,
                         0.9,
                         0.01,
                         1.0,
                         0.0002,
                         0.0004,
                         0.0002,
                         false,
                         false,
                         false,
                         markFixedNodes
                         );

        writer.write("GETMe_getme_seq_oneTet_after",gmds::N|gmds::R);

        minScaledJacobian = getMe.computeMinNormalizedScaledJacobian();
        std::cout<<"minNormalizedScaledJacobian "<<minScaledJacobian<<std::endl;

        mesh.unmarkAll<gmds::Node>(markFixedNodes);
        mesh.freeMark<gmds::Node>(markFixedNodes);

        EXPECT_LT(0.9,minScaledJacobian);
}
/*----------------------------------------------------------------------------*/
TEST_F(GETMeTest,getme_seq_onePyr) {

        MeshModel mod = DIM3|N|R|R2N;
        IGMesh mesh(mod);

        gmds::Node n0 = mesh.newNode( 0., 0., 0.);
        gmds::Node n1 = mesh.newNode( 1., 0., 0.);
        gmds::Node n2 = mesh.newNode( 1., 1., 0.);
        gmds::Node n3 = mesh.newNode( 0., 1., 0.);
	gmds::Node n4 = mesh.newNode( 0., 1., 0.5);
        //gmds::Node n4 = mesh.newNode( 0.5, 0.5,sqrt(0.5));

        int markFixedNodes = mesh.getNewMark<gmds::Node>();

        mesh.newPyramid(n0,n1,n2,n3,n4);

        gmds::VTKWriter<gmds::IGMesh> writer(mesh);
        writer.write("GETMe_getme_seq_onePyr_before",gmds::N|gmds::R);

        gmds::GETMe getMe(mesh);

        double minScaledJacobian = getMe.computeMinNormalizedScaledJacobian();
        std::cout<<"minNormalizedScaledJacobian "<<minScaledJacobian<<std::endl;

        //getMe.exec(100,0.8,false,markFixedNodes);
        getMe.execSeq(
                         10000,
                         0.9,
                         0.01,
                         1.0,
                         0.0002,
                         0.0004,
                         0.0002,
                         false,
                         false,
                         false,
                         markFixedNodes
                         );

        writer.write("GETMe_getme_seq_onePyr_after",gmds::N|gmds::R);

        minScaledJacobian = getMe.computeMinNormalizedScaledJacobian();
        std::cout<<"minNormalizedScaledJacobian "<<minScaledJacobian<<std::endl;

        mesh.unmarkAll<gmds::Node>(markFixedNodes);
        mesh.freeMark<gmds::Node>(markFixedNodes);

        EXPECT_LT(0.9,minScaledJacobian);
}
/*----------------------------------------------------------------------------*/
TEST_F(GETMeTest,getme_seq_onePrism3) {

        MeshModel mod = DIM3|N|R|R2N;
        IGMesh mesh(mod);

        gmds::Node n0 = mesh.newNode( 0., 0., 0.);
        gmds::Node n1 = mesh.newNode( 1., 0., 0.);
        gmds::Node n2 = mesh.newNode( 1., 1., 0.);
        gmds::Node n3 = mesh.newNode( 0., 0., 1.);
        gmds::Node n4 = mesh.newNode( 1., 0., 1.);
	gmds::Node n5 = mesh.newNode( 1., 1., 0.5);

        int markFixedNodes = mesh.getNewMark<gmds::Node>();

        mesh.newPrism3(n0,n1,n2,n3,n4,n5);

        gmds::VTKWriter<gmds::IGMesh> writer(mesh);
        writer.write("GETMe_getme_seq_onePrism3_before",gmds::N|gmds::R);

        gmds::GETMe getMe(mesh);

        double minScaledJacobian = getMe.computeMinNormalizedScaledJacobian();
        std::cout<<"minNormalizedScaledJacobian "<<minScaledJacobian<<std::endl;

        //getMe.exec(100,0.8,false,markFixedNodes);
        getMe.execSeq(
                         10000,
                         0.9,
                         0.01,
                         1.0,
                         0.0002,
                         0.0004,
                         0.0002,
                         false,
                         false,
                         false,
                         markFixedNodes
                         );

        writer.write("GETMe_getme_seq_onePrism3_after",gmds::N|gmds::R);

        minScaledJacobian = getMe.computeMinNormalizedScaledJacobian();
        std::cout<<"minNormalizedScaledJacobian "<<minScaledJacobian<<std::endl;

        mesh.unmarkAll<gmds::Node>(markFixedNodes);
        mesh.freeMark<gmds::Node>(markFixedNodes);

        EXPECT_LT(0.9,minScaledJacobian);
}
/*----------------------------------------------------------------------------*/
