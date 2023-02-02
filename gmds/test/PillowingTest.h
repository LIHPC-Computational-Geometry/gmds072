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
 * PillowingTest.h
 *
 *  Created on: 1 aug. 2016
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <GMDS/Algo/Pillowing.h>
#include <GMDS/IG/IGMeshDoctor.h>
#include <GMDS/IO/VTKWriter.h>
/*----------------------------------------------------------------------------*/
class PillowingTest: public ::testing::Test {
protected:
	PillowingTest(){;}
	virtual ~PillowingTest(){;}
};

using namespace gmds;

/*----------------------------------------------------------------------------*/
TEST_F(PillowingTest,pillow) {
	
	/*
                region OUT : [ ]
                region IN  : [X]
        
                [X] 
                
        */


	MeshModel mod = DIM3|R|F|N|R2F|R2N|F2N|N2R|F2R;
	IGMesh m(mod);

	m.initializeGeometryClassification();

	Node n1 = m.newNode(0,0,0);
	Node n2 = m.newNode(1,0,0);
	Node n3 = m.newNode(1,1,0);
	Node n4 = m.newNode(0,1,0);
	Node n5 = m.newNode(0,0,1);
	Node n6 = m.newNode(1,0,1);
	Node n7 = m.newNode(1,1,1);
	Node n8 = m.newNode(0,1,1);

	Region r = m.newHex(n1,n2,n3,n4,n5,n6,n7,n8);

	IGMeshDoctor doc(&m);
	doc.buildFacesAndR2F();
	doc.updateUpwardConnectivity();

	Pillowing op(m);

	std::vector<Region> toPillow, newRegions;
	toPillow.push_back(r);
	std::vector<Face> toPillowAlong = r.get<Face>();

	op.pillow(toPillowAlong, toPillow, newRegions, true);

	gmds::VTKWriter<gmds::IGMesh> writer(m);
        writer.write("PillowingTest_pillow",gmds::N|gmds::F|gmds::R);

	EXPECT_EQ(6,newRegions.size());
}
/*----------------------------------------------------------------------------*/
TEST_F(PillowingTest,pillow_plane) {

	/*
		region OUT : [ ]
		region IN  : [X]
	
		[ ][ ][ ]
		[X][X][X]
		[X][X][X] 
		
	*/

        MeshModel mod = DIM3|R|F|N|R2F|R2N|F2N|N2R|F2R;
        IGMesh mesh(mod);

	mesh.initializeGeometryClassification();

	const int nx = 3;
        const int ny = 1;
        const int nz = 3;

        gmds::Node nodes[nx+1][ny+1][nz+1];

	const double xmin = 0.;
	const double ymin = 0.;
	const double zmin = 0.; 

	const double dx = 1.;
        const double dy = 1.;
        const double dz = 1.;

	for(int i=0; i<=nx; i++) {
                for(int j=0; j<=ny; j++) {
                        for(int k=0; k<=nz; k++) {
                                nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
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

				nodes[i][j][k+1].add(r);
				nodes[i][j+1][k+1].add(r);
                              	nodes[i+1][j+1][k+1].add(r);
                                nodes[i+1][j][k+1].add(r);
                                nodes[i][j][k].add(r);
                                nodes[i][j+1][k].add(r);
                                nodes[i+1][j+1][k].add(r);
                                nodes[i+1][j][k].add(r);

				regions[i][j][k] = r;
                        }
                }
        }
	
	// faces creation
	// Create only the faces on plane k=2
	const int k = 2;

	std::vector<gmds::Face> toPillowAlong;
	
	for(int i=0; i<nx; i++) {
		for(int j=0; j<ny; j++) {
			gmds::Face f = mesh.newQuad(
					nodes[i][j][k],
					nodes[i+1][j][k],
					nodes[i+1][j+1][k],
					nodes[i][j+1][k]
			);

			regions[i][j][k-1].add(f);
			regions[i][j][k].add(f);
			f.add(regions[i][j][k-1]);
			f.add(regions[i][j][k]);

			toPillowAlong.push_back(f);
		}
	}	

	// update connectivities
	//IGMeshDoctor doc(&mesh);
	//doc.buildR2F(mod);
	//doc.buildN2R(mod);
	//doc.buildF2R(mod);
	
	Pillowing op(mesh);

        std::vector<Region> toPillow;

	for(int i=0; i<nx; i++) {
		for(int j=0; j<ny; j++) {
		  for(int k=0; k<2; k++) {
			toPillow.push_back(regions[i][j][k]);
		  }
		}
	}

	std::vector<Region> newRegions;

        op.pillow(toPillowAlong, toPillow, newRegions, true);
	
	gmds::VTKWriter<gmds::IGMesh> writer(mesh);
        writer.write("PillowingTest_pillow_plane",gmds::N|gmds::F|gmds::R);

	EXPECT_EQ(22,newRegions.size());
}
/*----------------------------------------------------------------------------*/
TEST_F(PillowingTest,pillow_corner1) {

	/*
                region OUT : [ ]
                region IN  : [X]
        
                [ ][ ][ ]
                [ ][X][X]
                [ ][X][X] 
                
        */

        MeshModel mod = DIM3|R|F|N|R2F|R2N|F2N|N2R|F2R;
        IGMesh mesh(mod);

	mesh.initializeGeometryClassification();

	const int nx = 3;
        const int ny = 1;
        const int nz = 3;

        gmds::Node nodes[nx+1][ny+1][nz+1];

	const double xmin = 0.;
	const double ymin = 0.;
	const double zmin = 0.; 

	const double dx = 1.;
        const double dy = 1.;
        const double dz = 1.;

	for(int i=0; i<=nx; i++) {
                for(int j=0; j<=ny; j++) {
                        for(int k=0; k<=nz; k++) {
                                nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
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

				nodes[i][j][k+1].add(r);
				nodes[i][j+1][k+1].add(r);
                              	nodes[i+1][j+1][k+1].add(r);
                                nodes[i+1][j][k+1].add(r);
                                nodes[i][j][k].add(r);
                                nodes[i][j+1][k].add(r);
                                nodes[i+1][j+1][k].add(r);
                                nodes[i+1][j][k].add(r);

				regions[i][j][k] = r;
                        }
                }
        }
	
	// faces creation
	std::vector<gmds::Face> toPillowAlong;
	{
		gmds::Face f = mesh.newQuad(
					nodes[1][1][0],
					nodes[1][0][0],
					nodes[1][0][1],
					nodes[1][1][1]
		);
	
		regions[0][0][0].add(f);
		regions[1][0][0].add(f);
		f.add(regions[0][0][0]);
		f.add(regions[1][0][0]);

		toPillowAlong.push_back(f);
	}
	{
                gmds::Face f = mesh.newQuad(
                                        nodes[1][1][1],
                                        nodes[1][0][1],
                                        nodes[1][0][2],
                                        nodes[1][1][2]
                );

                regions[0][0][1].add(f);
                regions[1][0][1].add(f);
                f.add(regions[0][0][1]);
                f.add(regions[1][0][1]);

                toPillowAlong.push_back(f);
        }
	{
                gmds::Face f = mesh.newQuad(
                                        nodes[1][0][2],
                                        nodes[2][0][2],
                                        nodes[2][1][2],
                                        nodes[1][1][2]
                );

                regions[1][0][1].add(f);
                regions[1][0][2].add(f);
                f.add(regions[1][0][1]);
                f.add(regions[1][0][2]);

                toPillowAlong.push_back(f);
        }
	{
                gmds::Face f = mesh.newQuad(
                                        nodes[2][0][2],
                                        nodes[3][0][2],
                                        nodes[3][1][2],
                                        nodes[2][1][2]
                );
        
                regions[2][0][1].add(f);
                regions[2][0][2].add(f);
                f.add(regions[2][0][1]);
                f.add(regions[2][0][2]);
                
                toPillowAlong.push_back(f);
        }

	// update connectivities
	//IGMeshDoctor doc(&mesh);
	//doc.buildR2F(mod);
	//doc.buildN2R(mod);
	//doc.buildF2R(mod);
	
	Pillowing op(mesh);

        std::vector<Region> toPillow;
	toPillow.push_back(regions[1][0][0]);
	toPillow.push_back(regions[2][0][0]);
	toPillow.push_back(regions[1][0][1]);
	toPillow.push_back(regions[2][0][1]);

	std::vector<Region> newRegions;

        op.pillow(toPillowAlong, toPillow,newRegions,true);

	gmds::VTKWriter<gmds::IGMesh> writer(mesh);
        writer.write("PillowingTest_pillow_corner1",gmds::N|gmds::F|gmds::R);
	
	EXPECT_EQ(16,newRegions.size());
}
/*----------------------------------------------------------------------------*/
TEST_F(PillowingTest,pillow_corner2) {

	/*
                region OUT : [ ]
                region IN  : [X]
        
                [ ][ ][X]
                [ ][ ][X]
                [X][X][X] 
                
        */

        MeshModel mod = DIM3|R|F|N|R2F|R2N|F2N|N2R|F2R;
        IGMesh mesh(mod);

	mesh.initializeGeometryClassification();

	const int nx = 3;
        const int ny = 1;
        const int nz = 3;

        gmds::Node nodes[nx+1][ny+1][nz+1];

	const double xmin = 0.;
	const double ymin = 0.;
	const double zmin = 0.; 

	const double dx = 1.;
        const double dy = 1.;
        const double dz = 1.;

	for(int i=0; i<=nx; i++) {
                for(int j=0; j<=ny; j++) {
                        for(int k=0; k<=nz; k++) {
                                nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
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

				nodes[i][j][k+1].add(r);
				nodes[i][j+1][k+1].add(r);
                              	nodes[i+1][j+1][k+1].add(r);
                                nodes[i+1][j][k+1].add(r);
                                nodes[i][j][k].add(r);
                                nodes[i][j+1][k].add(r);
                                nodes[i+1][j+1][k].add(r);
                                nodes[i+1][j][k].add(r);

				regions[i][j][k] = r;
                        }
                }
        }
	
	// faces creation
	std::vector<gmds::Face> toPillowAlong;
	{
		gmds::Face f = mesh.newQuad(
					nodes[0][0][1],
					nodes[1][0][1],
					nodes[1][1][1],
					nodes[0][1][1]
		);
	
		regions[0][0][0].add(f);
		regions[0][0][1].add(f);
		f.add(regions[0][0][0]);
		f.add(regions[0][0][1]);

		toPillowAlong.push_back(f);
	}
	{
                gmds::Face f = mesh.newQuad(
                                        nodes[1][0][1],
                                        nodes[2][0][1],
                                        nodes[2][1][1],
                                        nodes[1][1][1]
                );

                regions[1][0][0].add(f);
                regions[1][0][1].add(f);
                f.add(regions[1][0][0]);
                f.add(regions[1][0][1]);

                toPillowAlong.push_back(f);
        }
	{
                gmds::Face f = mesh.newQuad(
                                        nodes[2][1][1],
                                        nodes[2][0][1],
                                        nodes[2][0][2],
                                        nodes[2][1][2]
                );

                regions[1][0][1].add(f);
                regions[2][0][1].add(f);
                f.add(regions[1][0][1]);
                f.add(regions[2][0][1]);

                toPillowAlong.push_back(f);
        }
	{
                gmds::Face f = mesh.newQuad(
                                        nodes[2][1][2],
                                        nodes[2][0][2],
                                        nodes[2][0][3],
                                        nodes[2][1][3]
                );
        
                regions[1][0][2].add(f);
                regions[2][0][2].add(f);
                f.add(regions[1][0][2]);
                f.add(regions[2][0][2]);
                
                toPillowAlong.push_back(f);
        }

	// update connectivities
	//IGMeshDoctor doc(&mesh);
	//doc.buildR2F(mod);
	//doc.buildN2R(mod);
	//doc.buildF2R(mod);
	
	Pillowing op(mesh);

        std::vector<Region> toPillow;
	toPillow.push_back(regions[0][0][0]);
	toPillow.push_back(regions[1][0][0]);
	toPillow.push_back(regions[2][0][1]);
	toPillow.push_back(regions[2][0][2]);
	toPillow.push_back(regions[2][0][0]); // must add the region in the corner

	std::vector<Region> newRegions;

        op.pillow(toPillowAlong, toPillow,newRegions,true);

	gmds::VTKWriter<gmds::IGMesh> writer(mesh);
	writer.write("PillowingTest_pillow_corner2",gmds::N|gmds::F|gmds::R);

	EXPECT_EQ(22,newRegions.size());

}
/*----------------------------------------------------------------------------*/
TEST_F(PillowingTest,pillow_disjointSets) {

	/*
                region OUT : [ ]
                region IN  : [X]
        
                [X][ ][ ]
                [X][ ][ ]
                [ ][ ][X] 
                
        */

        MeshModel mod = DIM3|R|F|N|R2F|R2N|F2N|N2R|F2R;
        IGMesh mesh(mod);

	mesh.initializeGeometryClassification();

	const int nx = 3;
        const int ny = 1;
        const int nz = 3;

        gmds::Node nodes[nx+1][ny+1][nz+1];

	const double xmin = 0.;
	const double ymin = 0.;
	const double zmin = 0.; 

	const double dx = 1.;
        const double dy = 1.;
        const double dz = 1.;

	for(int i=0; i<=nx; i++) {
                for(int j=0; j<=ny; j++) {
                        for(int k=0; k<=nz; k++) {
                                nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
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

				nodes[i][j][k+1].add(r);
				nodes[i][j+1][k+1].add(r);
                              	nodes[i+1][j+1][k+1].add(r);
                                nodes[i+1][j][k+1].add(r);
                                nodes[i][j][k].add(r);
                                nodes[i][j+1][k].add(r);
                                nodes[i+1][j+1][k].add(r);
                                nodes[i+1][j][k].add(r);

				regions[i][j][k] = r;
                        }
                }
        }
	
	// faces creation
	std::vector<gmds::Face> toPillowAlong;
	/*
	{
		gmds::Face f = mesh.newQuad(
					nodes[1][1][0],
					nodes[1][0][0],
					nodes[1][0][1],
					nodes[1][1][1]
		);
	
		regions[0][0][0].add(f);
		regions[1][0][0].add(f);
		f.add(regions[0][0][0]);
		f.add(regions[1][0][0]);

		toPillowAlong.push_back(f);
	}
	{
                gmds::Face f = mesh.newQuad(
                                        nodes[1][1][1],
                                        nodes[1][0][1],
                                        nodes[1][0][2],
                                        nodes[1][1][2]
                );

                regions[0][0][1].add(f);
                regions[1][0][1].add(f);
                f.add(regions[0][0][1]);
                f.add(regions[1][0][1]);

                toPillowAlong.push_back(f);
        }
	{
                gmds::Face f = mesh.newQuad(
                                        nodes[1][0][2],
                                        nodes[2][0][2],
                                        nodes[2][1][2],
                                        nodes[1][1][2]
                );

                regions[1][0][1].add(f);
                regions[1][0][2].add(f);
                f.add(regions[1][0][1]);
                f.add(regions[1][0][2]);

                toPillowAlong.push_back(f);
        }
	{
                gmds::Face f = mesh.newQuad(
                                        nodes[2][0][2],
                                        nodes[3][0][2],
                                        nodes[3][1][2],
                                        nodes[2][1][2]
                );
        
                regions[2][0][1].add(f);
                regions[2][0][2].add(f);
                f.add(regions[2][0][1]);
                f.add(regions[2][0][2]);
                
                toPillowAlong.push_back(f);
        }
	*/

	// update connectivities
	//IGMeshDoctor doc(&mesh);
	//doc.buildR2F(mod);
	//doc.buildN2R(mod);
	//doc.buildF2R(mod);
	
	Pillowing op(mesh);

        std::vector<Region> toPillow;
	toPillow.push_back(regions[2][0][0]);
	toPillow.push_back(regions[0][0][1]);
	toPillow.push_back(regions[0][0][2]);

	std::vector<Region> newRegions;

        op.pillow(toPillowAlong, toPillow,newRegions,true);

	gmds::VTKWriter<gmds::IGMesh> writer(mesh);
        writer.write("PillowingTest_pillow_disjointSets",gmds::N|gmds::F|gmds::R);
	
	EXPECT_EQ(16,newRegions.size());
}
/*----------------------------------------------------------------------------*/
TEST_F(PillowingTest,pillow_nonmanifold) {

	/*
                region OUT : [ ]
                region IN  : [X]
        
                [X][ ][ ]
                [ ][X][X]
                [ ][X][X] 
                
        */

        MeshModel mod = DIM3|R|F|N|R2F|R2N|F2N|N2R|F2R;
        IGMesh mesh(mod);

	mesh.initializeGeometryClassification();

	const int nx = 3;
        const int ny = 1;
        const int nz = 3;

        gmds::Node nodes[nx+1][ny+1][nz+1];

	const double xmin = 0.;
	const double ymin = 0.;
	const double zmin = 0.; 

	const double dx = 1.;
        const double dy = 1.;
        const double dz = 1.;

	for(int i=0; i<=nx; i++) {
                for(int j=0; j<=ny; j++) {
                        for(int k=0; k<=nz; k++) {
                                nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
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

				nodes[i][j][k+1].add(r);
				nodes[i][j+1][k+1].add(r);
                              	nodes[i+1][j+1][k+1].add(r);
                                nodes[i+1][j][k+1].add(r);
                                nodes[i][j][k].add(r);
                                nodes[i][j+1][k].add(r);
                                nodes[i+1][j+1][k].add(r);
                                nodes[i+1][j][k].add(r);

				regions[i][j][k] = r;
                        }
                }
        }
	
	// faces creation
	std::vector<gmds::Face> toPillowAlong;
	/*
	{
		gmds::Face f = mesh.newQuad(
					nodes[1][1][0],
					nodes[1][0][0],
					nodes[1][0][1],
					nodes[1][1][1]
		);
	
		regions[0][0][0].add(f);
		regions[1][0][0].add(f);
		f.add(regions[0][0][0]);
		f.add(regions[1][0][0]);

		toPillowAlong.push_back(f);
	}
	{
                gmds::Face f = mesh.newQuad(
                                        nodes[1][1][1],
                                        nodes[1][0][1],
                                        nodes[1][0][2],
                                        nodes[1][1][2]
                );

                regions[0][0][1].add(f);
                regions[1][0][1].add(f);
                f.add(regions[0][0][1]);
                f.add(regions[1][0][1]);

                toPillowAlong.push_back(f);
        }
	{
                gmds::Face f = mesh.newQuad(
                                        nodes[1][0][2],
                                        nodes[2][0][2],
                                        nodes[2][1][2],
                                        nodes[1][1][2]
                );

                regions[1][0][1].add(f);
                regions[1][0][2].add(f);
                f.add(regions[1][0][1]);
                f.add(regions[1][0][2]);

                toPillowAlong.push_back(f);
        }
	{
                gmds::Face f = mesh.newQuad(
                                        nodes[2][0][2],
                                        nodes[3][0][2],
                                        nodes[3][1][2],
                                        nodes[2][1][2]
                );
        
                regions[2][0][1].add(f);
                regions[2][0][2].add(f);
                f.add(regions[2][0][1]);
                f.add(regions[2][0][2]);
                
                toPillowAlong.push_back(f);
        }
	*/

	// update connectivities
	//IGMeshDoctor doc(&mesh);
	//doc.buildR2F(mod);
	//doc.buildN2R(mod);
	//doc.buildF2R(mod);
	
	Pillowing op(mesh);

        std::vector<Region> toPillow;
	toPillow.push_back(regions[1][0][0]);
	toPillow.push_back(regions[2][0][0]);
	toPillow.push_back(regions[1][0][1]);
	toPillow.push_back(regions[2][0][1]);
	toPillow.push_back(regions[0][0][2]);

	std::vector<Region> newRegions;

        op.pillow(toPillowAlong, toPillow,newRegions,true);

	gmds::VTKWriter<gmds::IGMesh> writer(mesh);
        writer.write("PillowingTest_pillow_nonmanifold",gmds::N|gmds::F|gmds::R);
	
	EXPECT_EQ(22,newRegions.size());
}
/*----------------------------------------------------------------------------*/
