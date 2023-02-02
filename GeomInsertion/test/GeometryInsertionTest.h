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
 * GeometryInsertionTest.h
 *
 *  Created on: 25 mars 2015
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <GMDS/MeshInsertDetailInOut.h>
#include <GMDS/CAD/FACFacetedGeomReadAndWrite.h>
#include <GMDS/IO/VTKFacetedGeomReadAndWrite.h>
#include <GMDS/GeomMeshIntersectionService.h>
#include <GMDS/FacetedMeshIntersectionService.h>
#include <GMDS/LaplacianSmoothing.h>
/*----------------------------------------------------------------------------*/
class GeometryInsertionTest: public ::testing::Test {

  protected:
	GeometryInsertionTest(){;}
    virtual ~GeometryInsertionTest(){;}
};
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionTest,cas_test_cas0) {

	std::cout<<"=================================================="<<std::endl;
	std::cout<<"GeometryInsertionTest::cas_test_cas0"<<std::endl;

    	std::vector<gmds::geom::GeomVolume* > vols;
    	gmds::geom::FacetedGeomManager model_tmp;

    	gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    	gmds::IGMesh mesh(mod);

        gmds::geom::FacetedGeomManager model;
	gmds::VTKFacetedGeomReadAndWrite vtkReader;
	vtkReader.import(model,"/homePOYOP/travail/git_workingcopyGMDS_cmake/GeomInsertion/test/Samples/cas0.vtk");

	model.getVolumes(vols);
	gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
	vol->reorient();

	// generate initial grid
	{
	  int nx = 10;
	  int ny = 10;
	  int nz = 10;

	  gmds::Node nodes[nx+1][ny+1][nz+1];

	  double minXYZ[3], maxXYZ[3];
	  vol->computeBoundingBox(minXYZ,maxXYZ);
	  double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
	  double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
	  double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
	  double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
	  double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
	  double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);

	  double dx = (xmax - xmin) / nx;
	  double dy = (ymax - ymin) / ny;
	  double dz = (zmax - zmin) / nz;

	  for(int i=0; i<=nx; i++) {
    		for(int j=0; j<=ny; j++) {
    			for(int k=0; k<=nz; k++) {
    				nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
    			}
    		}
	  }

	  for(int i=0; i<nx; i++) {
    		for(int j=0; j<ny; j++) {
    			for(int k=0; k<nz; k++) {
    				mesh.newHex(
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

	}

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::geom::GeomMeshIntersectionService*  service =
                                new gmds::geom::FacetedMeshIntersectionService();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK("cas_test_cas0.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    insertion.project();
    insertion.exportMeshVTK("cas_test_cas0_project.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
    	smooth->exec();
    }
    insertion.exportMeshVTK("cas_test_cas0_smooth.mli",gmds::N|gmds::E|gmds::F|gmds::R);	

    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionTest,cas_test_cas1) {

	std::cout<<"=================================================="<<std::endl;
	std::cout<<"GeometryInsertionTest::cas_test_cas1"<<std::endl;

    	std::vector<gmds::geom::GeomVolume* > vols;
    	gmds::geom::FacetedGeomManager model_tmp;

    	gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    	gmds::IGMesh mesh(mod);

        gmds::geom::FacetedGeomManager model;
	gmds::VTKFacetedGeomReadAndWrite vtkReader;
	vtkReader.import(model,"/homePOYOP/travail/git_workingcopyGMDS_cmake/GeomInsertion/test/Samples/cas1.vtk");

	model.getVolumes(vols);
	gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
	vol->reorient();

	// generate initial grid
	{
	  int nx = 10;
	  int ny = 10;
	  int nz = 10;

	  gmds::Node nodes[nx+1][ny+1][nz+1];

	  double minXYZ[3], maxXYZ[3];
	  vol->computeBoundingBox(minXYZ,maxXYZ);
	  double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
	  double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
	  double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
	  double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
	  double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
	  double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);

	  double dx = (xmax - xmin) / nx;
	  double dy = (ymax - ymin) / ny;
	  double dz = (zmax - zmin) / nz;

	  for(int i=0; i<=nx; i++) {
    		for(int j=0; j<=ny; j++) {
    			for(int k=0; k<=nz; k++) {
    				nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
    			}
    		}
	  }

	  for(int i=0; i<nx; i++) {
    		for(int j=0; j<ny; j++) {
    			for(int k=0; k<nz; k++) {
    				mesh.newHex(
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

	}

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::geom::GeomMeshIntersectionService*  service =
                                new gmds::geom::FacetedMeshIntersectionService();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK("cas_test_cas1.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    insertion.project();
    insertion.exportMeshVTK("cas_test_cas1_project.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
    	smooth->exec();
    }
    insertion.exportMeshVTK("cas_test_cas1_smooth.mli",gmds::N|gmds::E|gmds::F|gmds::R);	

    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionTest,DISABLED_cas_test_cas2) {

	std::cout<<"=================================================="<<std::endl;
	std::cout<<"GeometryInsertionTest::cas_test_cas2"<<std::endl;

    	std::vector<gmds::geom::GeomVolume* > vols;
    	gmds::geom::FacetedGeomManager model_tmp;

    	gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    	gmds::IGMesh mesh(mod);

        gmds::geom::FacetedGeomManager model;
	gmds::VTKFacetedGeomReadAndWrite vtkReader;
	vtkReader.import(model,"/homePOYOP/travail/git_workingcopyGMDS_cmake/GeomInsertion/test/Samples/cas2.vtk");

	model.getVolumes(vols);
	gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
	vol->reorient();

	// generate initial grid
	{
	  int nx = 10;
	  int ny = 10;
	  int nz = 10;

	  gmds::Node nodes[nx+1][ny+1][nz+1];

	  double minXYZ[3], maxXYZ[3];
	  vol->computeBoundingBox(minXYZ,maxXYZ);
	  double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
	  double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
	  double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
	  double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
	  double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
	  double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);

	  double dx = (xmax - xmin) / nx;
	  double dy = (ymax - ymin) / ny;
	  double dz = (zmax - zmin) / nz;

	  for(int i=0; i<=nx; i++) {
    		for(int j=0; j<=ny; j++) {
    			for(int k=0; k<=nz; k++) {
    				nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
    			}
    		}
	  }

	  for(int i=0; i<nx; i++) {
    		for(int j=0; j<ny; j++) {
    			for(int k=0; k<nz; k++) {
    				mesh.newHex(
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

	}

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::geom::GeomMeshIntersectionService*  service =
                                new gmds::geom::FacetedMeshIntersectionService();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK("cas_test_cas2.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    insertion.project();
    insertion.exportMeshVTK("cas_test_cas2_project.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
    	smooth->exec();
    }
    insertion.exportMeshVTK("cas_test_cas2_smooth.mli",gmds::N|gmds::E|gmds::F|gmds::R);	

    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionTest,cas_test_cas3) {

	std::cout<<"=================================================="<<std::endl;
	std::cout<<"GeometryInsertionTest::cas_test_cas3"<<std::endl;

    	std::vector<gmds::geom::GeomVolume* > vols;
    	gmds::geom::FacetedGeomManager model_tmp;

    	gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    	gmds::IGMesh mesh(mod);

        gmds::geom::FacetedGeomManager model;
	gmds::VTKFacetedGeomReadAndWrite vtkReader;
	vtkReader.import(model,"/homePOYOP/travail/git_workingcopyGMDS_cmake/GeomInsertion/test/Samples/cas3.vtk");

	model.getVolumes(vols);
	gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
	vol->reorient();

	// generate initial grid
	{
	  int nx = 10;
	  int ny = 10;
	  int nz = 10;

	  gmds::Node nodes[nx+1][ny+1][nz+1];

	  double minXYZ[3], maxXYZ[3];
	  vol->computeBoundingBox(minXYZ,maxXYZ);
	  double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
	  double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
	  double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
	  double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
	  double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
	  double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);

	  double dx = (xmax - xmin) / nx;
	  double dy = (ymax - ymin) / ny;
	  double dz = (zmax - zmin) / nz;

	  for(int i=0; i<=nx; i++) {
    		for(int j=0; j<=ny; j++) {
    			for(int k=0; k<=nz; k++) {
    				nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
    			}
    		}
	  }

	  for(int i=0; i<nx; i++) {
    		for(int j=0; j<ny; j++) {
    			for(int k=0; k<nz; k++) {
    				mesh.newHex(
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

	}

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::geom::GeomMeshIntersectionService*  service =
                                new gmds::geom::FacetedMeshIntersectionService();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK("cas_test_cas3.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    insertion.project();
    insertion.exportMeshVTK("cas_test_cas3_project.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
    	smooth->exec();
    }
    insertion.exportMeshVTK("cas_test_cas3_smooth.mli",gmds::N|gmds::E|gmds::F|gmds::R);	

    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionTest,DISABLED_cas_test_cas4) {

	std::cout<<"=================================================="<<std::endl;
	std::cout<<"GeometryInsertionTest::cas_test_cas4"<<std::endl;

    	std::vector<gmds::geom::GeomVolume* > vols;
    	gmds::geom::FacetedGeomManager model_tmp;

    	gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    	gmds::IGMesh mesh(mod);

        gmds::geom::FacetedGeomManager model;
	gmds::VTKFacetedGeomReadAndWrite vtkReader;
	vtkReader.import(model,"/homePOYOP/travail/git_workingcopyGMDS_cmake/GeomInsertion/test/Samples/cas4.vtk");

	model.getVolumes(vols);
	gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
	vol->reorient();

	// generate initial grid
	{
	  int nx = 10;
	  int ny = 10;
	  int nz = 10;

	  gmds::Node nodes[nx+1][ny+1][nz+1];

	  double minXYZ[3], maxXYZ[3];
	  vol->computeBoundingBox(minXYZ,maxXYZ);
	  double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
	  double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
	  double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
	  double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
	  double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
	  double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);

	  double dx = (xmax - xmin) / nx;
	  double dy = (ymax - ymin) / ny;
	  double dz = (zmax - zmin) / nz;

	  for(int i=0; i<=nx; i++) {
    		for(int j=0; j<=ny; j++) {
    			for(int k=0; k<=nz; k++) {
    				nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
    			}
    		}
	  }

	  for(int i=0; i<nx; i++) {
    		for(int j=0; j<ny; j++) {
    			for(int k=0; k<nz; k++) {
    				mesh.newHex(
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

	}

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::geom::GeomMeshIntersectionService*  service =
                                new gmds::geom::FacetedMeshIntersectionService();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK("cas_test_cas4.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    insertion.project();
    insertion.exportMeshVTK("cas_test_cas4_project.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
    	smooth->exec();
    }
    insertion.exportMeshVTK("cas_test_cas4_smooth.mli",gmds::N|gmds::E|gmds::F|gmds::R);	

    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionTest,DISABLED_cas_test_cas5) {

	std::cout<<"=================================================="<<std::endl;
	std::cout<<"GeometryInsertionTest::cas_test_cas5"<<std::endl;

    	std::vector<gmds::geom::GeomVolume* > vols;
    	gmds::geom::FacetedGeomManager model_tmp;

    	gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    	gmds::IGMesh mesh(mod);

        gmds::geom::FacetedGeomManager model;
	gmds::VTKFacetedGeomReadAndWrite vtkReader;
	vtkReader.import(model,"/homePOYOP/travail/git_workingcopyGMDS_cmake/GeomInsertion/test/Samples/cas5.vtk");

	model.getVolumes(vols);
	gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
	vol->reorient();

	// generate initial grid
	{
	  int nx = 10;
	  int ny = 10;
	  int nz = 10;

	  gmds::Node nodes[nx+1][ny+1][nz+1];

	  double minXYZ[3], maxXYZ[3];
	  vol->computeBoundingBox(minXYZ,maxXYZ);
	  double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
	  double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
	  double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
	  double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
	  double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
	  double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);

	  double dx = (xmax - xmin) / nx;
	  double dy = (ymax - ymin) / ny;
	  double dz = (zmax - zmin) / nz;

	  for(int i=0; i<=nx; i++) {
    		for(int j=0; j<=ny; j++) {
    			for(int k=0; k<=nz; k++) {
    				nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
    			}
    		}
	  }

	  for(int i=0; i<nx; i++) {
    		for(int j=0; j<ny; j++) {
    			for(int k=0; k<nz; k++) {
    				mesh.newHex(
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

	}

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::geom::GeomMeshIntersectionService*  service =
                                new gmds::geom::FacetedMeshIntersectionService();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK("cas_test_cas5.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    insertion.project();
    insertion.exportMeshVTK("cas_test_cas5_project.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
    	smooth->exec();
    }
    insertion.exportMeshVTK("cas_test_cas5_smooth.mli",gmds::N|gmds::E|gmds::F|gmds::R);	

    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionTest,cas_test_cas6) {

	std::cout<<"=================================================="<<std::endl;
	std::cout<<"GeometryInsertionTest::cas_test_cas6"<<std::endl;

    	std::vector<gmds::geom::GeomVolume* > vols;
    	gmds::geom::FacetedGeomManager model_tmp;

    	gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    	gmds::IGMesh mesh(mod);

        gmds::geom::FacetedGeomManager model;
	gmds::VTKFacetedGeomReadAndWrite vtkReader;
	vtkReader.import(model,"/homePOYOP/travail/git_workingcopyGMDS_cmake/GeomInsertion/test/Samples/cas6.vtk");

	model.getVolumes(vols);
	gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
	vol->reorient();

	// generate initial grid
	{
	  int nx = 10;
	  int ny = 10;
	  int nz = 10;

	  gmds::Node nodes[nx+1][ny+1][nz+1];

	  double minXYZ[3], maxXYZ[3];
	  vol->computeBoundingBox(minXYZ,maxXYZ);
	  double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
	  double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
	  double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
	  double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
	  double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
	  double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);

	  double dx = (xmax - xmin) / nx;
	  double dy = (ymax - ymin) / ny;
	  double dz = (zmax - zmin) / nz;

	  for(int i=0; i<=nx; i++) {
    		for(int j=0; j<=ny; j++) {
    			for(int k=0; k<=nz; k++) {
    				nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
    			}
    		}
	  }

	  for(int i=0; i<nx; i++) {
    		for(int j=0; j<ny; j++) {
    			for(int k=0; k<nz; k++) {
    				mesh.newHex(
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

	}

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::geom::GeomMeshIntersectionService*  service =
                                new gmds::geom::FacetedMeshIntersectionService();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK("cas_test_cas6.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    insertion.project();
    insertion.exportMeshVTK("cas_test_cas6_project.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
    	smooth->exec();
    }
    insertion.exportMeshVTK("cas_test_cas6_smooth.mli",gmds::N|gmds::E|gmds::F|gmds::R);	

    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionTest,DISABLED_cas_test_cas7) {

	std::cout<<"=================================================="<<std::endl;
	std::cout<<"GeometryInsertionTest::cas_test_cas7"<<std::endl;

    	std::vector<gmds::geom::GeomVolume* > vols;
    	gmds::geom::FacetedGeomManager model_tmp;

    	gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    	gmds::IGMesh mesh(mod);

        gmds::geom::FacetedGeomManager model;
	gmds::VTKFacetedGeomReadAndWrite vtkReader;
	vtkReader.import(model,"/homePOYOP/travail/git_workingcopyGMDS_cmake/GeomInsertion/test/Samples/cas7.vtk");

	model.getVolumes(vols);
	gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
	vol->reorient();

	// generate initial grid
	{
	  int nx = 10;
	  int ny = 10;
	  int nz = 10;

	  gmds::Node nodes[nx+1][ny+1][nz+1];

	  double minXYZ[3], maxXYZ[3];
	  vol->computeBoundingBox(minXYZ,maxXYZ);
	  double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
	  double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
	  double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
	  double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
	  double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
	  double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);

	  double dx = (xmax - xmin) / nx;
	  double dy = (ymax - ymin) / ny;
	  double dz = (zmax - zmin) / nz;

	  for(int i=0; i<=nx; i++) {
    		for(int j=0; j<=ny; j++) {
    			for(int k=0; k<=nz; k++) {
    				nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
    			}
    		}
	  }

	  for(int i=0; i<nx; i++) {
    		for(int j=0; j<ny; j++) {
    			for(int k=0; k<nz; k++) {
    				mesh.newHex(
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

	}

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::geom::GeomMeshIntersectionService*  service =
                                new gmds::geom::FacetedMeshIntersectionService();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK("cas_test_cas7.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    insertion.project();
    insertion.exportMeshVTK("cas_test_cas7_project.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
    	smooth->exec();
    }
    insertion.exportMeshVTK("cas_test_cas7_smooth.mli",gmds::N|gmds::E|gmds::F|gmds::R);	

    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionTest,cas_test_cas8) {

	std::cout<<"=================================================="<<std::endl;
	std::cout<<"GeometryInsertionTest::cas_test_cas8"<<std::endl;

    	std::vector<gmds::geom::GeomVolume* > vols;
    	gmds::geom::FacetedGeomManager model_tmp;

    	gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    	gmds::IGMesh mesh(mod);

        gmds::geom::FacetedGeomManager model;
	gmds::VTKFacetedGeomReadAndWrite vtkReader;
	vtkReader.import(model,"/homePOYOP/travail/git_workingcopyGMDS_cmake/GeomInsertion/test/Samples/cas8.vtk");

	model.getVolumes(vols);
	gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
	vol->reorient();

	// generate initial grid
	{
	  int nx = 10;
	  int ny = 10;
	  int nz = 10;

	  gmds::Node nodes[nx+1][ny+1][nz+1];

	  double minXYZ[3], maxXYZ[3];
	  vol->computeBoundingBox(minXYZ,maxXYZ);
	  double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
	  double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
	  double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
	  double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
	  double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
	  double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);

	  double dx = (xmax - xmin) / nx;
	  double dy = (ymax - ymin) / ny;
	  double dz = (zmax - zmin) / nz;

	  for(int i=0; i<=nx; i++) {
    		for(int j=0; j<=ny; j++) {
    			for(int k=0; k<=nz; k++) {
    				nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
    			}
    		}
	  }

	  for(int i=0; i<nx; i++) {
    		for(int j=0; j<ny; j++) {
    			for(int k=0; k<nz; k++) {
    				mesh.newHex(
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

	}

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::geom::GeomMeshIntersectionService*  service =
                                new gmds::geom::FacetedMeshIntersectionService();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK("cas_test_cas8.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    insertion.project();
    insertion.exportMeshVTK("cas_test_cas8_project.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
    	smooth->exec();
    }
    insertion.exportMeshVTK("cas_test_cas8_smooth.mli",gmds::N|gmds::E|gmds::F|gmds::R);	

    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionTest,cas_test_cas10) {

	std::cout<<"=================================================="<<std::endl;
	std::cout<<"GeometryInsertionTest::cas_test_cas10"<<std::endl;

    	std::vector<gmds::geom::GeomVolume* > vols;
    	gmds::geom::FacetedGeomManager model_tmp;

    	gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    	gmds::IGMesh mesh(mod);

        gmds::geom::FacetedGeomManager model;
	gmds::VTKFacetedGeomReadAndWrite vtkReader;
	vtkReader.import(model,"/homePOYOP/travail/git_workingcopyGMDS_cmake/GeomInsertion/test/Samples/cas10.vtk");

	model.getVolumes(vols);
	gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
	vol->reorient();

	// generate initial grid
	{
	  int nx = 10;
	  int ny = 10;
	  int nz = 10;

	  gmds::Node nodes[nx+1][ny+1][nz+1];

	  double minXYZ[3], maxXYZ[3];
	  vol->computeBoundingBox(minXYZ,maxXYZ);
	  double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
	  double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
	  double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
	  double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
	  double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
	  double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);

	  double dx = (xmax - xmin) / nx;
	  double dy = (ymax - ymin) / ny;
	  double dz = (zmax - zmin) / nz;

	  for(int i=0; i<=nx; i++) {
    		for(int j=0; j<=ny; j++) {
    			for(int k=0; k<=nz; k++) {
    				nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
    			}
    		}
	  }

	  for(int i=0; i<nx; i++) {
    		for(int j=0; j<ny; j++) {
    			for(int k=0; k<nz; k++) {
    				mesh.newHex(
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

	}

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::geom::GeomMeshIntersectionService*  service =
                                new gmds::geom::FacetedMeshIntersectionService();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK("cas_test_cas10.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    insertion.project();
    insertion.exportMeshVTK("cas_test_cas10_project.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
    	smooth->exec();
    }
    insertion.exportMeshVTK("cas_test_cas10_smooth.mli",gmds::N|gmds::E|gmds::F|gmds::R);	

    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionTest,cas_test_cas11) {

	std::cout<<"=================================================="<<std::endl;
	std::cout<<"GeometryInsertionTest::cas_test_cas11"<<std::endl;

    	std::vector<gmds::geom::GeomVolume* > vols;
    	gmds::geom::FacetedGeomManager model_tmp;

    	gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    	gmds::IGMesh mesh(mod);

        gmds::geom::FacetedGeomManager model;
	gmds::VTKFacetedGeomReadAndWrite vtkReader;
	vtkReader.import(model,"/homePOYOP/travail/git_workingcopyGMDS_cmake/GeomInsertion/test/Samples/cas11.vtk");

	model.getVolumes(vols);
	gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
	vol->reorient();

	// generate initial grid
	{
	  int nx = 10;
	  int ny = 10;
	  int nz = 10;

	  gmds::Node nodes[nx+1][ny+1][nz+1];

	  double minXYZ[3], maxXYZ[3];
	  vol->computeBoundingBox(minXYZ,maxXYZ);
	  double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
	  double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
	  double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
	  double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
	  double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
	  double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);

	  double dx = (xmax - xmin) / nx;
	  double dy = (ymax - ymin) / ny;
	  double dz = (zmax - zmin) / nz;

	  for(int i=0; i<=nx; i++) {
    		for(int j=0; j<=ny; j++) {
    			for(int k=0; k<=nz; k++) {
    				nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
    			}
    		}
	  }

	  for(int i=0; i<nx; i++) {
    		for(int j=0; j<ny; j++) {
    			for(int k=0; k<nz; k++) {
    				mesh.newHex(
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

	}

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::geom::GeomMeshIntersectionService*  service =
                                new gmds::geom::FacetedMeshIntersectionService();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK("cas_test_cas11.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    insertion.project();
    insertion.exportMeshVTK("cas_test_cas11_project.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
    	smooth->exec();
    }
    insertion.exportMeshVTK("cas_test_cas11_smooth.mli",gmds::N|gmds::E|gmds::F|gmds::R);	

    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionTest,cas_test_cas12) {

	std::cout<<"=================================================="<<std::endl;
	std::cout<<"GeometryInsertionTest::cas_test_cas12"<<std::endl;

    	std::vector<gmds::geom::GeomVolume* > vols;
    	gmds::geom::FacetedGeomManager model_tmp;

    	gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    	gmds::IGMesh mesh(mod);

        gmds::geom::FacetedGeomManager model;
	gmds::VTKFacetedGeomReadAndWrite vtkReader;
	vtkReader.import(model,"/homePOYOP/travail/git_workingcopyGMDS_cmake/GeomInsertion/test/Samples/cas12.vtk");

	model.getVolumes(vols);
	gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
	vol->reorient();

	// generate initial grid
	{
	  int nx = 10;
	  int ny = 10;
	  int nz = 10;

	  gmds::Node nodes[nx+1][ny+1][nz+1];

	  double minXYZ[3], maxXYZ[3];
	  vol->computeBoundingBox(minXYZ,maxXYZ);
	  double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
	  double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
	  double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
	  double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
	  double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
	  double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);

	  double dx = (xmax - xmin) / nx;
	  double dy = (ymax - ymin) / ny;
	  double dz = (zmax - zmin) / nz;

	  for(int i=0; i<=nx; i++) {
    		for(int j=0; j<=ny; j++) {
    			for(int k=0; k<=nz; k++) {
    				nodes[i][j][k] = mesh.newNode(xmin+i*dx,ymin+j*dy,zmin+k*dz);
    			}
    		}
	  }

	  for(int i=0; i<nx; i++) {
    		for(int j=0; j<ny; j++) {
    			for(int k=0; k<nz; k++) {
    				mesh.newHex(
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

	}

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::geom::GeomMeshIntersectionService*  service =
                                new gmds::geom::FacetedMeshIntersectionService();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK("cas_test_cas12.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    insertion.project();
    insertion.exportMeshVTK("cas_test_cas12_project.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
    	smooth->exec();
    }
    insertion.exportMeshVTK("cas_test_cas12_smooth.mli",gmds::N|gmds::E|gmds::F|gmds::R);	

    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
