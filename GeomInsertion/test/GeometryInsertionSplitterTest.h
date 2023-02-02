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
 * GeometryInsertionSplitterTest.h
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
#include <GMDS/MeshSplitter3Refinement.h>
#include <GMDS/SmartLaplacianSmoothing.h>
/*----------------------------------------------------------------------------*/
class GeometryInsertionSplitterTest: public ::testing::Test {

  protected:
	GeometryInsertionSplitterTest(){;}
    virtual ~GeometryInsertionSplitterTest(){;}
};
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionSplitterTest,cas_test_cas2) {

    std::cout<<"=================================================="<<std::endl;
    std::cout<<"GeometryInsertionSplitterTest::cas_test_cas2"<<std::endl;

    std::vector<gmds::geom::GeomVolume* > vols;
    
    gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
      		gmds::E|gmds::E2N|gmds::F2E;
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

    gmds::geom::GeomMeshIntersectionService*  service =
                                new gmds::geom::FacetedMeshIntersectionService();

    std::vector<gmds::geom::GeomSurface* > surfs;
    model.getSurfaces(surfs);
    service->initialization(surfs);
//    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);

    gmds::MeshSplitter3Refinement splitter(mesh,model,*service);
    splitter.splitMesh();
    splitter.exportVTK("splitMesh.mli");

    // change model
    gmds::MeshModel mod_bis = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    mesh.changeModel(mod_bis,false);

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
//    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.exportModelVTK("cas_test_cas2_model.mli");

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK("cas_test_cas2.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    insertion.project();
    insertion.exportMeshVTK("cas_test_cas2_project.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::SmartLaplacianSmoothing* smooth = new gmds::SmartLaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
    	smooth->exec();
    }
    insertion.exportMeshVTK("cas_test_cas2_smartsmooth.mli",gmds::N|gmds::E|gmds::F|gmds::R);

    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
