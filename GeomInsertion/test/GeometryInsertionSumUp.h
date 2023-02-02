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
 * GeometryInsertionSumUp.h
 *
 *  Created on: 14 april 2015
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
#include <GMDS/SmartDampenedLaplacianSmoothing.h>
#include <GMDS/SmartLaplacianSmoothing.h>
#include <GMDS/OrderedSmartLaplacianSmoothing.h>
#include <GMDS/NormalLaplacianSmoothing.h>
#include <GMDS/MeshModelAlgo.h>

#include <GMDS/Algo/SheetOperator.h>
#include <GMDS/IO/VTKWriter.h>
#include <GMDS/IO/VTKReader.h>

#include <GMDS/Algo/GETMe.h>

#include <stdio.h>
/*----------------------------------------------------------------------------*/
class GeometryInsertionSumUp: public ::testing::Test {

  protected:
	GeometryInsertionSumUp(){;}
    virtual ~GeometryInsertionSumUp(){;}
};
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionSumUp,DISABLED_cas_0) {

    std::cout<<"=================================================="<<std::endl;
    std::cout<<"GeometryInsertionSumUp::cas_0"<<std::endl;

    std::string case_name("cas_0");

    std::vector<gmds::geom::GeomVolume* > vols;
    
    gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
      		gmds::E|gmds::E2N|gmds::F2E;
    gmds::IGMesh mesh(mod);

    gmds::geom::FacetedGeomManager model;
    gmds::VTKFacetedGeomReadAndWrite vtkReader;
    vtkReader.import(model,"/home/legoff/travail/GMDS/gmds_workingcopy_new/GeomInsertion/test/Samples/cas0.vtk");

    model.getVolumes(vols);
    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
    vol->reorient();

    // generate initial grid
    {
    	const int nx = 20;
    	const int ny = 20;
    	const int nz = 20;

    	gmds::Node nodes[nx+1][ny+1][nz+1];

        double minXYZ[3], maxXYZ[3];
        vol->computeBoundingBox(minXYZ,maxXYZ);
/*
    	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
    	double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);
*/
	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-9);
        double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-9);

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
    service->buildAABBSurfacesTriangulationTree(surfs);
//    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);

    gmds::MeshSplitter3Refinement splitter(mesh,model,*service);
    splitter.splitMesh();
    splitter.exportVTK(case_name+std::string("splitMesh"));

    // change model
    gmds::MeshModel mod_bis = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    mesh.changeModel(mod_bis,false);

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
//    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.exportModelVTK(case_name+std::string("_model"));

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK(case_name,gmds::N|gmds::E|gmds::F|gmds::R);
    insertion.exportCurvesEdgesVTK(case_name+std::string("_curves"));
    //insertion.project();
    insertion.exportCurvesEdgesVTK(case_name+std::string("_projectCurves"));
	insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_project"),gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    int markBoundaryNodes = mesh.getNewMark<gmds::Node>();
    meshModelAlgo.markBoundaryNodes(markBoundaryNodes);    

    //gmds::OrderedSmartLaplacianSmoothing* smooth = new gmds::OrderedSmartLaplacianSmoothing(mesh,model,*service);
    //gmds::SmartLaplacianSmoothing* smooth = new gmds::SmartLaplacianSmoothing(mesh,model,*service);
    gmds::SmartDampenedLaplacianSmoothing* smooth = new gmds::SmartDampenedLaplacianSmoothing(mesh,model,*service);
    //gmds::NormalLaplacianSmoothing* smooth = new gmds::NormalLaplacianSmoothing(mesh,model,*service);
    //gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
        //smooth->exec(markBoundaryNodes);
    }
    insertion.exportMeshVTK(case_name+std::string("_smartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);

    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    gmds::MeshModel mod_AAAAAA = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
                gmds::E|gmds::E2N|gmds::F2E;
    mesh.changeModel(mod_AAAAAA,false);
	std::cout<<"mesh.getNbEdges() "<<mesh.getNbEdges()<<std::endl;
    splitter.refineQuads2EdgesOnCurveImproved();
    //gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    meshModelAlgo.associateNodes();
    mesh.changeModel(mod_bis,false);
    doc.buildFacesAndR2F();  
    doc.updateUpwardConnectivity();
 
    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    insertion.exportMeshVTK(case_name+std::string("_projectrefine"),gmds::N|gmds::E|gmds::F|gmds::R);
    mesh.changeModel(mod_bis,false);

    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
       	//smooth->exec(markBoundaryNodes);
    }
	insertion.displayMeshQuality();
    insertion.exportMeshVTK(case_name+std::string("_projectrefinesmartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);

    // pillowing
    gmds::MeshModel mod_ter = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2R;
    mesh.changeModel(mod_ter,false);

    insertion.removeUnassociatedFaces();
    insertion.pillow(vols[0]);
	std::cout<<"pillowINT"<<std::endl;
    insertion.removeUnassociatedFaces();
    insertion.pillowExt(vols[0]);
	std::cout<<"pillowEXT"<<std::endl;
	insertion.displayMeshQuality();
    insertion.exportMeshVTK(case_name+std::string("_pillow"),gmds::N|gmds::E|gmds::F|gmds::R);
	
//	insertion.project();
	insertion.displayMeshQuality();
	insertion.exportMeshVTK(case_name+std::string("_poyop"),gmds::N|gmds::E|gmds::F|gmds::R);
    smooth->initNodesAdjacencies();

    for(int i=0; i<20; i++) {
	insertion.displayMeshQuality();
        smooth->exec(markBoundaryNodes);
	char s[10];
	sprintf(s,"%d",i);
	insertion.exportMeshVTK(case_name+std::string("_pillowsmooth")+std::string(s),gmds::N|gmds::F|gmds::R);
    }

	insertion.project();

    insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_pillowsmooth"),gmds::N|gmds::F|gmds::R);

	//mesh.changeModel(mod_AAAAAA,false); 
	//meshModelAlgo.associateNodes();
	gmds::NormalLaplacianSmoothing* smooth2 = new gmds::NormalLaplacianSmoothing(mesh,model,*service);
	smooth2->initNodesAdjacencies();
	for(int i=0; i<5; i++) {
		insertion.displayMeshQuality();
		smooth2->exec(markBoundaryNodes,vol);
	}
 
    insertion.displayMeshQuality();
    insertion.exportMeshVTK(case_name+std::string("_pillowsmoothNormal"),gmds::N|gmds::F|gmds::R);

	gmds::GETMe getme(mesh);
	getme.execSimult(     
                        100,
			0.2,
                        0.5,
                        0.,
                        1.,
                        1.,
                        true,
                        true,
			false,
                        markBoundaryNodes
        );	

	insertion.displayMeshQuality();
	insertion.exportMeshVTK(case_name+std::string("_getme"),gmds::N|gmds::F|gmds::R);

	insertion.project();
	insertion.displayMeshQuality();
        insertion.exportMeshVTK(case_name+std::string("_getmeproj"),gmds::N|gmds::F|gmds::R);	
 
    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionSumUp,DISABLED_cas_1) {

    std::cout<<"=================================================="<<std::endl;
    std::cout<<"GeometryInsertionSumUp::cas_1"<<std::endl;

    std::string case_name("cas_1");

    std::vector<gmds::geom::GeomVolume* > vols;
    
    gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
      		gmds::E|gmds::E2N|gmds::F2E;
    gmds::IGMesh mesh(mod);

    gmds::geom::FacetedGeomManager model;
    gmds::VTKFacetedGeomReadAndWrite vtkReader;
    vtkReader.import(model,"/home/legoff/travail/GMDS/gmds_workingcopy_new/GeomInsertion/test/Samples/cas1.vtk");

    model.getVolumes(vols);
    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
    vol->reorient();

    // generate initial grid
    {
    	const int nx = 20;
    	const int ny = 20;
    	const int nz = 20;

    	gmds::Node nodes[nx+1][ny+1][nz+1];

        double minXYZ[3], maxXYZ[3];
        vol->computeBoundingBox(minXYZ,maxXYZ);
/*
    	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
    	double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);
*/
	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-9);
        double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-9);

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
    service->buildAABBSurfacesTriangulationTree(surfs);
//    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);

    gmds::MeshSplitter3Refinement splitter(mesh,model,*service);
    splitter.splitMesh();
    splitter.exportVTK(case_name+std::string("splitMesh"));

    // change model
    gmds::MeshModel mod_bis = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    mesh.changeModel(mod_bis,false);

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
//    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.exportModelVTK(case_name+std::string("_model"));

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK(case_name,gmds::N|gmds::E|gmds::F|gmds::R);
    insertion.exportCurvesEdgesVTK(case_name+std::string("_curves"));
    //insertion.project();
    insertion.exportCurvesEdgesVTK(case_name+std::string("_projectCurves"));
	insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_project"),gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    int markBoundaryNodes = mesh.getNewMark<gmds::Node>();
    meshModelAlgo.markBoundaryNodes(markBoundaryNodes);    

    //gmds::OrderedSmartLaplacianSmoothing* smooth = new gmds::OrderedSmartLaplacianSmoothing(mesh,model,*service);
    //gmds::SmartLaplacianSmoothing* smooth = new gmds::SmartLaplacianSmoothing(mesh,model,*service);
    gmds::SmartDampenedLaplacianSmoothing* smooth = new gmds::SmartDampenedLaplacianSmoothing(mesh,model,*service);
    //gmds::NormalLaplacianSmoothing* smooth = new gmds::NormalLaplacianSmoothing(mesh,model,*service);
    //gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
        //smooth->exec(markBoundaryNodes);
    }
    insertion.exportMeshVTK(case_name+std::string("_smartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);

    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    gmds::MeshModel mod_AAAAAA = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
                gmds::E|gmds::E2N|gmds::F2E;
    mesh.changeModel(mod_AAAAAA,false);
	std::cout<<"mesh.getNbEdges() "<<mesh.getNbEdges()<<std::endl;
    splitter.refineQuads2EdgesOnCurveImproved();
    //gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    meshModelAlgo.associateNodes();
    mesh.changeModel(mod_bis,false);
    doc.buildFacesAndR2F();  
    doc.updateUpwardConnectivity();
 
    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    insertion.exportMeshVTK(case_name+std::string("_projectrefine"),gmds::N|gmds::E|gmds::F|gmds::R);
    mesh.changeModel(mod_bis,false);

    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
       	//smooth->exec(markBoundaryNodes);
    }
	insertion.displayMeshQuality();
    insertion.exportMeshVTK(case_name+std::string("_projectrefinesmartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);

    // pillowing
    gmds::MeshModel mod_ter = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2R;
    mesh.changeModel(mod_ter,false);

    insertion.removeUnassociatedFaces();
    insertion.pillow(vols[0]);
	std::cout<<"pillowINT"<<std::endl;
    insertion.removeUnassociatedFaces();
    insertion.pillowExt(vols[0]);
	std::cout<<"pillowEXT"<<std::endl;
	insertion.displayMeshQuality();
    insertion.exportMeshVTK(case_name+std::string("_pillow"),gmds::N|gmds::E|gmds::F|gmds::R);
	
//	insertion.project();
	insertion.displayMeshQuality();
	insertion.exportMeshVTK(case_name+std::string("_poyop"),gmds::N|gmds::E|gmds::F|gmds::R);
    smooth->initNodesAdjacencies();

    for(int i=0; i<20; i++) {
	insertion.displayMeshQuality();
        smooth->exec(markBoundaryNodes);
	char s[10];
	sprintf(s,"%d",i);
	insertion.exportMeshVTK(case_name+std::string("_pillowsmooth")+std::string(s),gmds::N|gmds::F|gmds::R);
    }

	insertion.project();

    insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_pillowsmooth"),gmds::N|gmds::F|gmds::R);

/*
	//mesh.changeModel(mod_AAAAAA,false); 
	//meshModelAlgo.associateNodes();
	gmds::NormalLaplacianSmoothing* smooth2 = new gmds::NormalLaplacianSmoothing(mesh,model,*service);
	smooth2->initNodesAdjacencies();
	for(int i=0; i<5; i++) {
		insertion.displayMeshQuality();
		smooth2->exec(markBoundaryNodes,vol);
	}
 */
    insertion.displayMeshQuality();
    insertion.exportMeshVTK(case_name+std::string("_pillowsmoothNormal"),gmds::N|gmds::F|gmds::R);

	gmds::GETMe getme(mesh);
	getme.execSimult(     
                        100,
			0.2,
                        0.5,
                        0.,
                        1.,
                        1.,
                        true,
                        true,
			false,
                        markBoundaryNodes
        );	

	insertion.displayMeshQuality();
	insertion.exportMeshVTK(case_name+std::string("_getme"),gmds::N|gmds::F|gmds::R);

	insertion.project();
	insertion.displayMeshQuality();
        insertion.exportMeshVTK(case_name+std::string("_getmeproj"),gmds::N|gmds::F|gmds::R);	
 
    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionSumUp,DISABLED_cas_2) {

    std::cout<<"=================================================="<<std::endl;
    std::cout<<"GeometryInsertionSumUp::cas_2"<<std::endl;

    std::string case_name("cas_2");

    std::vector<gmds::geom::GeomVolume* > vols;
    
    gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
      		gmds::E|gmds::E2N|gmds::F2E;
    gmds::IGMesh mesh(mod);

    gmds::geom::FacetedGeomManager model;
    gmds::VTKFacetedGeomReadAndWrite vtkReader;
    vtkReader.import(model,"/home/legoff/travail/GMDS/gmds_workingcopy_new/GeomInsertion/test/Samples/cas2.vtk");

    model.getVolumes(vols);
    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
    vol->reorient();

    // generate initial grid
    {
    	const int nx = 20;
    	const int ny = 20;
    	const int nz = 20;

    	gmds::Node nodes[nx+1][ny+1][nz+1];

        double minXYZ[3], maxXYZ[3];
        vol->computeBoundingBox(minXYZ,maxXYZ);
/*
    	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
    	double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);
*/
	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-9);
        double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-9);

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
    service->buildAABBSurfacesTriangulationTree(surfs);
//    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);

    gmds::MeshSplitter3Refinement splitter(mesh,model,*service);
    splitter.splitMesh();
    splitter.exportVTK(case_name+std::string("splitMesh"));

    // change model
    gmds::MeshModel mod_bis = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    mesh.changeModel(mod_bis,false);

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
//    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.exportModelVTK(case_name+std::string("_model"));

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK(case_name,gmds::N|gmds::E|gmds::F|gmds::R);
    insertion.exportCurvesEdgesVTK(case_name+std::string("_curves"));
    //insertion.project();
    insertion.exportCurvesEdgesVTK(case_name+std::string("_projectCurves"));
	insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_project"),gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    int markBoundaryNodes = mesh.getNewMark<gmds::Node>();
    meshModelAlgo.markBoundaryNodes(markBoundaryNodes);    

    //gmds::OrderedSmartLaplacianSmoothing* smooth = new gmds::OrderedSmartLaplacianSmoothing(mesh,model,*service);
    //gmds::SmartLaplacianSmoothing* smooth = new gmds::SmartLaplacianSmoothing(mesh,model,*service);
    gmds::SmartDampenedLaplacianSmoothing* smooth = new gmds::SmartDampenedLaplacianSmoothing(mesh,model,*service);
    //gmds::NormalLaplacianSmoothing* smooth = new gmds::NormalLaplacianSmoothing(mesh,model,*service);
    //gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
        //smooth->exec(markBoundaryNodes);
    }
    insertion.exportMeshVTK(case_name+std::string("_smartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);

    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    gmds::MeshModel mod_AAAAAA = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
                gmds::E|gmds::E2N|gmds::F2E;
    mesh.changeModel(mod_AAAAAA,false);
	std::cout<<"mesh.getNbEdges() "<<mesh.getNbEdges()<<std::endl;
    splitter.refineQuads2EdgesOnCurveImproved();
    //gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    meshModelAlgo.associateNodes();
    mesh.changeModel(mod_bis,false);
    doc.buildFacesAndR2F();  
    doc.updateUpwardConnectivity();
 
    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    insertion.exportMeshVTK(case_name+std::string("_projectrefine"),gmds::N|gmds::E|gmds::F|gmds::R);
    mesh.changeModel(mod_bis,false);

    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
       	//smooth->exec(markBoundaryNodes);
    }
	insertion.displayMeshQuality();
    insertion.exportMeshVTK(case_name+std::string("_projectrefinesmartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);

    // pillowing
    gmds::MeshModel mod_ter = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2R;
    mesh.changeModel(mod_ter,false);

    insertion.removeUnassociatedFaces();
    insertion.pillow(vols[0]);
	std::cout<<"pillowINT"<<std::endl;
    insertion.removeUnassociatedFaces();
    insertion.pillowExt(vols[0]);
	std::cout<<"pillowEXT"<<std::endl;
	insertion.displayMeshQuality();
    insertion.exportMeshVTK(case_name+std::string("_pillow"),gmds::N|gmds::E|gmds::F|gmds::R);
	
//	insertion.project();
	insertion.displayMeshQuality();
	insertion.exportMeshVTK(case_name+std::string("_poyop"),gmds::N|gmds::E|gmds::F|gmds::R);
    smooth->initNodesAdjacencies();

    for(int i=0; i<20; i++) {
	insertion.displayMeshQuality();
        smooth->exec(markBoundaryNodes);
	char s[10];
	sprintf(s,"%d",i);
	insertion.exportMeshVTK(case_name+std::string("_pillowsmooth")+std::string(s),gmds::N|gmds::F|gmds::R);
    }

	insertion.project();

    insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_pillowsmooth"),gmds::N|gmds::F|gmds::R);

	//mesh.changeModel(mod_AAAAAA,false); 
	//meshModelAlgo.associateNodes();
	gmds::NormalLaplacianSmoothing* smooth2 = new gmds::NormalLaplacianSmoothing(mesh,model,*service);
	smooth2->initNodesAdjacencies();
	for(int i=0; i<5; i++) {
		insertion.displayMeshQuality();
		smooth2->exec(markBoundaryNodes,vol);
	}
 
    insertion.displayMeshQuality();
    insertion.exportMeshVTK(case_name+std::string("_pillowsmoothNormal"),gmds::N|gmds::F|gmds::R);

	gmds::GETMe getme(mesh);
	getme.execSimult(     
                        100,
			0.2,
                        0.5,
                        0.,
                        1.,
                        1.,
                        true,
                        true,
			false,
                        markBoundaryNodes
        );	

	insertion.displayMeshQuality();
	insertion.exportMeshVTK(case_name+std::string("_getme"),gmds::N|gmds::F|gmds::R);

	insertion.project();
	insertion.displayMeshQuality();
        insertion.exportMeshVTK(case_name+std::string("_getmeproj"),gmds::N|gmds::F|gmds::R);	
 
    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionSumUp,DISABLED_cas_3) {

    std::cout<<"=================================================="<<std::endl;
    std::cout<<"GeometryInsertionSumUp::cas_3"<<std::endl;

    std::string case_name("cas_3");

    std::vector<gmds::geom::GeomVolume* > vols;
    
    gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
      		gmds::E|gmds::E2N|gmds::F2E;
    gmds::IGMesh mesh(mod);

    gmds::geom::FacetedGeomManager model;
    gmds::VTKFacetedGeomReadAndWrite vtkReader;
    vtkReader.import(model,"/home/legoff/travail/GMDS/gmds_workingcopy_new/GeomInsertion/test/Samples/cas3.vtk");

    model.getVolumes(vols);
    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
    vol->reorient();

    // generate initial grid
    {
    	const int nx = 20;
    	const int ny = 20;
    	const int nz = 20;

    	gmds::Node nodes[nx+1][ny+1][nz+1];

        double minXYZ[3], maxXYZ[3];
        vol->computeBoundingBox(minXYZ,maxXYZ);
/*
    	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
    	double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);
*/
	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-9);
        double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-9);

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
    service->buildAABBSurfacesTriangulationTree(surfs);
//    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);

    gmds::MeshSplitter3Refinement splitter(mesh,model,*service);
    splitter.splitMesh();
    splitter.exportVTK(case_name+std::string("splitMesh"));

    // change model
    gmds::MeshModel mod_bis = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    mesh.changeModel(mod_bis,false);

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
//    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.exportModelVTK(case_name+std::string("_model"));

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK(case_name,gmds::N|gmds::E|gmds::F|gmds::R);
    insertion.exportCurvesEdgesVTK(case_name+std::string("_curves"));
    //insertion.project();
    insertion.exportCurvesEdgesVTK(case_name+std::string("_projectCurves"));
	insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_project"),gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    int markBoundaryNodes = mesh.getNewMark<gmds::Node>();
    meshModelAlgo.markBoundaryNodes(markBoundaryNodes);    

    //gmds::OrderedSmartLaplacianSmoothing* smooth = new gmds::OrderedSmartLaplacianSmoothing(mesh,model,*service);
    //gmds::SmartLaplacianSmoothing* smooth = new gmds::SmartLaplacianSmoothing(mesh,model,*service);
    gmds::SmartDampenedLaplacianSmoothing* smooth = new gmds::SmartDampenedLaplacianSmoothing(mesh,model,*service);
    //gmds::NormalLaplacianSmoothing* smooth = new gmds::NormalLaplacianSmoothing(mesh,model,*service);
    //gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
        //smooth->exec(markBoundaryNodes);
    }
    insertion.exportMeshVTK(case_name+std::string("_smartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);

    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    gmds::MeshModel mod_AAAAAA = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
                gmds::E|gmds::E2N|gmds::F2E;
    mesh.changeModel(mod_AAAAAA,false);
	std::cout<<"mesh.getNbEdges() "<<mesh.getNbEdges()<<std::endl;
    splitter.refineQuads2EdgesOnCurveImproved();
    //gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    meshModelAlgo.associateNodes();
    mesh.changeModel(mod_bis,false);
    doc.buildFacesAndR2F();  
    doc.updateUpwardConnectivity();
 
    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    insertion.exportMeshVTK(case_name+std::string("_projectrefine"),gmds::N|gmds::E|gmds::F|gmds::R);
    mesh.changeModel(mod_bis,false);

    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
       	//smooth->exec(markBoundaryNodes);
    }
	insertion.displayMeshQuality();
    insertion.exportMeshVTK(case_name+std::string("_projectrefinesmartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);

    // pillowing
    gmds::MeshModel mod_ter = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2R;
    mesh.changeModel(mod_ter,false);

    insertion.removeUnassociatedFaces();
    insertion.pillow(vols[0]);
	std::cout<<"pillowINT"<<std::endl;
    insertion.removeUnassociatedFaces();
    insertion.pillowExt(vols[0]);
	std::cout<<"pillowEXT"<<std::endl;
	insertion.displayMeshQuality();
    insertion.exportMeshVTK(case_name+std::string("_pillow"),gmds::N|gmds::E|gmds::F|gmds::R);
	
//	insertion.project();
	insertion.displayMeshQuality();
	insertion.exportMeshVTK(case_name+std::string("_poyop"),gmds::N|gmds::E|gmds::F|gmds::R);
    smooth->initNodesAdjacencies();

    for(int i=0; i<20; i++) {
	insertion.displayMeshQuality();
        smooth->exec(markBoundaryNodes);
	char s[10];
	sprintf(s,"%d",i);
	insertion.exportMeshVTK(case_name+std::string("_pillowsmooth")+std::string(s),gmds::N|gmds::F|gmds::R);
    }

	insertion.project();

    insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_pillowsmooth"),gmds::N|gmds::F|gmds::R);

	//mesh.changeModel(mod_AAAAAA,false); 
	//meshModelAlgo.associateNodes();
	gmds::NormalLaplacianSmoothing* smooth2 = new gmds::NormalLaplacianSmoothing(mesh,model,*service);
	smooth2->initNodesAdjacencies();
	for(int i=0; i<5; i++) {
		insertion.displayMeshQuality();
		smooth2->exec(markBoundaryNodes,vol);
	}
 
    insertion.displayMeshQuality();
    insertion.exportMeshVTK(case_name+std::string("_pillowsmoothNormal"),gmds::N|gmds::F|gmds::R);

	gmds::GETMe getme(mesh);
	getme.execSimult(     
                        100,
			0.2,
                        0.5,
                        0.,
                        1.,
                        1.,
                        true,
                        true,
			false,
                        markBoundaryNodes
        );	

	insertion.displayMeshQuality();
	insertion.exportMeshVTK(case_name+std::string("_getme"),gmds::N|gmds::F|gmds::R);

	insertion.project();
	insertion.displayMeshQuality();
        insertion.exportMeshVTK(case_name+std::string("_getmeproj"),gmds::N|gmds::F|gmds::R);	
 
    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(GeometryInsertionSumUp,cas_pyramid) {

    std::cout<<"=================================================="<<std::endl;
    std::cout<<"GeometryInsertionSumUp::cas_pyramid"<<std::endl;

    std::string case_name("cas_pyramid");

    std::vector<gmds::geom::GeomVolume* > vols;
    
    gmds::MeshModel mod = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
      		gmds::E|gmds::E2N|gmds::F2E;
    gmds::IGMesh mesh(mod);

    gmds::geom::FacetedGeomManager model;
    gmds::VTKFacetedGeomReadAndWrite vtkReader;
    vtkReader.import(model,"/home/nicolas/travail/GMDS/gmds_workingcopy_new/GeomInsertion/test/Samples/pyramid.vtk");

    model.getVolumes(vols);
    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
    vol->reorient();

    // generate initial grid
    {
    	const int nx = 20;
    	const int ny = 20;
    	const int nz = 20;

    	gmds::Node nodes[nx+1][ny+1][nz+1];

        double minXYZ[3], maxXYZ[3];
        vol->computeBoundingBox(minXYZ,maxXYZ);
/*
    	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-3);
    	double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-3);
    	double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-3);
    	double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-3);
*/
	double xmin = minXYZ[0] - (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymin = minXYZ[1] - (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmin = minXYZ[2] - (maxXYZ[2] - minXYZ[2])/(nz-9);
        double xmax = maxXYZ[0] + (maxXYZ[0] - minXYZ[0])/(nx-9);
        double ymax = maxXYZ[1] + (maxXYZ[1] - minXYZ[1])/(ny-9);
        double zmax = maxXYZ[2] + (maxXYZ[2] - minXYZ[2])/(nz-9);

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
    service->buildAABBSurfacesTriangulationTree(surfs);
//    gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);

    gmds::MeshSplitter3Refinement splitter(mesh,model,*service);
    splitter.splitMesh();
    splitter.exportVTK(case_name+std::string("splitMesh"));

    // change model
    gmds::MeshModel mod_bis = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2F|gmds::N2R;
    mesh.changeModel(mod_bis,false);

    gmds::IGMeshDoctor doc(&mesh);
    doc.buildFacesAndR2F();
//    doc.buildEdgesAndX2E();
    doc.updateUpwardConnectivity();

    mesh.initializeGeometryClassification();

    gmds::MeshInsertDetailInOut insertion(mesh,model,*service);

    insertion.exportModelVTK(case_name+std::string("_model"));

    insertion.autoInsert(vols[0]);

    insertion.exportMeshVTK(case_name,gmds::N|gmds::E|gmds::F|gmds::R);
    insertion.exportCurvesEdgesVTK(case_name+std::string("_curves"));
    //insertion.project();
    insertion.exportCurvesEdgesVTK(case_name+std::string("_projectCurves"));
	insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_project"),gmds::N|gmds::E|gmds::F|gmds::R);

    gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    int markBoundaryNodes = mesh.getNewMark<gmds::Node>();
    meshModelAlgo.markBoundaryNodes(markBoundaryNodes);    

    //gmds::OrderedSmartLaplacianSmoothing* smooth = new gmds::OrderedSmartLaplacianSmoothing(mesh,model,*service);
    //gmds::SmartLaplacianSmoothing* smooth = new gmds::SmartLaplacianSmoothing(mesh,model,*service);
    gmds::SmartDampenedLaplacianSmoothing* smooth = new gmds::SmartDampenedLaplacianSmoothing(mesh,model,*service);
    //gmds::NormalLaplacianSmoothing* smooth = new gmds::NormalLaplacianSmoothing(mesh,model,*service);
    //gmds::LaplacianSmoothing* smooth = new gmds::LaplacianSmoothing(mesh,model,*service);
    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
        //smooth->exec(markBoundaryNodes);
    }
    insertion.exportMeshVTK(case_name+std::string("_smartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);

    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    gmds::MeshModel mod_AAAAAA = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|
                gmds::E|gmds::E2N|gmds::F2E;
    mesh.changeModel(mod_AAAAAA,false);
	std::cout<<"mesh.getNbEdges() "<<mesh.getNbEdges()<<std::endl;
    splitter.refineQuads2EdgesOnCurveImproved();
    //gmds::MeshModelAlgo meshModelAlgo(mesh,model);
    meshModelAlgo.associateNodes();
    mesh.changeModel(mod_bis,false);
    doc.buildFacesAndR2F();  
    doc.updateUpwardConnectivity();
 
    std::cout<<"mesh.getNbFaces() "<<mesh.getNbFaces()<<std::endl;
    insertion.exportMeshVTK(case_name+std::string("_projectrefine"),gmds::N|gmds::E|gmds::F|gmds::R);
    mesh.changeModel(mod_bis,false);

    smooth->initNodesAdjacencies();
    for(int i=0; i<10; i++) {
       	//smooth->exec(markBoundaryNodes);
    }
	insertion.displayMeshQuality();
    insertion.exportMeshVTK(case_name+std::string("_projectrefinesmartsmooth"),gmds::N|gmds::E|gmds::F|gmds::R);

    // pillowing
    gmds::MeshModel mod_ter = gmds::DIM3|gmds::N|gmds::F|gmds::R|gmds::F2N|gmds::R2N|gmds::R2F|gmds::F2R|gmds::N2R;
    mesh.changeModel(mod_ter,false);

    insertion.removeUnassociatedFaces();
    insertion.pillow(vols[0]);
	std::cout<<"pillowINT"<<std::endl;
    insertion.removeUnassociatedFaces();
    insertion.pillowExt(vols[0]);
	std::cout<<"pillowEXT"<<std::endl;
	insertion.displayMeshQuality();
    insertion.exportMeshVTK(case_name+std::string("_pillow"),gmds::N|gmds::E|gmds::F|gmds::R);
	
//	insertion.project();
	insertion.displayMeshQuality();
	insertion.exportMeshVTK(case_name+std::string("_poyop"),gmds::N|gmds::E|gmds::F|gmds::R);
    smooth->initNodesAdjacencies();
/*
    for(int i=0; i<20; i++) {
	insertion.displayMeshQuality();
        smooth->exec(markBoundaryNodes);
	char s[10];
	sprintf(s,"%d",i);
	insertion.exportMeshVTK(case_name+std::string("_pillowsmooth")+std::string(s),gmds::N|gmds::F|gmds::R);
    }

	insertion.project();

    insertion.displayMeshQuality();

    insertion.exportMeshVTK(case_name+std::string("_pillowsmooth"),gmds::N|gmds::F|gmds::R);

	//mesh.changeModel(mod_AAAAAA,false); 
	//meshModelAlgo.associateNodes();
	gmds::NormalLaplacianSmoothing* smooth2 = new gmds::NormalLaplacianSmoothing(mesh,model,*service);
	smooth2->initNodesAdjacencies();
	for(int i=0; i<5; i++) {
		insertion.displayMeshQuality();
		smooth2->exec(markBoundaryNodes,vol);
	}
 */
    insertion.displayMeshQuality();
    insertion.exportMeshVTK(case_name+std::string("_pillowsmoothNormal"),gmds::N|gmds::F|gmds::R);

	gmds::GETMe getme(mesh);
	getme.execSimult(     
                        100,
			0.2,
                        0.5,
                        0.,
                        1.,
                        1.,
                        true,
                        true,
			false,
                        markBoundaryNodes
        );	

	insertion.displayMeshQuality();
	insertion.exportMeshVTK(case_name+std::string("_getme"),gmds::N|gmds::F|gmds::R);

	insertion.project();
	insertion.displayMeshQuality();
        insertion.exportMeshVTK(case_name+std::string("_getmeproj"),gmds::N|gmds::F|gmds::R);	
 
    EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
