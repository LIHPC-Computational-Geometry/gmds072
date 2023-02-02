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
 * SubMappingTestSuite.h
 *
 *  Created on: 08 june 2015
 *      Author: legoff
 */
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <GMDS/SubMapping/DummyClassSubMapping.h>
#include <GMDS/SubMapping/SubMapping.h>

#include <GMDS/IO/VTKFacetedGeomReadAndWrite.h>
#include <GMDS/CAD/FacetedGeomManager.h>
#include <GMDS/CAD/GeomVolume.h>
#include <GMDS/CAD/FacetedVolume.h>

/*----------------------------------------------------------------------------*/
class SubMappingTestSuite: public ::testing::Test {

  protected:
	SubMappingTestSuite(){;}
    virtual ~SubMappingTestSuite(){;}
};
/*----------------------------------------------------------------------------*/
TEST_F(SubMappingTestSuite,cas_test_submapping) {

        std::cout<<"=================================================="<<std::endl;
        std::cout<<"SubMappingTestSuite::cas_test_submapping"<<std::endl;


	gmds::geom::FacetedGeomManager model;
	gmds::VTKFacetedGeomReadAndWrite vtkReader;
	vtkReader.import(model,"/home/legoff/travail/GMDS/gmds_workingcopy/GeomInsertion/test/Samples/cas2.vtk");

	std::vector<gmds::geom::GeomVolume* > vols;
	model.getVolumes(vols);
	gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
	vol->reorient();

	gmds::submapping::SubMapping sub;
	sub.exec(vol);

	std::vector<gmds::geom::GeomSurface*> surfaces;
        vol->get(surfaces);	
	sub.printVerticesClassification(surfaces[9]);
	sub.printSubmappingInfo(surfaces[9]);

	exit(-1);

	EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(SubMappingTestSuite,cas_test_submapping_surf) {

        std::cout<<"=================================================="<<std::endl;
        std::cout<<"SubMappingTestSuite::cas_test_submapping_surf"<<std::endl;

        gmds::geom::FacetedGeomManager model;
        gmds::VTKFacetedGeomReadAndWrite vtkReader;
        vtkReader.import(model,"/home/legoff/travail/GMDS/gmds_workingcopy/GeomInsertion/test/Samples/cas2.vtk");

        std::vector<gmds::geom::GeomVolume* > vols;
        model.getVolumes(vols);
        gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
        vol->reorient();

        std::vector<gmds::geom::GeomSurface*> surfaces;
        vol->get(surfaces);

	gmds::submapping::SubMapping sub;
        sub.exec(surfaces[9]);
        
	sub.printVerticesClassification(surfaces[9]);
        sub.printSubmappingInfo(surfaces[9]);

	exit(-1);

        EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(SubMappingTestSuite,cas_test_submapping_cas4) {

        std::cout<<"=================================================="<<std::endl;
        std::cout<<"SubMappingTestSuite::cas_test_submapping_cas4"<<std::endl;


        gmds::geom::FacetedGeomManager model;
        gmds::VTKFacetedGeomReadAndWrite vtkReader;
        vtkReader.import(model,"/home/legoff/travail/GMDS/gmds_workingcopy/GeomInsertion/test/Samples/cas4.vtk");

        std::vector<gmds::geom::GeomVolume* > vols;
        model.getVolumes(vols);
        gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
        vol->reorient();

        gmds::submapping::SubMapping sub;
        sub.exec(vol);

        EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(SubMappingTestSuite,cas_test_submapping_cas5) {

        std::cout<<"=================================================="<<std::endl;
        std::cout<<"SubMappingTestSuite::cas_test_submapping_cas5"<<std::endl;


        gmds::geom::FacetedGeomManager model;
        gmds::VTKFacetedGeomReadAndWrite vtkReader;
        vtkReader.import(model,"/home/legoff/travail/GMDS/gmds_workingcopy/GeomInsertion/test/Samples/cas5.vtk");

        std::vector<gmds::geom::GeomVolume* > vols;
        model.getVolumes(vols);
        gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
        vol->reorient();

        gmds::submapping::SubMapping sub;
        sub.exec(vol);

        EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(SubMappingTestSuite,cas_test_submapping_cas8) {

        std::cout<<"=================================================="<<std::endl;
        std::cout<<"SubMappingTestSuite::cas_test_submapping_cas8"<<std::endl;
        

        gmds::geom::FacetedGeomManager model;
        gmds::VTKFacetedGeomReadAndWrite vtkReader;
        vtkReader.import(model,"/home/legoff/travail/GMDS/gmds_workingcopy/GeomInsertion/test/Samples/cas8.vtk");
        
        std::vector<gmds::geom::GeomVolume* > vols;
        model.getVolumes(vols);
        gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
        vol->reorient();
        
        gmds::submapping::SubMapping sub;
        sub.exec(vol);

	sub.exportVTKSubmappingInfo(vol,"poyop");
        
        EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(SubMappingTestSuite,cas_test_submapping_easy1) {

        std::cout<<"=================================================="<<std::endl;
        std::cout<<"SubMappingTestSuite::cas_test_submapping_easy1"<<std::endl;


        gmds::geom::FacetedGeomManager model;
        gmds::VTKFacetedGeomReadAndWrite vtkReader;
        vtkReader.import(model,"/home/legoff/travail/GMDS/gmds_workingcopy/GeomInsertion/test/Samples/easy1.vtk");

        std::vector<gmds::geom::GeomVolume* > vols;
        model.getVolumes(vols);
        gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
        vol->reorient();

        gmds::submapping::SubMapping sub;
        sub.exec(vol);

        EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(SubMappingTestSuite,cas_test_submapping_easy2) {

        std::cout<<"=================================================="<<std::endl;
        std::cout<<"SubMappingTestSuite::cas_test_submapping_easy2"<<std::endl;


        gmds::geom::FacetedGeomManager model;
        gmds::VTKFacetedGeomReadAndWrite vtkReader;
        vtkReader.import(model,"/home/legoff/travail/GMDS/gmds_workingcopy/GeomInsertion/test/Samples/easy2.vtk");

        std::vector<gmds::geom::GeomVolume* > vols;
        model.getVolumes(vols);
        gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
        vol->reorient();

        gmds::submapping::SubMapping sub;
        sub.exec(vol);

        EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(SubMappingTestSuite,cas_test_submapping_medium2) {

        std::cout<<"=================================================="<<std::endl;
        std::cout<<"SubMappingTestSuite::cas_test_submapping_medium2"<<std::endl;


        gmds::geom::FacetedGeomManager model;
        gmds::VTKFacetedGeomReadAndWrite vtkReader;
        vtkReader.import(model,"/home/legoff/travail/GMDS/gmds_workingcopy/GeomInsertion/test/Samples/medium2.vtk");

        std::vector<gmds::geom::GeomVolume* > vols;
        model.getVolumes(vols);
        gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
        vol->reorient();

        gmds::submapping::SubMapping sub;
        sub.exec(vol);

        EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/
TEST_F(SubMappingTestSuite,cas_test_submapping_fan) {
        
        std::cout<<"=================================================="<<std::endl;
        std::cout<<"SubMappingTestSuite::cas_test_submapping_fan"<<std::endl;

        
        gmds::geom::FacetedGeomManager model;
        gmds::VTKFacetedGeomReadAndWrite vtkReader;
        vtkReader.import(model,"/home/legoff/travail/GMDS/gmds_workingcopy/GeomInsertion/test/Samples/fan.vtk");
        
        std::vector<gmds::geom::GeomVolume* > vols;
        model.getVolumes(vols);
        gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);
        vol->reorient();
        
        gmds::submapping::SubMapping sub;
        sub.exec(vol);

        EXPECT_EQ(0.0, 1.0);
}
/*----------------------------------------------------------------------------*/ 
