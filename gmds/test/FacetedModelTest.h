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
#include <GMDS/CAD/FacetedGeomManager.h>
#include <GMDS/CAD/FacetedSurface.h>
#include <GMDS/CAD/FacetedVolume.h>
#include <GMDS/CAD/FACFacetedGeomReadAndWrite.h>
#include <GMDS/IO/VTKFacetedGeomReadAndWrite.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
class FacetedModelTest: public ::testing::Test {

  protected:
	FacetedModelTest(){;}
    virtual ~FacetedModelTest(){;}
};
/*----------------------------------------------------------------------------*/
TEST_F(FacetedModelTest,reorient_cas0) {

	gmds::geom::FacetedGeomManager model;
    	gmds::VTKFacetedGeomReadAndWrite vtkReader;
    	vtkReader.import(model,"Samples/cas0.vtk");

	std::vector<gmds::geom::GeomVolume* > vols;
    	model.getVolumes(vols);
    	gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);

	std::vector<gmds::geom::GeomSurface*> surfaces;
	vol->get(surfaces);

	gmds::geom::FacetedSurface* surf0 = dynamic_cast<gmds::geom::FacetedSurface*> (surfaces[0]);
	gmds::geom::FacetedSurface* surf1 = dynamic_cast<gmds::geom::FacetedSurface*> (surfaces[1]);
	gmds::geom::FacetedSurface* surf2 = dynamic_cast<gmds::geom::FacetedSurface*> (surfaces[2]);

	bool isOutward0 = surf0->isOutwardDirection(vol);
	bool isOutward1 = surf1->isOutwardDirection(vol);
	bool isOutward2 = surf2->isOutwardDirection(vol);

	EXPECT_FALSE(isOutward0);
        EXPECT_FALSE(isOutward1);
        EXPECT_FALSE(isOutward2);

	vol->reorient();

	isOutward0 = surf0->isOutwardDirection(vol);
        isOutward1 = surf1->isOutwardDirection(vol);
        isOutward2 = surf2->isOutwardDirection(vol);

	EXPECT_TRUE(isOutward0);
	EXPECT_TRUE(isOutward1);
	EXPECT_TRUE(isOutward2);
}
/*----------------------------------------------------------------------------*/
TEST_F(FacetedModelTest,reorient_pyramid) {

        gmds::geom::FacetedGeomManager model;
        gmds::VTKFacetedGeomReadAndWrite vtkReader;
        vtkReader.import(model,"Samples/pyramid.vtk");

        std::vector<gmds::geom::GeomVolume* > vols;
        model.getVolumes(vols);
        gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);

        std::vector<gmds::geom::GeomSurface*> surfaces;
        vol->get(surfaces);

        gmds::geom::FacetedSurface* surf0 = dynamic_cast<gmds::geom::FacetedSurface*> (surfaces[0]);
        gmds::geom::FacetedSurface* surf1 = dynamic_cast<gmds::geom::FacetedSurface*> (surfaces[1]);
        gmds::geom::FacetedSurface* surf2 = dynamic_cast<gmds::geom::FacetedSurface*> (surfaces[2]);
	gmds::geom::FacetedSurface* surf3 = dynamic_cast<gmds::geom::FacetedSurface*> (surfaces[3]);
	gmds::geom::FacetedSurface* surf4 = dynamic_cast<gmds::geom::FacetedSurface*> (surfaces[4]);

        bool isOutward0 = surf0->isOutwardDirection(vol);
        bool isOutward1 = surf1->isOutwardDirection(vol);
        bool isOutward2 = surf2->isOutwardDirection(vol);
	bool isOutward3 = surf3->isOutwardDirection(vol);
	bool isOutward4 = surf4->isOutwardDirection(vol);

        EXPECT_FALSE(isOutward0);
        EXPECT_TRUE(isOutward1);
        EXPECT_TRUE(isOutward2);
	EXPECT_TRUE(isOutward3);
	EXPECT_TRUE(isOutward4);
}
/*----------------------------------------------------------------------------*/
TEST_F(FacetedModelTest,reorient_cas1) {

        gmds::geom::FacetedGeomManager model;
        gmds::VTKFacetedGeomReadAndWrite vtkReader;
        vtkReader.import(model,"Samples/cas1.vtk");

        std::vector<gmds::geom::GeomVolume* > vols;
        model.getVolumes(vols);
        gmds::geom::FacetedVolume* vol = dynamic_cast<gmds::geom::FacetedVolume*> (vols[0]);

        std::vector<gmds::geom::GeomSurface*> surfaces;
        vol->get(surfaces);

        gmds::geom::FacetedSurface* surf0 = dynamic_cast<gmds::geom::FacetedSurface*> (surfaces[0]);
	gmds::geom::FacetedSurface* surf1 = dynamic_cast<gmds::geom::FacetedSurface*> (surfaces[1]);
	gmds::geom::FacetedSurface* surf2 = dynamic_cast<gmds::geom::FacetedSurface*> (surfaces[2]);
	gmds::geom::FacetedSurface* surf3 = dynamic_cast<gmds::geom::FacetedSurface*> (surfaces[3]);
	gmds::geom::FacetedSurface* surf4 = dynamic_cast<gmds::geom::FacetedSurface*> (surfaces[4]);
	gmds::geom::FacetedSurface* surf5 = dynamic_cast<gmds::geom::FacetedSurface*> (surfaces[5]);
	gmds::geom::FacetedSurface* surf6 = dynamic_cast<gmds::geom::FacetedSurface*> (surfaces[6]);
	gmds::geom::FacetedSurface* surf7 = dynamic_cast<gmds::geom::FacetedSurface*> (surfaces[7]);
	gmds::geom::FacetedSurface* surf8 = dynamic_cast<gmds::geom::FacetedSurface*> (surfaces[8]);
	gmds::geom::FacetedSurface* surf9 = dynamic_cast<gmds::geom::FacetedSurface*> (surfaces[9]);

        bool isOutward0 = surf0->isOutwardDirection(vol);
	bool isOutward1 = surf1->isOutwardDirection(vol);
	bool isOutward2 = surf2->isOutwardDirection(vol);
	bool isOutward3 = surf3->isOutwardDirection(vol);
	bool isOutward4 = surf4->isOutwardDirection(vol);
	bool isOutward5 = surf5->isOutwardDirection(vol);
	bool isOutward6 = surf6->isOutwardDirection(vol);
	bool isOutward7 = surf7->isOutwardDirection(vol);
	bool isOutward8 = surf8->isOutwardDirection(vol);
	bool isOutward9 = surf9->isOutwardDirection(vol);

        EXPECT_TRUE(isOutward0);
	EXPECT_FALSE(isOutward1);
	EXPECT_FALSE(isOutward2);
	EXPECT_TRUE(isOutward3);
	EXPECT_TRUE(isOutward4);
	EXPECT_FALSE(isOutward5);
	EXPECT_FALSE(isOutward6);
	EXPECT_FALSE(isOutward7);
	EXPECT_TRUE(isOutward8);
	EXPECT_FALSE(isOutward9);
}
/*----------------------------------------------------------------------------*/
