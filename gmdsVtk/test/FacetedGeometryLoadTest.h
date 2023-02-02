/*
 * FacetedGeometryTest.h
 *
 *  Created on: 6 juin 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <gtest/gtest.h>
/*----------------------------------------------------------------------------*/
#include <GMDS/CAD/FacetedGeomManager.h>
#include <GMDS/IOVTK/VTKFacetedGeomReadAndWrite.h>


/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
class FacetedGeometryLoadTest: public ::testing::Test {

  protected:
	FacetedGeometryLoadTest(){;}
    virtual ~FacetedGeometryLoadTest(){;}
};/*----------------------------------------------------------------------------*/
TEST_F(FacetedGeometryLoadTest,loadVTK_Box) {

	geom::FacetedGeomManager m;
	VTKFacetedGeomReadAndWrite rw;
	rw.import(m,"Samples/box_importVTK.vtp");

	EXPECT_TRUE(m.getNbSurfaces()!=0);
	EXPECT_TRUE(m.getNbCurves()!=0);
	EXPECT_TRUE(m.getNbPoints()!=0);

	std::vector<geom::GeomSurface* > surfs;
	m.getSurfaces(surfs);
	for(unsigned int i=0;i<surfs.size();i++)
	{

		//DISTANCE + PROJECTION TESTS
		geom::GeomSurface* s = surfs[i];

		math::Point p(0.5,0.5,0.5), proj=p;
		s->project(proj);
		TCoord dist = p.distance(proj);
		EXPECT_DOUBLE_EQ(0.5,dist);

		math::Point p2(0.5,0.5,0.0), proj2=p2;
		s->project(proj2);
		TCoord dist2 = p2.distance(proj2);
		EXPECT_TRUE(isZero(dist2) || isZero(dist2-1.0)  || isZero(dist2-0.5));

		//COMPUTATION AREA
		EXPECT_NEAR(1, s->computeArea(), TCoord_Epsilon);

		//ADJACENCY RELATIONS
		std::vector<geom::GeomCurve* > local_curves;
		s->get(local_curves);
		EXPECT_EQ(4, local_curves.size());

		std::vector<geom::GeomPoint* > local_vertices;
		s->get(local_vertices);
		EXPECT_EQ(4, local_vertices.size());


	}
	rw.exportVTK(m,"Samples/OUT/box_importVTK_exportVTK");
}
/*----------------------------------------------------------------------------*/
TEST_F(FacetedGeometryLoadTest,loadVTK_Cylinder) {

	geom::FacetedGeomManager m;
	VTKFacetedGeomReadAndWrite rw;
	rw.import(m,"Samples/cylinder_importVTK.vtp");


	EXPECT_TRUE(m.getNbSurfaces()!=0);
	EXPECT_TRUE(m.getNbCurves()!=0);
	EXPECT_TRUE(m.getNbPoints()!=0);
	rw.exportVTK(m,"Samples/OUT/cylinder_importVTK_exportVTK");
}
/*----------------------------------------------------------------------------*/
TEST_F(FacetedGeometryLoadTest,loadVTK_step) {

	geom::FacetedGeomManager m;
	VTKFacetedGeomReadAndWrite rw;
	rw.import(m,"Samples/step_importVTK.vtp");

	EXPECT_TRUE(m.getNbSurfaces()!=0);
	EXPECT_TRUE(m.getNbCurves()!=0);
	EXPECT_TRUE(m.getNbPoints()!=0);
	rw.exportVTK(m,"Samples/OUT/step_importVTK_exportVTK");
}
/*----------------------------------------------------------------------------*/
