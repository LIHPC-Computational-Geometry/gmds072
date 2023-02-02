/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux and N. Le Goff (2015)
 *
 * franck.ledoux@cea.fr
 * nicolas.le-goff@cea.fr
 *
 * This software is a computer program whose purpose is to provide a set of
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
 * data to be ensured and, more generally, to use and operate it in the 
 * same conditions as regards security. 
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*
 * FacetedMeshIntersectionService.t.h
 *
 *  Created on: 1 juil. 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <iostream>
#include <fstream>
#include <glib.h>
/*----------------------------------------------------------------------------*/
#include <GMDS/Math/Ray.h>
#include <GMDS/Math/Plane.h>
#include <GMDS/Math/Tetrahedron.h>
#include <GMDS/Math/Pyramid.h>
#include <GMDS/Math/Prism3.h>
/*----------------------------------------------------------------------------*/
#include <GMDS/FacetedMeshIntersectionService.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace geom{
/*----------------------------------------------------------------------------*/
FacetedMeshIntersectionService::FacetedMeshIntersectionService()
{

}
/*----------------------------------------------------------------------------*/
FacetedMeshIntersectionService::~FacetedMeshIntersectionService()
{

}
/*----------------------------------------------------------------------------*/
//template<typename TBase >
//void FacetedMeshIntersectionService<TBase>::
//computeBoundingBox(std::vector<GeomSurface<TBase>* >& ASurfaces, TBase minXYZ[3], TBase maxXYZ[3]) const
//{
//	if(ASurfaces.size() == 0)
//	{
//		minXYZ[0] = minXYZ[1] = minXYZ[2] = maxXYZ[0] = maxXYZ[1] = maxXYZ[2] = 0;
//		return;
//	}
//
//	FacetedSurface<TBase>* s = dynamic_cast<FacetedSurface<TBase>*>(ASurfaces[0]);
//	if(s==0)
//		throw GMDSException("Wrong geometric surface type");
//	s->computeBoundingBox(minXYZ,maxXYZ);
//
//	TBase minXYZ_tmp[3], maxXYZ_tmp[3];
//
//	for(int iSurface=1; iSurface<ASurfaces.size(); iSurface++)
//	{
//		s = dynamic_cast<FacetedSurface<TBase>*>(ASurfaces[iSurface]);
//
//		if(s==0)
//			throw GMDSException("Wrong geometric surface type");
//		else
//		{
//			s->computeBoundingBox(minXYZ_tmp,maxXYZ_tmp);
//
//			if(minXYZ_tmp[0]<minXYZ[0])
//				minXYZ[0] = minXYZ_tmp[0];
//			if(minXYZ_tmp[1]<minXYZ[1])
//				minXYZ[1] = minXYZ_tmp[1];
//			if(minXYZ_tmp[2]<minXYZ[2])
//				minXYZ[2] = minXYZ_tmp[2];
//			if(maxXYZ_tmp[0]>maxXYZ[0])
//				maxXYZ[0] = maxXYZ_tmp[0];
//			if(maxXYZ_tmp[1]>maxXYZ[1])
//				maxXYZ[1] = maxXYZ_tmp[1];
//			if(maxXYZ_tmp[2]>maxXYZ[2])
//				maxXYZ[2] = maxXYZ_tmp[2];
//		}
//	}
//
//	return;
//}
//
/*----------------------------------------------------------------------------*/
void
FacetedMeshIntersectionService::initialization(std::vector<GeomSurface*>& ASurfaces)
{
	for(size_t iSurface=0; iSurface<ASurfaces.size(); iSurface++)
	{
		std::vector<gmds::math::Triangle> triangles;
		ASurfaces[iSurface]->getTriangulation(triangles);
		m_surfacesTriangulation[ASurfaces[iSurface]] = triangles;
	}

	for(size_t iSurface=0; iSurface<ASurfaces.size(); iSurface++)
	{
		double coordMin[3];
		double coordMax[3];
		ASurfaces[iSurface]->computeBoundingBox(coordMin,coordMax);

		std::vector<gmds::math::Point> points;
		points.push_back(gmds::math::Point(coordMin[0],coordMin[1],coordMin[2]));
		points.push_back(gmds::math::Point(coordMax[0],coordMax[1],coordMax[2]));
		m_surfacesBoundingBox[ASurfaces[iSurface]] = points;

		// valid in c++0x
//		m_surfacesTriangulation[ASurfaces[iSurface]].push_back
//				gmds::math::Point(coordMin[0],coordMin[1],coordMin[2]),
//				gmds::math::Point(coordMax[0],coordMax[1],coordMax[2])
//		};
	}

}
/*----------------------------------------------------------------------------*/
bool FacetedMeshIntersectionService::
intersects(std::vector<GeomSurface* >& ASurfaces, gmds::math::Hexahedron& AHex) const
{

	for(size_t iSurface=0; iSurface<ASurfaces.size(); iSurface++)
	{
		// first discriminate against the bounding box
		gmds::math::Point points[8];
		points[0] = AHex.getPoint(0);
		points[1] = AHex.getPoint(1);
		points[2] = AHex.getPoint(2);
		points[3] = AHex.getPoint(3);
		points[4] = AHex.getPoint(4);
		points[5] = AHex.getPoint(5);
		points[6] = AHex.getPoint(6);
		points[7] = AHex.getPoint(7);

		gmds::math::Point pmin = m_surfacesBoundingBox[ASurfaces[iSurface]][0];
		gmds::math::Point pmax = m_surfacesBoundingBox[ASurfaces[iSurface]][1];

		if((points[0].X()<pmin.X() && points[1].X()<pmin.X() && points[2].X()<pmin.X() && points[3].X()<pmin.X() &&
		    points[4].X()<pmin.X() && points[5].X()<pmin.X() && points[6].X()<pmin.X() && points[7].X()<pmin.X()) ||
		   (points[0].Y()<pmin.Y() && points[1].Y()<pmin.Y() && points[2].Y()<pmin.Y() && points[3].Y()<pmin.Y() &&
			points[4].Y()<pmin.Y() && points[5].Y()<pmin.Y() && points[6].Y()<pmin.Y() && points[7].Y()<pmin.Y()) ||
		   (points[0].Z()<pmin.Z() && points[1].Z()<pmin.Z() && points[2].Z()<pmin.Z() && points[3].Z()<pmin.Z() &&
			points[4].Z()<pmin.Z() && points[5].Z()<pmin.Z() && points[6].Z()<pmin.Z() && points[7].Z()<pmin.Z()) ||
		   (points[0].X()>pmax.X() && points[1].X()>pmax.X() && points[2].X()>pmax.X() && points[3].X()>pmax.X() &&
			points[4].X()>pmax.X() && points[5].X()>pmax.X() && points[6].X()>pmax.X() && points[7].X()>pmax.X()) ||
		   (points[0].Y()>pmax.Y() && points[1].Y()>pmax.Y() && points[2].Y()>pmax.Y() && points[3].Y()>pmax.Y() &&
			points[4].Y()>pmax.Y() && points[5].Y()>pmax.Y() && points[6].Y()>pmax.Y() && points[7].Y()>pmax.Y()) ||
		   (points[0].Z()>pmax.Z() && points[1].Z()>pmax.Z() && points[2].Z()>pmax.Z() && points[3].Z()>pmax.Z() &&
			points[4].Z()>pmax.Z() && points[5].Z()>pmax.Z() && points[6].Z()>pmax.Z() && points[7].Z()>pmax.Z())
				   ) {
			// outside of this surface bounding box, so we check the next surface
			continue;
		}

		for(size_t iTriangle=0; iTriangle<m_surfacesTriangulation[ASurfaces[iSurface]].size(); iTriangle++)
		{
			bool result = AHex.intersect(m_surfacesTriangulation[ASurfaces[iSurface]][iTriangle]);
			if(result) {
				return true;
			}
		}
	}

	return false;
}
/*----------------------------------------------------------------------------*/
bool FacetedMeshIntersectionService::
intersects(std::vector<GeomSurface* >& ASurfaces, gmds::Region& ARegion) const
{
	TCoord minXYZ[3];
	TCoord maxXYZ[3];

	ARegion.computeBoundingBox(minXYZ,maxXYZ);

	for(size_t iSurface=0; iSurface<ASurfaces.size(); iSurface++)
	{
		// first discriminate against the bounding box
		gmds::math::Point pmin = m_surfacesBoundingBox[ASurfaces[iSurface]][0];
		gmds::math::Point pmax = m_surfacesBoundingBox[ASurfaces[iSurface]][1];

		if((maxXYZ[0]<pmin.X()) || (maxXYZ[1]<pmin.Y()) || (maxXYZ[2]<pmin.Z()) ||
		   (minXYZ[0]>pmax.X()) || (minXYZ[1]>pmax.Y()) || (minXYZ[2]>pmax.Z())) {
			// outside of this surface bounding box, so we check the next surface
			continue;
		}

		std::vector<gmds::Node> nodes = ARegion.get<gmds::Node>();

		switch(ARegion.getType())
		{
		case gmds::GMDS_HEX :
			{
				gmds::math::Hexahedron hexa(nodes[4].getPoint(),nodes[7].getPoint(),nodes[6].getPoint(),nodes[5].getPoint(),
					nodes[0].getPoint(),nodes[3].getPoint(),nodes[2].getPoint(),nodes[1].getPoint());

				for(size_t iTriangle=0; iTriangle<m_surfacesTriangulation[ASurfaces[iSurface]].size(); iTriangle++)
				{
					bool result = hexa.intersect(m_surfacesTriangulation[ASurfaces[iSurface]][iTriangle]);
					if(result) {
						return true;
					}
				}
				break;
			}
		case gmds::GMDS_TETRA :
			{
				gmds::math::Tetrahedron tet(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),nodes[3].getPoint());

				for(size_t iTriangle=0; iTriangle<m_surfacesTriangulation[ASurfaces[iSurface]].size(); iTriangle++)
				{
					bool result = tet.intersect(m_surfacesTriangulation[ASurfaces[iSurface]][iTriangle]);
					if(result) {
						return true;
					}
				}
				break;
			}
		case gmds::GMDS_PYRAMID :
			{
				gmds::math::Pyramid pyr(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),nodes[3].getPoint(),nodes[4].getPoint());

				for(size_t iTriangle=0; iTriangle<m_surfacesTriangulation[ASurfaces[iSurface]].size(); iTriangle++)
				{
					bool result = pyr.intersect(m_surfacesTriangulation[ASurfaces[iSurface]][iTriangle]);
					if(result) {
						return true;
					}
				}
				break;
			}
		case gmds::GMDS_PRISM3 :
			{
				gmds::math::Prism3 prsm(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),
						nodes[3].getPoint(),nodes[4].getPoint(),nodes[5].getPoint());

				for(size_t iTriangle=0; iTriangle<m_surfacesTriangulation[ASurfaces[iSurface]].size(); iTriangle++)
				{
					bool result = prsm.intersect(m_surfacesTriangulation[ASurfaces[iSurface]][iTriangle]);
					if(result) {
						return true;
					}
				}
				break;
			}
		default:
			throw GMDSException("FacetedMeshIntersectionService::intersects cell type not treated yet");
		}
	}

	return false;
}
/*----------------------------------------------------------------------------*/
bool
FacetedMeshIntersectionService::intersects(
		std::vector<GeomSurface* >& ASurfaces,
		gmds::Region& ARegion,
		std::vector<gmds::math::Triangle>& ATriangles) const
{
	ATriangles.clear();

	TCoord minXYZ[3];
	TCoord maxXYZ[3];

	ARegion.computeBoundingBox(minXYZ,maxXYZ);

	for(size_t iSurface=0; iSurface<ASurfaces.size(); iSurface++)
	{
		// first discriminate against the bounding box
		gmds::math::Point pmin = m_surfacesBoundingBox[ASurfaces[iSurface]][0];
		gmds::math::Point pmax = m_surfacesBoundingBox[ASurfaces[iSurface]][1];

		if((maxXYZ[0]<pmin.X()) || (maxXYZ[1]<pmin.Y()) || (maxXYZ[2]<pmin.Z()) ||
		   (minXYZ[0]>pmax.X()) || (minXYZ[1]>pmax.Y()) || (minXYZ[2]>pmax.Z())) {
			// outside of this surface bounding box, so we check the next surface
			continue;
		}

		std::vector<gmds::math::Triangle>& surfaceTriangulation = m_surfacesTriangulation[ASurfaces[iSurface]];

		std::vector<gmds::Node> nodes = ARegion.get<gmds::Node>();

		switch(ARegion.getType())
		{
		case gmds::GMDS_HEX :
			{
				gmds::math::Hexahedron hexa(nodes[4].getPoint(),nodes[7].getPoint(),nodes[6].getPoint(),nodes[5].getPoint(),
					nodes[0].getPoint(),nodes[3].getPoint(),nodes[2].getPoint(),nodes[1].getPoint());

				for(size_t iTriangle=0; iTriangle<surfaceTriangulation.size(); iTriangle++)
				{
					bool result = hexa.intersect(surfaceTriangulation[iTriangle]);
					if(result) {
						ATriangles.push_back(surfaceTriangulation[iTriangle]);
					}
				}
				break;
			}
		case gmds::GMDS_TETRA :
			{
				gmds::math::Tetrahedron tet(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),nodes[3].getPoint());

				for(size_t iTriangle=0; iTriangle<surfaceTriangulation.size(); iTriangle++)
				{
					bool result = tet.intersect(surfaceTriangulation[iTriangle]);
					if(result) {
						ATriangles.push_back(surfaceTriangulation[iTriangle]);
					}
				}
				break;
			}
		case gmds::GMDS_PYRAMID :
			{
				gmds::math::Pyramid pyr(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),nodes[3].getPoint(),nodes[4].getPoint());

				for(size_t iTriangle=0; iTriangle<surfaceTriangulation.size(); iTriangle++)
				{
					bool result = pyr.intersect(surfaceTriangulation[iTriangle]);
					if(result) {
						ATriangles.push_back(surfaceTriangulation[iTriangle]);
					}
				}
				break;
			}
		case gmds::GMDS_PRISM3 :
			{
				gmds::math::Prism3 prsm(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),
						nodes[3].getPoint(),nodes[4].getPoint(),nodes[5].getPoint());

				for(size_t iTriangle=0; iTriangle<surfaceTriangulation.size(); iTriangle++)
				{
					bool result = prsm.intersect(surfaceTriangulation[iTriangle]);
					if(result) {
						ATriangles.push_back(surfaceTriangulation[iTriangle]);
					}
				}
				break;
			}
		default:
			throw GMDSException("FacetedMeshIntersectionService::intersects cell type not treated yet");
		}
	}

	return (!ATriangles.empty());
}
/*----------------------------------------------------------------------------*/
bool
FacetedMeshIntersectionService::intersects(
		gmds::Region& ARegion,
		std::vector<gmds::math::Triangle>& ATriangles) const
{
	ATriangles.clear();

	TCoord minXYZ[3];
	TCoord maxXYZ[3];

	ARegion.computeBoundingBox(minXYZ,maxXYZ);

	std::map<gmds::geom::GeomSurface*,std::vector<gmds::math::Triangle> >::iterator itSurf = m_surfacesTriangulation.begin();

	for(; itSurf != m_surfacesTriangulation.end(); itSurf++) {

		// first discriminate against the bounding box
		gmds::math::Point pmin = m_surfacesBoundingBox[itSurf->first][0];
		gmds::math::Point pmax = m_surfacesBoundingBox[itSurf->first][1];

		if((maxXYZ[0]<pmin.X()) || (maxXYZ[1]<pmin.Y()) || (maxXYZ[2]<pmin.Z()) ||
		   (minXYZ[0]>pmax.X()) || (minXYZ[1]>pmax.Y()) || (minXYZ[2]>pmax.Z())) {
			// outside of this surface bounding box, so we check the next surface
			continue;
		}

		std::vector<gmds::math::Triangle>& surfaceTriangulation = itSurf->second;

		std::vector<gmds::Node> nodes = ARegion.get<gmds::Node>();

		switch(ARegion.getType())
		{
		case gmds::GMDS_HEX :
			{
				gmds::math::Hexahedron hexa(nodes[4].getPoint(),nodes[7].getPoint(),nodes[6].getPoint(),nodes[5].getPoint(),
					nodes[0].getPoint(),nodes[3].getPoint(),nodes[2].getPoint(),nodes[1].getPoint());

				for(size_t iTriangle=0; iTriangle<surfaceTriangulation.size(); iTriangle++)
				{
					bool result = hexa.intersect(surfaceTriangulation[iTriangle]);
					if(result) {
						ATriangles.push_back(surfaceTriangulation[iTriangle]);
					}
				}
				break;
			}
		case gmds::GMDS_TETRA :
			{
				gmds::math::Tetrahedron tet(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),nodes[3].getPoint());

				for(size_t iTriangle=0; iTriangle<surfaceTriangulation.size(); iTriangle++)
				{
					bool result = tet.intersect(surfaceTriangulation[iTriangle]);
					if(result) {
						ATriangles.push_back(surfaceTriangulation[iTriangle]);
					}
				}
				break;
			}
		case gmds::GMDS_PYRAMID :
			{
				gmds::math::Pyramid pyr(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),nodes[3].getPoint(),nodes[4].getPoint());

				for(size_t iTriangle=0; iTriangle<surfaceTriangulation.size(); iTriangle++)
				{
					bool result = pyr.intersect(surfaceTriangulation[iTriangle]);
					if(result) {
						ATriangles.push_back(surfaceTriangulation[iTriangle]);
					}
				}
				break;
			}
		case gmds::GMDS_PRISM3 :
			{
				gmds::math::Prism3 prsm(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),
						nodes[3].getPoint(),nodes[4].getPoint(),nodes[5].getPoint());

				for(size_t iTriangle=0; iTriangle<surfaceTriangulation.size(); iTriangle++)
				{
					bool result = prsm.intersect(surfaceTriangulation[iTriangle]);
					if(result) {
						ATriangles.push_back(surfaceTriangulation[iTriangle]);
					}
				}
				break;
			}
		default:
			throw GMDSException("FacetedMeshIntersectionService::intersects cell type not treated yet");
		}
	}

	return (!ATriangles.empty());
}
/*----------------------------------------------------------------------------*/
bool
FacetedMeshIntersectionService::intersects(
		GNode* ATree,
		std::vector<GeomSurface* >& ASurfaces,
		gmds::Region& ARegion,
		std::vector<gmds::math::Triangle>& ATriangles) const
{
	ATriangles.clear();

	TCoord minXYZ[3];
	TCoord maxXYZ[3];
	ARegion.computeBoundingBox(minXYZ,maxXYZ);

	gpointer pointer = NULL;
	GtsBBox* bbox = gts_bbox_new(
							gts_bbox_class (),
							pointer,
							minXYZ[0],minXYZ[1],minXYZ[2],
							maxXYZ[0],maxXYZ[1],maxXYZ[2]);

	GSList* boxList = gts_bb_tree_overlap(ATree,bbox);
	if(boxList == NULL) {
		return false;
	}

	std::vector<gmds::Node> nodes = ARegion.get<gmds::Node>();

	switch(ARegion.getType())
	{
	case gmds::GMDS_HEX :
		{
			gmds::math::Hexahedron hexa(nodes[4].getPoint(),nodes[7].getPoint(),nodes[6].getPoint(),nodes[5].getPoint(),
				nodes[0].getPoint(),nodes[3].getPoint(),nodes[2].getPoint(),nodes[1].getPoint());

			while(boxList != NULL) {
				GtsBBox* box = (GtsBBox*)(boxList->data);
				gmds::math::Triangle* triangle = (gmds::math::Triangle*)(box->bounded);
				if(triangle == NULL) {
					throw GMDSException("FacetedMeshIntersectionService::intersects error dynamic_cast");
				}

				bool result = hexa.intersect(*triangle);
				if(result) {
					ATriangles.push_back(*triangle);
				}

				boxList = g_slist_next(boxList);
			}
			break;
		}
		case gmds::GMDS_TETRA :
			{
				gmds::math::Tetrahedron tet(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),nodes[3].getPoint());

				while(boxList != NULL) {

					gmds::math::Triangle* triangle = (gmds::math::Triangle*)(boxList->data);
					if(triangle == NULL) {
						throw GMDSException("FacetedMeshIntersectionService::intersects error dynamic_cast");
					}

					bool result = tet.intersect(*triangle);
					if(result) {
						ATriangles.push_back(*triangle);
					}

					boxList = g_slist_next(boxList);
				}
				break;
			}
		case gmds::GMDS_PYRAMID :
			{
				gmds::math::Pyramid pyr(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),nodes[3].getPoint(),nodes[4].getPoint());

				while(boxList != NULL) {

					gmds::math::Triangle* triangle = (gmds::math::Triangle*)(boxList->data);
					if(triangle == NULL) {
						throw GMDSException("FacetedMeshIntersectionService::intersects error dynamic_cast");
					}

					bool result = pyr.intersect(*triangle);
					if(result) {
						ATriangles.push_back(*triangle);
					}

					boxList = g_slist_next(boxList);
				}
				break;
			}
		case gmds::GMDS_PRISM3 :
			{
				gmds::math::Prism3 prsm(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),
						nodes[3].getPoint(),nodes[4].getPoint(),nodes[5].getPoint());

				while(boxList != NULL) {

					gmds::math::Triangle* triangle = (gmds::math::Triangle*)(boxList->data);
					if(triangle == NULL) {
						throw GMDSException("FacetedMeshIntersectionService::intersects error dynamic_cast");
					}

					bool result = prsm.intersect(*triangle);
					if(result) {
						ATriangles.push_back(*triangle);
					}

					boxList = g_slist_next(boxList);
				}
				break;
			}
		default:
			throw GMDSException("FacetedMeshIntersectionService::intersects cell type not treated yet");
		}

	return (!ATriangles.empty());
}
/*----------------------------------------------------------------------------*/
void
FacetedMeshIntersectionService::intersects(
		GNode* ATree,
		std::map<gmds::TCellID, std::vector<gmds::math::Triangle*> >& AInOutTriangles) const
{

	std::map<gmds::geom::GeomSurface*,std::vector<gmds::math::Triangle> >::iterator itSurf = m_surfacesTriangulation.begin();

	for(; itSurf != m_surfacesTriangulation.end(); itSurf++) {

		std::vector<gmds::math::Triangle>& surfaceTriangulation = itSurf->second;

		for(unsigned int iTriangle=0; iTriangle<surfaceTriangulation.size(); iTriangle++) {

			TCoord minXYZ[3];
			TCoord maxXYZ[3];
			surfaceTriangulation[iTriangle].computeBoundingBox(minXYZ,maxXYZ);

			gpointer pointer = &(surfaceTriangulation[iTriangle]);
			GtsBBox* bbox = gts_bbox_new(
									gts_bbox_class (),
									pointer,
									minXYZ[0],minXYZ[1],minXYZ[2],
									maxXYZ[0],maxXYZ[1],maxXYZ[2]);

			GSList* boxList = gts_bb_tree_overlap(ATree,bbox);
			while(boxList != NULL) {

				GtsBBox* box = (GtsBBox*)(boxList->data);
				AInOutTriangles[GPOINTER_TO_INT(box->bounded)].push_back(&(surfaceTriangulation[iTriangle]));

				boxList = g_slist_next(boxList);
			}

		}
	}

}
/*----------------------------------------------------------------------------*/
void
FacetedMeshIntersectionService::intersects(
		GNode* ATree,
		std::map<gmds::TCellID, std::map<gmds::geom::GeomSurface*,std::vector<gmds::math::Triangle*> > >& AInOutSurfTriangles) const
{

	std::map<gmds::geom::GeomSurface*,std::vector<gmds::math::Triangle> >::iterator itSurf = m_surfacesTriangulation.begin();

	for(; itSurf != m_surfacesTriangulation.end(); itSurf++) {

		std::vector<gmds::math::Triangle>& surfaceTriangulation = itSurf->second;

		for(unsigned int iTriangle=0; iTriangle<surfaceTriangulation.size(); iTriangle++) {

			TCoord minXYZ[3];
			TCoord maxXYZ[3];
			surfaceTriangulation[iTriangle].computeBoundingBox(minXYZ,maxXYZ);

			gpointer pointer = &(surfaceTriangulation[iTriangle]);
			GtsBBox* bbox = gts_bbox_new(
									gts_bbox_class (),
									pointer,
									minXYZ[0],minXYZ[1],minXYZ[2],
									maxXYZ[0],maxXYZ[1],maxXYZ[2]);

			GSList* boxList = gts_bb_tree_overlap(ATree,bbox);
			while(boxList != NULL) {

				GtsBBox* box = (GtsBBox*)(boxList->data);

				AInOutSurfTriangles[GPOINTER_TO_INT(box->bounded)][itSurf->first].push_back(&(surfaceTriangulation[iTriangle]));

				boxList = g_slist_next(boxList);
			}

		}
	}

}
/*----------------------------------------------------------------------------*/
bool
FacetedMeshIntersectionService::intersects(
		gmds::math::Triangle& ATri, gmds::Region& ARegion) const
{
	std::vector<gmds::Node> nodes = ARegion.get<gmds::Node>();

	switch(ARegion.getType())
	{
	case gmds::GMDS_HEX :
	{
		gmds::math::Hexahedron hexa(nodes[4].getPoint(),nodes[7].getPoint(),nodes[6].getPoint(),nodes[5].getPoint(),
				nodes[0].getPoint(),nodes[3].getPoint(),nodes[2].getPoint(),nodes[1].getPoint());

		bool result = hexa.intersect(ATri);
		if(result) {
			return true;
		}
		break;
	}
	case gmds::GMDS_TETRA :
	{
		gmds::math::Tetrahedron tet(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),nodes[3].getPoint());

		bool result = tet.intersect(ATri);
		if(result) {
			return true;
		}
		break;
	}
	case gmds::GMDS_PYRAMID :
	{
		gmds::math::Pyramid pyr(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),nodes[3].getPoint(),nodes[4].getPoint());

		bool result = pyr.intersect(ATri);
		if(result) {
			return true;
		}
		break;
	}
	case gmds::GMDS_PRISM3 :
	{
		gmds::math::Prism3 prsm(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint(),
				nodes[3].getPoint(),nodes[4].getPoint(),nodes[5].getPoint());

		bool result = prsm.intersect(ATri);
		if(result) {
			return true;
		}
		break;
	}
	default:
		throw GMDSException("FacetedMeshIntersectionService::intersects cell type not treated yet");
	}

	return false;
}
/*----------------------------------------------------------------------------*/
bool FacetedMeshIntersectionService::
intersectsTwiceTheModel(std::vector<GeomSurface*>& ASurfaces,
		std::map<gmds::geom::GeomSurface*, std::set<gmds::geom::GeomSurface*> > ASurfacesNeighbors,
		gmds::Region& ARegion) const
{
	std::vector<std::vector<gmds::math::Triangle> > triangles;

	for(size_t iSurface=0; iSurface<ASurfaces.size(); iSurface++)
	{
		std::vector<gmds::math::Triangle> triangles_tmp;
		ASurfaces[iSurface]->getTriangulation(triangles_tmp);

		triangles.push_back(triangles_tmp);
	}

	std::vector<unsigned int> iSurfaces_found;

	for(size_t iSurface=0; iSurface<ASurfaces.size(); iSurface++)
	{
		TCoord minXYZ[3];
		TCoord maxXYZ[3];

		ARegion.computeBoundingBox(minXYZ,maxXYZ);

		// first discriminate against the bounding box
		gmds::math::Point pmin = m_surfacesBoundingBox[ASurfaces[iSurface]][0];
		gmds::math::Point pmax = m_surfacesBoundingBox[ASurfaces[iSurface]][1];

		if((maxXYZ[0]<pmin.X()) || (maxXYZ[1]<pmin.Y()) || (maxXYZ[2]<pmin.Z()) ||
		   (minXYZ[0]>pmax.X()) || (minXYZ[1]>pmax.Y()) || (minXYZ[2]>pmax.Z())) {
			// outside of this surface bounding box, so we check the next surface
			continue;
		}

		std::vector<gmds::Node> nodes = ARegion.get<gmds::Node>();

		gmds::math::Hexahedron hexa(nodes[4].getPoint(),nodes[7].getPoint(),nodes[6].getPoint(),nodes[5].getPoint(),
				nodes[0].getPoint(),nodes[3].getPoint(),nodes[2].getPoint(),nodes[1].getPoint());

		for(size_t iTriangle=0; iTriangle<triangles[iSurface].size(); iTriangle++)
		{
			if(hexa.intersect(triangles[iSurface][iTriangle])) {
				iSurfaces_found.push_back(iSurface);
				break;
			}
		}
	}

	for(unsigned int iSurface_found=0; iSurface_found<iSurfaces_found.size(); iSurface_found++) {
		for(unsigned int iSurface_found2=iSurface_found+1; iSurface_found2<iSurfaces_found.size(); iSurface_found2++) {
			if((ASurfacesNeighbors[ASurfaces[iSurfaces_found[iSurface_found]]]).find(ASurfaces[iSurfaces_found[iSurface_found2]]) == (ASurfacesNeighbors[ASurfaces[iSurfaces_found[iSurface_found]]]).end()) {
				return true;
			}
		}
	}

	return false;
}
///*----------------------------------------------------------------------------*/
//template<typename TBase >
//bool FacetedMeshIntersectionService<TBase>::
//intersects(GeomCurve<TBase>& ACurve, GEPETO::Hexahedron<TBase>& AHex) const
//{
//	// epsilon in order to switch between TBase and  with epsilon precision
//	const int Epsilon = FacetedMeshIntersectionService_intersect_curve;
//
//	GEPETO::Point<3, GEPETO::NumericEps<double,Epsilon> > p[8];
//	for(int iPoint=0; iPoint<8; iPoint++)
//	{
//		p[iPoint] = AHex.getPoint(iPoint);
//	}
//
//	GEPETO::Hexahedron<GEPETO::NumericEps<double,Epsilon> > hex_eps(p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7]);
//
//	FacetedCurve<TBase>* curve = dynamic_cast<FacetedCurve<TBase>*>(&ACurve);
//
//	for(unsigned int iSegment=0; iSegment<curve->getNbSegments(); iSegment++) {
//
//		GEPETO::Segment<3, TBase> segment = curve->getSegment(iSegment);
//		GEPETO::Point<3, GEPETO::NumericEps<double,Epsilon> > point1_eps = segment.getPoint(0);
//		GEPETO::Point<3, GEPETO::NumericEps<double,Epsilon> > point2_eps = segment.getPoint(1);
//		GEPETO::Segment<3, GEPETO::NumericEps<double,Epsilon> > segment_eps(point1_eps,point2_eps);
//
//		GEPETO::EGeomPredicateResult result = hex_eps.intersect(segment_eps);
//
//		if((result == GEPETO::GEOM_YES) || (result == GEPETO::GEOM_UNDEF)) {
//			return true;
//		}
//	}
//
//	return false;
//}
///*----------------------------------------------------------------------------*/
//template<typename TBase >
//bool FacetedMeshIntersectionService<TBase>::
//isInside(std::vector<GeomSurface<TBase>* >& ASurfaces, GEPETO::Hexahedron<TBase>& AHex) const
//{
////	// count number of intersections of
////	GEPETO::Line<3,TBase> l1(AHex.p[0],AHex.p[1]);
////	GEPETO::Line<3,TBase> l2(AHex.p[0],AHex.p[2]);
////	GEPETO::Line<3,TBase> l3(AHex.p[0],AHex.p[3]);
////	GEPETO::Line<3,TBase> l4(AHex.p[0],AHex.p[4]);
////	GEPETO::Line<3,TBase> l5(AHex.p[0],AHex.p[5]);
////	GEPETO::Line<3,TBase> l6(AHex.p[0],AHex.p[6]);
////	GEPETO::Line<3,TBase> l7(AHex.p[0],AHex.p[7]);
////
////	size_t nbIntersections[7] = {0,0,0,0,0,0,0};
////
////	for(size_t iSurface=0; iSurface<ASurfaces.size(); iSurface++)
////	{
////		std::vector<GEPETO::Triangle<3,TBase> > triangles;
////
////		FacetedSurface<TBase>* s = dynamic_cast<FacetedSurface<TBase>*>(ASurfaces[iSurface]);
////		if(s==0)
////			throw GMDSException("Wrong geometric type");
////		else
////			s->getTriangulation(triangles);
////
////		for(size_t iTriangle=0; iTriangle<triangles.size(); iTriangle++)
////		{
////			if(GEPETO::GEOM_YES == GEPETO::GeomToolkit<3,TBase>::intersect(l1,triangles[iTriangle]))
////				nbIntersections[0]++;
////			if(GEPETO::GEOM_YES == GEPETO::GeomToolkit<3,TBase>::intersect(l2,triangles[iTriangle]))
////				nbIntersections[1]++;
////			if(GEPETO::GEOM_YES == GEPETO::GeomToolkit<3,TBase>::intersect(l3,triangles[iTriangle]))
////				nbIntersections[2]++;
////			if(GEPETO::GEOM_YES == GEPETO::GeomToolkit<3,TBase>::intersect(l4,triangles[iTriangle]))
////				nbIntersections[3]++;
////			if(GEPETO::GEOM_YES == GEPETO::GeomToolkit<3,TBase>::intersect(l5,triangles[iTriangle]))
////				nbIntersections[4]++;
////			if(GEPETO::GEOM_YES == GEPETO::GeomToolkit<3,TBase>::intersect(l6,triangles[iTriangle]))
////				nbIntersections[5]++;
////			if(GEPETO::GEOM_YES == GEPETO::GeomToolkit<3,TBase>::intersect(l7,triangles[iTriangle]))
////				nbIntersections[6]++;
////
////		}
////	}
////
////	bool isInner[7];
////
////	for(size_t iLine=0; iLine<7; iLine++)
////	{
////		isInner[iLine] = ((nbIntersections[iLine]%2) == 0);
////	}
////
////	if(!(isInner[0] == isInner[1] == isInner[2] == isInner[3] == isInner[4] == isInner[5] == isInner[6]))
////		throw GMDSException("Problem in isInside detection!");
////
////	return isInner[0];
//
//	// count number of intersections of
//	const int Epsilon = 14;
//
//	GEPETO::Point<3, GEPETO::NumericEps<double,Epsilon> > p[8];
//	for(int iPoint=0; iPoint<8; iPoint++)
//	{
//		p[iPoint] = AHex.getPoint(iPoint);
//	}
//
//	GEPETO::Ray<3,GEPETO::NumericEps<double,Epsilon> > r1(p[0],p[1]);
//	GEPETO::Ray<3,GEPETO::NumericEps<double,Epsilon> > r2(p[0],p[2]);
//	GEPETO::Ray<3,GEPETO::NumericEps<double,Epsilon> > r3(p[0],p[3]);
//	GEPETO::Ray<3,GEPETO::NumericEps<double,Epsilon> > r4(p[0],p[4]);
//	GEPETO::Ray<3,GEPETO::NumericEps<double,Epsilon> > r5(p[0],p[5]);
//	GEPETO::Ray<3,GEPETO::NumericEps<double,Epsilon> > r6(p[0],p[6]);
//	GEPETO::Ray<3,GEPETO::NumericEps<double,Epsilon> > r7(p[0],p[7]);
//
//	size_t nbIntersections[7] = {0,0,0,0,0,0,0};
//
//	for(size_t iSurface=0; iSurface<ASurfaces.size(); iSurface++)
//	{
//		std::vector<GEPETO::Triangle<3,TBase> > triangles;
//
//		FacetedSurface<TBase>* s = dynamic_cast<FacetedSurface<TBase>*>(ASurfaces[iSurface]);
//		if(s==0)
//			throw GMDSException("Wrong geometric type");
//		else
//			s->getTriangulation(triangles);
//
//		for(size_t iTriangle=0; iTriangle<triangles.size(); iTriangle++)
//		{
//			GEPETO::Point<3, GEPETO::NumericEps<double,Epsilon> > pt[3];
//
//			pt[0] = triangles[iTriangle].getPoint(0);
//			pt[1] = triangles[iTriangle].getPoint(1);
//			pt[2] = triangles[iTriangle].getPoint(2);
//
//			GEPETO::Triangle<3,GEPETO::NumericEps<double,Epsilon> > triangle_tmp(pt[0],pt[1],pt[2]);
//
//			if(GEPETO::GEOM_YES == triangle_tmp.intersect(r1))
//				nbIntersections[0]++;
//			if(GEPETO::GEOM_YES == triangle_tmp.intersect(r2))
//				nbIntersections[1]++;
//			if(GEPETO::GEOM_YES == triangle_tmp.intersect(r3))
//				nbIntersections[2]++;
//			if(GEPETO::GEOM_YES == triangle_tmp.intersect(r4))
//				nbIntersections[3]++;
//			if(GEPETO::GEOM_YES == triangle_tmp.intersect(r5))
//				nbIntersections[4]++;
//			if(GEPETO::GEOM_YES == triangle_tmp.intersect(r6))
//				nbIntersections[5]++;
//			if(GEPETO::GEOM_YES == triangle_tmp.intersect(r7))
//				nbIntersections[6]++;
//
//		}
//	}
//
//	bool isInner[7];
//
//	for(size_t iLine=0; iLine<7; iLine++)
//	{
//		isInner[iLine] = ((nbIntersections[iLine]%2) == 1);
//	}
//
////	if((isInner[0] != isInner[1]) ||
////	   (isInner[0] != isInner[2]) ||
////	   (isInner[0] != isInner[3]) ||
////	   (isInner[0] != isInner[4]) ||
////	   (isInner[0] != isInner[5]) ||
////	   (isInner[0] != isInner[6]))
////		throw GMDSException("Problem in isInside detection!");
////	return isInner[0];
//
//	int nbIn = 0;
//	for(int i=0; i<7; i++)
//	{
//		if(isInner[i])
//			nbIn++;
//	}
//	return (nbIn<3)?false:true;
//
//}
/*----------------------------------------------------------------------------*/
bool FacetedMeshIntersectionService::
isMainlyInsideMC(
		std::vector<gmds::math::Triangle>& ATriangles,
		gmds::math::Hexahedron& AHex,
		double& ARatio) const
{
	// we get the 8 corners of the hex
	gmds::math::Point p[8];
	for(int iPoint=0; iPoint<8; iPoint++)
	{
		p[iPoint] = AHex.getPoint(iPoint);
	}

	// Now we build an inner grid of  27 points
	gmds::math::Point g[27];

	// we build intermediate nodes to create our 27 inner points

//	// middle of edge
	gmds::math::Point p01, p12, p23, p03,
		p45, p56, p67, p47, p04, p15, p26, p37;
//
	p01 = (p[0]+p[1])*0.5;
	p12 = (p[1]+p[2])*0.5;
	p23 = (p[2]+p[3])*0.5;
	p03 = (p[0]+p[3])*0.5;

	p45 = (p[4]+p[5])*0.5;
	p56 = (p[5]+p[6])*0.5;
	p67 = (p[6]+p[7])*0.5;
	p47 = (p[4]+p[7])*0.5;

	p04 = (p[0]+p[4])*0.5;
	p15 = (p[1]+p[5])*0.5;
	p26 = (p[2]+p[6])*0.5;
	p37 = (p[3]+p[7])*0.5;

	//face center
	gmds::math::Point p0123, p4567, p0154, p1265, p2376, p0374;

	p0123 = (p[0]+p[1]+p[2]+p[3])*0.25;
	p4567 = (p[4]+p[5]+p[6]+p[7])*0.25;
	p0154 = (p[0]+p[1]+p[5]+p[4])*0.25;
	p1265 = (p[1]+p[2]+p[6]+p[5])*0.25;
	p2376 = (p[2]+p[3]+p[7]+p[6])*0.25;
	p0374 = (p[0]+p[3]+p[7]+p[4])*0.25;

	//hex center
	gmds::math::Point bary = (p[0]+p[1]+p[2]+p[3]+p[4]+p[5]+p[6]+p[7])*0.125;

	g[0] = bary;

	g[1] = (bary+p[0])*0.5;
	g[2] = (bary+p[1])*0.5;
	g[3] = (bary+p[2])*0.5;
	g[4] = (bary+p[3])*0.5;
	g[5] = (bary+p[4])*0.5;
	g[6] = (bary+p[5])*0.5;
	g[7] = (bary+p[6])*0.5;
	g[8] = (bary+p[7])*0.5;

	g[9] = (bary+p01)*0.5;
	g[10]= (bary+p12)*0.5;
	g[11]= (bary+p23)*0.5;
	g[12]= (bary+p03)*0.5;

	g[13]= (bary+p45)*0.5;
	g[14]= (bary+p56)*0.5;
	g[15]= (bary+p67)*0.5;
	g[16]= (bary+p47)*0.5;

	g[17]= (bary+p04)*0.5;
	g[18]= (bary+p15)*0.5;
	g[19]= (bary+p26)*0.5;
	g[20]= (bary+p37)*0.5;


	g[21]= (bary+p0123)*0.5;
	g[22]= (bary+p4567)*0.5;
	g[23]= (bary+p0154)*0.5;
	g[24]= (bary+p1265)*0.5;
	g[25]= (bary+p2376)*0.5;
	g[26]= (bary+p0374)*0.5;

	int nbIn = 0;
	for(int iPoint=0; iPoint<27;iPoint++){
		gmds::math::Point currentPoint = g[iPoint];

		TCoord distance_max = HUGE_VALF;
		gmds::math::Point p_max;
		size_t iTriangle_max;

		for(size_t iTriangle=0; iTriangle<ATriangles.size(); iTriangle++)
		{
			TCoord distance2 = ATriangles[iTriangle].distance2(currentPoint);
			if(distance2<distance_max) {
				distance_max = distance2;
				p_max = ATriangles[iTriangle].project(currentPoint);
				iTriangle_max = iTriangle;
			}
		}

		gmds::math::Vector triNormal = ATriangles[iTriangle_max].getNormal();
		gmds::math::Vector vectProject(p_max,currentPoint);
		TCoord dotproduct = triNormal.dot(vectProject);

		if(dotproduct<=0.){
			nbIn++;
		}

	}
//	std::cout<<"NB IN = "<<nbIn<<std::endl;
	ARatio = nbIn/27.;
	return (nbIn>15);
}
/*----------------------------------------------------------------------------*/
bool FacetedMeshIntersectionService::
isMainlyInsideMC(
		std::vector<gmds::math::Triangle>& ATriangles,
		gmds::Region& ARegion,
		double& ARatio) const
{
	// we compute interpolation points that are inside the region
	std::vector<gmds::math::Point> ngllPoints = ARegion.computeNGLLPoints(0);

	int nbIn = 0;
	for(int iPoint=0; iPoint<ngllPoints.size(); iPoint++){
		gmds::math::Point currentPoint = ngllPoints[iPoint];

		TCoord distance_max = HUGE_VALF;
		gmds::math::Point p_max;
		size_t iTriangle_max;

		for(size_t iTriangle=0; iTriangle<ATriangles.size(); iTriangle++)
		{
			TCoord distance2 = ATriangles[iTriangle].distance2(currentPoint);
			if(distance2<distance_max) {
				distance_max = distance2;
				p_max = ATriangles[iTriangle].project(currentPoint);
				iTriangle_max = iTriangle;
			}
		}

		gmds::math::Vector triNormal = ATriangles[iTriangle_max].getNormal();
		gmds::math::Vector vectProject(p_max,currentPoint);
		TCoord dotproduct = triNormal.dot(vectProject);

		if(dotproduct<=0.){
			nbIn++;
		}

	}
//	std::cout<<"NB IN = "<<nbIn<<std::endl;
	ARatio = (double) nbIn / (double) ngllPoints.size();
	return (ARatio > 0.5);
}
/*----------------------------------------------------------------------------*/
bool FacetedMeshIntersectionService::
isInside(
		std::vector<gmds::math::Triangle>& ATriangles,
		gmds::math::Point& APoint) const
{

	TCoord distance_max = HUGE_VALF;
	gmds::math::Point p_max;
	size_t iTriangle_max;

	for(size_t iTriangle=0; iTriangle<ATriangles.size(); iTriangle++)
	{
		TCoord distance2 = ATriangles[iTriangle].distance2(APoint);
		if(distance2<distance_max) {
			distance_max = distance2;
			p_max = ATriangles[iTriangle].project(APoint);
			iTriangle_max = iTriangle;
		}
	}

	gmds::math::Vector triNormal = ATriangles[iTriangle_max].getNormal();
	gmds::math::Vector vectProject(p_max,APoint);
	TCoord dotproduct = triNormal.dot(vectProject);

	if(dotproduct<=0.) {
		return true;
	} else {
		return false;
	}

}
///*----------------------------------------------------------------------------*/
//template<typename TBase >
//bool FacetedMeshIntersectionService<TBase>::
//isInsideMC(
//		std::vector<GeomSurface<TBase>* >& ASurfaces,
//		std::vector<GEPETO::Vector<3,TBase > >& AVec,
//		GEPETO::Point<3, TBase >  & APoint) const
//{
//	const int Epsilon = 12;
//	GEPETO::Point<3,GEPETO::NumericEps<double,Epsilon>  > pnt(APoint.getX(),APoint.getY(),APoint.getZ());
//	std::vector<GEPETO::Vector<3,GEPETO::NumericEps<double,Epsilon> > > vec;
//	for(unsigned int i=0;i<AVec.size();i++){
//		GEPETO::Vector<3,TBase > vinit = AVec[i];
//		vec.push_back(GEPETO::Vector<3,GEPETO::NumericEps<double,Epsilon> >(vinit.get(0),vinit.get(1),vinit.get(2)));
//	}
//	size_t nbInner = 0;
//	int nbRays = AVec.size();
//
//	std::vector<GEPETO::Ray<3,GEPETO::NumericEps<double,Epsilon> > > rays;
//	size_t nbIntersections[nbRays];
//	for(size_t iRay=0; iRay<nbRays; iRay++)
//	{
//		rays.push_back(GEPETO::Ray<3,GEPETO::NumericEps<double,Epsilon>  >(pnt,vec[iRay]));
//		nbIntersections[iRay] = 0;
//	}
//
//	for(size_t iSurface=0; iSurface<ASurfaces.size(); iSurface++)
//	{
//		std::vector<GEPETO::Triangle<3,TBase> > triangles;
//
////		FacetedSurface<TBase>* s = dynamic_cast<FacetedSurface<TBase>*>(ASurfaces[iSurface]);
////		if(s==0)
////			throw GMDSException("Wrong geometric type");
////		else
//		ASurfaces[iSurface]->getTriangulation(triangles);
//
//		for(size_t iTriangle=0; iTriangle<triangles.size(); iTriangle++)
//		{
//			GEPETO::Point<3, TBase > pt[3];
//
//			pt[0] = triangles[iTriangle].getPoint(0);
//			pt[1] = triangles[iTriangle].getPoint(1);
//			pt[2] = triangles[iTriangle].getPoint(2);
//
//			GEPETO::Triangle<3,GEPETO::NumericEps<double,Epsilon> > triangle_tmp(pt[0],pt[1],pt[2]);
//
//			for(size_t iRay=0; iRay<rays.size(); iRay++)
//			{
//				if(GEPETO::GEOM_YES == triangle_tmp.intersect(rays[iRay]) ||
//						GEPETO::GEOM_UNDEF == triangle_tmp.intersect(rays[iRay]))
//					nbIntersections[iRay]++;
//			}
//		}
//	}
//
//	for(size_t iRay=0; iRay<rays.size(); iRay++)
//	{
//		if(((nbIntersections[iRay]%2) == 1))
//		{
//			nbInner++;
//		}
//	}
//
//	//td::cout<<"number of inside rays "<<nbInner<<" of "<<rays.size()<<std::endl;
//
//	return (nbInner<(nbRays/2))?false:true;
//}
/*----------------------------------------------------------------------------*/
bool FacetedMeshIntersectionService::
isInsideMC(std::vector<GeomSurface* >& ASurfaces, gmds::math::Hexahedron& AHex) const
{
	// number of rays emitted : the more the better
	const int nbRays = FacetedMeshIntersectionService_isInsideMC_NBDRAWS;

	// get the center of the hexahedron
	gmds::math::Point p0(AHex.getCenter());

	size_t nbInner = 0;

	// create several rays at random starting from the center of the hexahedron
	// and check the number of intersections with the surface

	this->random_init();

	std::vector<gmds::math::Ray> rays;

	for(size_t iRay=0; iRay<nbRays; iRay++)
	{
		// create random vector
		TCoord xyz_rand[3];
		xyz_rand[0] = this->random_value();
		xyz_rand[1] = this->random_value();
		xyz_rand[2] = this->random_value();
		gmds::math::Vector random_vector(xyz_rand[0],xyz_rand[1],xyz_rand[2]);

		// create ray
		if(!random_vector.isZero())
		{
			rays.push_back(gmds::math::Ray(p0,random_vector));
		}
	}

	// compute number of intersections
	size_t nbIntersections[nbRays];
	for(size_t iRay=0; iRay<nbRays; iRay++)
	{
		nbIntersections[iRay] = 0;
	}

	for(size_t iSurface=0; iSurface<ASurfaces.size(); iSurface++)
	{
		for(size_t iTriangle=0; iTriangle<m_surfacesTriangulation[ASurfaces[iSurface]].size(); iTriangle++)
		{
			for(size_t iRay=0; iRay<rays.size(); iRay++)
			{
				if(m_surfacesTriangulation[ASurfaces[iSurface]][iTriangle].intersect(rays[iRay]))
					nbIntersections[iRay]++;
			}
		}
	}

	for(size_t iRay=0; iRay<rays.size(); iRay++)
	{
		if(((nbIntersections[iRay]%2) == 1))
		{
			nbInner++;
		}
	}

	std::cout<<"number of inside rays "<<nbInner<<" of "<<rays.size()<<std::endl;

	return (nbInner<(nbRays/2))?false:true;

}
/*----------------------------------------------------------------------------*/
bool FacetedMeshIntersectionService::
isInsideMC(std::vector<GeomSurface* >& ASurfaces, gmds::Region& ARegion) const
{
	// number of rays emitted : the more the better
	const int nbRays = FacetedMeshIntersectionService_isInsideMC_NBDRAWS;

	// get the center of the hexahedron
	gmds::math::Point p0(ARegion.center());

	size_t nbInner = 0;

	// create several rays at random starting from the center of the hexahedron
	// and check the number of intersections with the surface

	this->random_init();

	std::vector<gmds::math::Ray> rays;

	for(size_t iRay=0; iRay<nbRays; iRay++)
	{
		// create random vector
		TCoord xyz_rand[3];
		xyz_rand[0] = this->random_value();
		xyz_rand[1] = this->random_value();
		xyz_rand[2] = this->random_value();
		gmds::math::Vector random_vector(xyz_rand[0],xyz_rand[1],xyz_rand[2]);

		// create ray
		if(!random_vector.isZero())
		{
			rays.push_back(gmds::math::Ray(p0,random_vector));
		}
	}

	// compute number of intersections
	size_t nbIntersections[nbRays];
	for(size_t iRay=0; iRay<nbRays; iRay++)
	{
		nbIntersections[iRay] = 0;
	}

	for(size_t iSurface=0; iSurface<ASurfaces.size(); iSurface++)
	{
		for(size_t iTriangle=0; iTriangle<m_surfacesTriangulation[ASurfaces[iSurface]].size(); iTriangle++)
		{
			for(size_t iRay=0; iRay<rays.size(); iRay++)
			{
				if(m_surfacesTriangulation[ASurfaces[iSurface]][iTriangle].intersect(rays[iRay]))
					nbIntersections[iRay]++;
			}
		}
	}

	for(size_t iRay=0; iRay<rays.size(); iRay++)
	{
		if(((nbIntersections[iRay]%2) == 1))
		{
			nbInner++;
		}
	}

	std::cout<<"number of inside rays "<<nbInner<<" of "<<rays.size()<<std::endl;

	return (nbInner<(nbRays/2))?false:true;

}
/*----------------------------------------------------------------------------*/
//template<typename TBase >
//GEPETO::Point<3,TBase>
//FacetedMeshIntersectionService<TBase>::
//project(const std::vector<GeomSurface<TBase>* >& ASurfaces,
//		const GEPETO::Point<3,TBase>& APoint) const
//{
//	if(ASurfaces.size() == 0)
//	{
//		throw GMDSException("Cannot project when there are no surfaces");
//	}
//
//	FacetedSurface<TBase>* s = dynamic_cast<FacetedSurface<TBase>*>(ASurfaces[0]);
//	if(s==0)
//		throw GMDSException("Wrong geometric surface type");
//
//	std::vector<GEPETO::Triangle<3,TBase> > triangles;
//	s->getTriangulation(triangles);
//
//	if(triangles.size() == 0)
//		throw GMDSException("No triangle in surface.");
//
//	GEPETO::Point<3,TBase> P = triangles[0].project(APoint);
//
//	// we compute projection on all the triangles of the surface and keep the
//	// nearest
//	for(size_t iSurface=0; iSurface<ASurfaces.size(); iSurface++)
//	{
//		std::vector<GEPETO::Triangle<3,TBase> > triangles;
//
//		FacetedSurface<TBase>* s = dynamic_cast<FacetedSurface<TBase>*>(ASurfaces[iSurface]);
//		if(s==0)
//			throw GMDSException("Wrong geometric type");
//		else
//			s->getTriangulation(triangles);
//
//		for(size_t iTriangle=0; iTriangle<triangles.size(); iTriangle++)
//		{
//			GEPETO::Point<3,TBase> P_tmp = triangles[iTriangle].project(APoint);
//
//			if(APoint.distance(P) > APoint.distance(P_tmp))
//			{
//				P = P_tmp;
//			}
//		}
//	}
//
//	return P;
//}
//
///*----------------------------------------------------------------------------*/
//template<typename TBase >
//GEPETO::Point<3,TBase>
//FacetedMeshIntersectionService<TBase>::
//project(const std::vector<GeomSurface<TBase>* >& ASurfaces,
//		const GEPETO::Point<3,TBase>& APoint,
//		int& ASurfaceDestID) const
//{
//	if(ASurfaces.size() == 0)
//	{
//		throw GMDSException("Cannot project when there are no surfaces");
//	}
//
//	FacetedSurface<TBase>* s = dynamic_cast<FacetedSurface<TBase>*>(ASurfaces[0]);
//	if(s==0)
//		throw GMDSException("Wrong geometric surface type");
//
//	std::vector<GEPETO::Triangle<3,TBase> > triangles;
//	s->getTriangulation(triangles);
//
//	if(triangles.size() == 0)
//		throw GMDSException("No triangle in surface.");
//
//	GEPETO::Point<3,TBase> P = triangles[0].project(APoint);
//	ASurfaceDestID = s->getId();
//
//	// we compute projection on all the triangles of the surface and keep the
//	// nearest
//	for(size_t iSurface=0; iSurface<ASurfaces.size(); iSurface++)
//	{
//		std::vector<GEPETO::Triangle<3,TBase> > triangles;
//
//		FacetedSurface<TBase>* s = dynamic_cast<FacetedSurface<TBase>*>(ASurfaces[iSurface]);
//		if(s==0)
//			throw GMDSException("Wrong geometric type");
//		else
//			s->getTriangulation(triangles);
//
//		for(size_t iTriangle=0; iTriangle<triangles.size(); iTriangle++)
//		{
//			GEPETO::Point<3,TBase> P_tmp = triangles[iTriangle].project(APoint);
//
//			if(APoint.distance(P) > APoint.distance(P_tmp))
//			{
//				P = P_tmp;
//				ASurfaceDestID = s->getId();
//			}
//		}
//	}
//
//	return P;
//}
//
///*----------------------------------------------------------------------------*/
//template<typename TBase >
//GEPETO::Point<3,TBase>
//FacetedMeshIntersectionService<TBase>::
//project(GeomCurve<TBase>& ACurve,
//		const GEPETO::Point<3,TBase>& APoint) const
//{
//	FacetedCurve<TBase>* c = dynamic_cast<FacetedCurve<TBase>* >(&ACurve);
//	if(c==0)
//		throw GMDSException("Wrong geometric curve type");
//
//
//	GEPETO::Point<3,TBase> P(APoint);
//
//	c->project(P);
//
//	return P;
//}
//
///*----------------------------------------------------------------------------*/
//template<typename TBase >
//GEPETO::Point<3,TBase>
//FacetedMeshIntersectionService<TBase>::
//project(const std::vector<GeomSurface<TBase>* >& ASurfaces,
//		const GEPETO::Point<3,TBase>& APoint,
//		const GEPETO::Vector<3,TBase>& AVector) const
//{
//	if(ASurfaces.size() == 0)
//	{
//		throw GMDSException("Cannot project when there are no surfaces");
//	}
//
//	FacetedSurface<TBase>* s = dynamic_cast<FacetedSurface<TBase>*>(ASurfaces[0]);
//	if(s==0)
//		throw GMDSException("Wrong geometric surface type");
//
//	std::vector<GEPETO::Triangle<3,TBase> > triangles;
//	s->getTriangulation(triangles);
//
//	if(triangles.size() == 0)
//		throw GMDSException("No triangle in surface.");
//
//	GEPETO::Point<3,TBase> P(0.,0.,0.);
//	bool found = false;
//
////	// we compute projection on all the triangles of the surface and keep the
////	// nearest
////	for(size_t iSurface=0; iSurface<ASurfaces.size(); iSurface++)
////	{
////		std::vector<GEPETO::Triangle<3,TBase> > triangles;
////
////		FacetedSurface<TBase>* s = dynamic_cast<FacetedSurface<TBase>*>(ASurfaces[iSurface]);
////		if(s==0)
////			throw GMDSException("Wrong geometric type");
////		else
////			s->getTriangulation(triangles);
////
////		for(size_t iTriangle=0; iTriangle<triangles.size(); iTriangle++)
////		{
////			GEPETO::Ray<3,TBase> ray(APoint,AVector);
////			GEPETO::EGeomPredicateResult result = GEPETO::GeomToolkit<3,TBase>::intersect(ray,triangles[iTriangle]);
////			if(result != GEPETO::GEOM_NO)
////			{
////				if(!found)
////				{
////					P = GEPETO::GeomToolkit<3,TBase>::computeIntersection(ray,triangles[iTriangle]);
////					found = true;
////				}
////
////				GEPETO::Point<3,TBase> P_tmp = GEPETO::GeomToolkit<3,TBase>::computeIntersection(ray,triangles[iTriangle]);
////
////				if(GEPETO::GeomToolkit<3,TBase>::distance(P,APoint) > GEPETO::GeomToolkit<3,TBase>::distance(P_tmp,APoint))
////				{
////					P = P_tmp;
////				}
////			}
////		}
////	}
//
//	const int Epsilon = FacetedMeshIntersectionService_project_EPS;
//
//	GEPETO::Point<3, GEPETO::NumericEps<double,Epsilon> > APoint_prime(APoint.getX().toDouble(),APoint.getY().toDouble(),APoint.getZ().toDouble());
//	GEPETO::Vector<3, GEPETO::NumericEps<double,Epsilon> > AVector_prime(AVector.get(0).toDouble(),AVector.get(1).toDouble(),AVector.get(2).toDouble());
//	GEPETO::Point<3, GEPETO::NumericEps<double,Epsilon> > P_prime(0.,0.,0.);
//
//	// we compute projection on all the triangles of the surface and keep the
//	// nearest
//	for(size_t iSurface=0; iSurface<ASurfaces.size(); iSurface++)
//	{
//		std::vector<GEPETO::Triangle<3,TBase> > triangles;
//
//		FacetedSurface<TBase>* s = dynamic_cast<FacetedSurface<TBase>*>(ASurfaces[iSurface]);
//		if(s==0)
//			throw GMDSException("Wrong geometric type");
//		else
//			s->getTriangulation(triangles);
//
//		for(size_t iTriangle=0; iTriangle<triangles.size(); iTriangle++)
//		{
//			GEPETO::Point<3, GEPETO::NumericEps<double,Epsilon> > pt[3];
//
//			pt[0] = triangles[iTriangle].getPoint(0);
//			pt[1] = triangles[iTriangle].getPoint(1);
//			pt[2] = triangles[iTriangle].getPoint(2);
//
//			GEPETO::Triangle<3,GEPETO::NumericEps<double,Epsilon> > triangle_tmp(pt[0],pt[1],pt[2]);
//
//			GEPETO::Ray<3,GEPETO::NumericEps<double,Epsilon> > ray(APoint_prime,AVector_prime);
//
//			GEPETO::EGeomPredicateResult result = triangle_tmp.intersect(ray);
//			if(result != GEPETO::GEOM_NO)
//			{
//				if(!found)
//				{
//					triangle_tmp.intersect(ray,P_prime);
//					found = true;
//				}
//
//				GEPETO::Point<3,GEPETO::NumericEps<double,Epsilon> > P_tmp;
//				triangle_tmp.intersect(ray,P_tmp);
//
//				if(APoint_prime.distance(P_prime) > APoint_prime.distance(P_tmp))
//				{
//					P_prime = P_tmp;
//				}
//			}
//		}
//	}
//
//	if(!found)
//	{
//		throw GMDSException("the projection along a vector could not be done");
//	}
//
//	P.set(0,P_prime.getX().toDouble());
//	P.set(1,P_prime.getY().toDouble());
//	P.set(2,P_prime.getZ().toDouble());
//
//	return P;
//}
//
///*----------------------------------------------------------------------------*/
//template<typename TBase >
//GEPETO::Point<3,TBase>
//FacetedMeshIntersectionService<TBase>::
//project(const std::vector<GeomSurface<TBase>* >& ASurfaces,
//		const GEPETO::Point<3,TBase>& APoint,
//		const GEPETO::Vector<3,TBase>& AVector,
//		int& ASurfaceDestID) const
//{
//	if(ASurfaces.size() == 0)
//	{
//		throw GMDSException("Cannot project when there are no surfaces");
//	}
//
//	FacetedSurface<TBase>* s = dynamic_cast<FacetedSurface<TBase>*>(ASurfaces[0]);
//	if(s==0)
//		throw GMDSException("Wrong geometric surface type");
//
//	std::vector<GEPETO::Triangle<3,TBase> > triangles;
//	s->getTriangulation(triangles);
//
//	if(triangles.size() == 0)
//		throw GMDSException("No triangle in surface.");
//
//	GEPETO::Point<3,TBase> P(0.,0.,0.);
//	bool found = false;
//
////	// we compute projection on all the triangles of the surface and keep the
////	// nearest
////	for(size_t iSurface=0; iSurface<ASurfaces.size(); iSurface++)
////	{
////		std::vector<GEPETO::Triangle<3,TBase> > triangles;
////
////		FacetedSurface<TBase>* s = dynamic_cast<FacetedSurface<TBase>*>(ASurfaces[iSurface]);
////		if(s==0)
////			throw GMDSException("Wrong geometric type");
////		else
////			s->getTriangulation(triangles);
////
////		for(size_t iTriangle=0; iTriangle<triangles.size(); iTriangle++)
////		{
////			GEPETO::Ray<3,TBase> ray(APoint,AVector);
////			GEPETO::EGeomPredicateResult result = GEPETO::GeomToolkit<3,TBase>::intersect(ray,triangles[iTriangle]);
////			if(result != GEPETO::GEOM_NO)
////			{
////				if(!found)
////				{
////					P = GEPETO::GeomToolkit<3,TBase>::computeIntersection(ray,triangles[iTriangle]);
////					found = true;
////				}
////
////				GEPETO::Point<3,TBase> P_tmp = GEPETO::GeomToolkit<3,TBase>::computeIntersection(ray,triangles[iTriangle]);
////
////				if(GEPETO::GeomToolkit<3,TBase>::distance(P,APoint) > GEPETO::GeomToolkit<3,TBase>::distance(P_tmp,APoint))
////				{
////					P = P_tmp;
////				}
////			}
////		}
////	}
//
//	const int Epsilon = FacetedMeshIntersectionService_project_EPS;
//
//	GEPETO::Point<3, GEPETO::NumericEps<double,Epsilon> > APoint_prime(APoint.getX().toDouble(),APoint.getY().toDouble(),APoint.getZ().toDouble());
//	GEPETO::Vector<3, GEPETO::NumericEps<double,Epsilon> > AVector_prime(AVector.get(0).toDouble(),AVector.get(1).toDouble(),AVector.get(2).toDouble());
//	GEPETO::Point<3, GEPETO::NumericEps<double,Epsilon> > P_prime(0.,0.,0.);
//
//	//AAAAAAAAAAAA
//	std::cout<<"ray"<<std::endl;
//	std::cout<<APoint_prime<<std::endl;
//	std::cout<<AVector_prime<<std::endl;
//	std::cout<<APoint_prime.getX()<<" "<<APoint_prime.getY()<<" "<<APoint_prime.getZ()<<std::endl;
//	std::cout<<(APoint_prime+3.*AVector_prime).getX()<<" "<<(APoint_prime+3.*AVector_prime).getY()<<" "<<(APoint_prime+3.*AVector_prime).getZ()<<std::endl;
//
//	//AAAAAAAAAAAA
//
//	// we compute projection on all the triangles of the surface and keep the
//	// nearest
//	for(size_t iSurface=0; iSurface<ASurfaces.size(); iSurface++)
//	{
//		std::vector<GEPETO::Triangle<3,TBase> > triangles;
//
//		FacetedSurface<TBase>* s = dynamic_cast<FacetedSurface<TBase>*>(ASurfaces[iSurface]);
//		if(s==0)
//			throw GMDSException("Wrong geometric type");
//		else
//			s->getTriangulation(triangles);
//
//		for(size_t iTriangle=0; iTriangle<triangles.size(); iTriangle++)
//		{
//			GEPETO::Point<3, GEPETO::NumericEps<double,Epsilon> > pt[3];
//
//			pt[0] = triangles[iTriangle].getPoint(0);
//			pt[1] = triangles[iTriangle].getPoint(1);
//			pt[2] = triangles[iTriangle].getPoint(2);
//
//			GEPETO::Triangle<3,GEPETO::NumericEps<double,Epsilon> > triangle_tmp(pt[0],pt[1],pt[2]);
//
//			GEPETO::Ray<3,GEPETO::NumericEps<double,Epsilon> > ray(APoint_prime,AVector_prime);
//
//			GEPETO::EGeomPredicateResult result;
//			try{
//				result = triangle_tmp.intersect(ray);
//			}
//			catch(GEPETO::GepetoException& e){
//				result = GEPETO::GEOM_NO;
//			}
//
//			if(result != GEPETO::GEOM_NO)
//			{
//				if(!found)
//				{
//					triangle_tmp.intersect(ray,P_prime);
//					found = true;
//					ASurfaceDestID = s->getId();
//				}
//
//				GEPETO::Point<3,GEPETO::NumericEps<double,Epsilon> > P_tmp;
//				triangle_tmp.intersect(ray,P_tmp);
//
//				if(APoint_prime.distance(P_prime) > APoint_prime.distance(P_tmp))
//				{
//					P_prime = P_tmp;
//					ASurfaceDestID = s->getId();
//				}
//			}
//
//			// AAAAAAAAA
//			//std::cout<<"triangle"<<std::endl;
//			if(result == GEPETO::GEOM_NO)
//			{
//				if(triangle_tmp.getSphereIncluding().intersect(ray,GEPETO::EGPP_IN_PEEL))
//				{
//					//std::cout<<triangle_tmp<<std::endl;
//					//std::cout<<triangle_tmp.getPoint(0)<<std::endl;
//					//std::cout<<triangle_tmp.getPoint(1)<<std::endl;
//					//std::cout<<triangle_tmp.getPoint(2)<<std::endl;
////					std::cout<<triangle_tmp.getPoint(0).getX()<<" "<<triangle_tmp.getPoint(0).getY()<<" "<<triangle_tmp.getPoint(0).getZ()<<std::endl;
////					std::cout<<triangle_tmp.getPoint(1).getX()<<" "<<triangle_tmp.getPoint(1).getY()<<" "<<triangle_tmp.getPoint(1).getZ()<<std::endl;
////					std::cout<<triangle_tmp.getPoint(2).getX()<<" "<<triangle_tmp.getPoint(2).getY()<<" "<<triangle_tmp.getPoint(2).getZ()<<std::endl;
////					std::cout<<triangle_tmp.getPoint(0).getX()<<" "<<triangle_tmp.getPoint(0).getY()<<" "<<triangle_tmp.getPoint(0).getZ()<<std::endl;
////					std::cout<<std::endl;
////					std::cout<<std::endl;
//				}
//			}
//			// AAAAAAAAA
//
//		}
//	}
//
//	if(!found)
//	{
//		//throw GMDSException("the projection along a vector could not be done");
//	}
//
//	P.set(0,P_prime.getX().toDouble());
//	P.set(1,P_prime.getY().toDouble());
//	P.set(2,P_prime.getZ().toDouble());
//
//	return P;
//}
/*----------------------------------------------------------------------------*/
bool
FacetedMeshIntersectionService::intersect(gmds::math::Hexahedron& AHex, gmds::math::Triangle& ATri) const
{
	// first check that the hexahedron intersects the plane of the triangle
	//gmds::math::Plane plane(ATri);

	throw GMDSException("FacetedMeshIntersectionService::intersect(gmds::math::Hexahedron& AHex, gmds::math::Triangle& ATri)" 
			"not implemented yet.");	
};
/*----------------------------------------------------------------------------*/
GNode*
FacetedMeshIntersectionService::buildAABBSurfacesTriangulationTree(
		std::vector<GeomSurface*>& ASurfaces,
		std::map<GeomSurface*,GNode*>& ASurfacesTriangulationTrees)
{
	GSList* boxList = NULL;

	for(unsigned int iSurface=0; iSurface<ASurfaces.size(); iSurface++) {

		std::vector<gmds::math::Triangle>& surfaceTriangulation = m_surfacesTriangulation[ASurfaces[iSurface]];

		GSList* boxList_local = NULL;

		for(unsigned int iTriangle=0; iTriangle<surfaceTriangulation.size(); iTriangle++) {
			double minXYZ[3];
			double maxXYZ[3];

			surfaceTriangulation[iTriangle].computeBoundingBox(minXYZ,maxXYZ);

			gpointer pointer = &(m_surfacesTriangulation[ASurfaces[iSurface]][iTriangle]);
			GtsBBox* bbox = gts_bbox_new(
									gts_bbox_class (),
									pointer,
									minXYZ[0],minXYZ[1],minXYZ[2],
									maxXYZ[0],maxXYZ[1],maxXYZ[2]);

			boxList = g_slist_prepend(boxList,bbox);
			boxList_local = g_slist_prepend(boxList_local,bbox);
		}
		GNode* boxTree_local = gts_bb_tree_new(boxList_local);
		if(boxTree_local == NULL) {
			throw GMDSException("FacetedMeshIntersectionService::buildAABBSurfacesTriangulationTree : failed to build local tree");
		}
		ASurfacesTriangulationTrees[ASurfaces[iSurface]] = boxTree_local;
	}


	std::cout<<"boxList "<<g_slist_length(boxList)<<std::endl;
	GNode* boxTree = gts_bb_tree_new(boxList);

	if(boxTree == NULL) {
		throw GMDSException("FacetedMeshIntersectionService::buildAABBSurfacesTriangulationTree : failed to build tree");
	}

	return boxTree;
}
/*----------------------------------------------------------------------------*/
void
FacetedMeshIntersectionService::buildAABBSurfacesTriangulationTree(
                std::vector<GeomSurface*>& ASurfaces)
{
        GSList* boxList = NULL;

        for(unsigned int iSurface=0; iSurface<ASurfaces.size(); iSurface++) {

                std::vector<gmds::math::Triangle>& surfaceTriangulation = m_surfacesTriangulation[ASurfaces[iSurface]];

                GSList* boxList_local = NULL;

                for(unsigned int iTriangle=0; iTriangle<surfaceTriangulation.size(); iTriangle++) {
                        double minXYZ[3];
                        double maxXYZ[3];

                        surfaceTriangulation[iTriangle].computeBoundingBox(minXYZ,maxXYZ);

                        gpointer pointer = &(m_surfacesTriangulation[ASurfaces[iSurface]][iTriangle]);
                        GtsBBox* bbox = gts_bbox_new(
                                                                        gts_bbox_class (),
                                                                        pointer,
                                                                        minXYZ[0],minXYZ[1],minXYZ[2],
                                                                        maxXYZ[0],maxXYZ[1],maxXYZ[2]);

                        boxList = g_slist_prepend(boxList,bbox);
                        boxList_local = g_slist_prepend(boxList_local,bbox);
                }
		GNode* boxTree_local = gts_bb_tree_new(boxList_local);
                if(boxTree_local == NULL) {
                        throw GMDSException("FacetedMeshIntersectionService::buildAABBSurfacesTriangulationTree : failed to build local tree");
                }
                this->aabbSurfacesTrianglesTrees_[ASurfaces[iSurface]] = boxTree_local;
        }


        std::cout<<"boxList "<<g_slist_length(boxList)<<std::endl;
        GNode* boxTree = gts_bb_tree_new(boxList);

        if(boxTree == NULL) {
                throw GMDSException("FacetedMeshIntersectionService::buildAABBSurfacesTriangulationTree : failed to build tree");
        }

        this->aabbSurfacesTrianglesTree_ = boxTree;
}
/*----------------------------------------------------------------------------*/
void
FacetedMeshIntersectionService::project(gmds::math::Point& APoint, GNode* ATree) const
{
	GtsPoint* p = gts_point_new(gts_point_class (),
				APoint.X(),APoint.Y(),APoint.Z());

	gdouble* distance = NULL;
	GtsPoint* newP = gts_bb_tree_point_closest(
					ATree,
                                        p,
                                        FacetedMeshIntersectionService_triangle_project,
                                        distance);
	
	APoint.setXYZ(newP->x,newP->y,newP->z);
}
/*----------------------------------------------------------------------------*/
void
FacetedMeshIntersectionService::project(
			gmds::geom::GeomSurface* ASurf,
			gmds::math::Point& APoint) const
{
	if((this->aabbSurfacesTrianglesTrees_).find(ASurf) == (this->aabbSurfacesTrianglesTrees_).end()) {
		throw GMDSException("FacetedMeshIntersectionService::project : surface not found in aabbSurfacesTrianglesTrees_. Probably missing init.");
	}

	this->project(APoint,(this->aabbSurfacesTrianglesTrees_).find(ASurf)->second);
}
/*----------------------------------------------------------------------------*/
void 
FacetedMeshIntersectionService::computeNormal(
                        const gmds::math::Point APoint,
                        gmds::geom::GeomEntity* AGeomEntity,
                        const bool AOutward,
                        gmds::math::Vector& ANormal) const
{
	ANormal = gmds::math::Vector(0.,0.,0.);

	if(AGeomEntity->getDim() == 3) {
		throw GMDSException("FacetedMeshIntersectionService::computeNormal cannot compute a normal to a volume");
	}

	if(AGeomEntity->getDim() == 2) {
		gmds::geom::GeomSurface* surface = dynamic_cast<gmds::geom::GeomSurface*> (AGeomEntity);
		surface->computeNormal(APoint,ANormal);

	} else if(AGeomEntity->getDim() == 1) {
		gmds::geom::GeomCurve* curve = dynamic_cast<gmds::geom::GeomCurve*> (AGeomEntity);
		std::vector<gmds::geom::GeomSurface*> surfaces;
		curve->get(surfaces);

		for(int iSurf=0; iSurf<surfaces.size(); iSurf++) {
			gmds::math::Vector normal;
			surfaces[iSurf]->computeNormal(APoint,normal);
			ANormal = ANormal + normal;
		}
	} else {
		gmds::geom::GeomPoint* point = dynamic_cast<gmds::geom::GeomPoint*> (AGeomEntity);
                std::vector<gmds::geom::GeomSurface*> surfaces; 
                point->get(surfaces);

                for(int iSurf=0; iSurf<surfaces.size(); iSurf++) {
                        gmds::math::Vector normal;
                        surfaces[iSurf]->computeNormal(APoint,normal);
                        ANormal = ANormal + normal;
                }
	}

	if(!AOutward) {
		ANormal = -1. * ANormal;
	}
	ANormal.normalize();
}
/*----------------------------------------------------------------------------*/
void
FacetedMeshIntersectionService::
random_init() const
{
  unsigned int seed;

  std::string input_file_name;
  input_file_name.append("/dev/random");
  std::ifstream input_file_stream(input_file_name.c_str(),std::ios::in);
  input_file_stream >> seed;
  input_file_stream.close();

  //seed = 2;
  srand(time(NULL) + seed);
};
/*----------------------------------------------------------------------------*/
double
FacetedMeshIntersectionService::
random_value() const
{
  return ((double)rand() / (double) RAND_MAX);
};
/*----------------------------------------------------------------------------*/
} // namespace geom
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
