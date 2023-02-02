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
/** \file    FacetedSurface.t.h
 *  \author  F. LEDOUX
 *  \date    30/05/2011
 */
/*----------------------------------------------------------------------------*/
#include <map>
#include <set>
/*----------------------------------------------------------------------------*/
#include "GMDS/CAD/FacetedSurface.h"
/*----------------------------------------------------------------------------*/
#include <GMDS/CAD/FacetedGeomManager.h>
#include "GMDS/CAD/FacetedVolume.h"
#include "GMDS/CAD/FacetedCurve.h"
#include "GMDS/CAD/FacetedPoint.h"
#include "GMDS/Math/Ray.h"
#include "GMDS/Math/Triangle.h"
#include "GMDS/Utils/RandomGenerator.h"
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace geom{
/*----------------------------------------------------------------------------*/
int FacetedSurface::next_id_=1;
/*----------------------------------------------------------------------------*/
FacetedSurface::FacetedSurface()
{}
/*----------------------------------------------------------------------------*/

FacetedSurface::
FacetedSurface(std::vector<FacetedPoint* >& AP,
					std::vector<FacetedCurve* >& AC,
					std::vector<Face>& ADiscret,
					const std::string& AName)
:GeomSurface(AName),id_(next_id_++),points_(AP),curves_(AC),mesh_representation_(ADiscret)
{
	for(unsigned int i=0;i<mesh_representation_.size();i++){
		std::vector<Node> nodes_fi;
		mesh_representation_[i].get<Node>(nodes_fi);
		if(nodes_fi.size()!=3)
			throw GMDSException("An error of surface discretisation occured in the FacetedModel generation process");

		math::Point p1(nodes_fi[0].getPoint());
		math::Point p2(nodes_fi[1].getPoint());
		math::Point p3(nodes_fi[2].getPoint());

//		discrete_representation_.push_back(GEPETO::Triangle<3,TBase>(p1,p2,p3));
	}
}
/*----------------------------------------------------------------------------*/

FacetedSurface::
~FacetedSurface()
{}
/*----------------------------------------------------------------------------*/

void FacetedSurface::add(FacetedVolume* AVol)
{
	volumes_.push_back(AVol);
}
/*----------------------------------------------------------------------------*/

void FacetedSurface::replace(std::vector<FacetedPoint* > APoints)
{
	points_.clear();
	points_ = APoints;
}
/*----------------------------------------------------------------------------*/

void FacetedSurface::replace(std::vector<FacetedCurve* > ACurves)
{
	curves_.clear();
	curves_ = ACurves;
}
/*----------------------------------------------------------------------------*/

void FacetedSurface::
get(std::vector<GeomPoint*>& points) const
{
	points.clear();
	points.resize(points_.size());
	for(unsigned int i=0;i<points_.size();i++){
		points[i]=points_[i];
	}
}
/*----------------------------------------------------------------------------*/

void FacetedSurface::
get(std::vector<GeomCurve*>& curves) const
{
	curves.clear();
	curves.resize(curves_.size());
	for(unsigned int i=0;i<curves_.size();i++){
		curves[i]=curves_[i];
	}
}
/*----------------------------------------------------------------------------*/

void FacetedSurface::
get(std::vector<GeomVolume*>& AVol) const
{
	AVol.clear();
	AVol.resize(curves_.size());
	for(unsigned int i=0;i<volumes_.size();i++){
		AVol[i]=volumes_[i];
	}
}
/*----------------------------------------------------------------------------*/
TCoord FacetedSurface::computeArea() const
{
	TCoord totalArea = 0.;

	for(unsigned int i=0;i<mesh_representation_.size();i++){

		totalArea +=  mesh_representation_[i].area();
	}
	return totalArea;
}
/*----------------------------------------------------------------------------*/

void FacetedSurface::
computeBoundingBox(TCoord minXYZ[3], TCoord maxXYZ[3]) const
{
	std::vector<math::Point > pnts;
	std::vector<math::Triangle > triangles;
	getTriangulation(triangles);
	for(unsigned int i=0; i<triangles.size(); i++) {
		pnts.push_back(triangles[i].getPoint(0));
		pnts.push_back(triangles[i].getPoint(1));
		pnts.push_back(triangles[i].getPoint(2));
	}

	// too many comparisons (3 factor)
	math::Point pi = pnts[0];
	minXYZ[0]=pi.X();
	maxXYZ[0]=pi.X();
	minXYZ[1]=pi.Y();
	maxXYZ[1]=pi.Y();
	minXYZ[2]=pi.Z();
	maxXYZ[2]=pi.Z();
	for(unsigned int i=1;i<pnts.size();i++){
		 pi = pnts[i];
		if (pi.X()<minXYZ[0])
			minXYZ[0]=pi.X();
		else if (pi.X()>maxXYZ[0])
			maxXYZ[0]=pi.X();

		if (pi.Y()<minXYZ[1])
			minXYZ[1]=pi.Y();
		else if (pi.Y()>maxXYZ[1])
			maxXYZ[1]=pi.Y();

		if (pi.Z()<minXYZ[2])
			minXYZ[2]=pi.Z();
		else if (pi.Z()>maxXYZ[2])
			maxXYZ[2]=pi.Z();
	}
}
/*----------------------------------------------------------------------------*/

void FacetedSurface::
computeNormal( const math::Point& AP,
			   math::Vector& AV) const
{
	Face f0 = mesh_representation_[0];
	TCoord min_dist = f0.distance(AP);
	TCellID index = 0;
	for(unsigned int i=1;i<mesh_representation_.size();i++){
		TCoord dist = mesh_representation_[i].distance(AP);
		if(dist<min_dist){
			min_dist=dist;
			index = i;
		}
	}
	AV = mesh_representation_[index].normal();
}
/*----------------------------------------------------------------------------*/
math::Point FacetedSurface::
closestPoint(const math::Point& AP) const
{
	Face f0 = mesh_representation_[0];
	TCoord min_dist = f0.distance(AP);
	TCellID index = 0;
	for(unsigned int i=1;i<mesh_representation_.size();i++){
		TCoord dist = mesh_representation_[i].distance(AP);
		if(dist<min_dist){
			min_dist=dist;
			index = i;
		}
	}
	return mesh_representation_[index].project(AP);

}
/*----------------------------------------------------------------------------*/
math::Point FacetedSurface::getCenter() const
{
	throw GMDSException("Not yet implemented");
//	TBase totalArea = 0.;
//	GEPETO::Point<3,TBase> center(0.,0.,0.);
//
//	for(unsigned int i=0;i<discrete_representation_.size();i++){
//		GEPETO::Point<3,TBase> triCenter = discrete_representation_[i].getCenter();
//		TBase area = discrete_representation_[i].computeArea();
//		totalArea += area;
//		center += area*triCenter;
//	}
//
//	center /= totalArea;
//	return center;

}
/*----------------------------------------------------------------------------*/
void FacetedSurface::
getMeshFaces(std::vector<Face>& AFaces) const
{
	AFaces.clear();
	AFaces = mesh_representation_;
}
/*----------------------------------------------------------------------------*/
void FacetedSurface::getTriangulation(std::vector<math::Triangle >& ATri) const
{
	ATri.clear();
	ATri.resize(mesh_representation_.size());
	for(unsigned int i=0;i<mesh_representation_.size();i++)
	{
		Face current = mesh_representation_[i];
		std::vector<Node> nodes = current.get<Node>();
		ATri[i]=math::Triangle(nodes[0].getPoint(),nodes[1].getPoint(),nodes[2].getPoint());
	}
}
/*----------------------------------------------------------------------------*/
void
FacetedSurface::reorient(FacetedVolume* AVol)
{
	std::vector<Face> faces;
	getMeshFaces(faces);
	
	if(0 == faces.size()) {
		return;
	}

	// we build a submesh containing only the faces of this surface
	MeshModel mod = DIM3|N|F|F2N|F2F;
	IGMesh mesh(mod);

	std::map<TCellID,Node> mapNodes;
	for(int iFace=0; iFace<faces.size(); iFace++) {
		std::vector<Node> nodes = faces[iFace].get<Node>();
		for(int iNode=0; iNode<nodes.size(); iNode++) {
			mapNodes[nodes[iNode].getID()] = nodes[iNode];
		}
	}

	std::map<TCellID,Node> oldID2newNodes;
	std::map<TCellID,Node>::iterator it = mapNodes.begin();
	for(;it!=mapNodes.end(); it++) {
		Node newNode = mesh.newNode(it->second.getPoint());
		oldID2newNodes[it->first] = newNode;
	}
	std::map<TCellID,Face> newID2oldFaces;
	for(int iFace=0; iFace<faces.size(); iFace++) {
		std::vector<Node> nodes = faces[iFace].get<Node>();
		Face newFace = mesh.newTriangle(oldID2newNodes[nodes[0].getID()],
				oldID2newNodes[nodes[1].getID()],
				oldID2newNodes[nodes[2].getID()]);
		newID2oldFaces[newFace.getID()] = faces[iFace];
	}

	gmds::IGMeshDoctor doc(&mesh);
	//doc.buildFacesAndR2F();
	doc.updateUpwardConnectivity();	

	// first reorient all the faces in the same direction
	int markTreatedFaces = mesh.getNewMark<Face>();
	int markFacesInverted = mesh.getNewMark<Face>();
	Face current_face = mesh.faces_begin().value();

	propagateOrient(current_face,markTreatedFaces,markFacesInverted,&mesh);

	IGMesh::face_iterator itf = mesh.faces_begin();
	for(; !itf.isDone(); itf.next()) {
		Face f = itf.value();
		if(mesh.isMarked(f,markFacesInverted)) {
			invertFace(newID2oldFaces[f.getID()]);
		}
	}

	// then check the direction; inward or outward compared to the volume
	bool isOutward = isOutwardDirection(AVol);	

	// if the orientation is inward, inverse all the faces
	if(!isOutward) {
                invertAllFaces();
	}

	mesh.unmarkAll<Face>(markTreatedFaces);
	mesh.unmarkAll<Face>(markFacesInverted);
	mesh.freeMark<Face>(markTreatedFaces);
	mesh.freeMark<Face>(markFacesInverted);	
}
/*----------------------------------------------------------------------------*/
void
FacetedSurface::propagateOrient(Face AFace, int AMarkTreatedFaces, int AMarkFacesInverted, IGMesh* AMesh)
{
	std::vector<Face> faces = AFace.get<Face>();

	for(int iFace=0; iFace<faces.size(); iFace++) {
		if(!AMesh->isMarked(faces[iFace],AMarkTreatedFaces)) {
			bool faceIsSameOrientation = checkSameOrientFace(AFace,faces[iFace]);
			
			if(!faceIsSameOrientation) {
				invertFace(faces[iFace]);
				AMesh->mark(faces[iFace],AMarkFacesInverted);	
			}

			AMesh->mark(faces[iFace],AMarkTreatedFaces);
			propagateOrient(faces[iFace],AMarkTreatedFaces,AMarkFacesInverted,AMesh);
		}
	}
}
/*----------------------------------------------------------------------------*/
bool
FacetedSurface::checkSameOrientFace(Face AFaceRef, Face AFaceCheck)
{
	std::vector<Node> nodes = AFaceRef.get<Node>();
	std::vector<Node> nodesBis = AFaceCheck.get<Node>();

	for(int iNode=0; iNode<nodes.size(); iNode++) {
		for(int iNodeBis=0; iNodeBis<nodesBis.size(); iNodeBis++) {
			if(nodesBis[iNodeBis] == nodes[iNode]) {
				if(nodesBis[(iNodeBis+1)%nodesBis.size()] == nodes[(iNode+1)%nodes.size()]) {
					return true;
				} else {
					return false;
				}
			}
		}
	}

	throw GMDSException("FacetedSurface::checkOrientFace we should not be in this part of the code.");
}
/*----------------------------------------------------------------------------*/
bool
FacetedSurface::isOutwardDirection(FacetedVolume* AVol)
{
	std::vector<GeomSurface*> surfaces;
	AVol->get(surfaces);

	// we get the triangulation of this surface
	std::vector<math::Triangle> triangles;
	this->getTriangulation(triangles);

	int nbTrianglesFromThisSurface = triangles.size();

	// get a subset of the triangles;
	// we will check the normal ray intersection evenness only on those, monte-carlo style
	const int FACETED_SURFACE_ORIENT_TRIANGLE_SAMPLING = 101;
	std::vector<int> sampledTrianglesIndex;
	if(nbTrianglesFromThisSurface < FACETED_SURFACE_ORIENT_TRIANGLE_SAMPLING) {
		for(int iTriangle=0; iTriangle<nbTrianglesFromThisSurface; iTriangle++) {
			sampledTrianglesIndex.push_back(iTriangle);
		}
	} else {
		gmds::RandomGenerator rgen;
		rgen.init();
		for(int i=0; i<FACETED_SURFACE_ORIENT_TRIANGLE_SAMPLING; i++) {
			int index = nearbyint(rgen.value() * FACETED_SURFACE_ORIENT_TRIANGLE_SAMPLING);
			sampledTrianglesIndex.push_back(index);
		}
	}
 
	// we add all the triangles of the surfaces except this one
	for(int iSurf=0; iSurf<surfaces.size(); iSurf++) {
		if(surfaces[iSurf] != this) {
			std::vector<math::Triangle> triangles_tmp;
			surfaces[iSurf]->getTriangulation(triangles_tmp);

			for(int iTriangle=0; iTriangle<triangles_tmp.size(); iTriangle++) {
				triangles.push_back(triangles_tmp[iTriangle]);
			}	
		}
	}

	// for every triangle of this surface, draw a ray 
	// and check its number of intersections
	int nbInside = 0;
	
	for(int index=0; index<sampledTrianglesIndex.size(); index++) {
		int iTriangle = sampledTrianglesIndex[index];
		math::Point pt(triangles[iTriangle].getCenter());
		math::Vector dir(triangles[iTriangle].getNormal());
		math::Ray ray(pt,dir);
		
		int nbIntersect = 0;

		for(int iTriangleBis=0; iTriangleBis<triangles.size(); iTriangleBis++) {
			if(iTriangleBis != iTriangle) {
				bool doesIntersect = triangles[iTriangleBis].intersect(ray);
				if(doesIntersect) {
					nbIntersect++;
				}
			}
		}
		if(nbIntersect%2 == 1) {
			nbInside++;
		}
	}

	std::cout<<"nbInside "<<nbInside<<" of "<<sampledTrianglesIndex.size()<<std::endl;
	if(nbInside>(sampledTrianglesIndex.size()/2)) {
		return false;
	} else {
		return true;
	}
}
/*----------------------------------------------------------------------------*/
void
FacetedSurface::invertAllFaces()
{
	std::vector<Face> faces;
	getMeshFaces(faces);
	for(int iFace=0; iFace<faces.size(); iFace++) {
		invertFace(faces[iFace]);
	}
}
/*----------------------------------------------------------------------------*/
void
FacetedSurface::invertFace(Face AFace)
{
	std::vector<Node> nodes = AFace.get<Node>();
	std::vector<Node> nodestmp(nodes.rbegin(),nodes.rend());

	AFace.set<Node>(nodestmp);	
}
/*----------------------------------------------------------------------------*/
} // namespace geom
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
