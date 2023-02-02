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
/** \file    FacetedCurve.cpp
 *  \author  F. LEDOUX
 *  \date    30/05/2011
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/CAD/FacetedCurve.h>
#include <GMDS/CAD/FacetedSurface.h>
#include <GMDS/CAD/FacetedPoint.h>
#include <GMDS/CAD/FacetedVolume.h>
#include <GMDS/Math/Constants.h>
#include <GMDS/Math/Plane.h>
#include <set>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace geom{
/*----------------------------------------------------------------------------*/
int FacetedCurve::next_id_=1;
/*----------------------------------------------------------------------------*/
FacetedCurve::FacetedCurve()
:id_(next_id_++),p1_(0),p2_(0)
{}
/*----------------------------------------------------------------------------*/
FacetedCurve::FacetedCurve(
		FacetedPoint* AP1,
		FacetedPoint* AP2,
		std::vector<Node>& APoints,
		std::vector<Edge>& AEdges,
		const std::string& AName)
:GeomCurve(AName),id_(next_id_++),p1_(AP1),p2_(AP2),
mesh_representation_(APoints),mesh_representation_edges_(AEdges)
{

	for(unsigned int i=0;i<APoints.size()-1;i++){
		Node n1 = APoints[i];
		Node n2 = APoints[i+1];
		math::Point p1(n1.X(),n1.Y(),n1.Z());
		math::Point p2(n2.X(),n2.Y(),n2.Z());
		geom_representation_.push_back(math::Segment(p1,p2));
	}

}
/*----------------------------------------------------------------------------*/
FacetedCurve::
~FacetedCurve()
{}
/*----------------------------------------------------------------------------*/
void FacetedCurve::add(FacetedSurface* ASurf)
{
	surfaces_.push_back(ASurf);
}
/*----------------------------------------------------------------------------*/
void FacetedCurve::get(std::vector<GeomPoint*>& points) const
{
	points.clear();
	points.push_back(p1_);
	points.push_back(p2_);
}
/*----------------------------------------------------------------------------*/
void FacetedCurve::get(std::vector<GeomSurface*>& ASurf) const
{
	ASurf.clear();
	ASurf.resize(surfaces_.size());
	for(unsigned int i=0;i<surfaces_.size();i++)
		ASurf[i]=surfaces_[i];
}
/*----------------------------------------------------------------------------*/
void FacetedCurve::get(std::vector<GeomVolume*>& AVol) const
{
	std::set<GeomVolume* > vol_set;
	for(unsigned int i=0;i<surfaces_.size();i++)
	{
		std::vector<GeomVolume* > vol_i;
		surfaces_[i]->get(vol_i);
		vol_set.insert(vol_i.begin(),vol_i.end());
	}
	AVol.clear();
	AVol.insert(AVol.begin(),vol_set.begin(),vol_set.end());
}
/*----------------------------------------------------------------------------*/
double FacetedCurve::length() const
{
	double length = 0.;
	for(int i=0; i<geom_representation_.size(); i++) {
		length += geom_representation_[i].computeLength();
	}

	return length;
}
/*----------------------------------------------------------------------------*/
void FacetedCurve::project(math::Point& AP) const
{
	TCoord min_dist = geom_representation_[0].distance(AP);
	int index = 0;
	for(unsigned int i=1;i<geom_representation_.size();i++){
		TCoord dist = geom_representation_[i].distance(AP);
		if(dist<min_dist){
			min_dist=dist;
			index = i;
		}
	}
	AP = geom_representation_[index].project(AP);
}
/*----------------------------------------------------------------------------*/
void FacetedCurve::project(math::Point& AP, math::Vector& AV) const
{
	TCoord min_dist = geom_representation_[0].distance(AP);
	int index = 0;
	for(unsigned int i=1;i<geom_representation_.size();i++){
		TCoord dist = geom_representation_[i].distance(AP);
		if(dist<min_dist){
			min_dist=dist;
			index = i;
		}
	}
	AP = geom_representation_[index].project(AP);

	math::Vector vector_tmp( geom_representation_[index].getPoint(0),
										geom_representation_[index].getPoint(1));
	AV = vector_tmp;
}
/*----------------------------------------------------------------------------*/
TCoord FacetedCurve::computeArea() const
{
	throw GMDSException("Not yet implemented!");
}
/*----------------------------------------------------------------------------*/
void FacetedCurve::computeBoundingBox(TCoord minXYZ[3], TCoord maxXYZ[3]) const
{
	std::vector<math::Point > pnts;
	for(unsigned int i=0;i<geom_representation_.size();i++){
		pnts.push_back(geom_representation_[i].getPoint(0));
		pnts.push_back(geom_representation_[i].getPoint(1));
	}
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
bool FacetedCurve::isALoop()const
{
	std::vector<GeomPoint*> points;
	this->get(points);

	if(2 != points.size()) {
		throw GMDSException("FacetedCurve::isALoop this curve has only one point.");
	}

	// in case of a loop, both points are NULL
	if(points[0] != points[1]) {
		return false;
	} else {
		return true;
	}
}
/*----------------------------------------------------------------------------*/
TCoord FacetedCurve::computeDihedralAngle() const
{
	TCoord anglemax = 0.;

	// we go only to size-1 because we are interested in edges (segments)
	for(unsigned int iNode=0; iNode<mesh_representation_.size()-1; iNode++)
	{
		std::vector<Edge> edges_possibilities;
		mesh_representation_[iNode].get<Edge>(edges_possibilities);

		// we discard edges that have only one point on the curve
		Edge current_edge;

		for(unsigned int iEdge=0; iEdge<edges_possibilities.size(); iEdge++)
		{
			std::vector<Node> nodes;
			edges_possibilities[iEdge].get<Node>(nodes);

			if(nodes[0] == mesh_representation_[iNode])
			{
				if(nodes[1] == mesh_representation_[iNode+1])
				{
					current_edge = edges_possibilities[iEdge];
					break;
				}
			}
			else
			{
				if(nodes[0] == mesh_representation_[iNode+1])
				{
					current_edge = edges_possibilities[iEdge];
					break;
				}
			}
		}

		// cannot use MeshQuery because we do not have mask?
//		{
//			const int mask = DIM3|N|E|F|F2N|E2N|E2F|F2E|N2E;
//
//			MeshQuery<mask> query;
//			TCoord angle = query.computeDihedralAngle(current_edge);
//
//			if(angle > anglemax)
//				anglemax  = angle;
//		}

		{
			std::vector<Face> faces;
			current_edge.get<Face>(faces);
			if(faces.size()!=2)
				throw GMDSException("Dihedral angles can only be computed for edges having 2 adjacent faces");

			math::Vector n0 = faces[0].normal();
			math::Vector n1 = faces[1].normal();
			//norm are equals to 1

			// first get angle between [0,PI]
			TCoord angle = acos(-(n0.dot(n1)))*math::Constants::INVPIDIV180;

			// then determine the sign of the sin
			// it works because the degenerate case here would be around Pi, and if the angle is Pi the curve is
			// not a sharp curve.
			math::Point p0 = faces[0].center();
			math::Point p1 = faces[1].center();
			math::Point p = (math::Plane (p1,n1)).project(p0);

			TCoord sign = (math::Vector (p0,p)).dot(n1);

			if(sign < 0.)
				angle = 360. - angle;

			if(angle > anglemax)
				anglemax  = angle;
		}

//		{
//			std::vector<Face*> faces;
//			current_edge->getFaces(faces);
//			if(faces.size()!=2)
//				throw GMDSException("Dihedral angles can only be computed for edges having 2 adjacent faces");
//
//			math::Vector n0 = faces[0]->normal();
//			math::Vector n1 = faces[1]->normal();
//
//			n0 = n0 * (-1.);
//
//			math::Vector n2 = n0.cross(n1);
//
//			TCoord cosPhi = n0.dot(n1);
//
//			math::Vector sinPhin2crossn0 = n1 - cosPhi*n0; //- (1. - cosPhi)*(n0.dot(n2))*n2;
//
//
//			TCoord angle = TCoord::acos(cosPhi)*math::NumericConstants<TCoord>::INVPIDIV180;
//
//			if(((n2.cross(n0)).dot(sinPhin2crossn0)) < 0.)
//				angle = angle+180;
//
//			if(angle > anglemax)
//				anglemax  = angle;
//		}

	}

	return anglemax;
}
/*----------------------------------------------------------------------------*/
void FacetedCurve::
computeVector(const GeomPoint& AP, math::Vector& AV) const
{

	math::Vector tangent;

	if(&AP == p1_)
	{
		unsigned int index = 0;

		tangent = math::Vector (geom_representation_[index].getPoint(0),geom_representation_[index].getPoint(1));
	}
	else
	{
		unsigned int index = geom_representation_.size() -1;

		tangent = math::Vector (geom_representation_[index].getPoint(1),geom_representation_[index].getPoint(0));
	}

	AV = tangent;
	AV.normalize();

}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
GeomSurface* FacetedCurve::
getLeftSurface() const
{
//	throw GMDSException("Not yet implemented");
	std::vector<GeomSurface* > surfaces;
	get(surfaces);

	gmds::geom::FacetedSurface* surface1 = (dynamic_cast<gmds::geom::FacetedSurface* >(surfaces[0]));
	if(surface1 == NULL) {
		throw GMDSException("FacetedCurve::getLeftSurface dynamic_cast failed for surface1.");
	}
	gmds::geom::FacetedSurface* surface2 = (dynamic_cast<gmds::geom::FacetedSurface* >(surfaces[1]));
	if(surface2 == NULL) {
                throw GMDSException("FacetedCurve::getLeftSurface dynamic_cast failed for surface2.");
        }

	Node node1 = mesh_representation_[0];
	Node node2 = mesh_representation_[1];

	std::vector<math::Triangle > triangles;
	surface1->getTriangulation(triangles);

	math::Point point1 = node1.getPoint();
	math::Point point2 = node2.getPoint();

	for(unsigned int iTriangle=0; iTriangle<triangles.size(); iTriangle++) {
		std::vector<math::Point > points(triangles[iTriangle].getNbPoints());
		points[0] = triangles[iTriangle].getPoint(0);
		points[1] = triangles[iTriangle].getPoint(1);
		points[2] = triangles[iTriangle].getPoint(2);

		for(unsigned int iPoint=0; iPoint<points.size(); iPoint++) {
			if(points[iPoint] == point1) {
				if(points[(iPoint+1)%triangles[iTriangle].getNbPoints()] == point2) {
					return surface1;
				} else {
					if(points[(iPoint+(triangles[iTriangle].getNbPoints()-1))%triangles[iTriangle].getNbPoints()] == point2) {
						return surface2;
					}
				}
			}
		}
	}

	throw GMDSException("Did not find the left surface");
}
/*----------------------------------------------------------------------------*/

GeomSurface* FacetedCurve::
getRightSurface() const
{
	throw GMDSException("Not yet implemented");
//	std::vector<GeomSurface* > surfaces;
//	get(surfaces);
//
//	gmds::geom::FacetedSurface* surface1 = (dynamic_cast<gmds::geom::FacetedSurface* >(surfaces[0]));
//	gmds::geom::FacetedSurface* surface2 = (dynamic_cast<gmds::geom::FacetedSurface* >(surfaces[1]));
//
//	Node* node1 = mesh_representation_[0];
//	Node* node2 = mesh_representation_[1];
//
//	std::vector<math::Triangle > triangles;
//	surface1->getTriangulation(triangles);
//
//	math::Point point1 = node1->getPoint();
//	math::Point point2 = node2->getPoint();
//
//	for(unsigned int iTriangle=0; iTriangle<triangles.size(); iTriangle++) {
//		std::vector<math::Point > points(triangles[iTriangle].getNbPoints());
//		points[0] = triangles[iTriangle].getPoint(0);
//		points[1] = triangles[iTriangle].getPoint(1);
//		points[2] = triangles[iTriangle].getPoint(2);
//
//		for(unsigned int iPoint=0; iPoint<points.size(); iPoint++) {
//			if(points[iPoint] == point1) {
//				if(points[(iPoint+1)%triangles[iTriangle].getNbPoints()] == point2) {
//					return surface2;
//				} else {
//					if(points[(iPoint+(triangles[iTriangle].getNbPoints()-1))%triangles[iTriangle].getNbPoints()] == point2) {
//						return surface1;
//					}
//				}
//			}
//		}
//	}
//
//	throw GMDSException("Did not find the right surface");
}
/*----------------------------------------------------------------------------*/

unsigned int FacetedCurve::
getNbSegments() const
{
	return this->geom_representation_.size();
}
/*----------------------------------------------------------------------------*/

math::Segment FacetedCurve::
getSegment(const unsigned int& AISegment) const
{
	return this->geom_representation_[AISegment];
}
/*----------------------------------------------------------------------------*/

void FacetedCurve::
getMultiplePoints(
		const int& ANbPoints,
		math::Point* APoints) const
{
	// first compute length of the curve
	TCoord length = 0.;

	for(unsigned int iSegment=0; iSegment<geom_representation_.size(); iSegment++) {
		length += geom_representation_[iSegment].computeLength();
	}

	TCoord step_length = length / ANbPoints;

	// create all points
	APoints[0] = mesh_representation_[0].getPoint();

	unsigned int iPoint = 1;
	unsigned int iSegment = 0;
	bool isNewInSegment = true;
	math::Point current_point = APoints[0];


	length = 0.;

	while(iPoint<ANbPoints) {
		if((length < step_length*iPoint)
		&& (length + current_point.distance(geom_representation_[iSegment].getPoint(1)) >= step_length*iPoint)) {

			math::Point point;
			if(isNewInSegment) {
				point  = geom_representation_[iSegment].getPoint(0);
			} else {
				point = APoints[iPoint-1];
			}

			math::Vector vect(geom_representation_[iSegment].getPoint(0),geom_representation_[iSegment].getPoint(1));
			vect.normalize();
			point = point + vect*(step_length*iPoint - length);

			APoints[iPoint] = point;
			current_point = APoints[iPoint];
			length = step_length*iPoint;
			iPoint++;
			isNewInSegment = false;

		} else {
			if(isNewInSegment) {
				length += geom_representation_[iSegment].computeLength();
			} else {
				length += APoints[iPoint-1].distance(geom_representation_[iSegment].getPoint(1));
			}
			current_point = geom_representation_[iSegment].getPoint(1);
			isNewInSegment = true;
			iSegment++;
		}

	} // while(iPoint<ANbPoints) {

//	this->project(starting_point);


}
/*----------------------------------------------------------------------------*/

TCoord FacetedCurve::
computeDistanceHaussdorf(
		const math::Segment& ASegment) const
{
/*
	TCoord dist = HUGE_VALF;

	for(unsigned int iSegment=0; iSegment<geom_representation_.size(); iSegment++) {
//		TCoord dist_tmp = ASegment.distanceInf(geom_representation_[iSegment]);
		TCoord dist_tmp = geom_representation_[iSegment].distanceInf(ASegment);
		if(dist_tmp<dist)
			dist = dist_tmp;
	}

	return dist;
*/
	TCoord dist0 = HUGE_VALF;
	TCoord dist1 = HUGE_VALF;

	math::Point pt0 = ASegment.getPoint(0);
	math::Point pt1 = ASegment.getPoint(1);

	for(unsigned int iSegment=0; iSegment<geom_representation_.size(); iSegment++) {
		TCoord dist0_tmp = (geom_representation_[iSegment]).distance(pt0);
		TCoord dist1_tmp = (geom_representation_[iSegment]).distance(pt1);
		
		if(dist0_tmp<dist0) {
                        dist0 = dist0_tmp;
		}
		if(dist1_tmp<dist1) {
                        dist1 = dist1_tmp;
                }

	}

	return std::max(dist0,dist1);
}
/*----------------------------------------------------------------------------*/

TCoord FacetedCurve::
computeDistanceHaussdorfSubCurve(const int& ANbPoints, const int& AIPoint, const math::Segment& ASegment) const
{
	// first compute length of the curve
	TCoord length = 0.;

	for(unsigned int iSegment=0; iSegment<geom_representation_.size(); iSegment++) {
		length += geom_representation_[iSegment].computeLength();
	}

	TCoord step_length = length / ANbPoints;

	// create all points
	std::vector<math::Point > points;
	points.resize(ANbPoints);
	points[0] = mesh_representation_[0].getPoint();

	unsigned int iPoint = 1;
	unsigned int iSegment = 0;
	bool isNewInSegment = true;
	math::Point current_point = points[0];

	std::vector<math::Point > points_tmp;

	length = 0.;

	while(iPoint<ANbPoints+1) {

		if((iPoint == AIPoint+1) && (points_tmp.size() == 0)) {
			points_tmp.push_back(current_point);
		}
		if((iPoint == AIPoint+1) && (AIPoint == ANbPoints-1)) {
			for(unsigned int iSegment_tmp=iSegment; iSegment_tmp<geom_representation_.size(); iSegment_tmp++) {
				math::Point point = geom_representation_[iSegment_tmp].getPoint(1);
				points_tmp.push_back(point);
			}
			break;
		}


		if((length < step_length*iPoint)
		&& (length + current_point.distance(geom_representation_[iSegment].getPoint(1)) >= step_length*iPoint)) {

			math::Point point;
			if(isNewInSegment) {
				point  = geom_representation_[iSegment].getPoint(0);
			} else {
				point = points[iPoint-1];
			}

			math::Vector vect(geom_representation_[iSegment].getPoint(0),geom_representation_[iSegment].getPoint(1));
			vect.normalize();
			point = point + vect*(step_length*iPoint - length);

			points[iPoint] = point;
			current_point = points[iPoint];
			length = step_length*iPoint;
			iPoint++;
			isNewInSegment = false;

			if(iPoint == AIPoint+2) {
				points_tmp.push_back(point);
				break;
			}

		} else {
			if(isNewInSegment) {
				length += geom_representation_[iSegment].computeLength();
			} else {
				length += points[iPoint-1].distance(geom_representation_[iSegment].getPoint(1));
			}
			current_point = geom_representation_[iSegment].getPoint(1);
			isNewInSegment = true;
			iSegment++;

			points_tmp.push_back(current_point);
		}

	} // while(iPoint<ANbPoints) {

	// compute haussdorf distance
	TCoord dist = HUGE_VALF;

	for(unsigned int iSegment=0; iSegment<points_tmp.size()-1; iSegment++) {

		math::Segment segment(points_tmp[iSegment],points_tmp[iSegment+1]);

		TCoord dist_tmp = segment.distanceInf(ASegment);
		if(dist_tmp<dist)
			dist = dist_tmp;
	}

	return dist;

}
/*----------------------------------------------------------------------------*/
void FacetedCurve::getMeshNodes(std::vector<Node>& ANodes) const
{
	ANodes.clear();
	ANodes = mesh_representation_;
}
/*----------------------------------------------------------------------------*/
void FacetedCurve::getMeshEdges(std::vector<Edge>& AEdges) const
{
	AEdges.clear();
	AEdges = mesh_representation_edges_;
}
/*----------------------------------------------------------------------------*/
} // namespace geom
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
