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
 * FacetedGeomManager.t.h
 *
 *  Created on: 1 juil. 2011
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/CAD/FacetedGeomManager.h>
/*----------------------------------------------------------------------------*/
#include <set>
/*----------------------------------------------------------------------------*/
namespace gmds{
/*----------------------------------------------------------------------------*/
namespace geom{
/*----------------------------------------------------------------------------*/
FacetedGeomManager::FacetedGeomManager()
:mesh_(DIM3|N|E|F|F2N|E2N|E2F|F2E|N2E)
{}
/*----------------------------------------------------------------------------*/
FacetedGeomManager::~FacetedGeomManager()
{
	for(unsigned int i=0;i<volumes_.size();i++)
		if(volumes_[i]!=0)
			delete volumes_[i];

	for(unsigned int i=0;i<surfaces_.size();i++)
		if(surfaces_[i]!=0)
			delete surfaces_[i];

	for(unsigned int i=0;i<curves_.size();i++)
		if(curves_[i]!=0)
			delete curves_[i];

	for(unsigned int i=0;i<points_.size();i++)
		if(points_[i]!=0)
			delete points_[i];


	for(unsigned int i=0;i<services_.size();i++)
		if(services_[i]!=0)
			delete services_[i];
}
/*----------------------------------------------------------------------------*/
GeomVolume* FacetedGeomManager::newVolume()
{
	FacetedVolume* v = new FacetedVolume();
	volumes_.push_back(v);
	return v;
}
/*----------------------------------------------------------------------------*/
GeomSurface* FacetedGeomManager::newSurface()
{
	FacetedSurface* s = new FacetedSurface();
	surfaces_.push_back(s);
	return s;
}
/*----------------------------------------------------------------------------*/
GeomCurve* FacetedGeomManager::newCurve()
{
	FacetedCurve* c = new FacetedCurve();
	curves_.push_back(c);
	return c;
}
/*----------------------------------------------------------------------------*/
GeomPoint* FacetedGeomManager::newPoint()
{
	FacetedPoint* p = new FacetedPoint();
	points_.push_back(p);
	return p;
}
/*----------------------------------------------------------------------------*/
GeomTriangulationService*
FacetedGeomManager::newGeomTriangulationService(){
	GeomTriangulationService* s = new FacetedTriangulationService();
	services_.push_back(s);
	return s;
}
/*----------------------------------------------------------------------------*/
IGMesh& FacetedGeomManager::getMeshView(){
	return mesh_;
}
/*----------------------------------------------------------------------------*/
TInt FacetedGeomManager::getNbPoints() const
{
	return points_.size();
}
/*----------------------------------------------------------------------------*/
TInt FacetedGeomManager::getNbCurves() const
{
	return curves_.size();
}
/*----------------------------------------------------------------------------*/
TInt FacetedGeomManager::getNbSurfaces() const
{
	return surfaces_.size();
}
/*----------------------------------------------------------------------------*/
TInt FacetedGeomManager::getNbVolumes() const
{
	return volumes_.size();
}
/*----------------------------------------------------------------------------*/
void FacetedGeomManager::getVolumes(std::vector<GeomVolume*>& volumes) const
{
	volumes.clear();
	volumes.resize(volumes_.size());
	for(unsigned int i=0;i<volumes_.size();i++)
		volumes[i]=volumes_[i];
}
/*----------------------------------------------------------------------------*/
void FacetedGeomManager::getSurfaces(std::vector<GeomSurface*>& surfaces) const
{
	surfaces.clear();
	surfaces.resize(surfaces_.size());
	for(unsigned int i=0;i<surfaces_.size();i++)
		surfaces[i]=surfaces_[i];
}
/*----------------------------------------------------------------------------*/
void FacetedGeomManager::getCurves(std::vector<GeomCurve*>& curves) const
{
	curves.clear();
	curves.resize(curves_.size());
	for(unsigned int i=0;i<curves_.size();i++)
		curves[i]=curves_[i];
}
/*----------------------------------------------------------------------------*/
void FacetedGeomManager::
getPoints(std::vector<GeomPoint*>& points) const
{
	points.clear();
	points.resize(points_.size());
	for(unsigned int i=0;i<points_.size();i++)
		points[i]=points_[i];
}
/*----------------------------------------------------------------------------*/
void FacetedGeomManager::reorient()
{
	throw GMDSException("not yet implemented");
//	TInt nbFacesDone = 0;
//	int mark_done = mesh_.getNewMark();
//
//	// we proceed surface by surface, selecting a face and propagating.
//	// if there are still faces not treated, it means the mesh is divided among
//	// several disconnected surfaces, and we locate and treat them.
//	while(nbFacesDone < mesh_.getNbFaces()){
//
//		// Initialization
//		typename Mesh<mask>::faces_iterator it_faces = mesh_.faces_begin();
//
//		// locate a face pertaining to a surface not yet treated
//		while(mesh_.isMarked(it_faces->currentItem(),mark_done)){
//			it_faces->next();
//		}
//		Face* f = it_faces->currentItem();
//		OutwardNormal(f);
//		int mark_wait = mesh_.getNewMark();
//
//		// Stack of the Treated Faces
//		mesh_.mark(f,mark_wait);
//		mesh_.mark(f,mark_done);
//		nbFacesDone++;
//
//		// Stack of waiting Faces / Waiting Nodes of Edge
//		std::stack<Face*> stackWaitF;
//		std::stack<Node*> stackWaitN1;
//		std::stack<Node*> stackWaitN2;
//
//		// Initialization for the first triangle
//		std::vector<Node*> NodesOf_f = f->getNodes();
//		std::vector<Edge*> vectEdge = f->getEdges();
//		for (unsigned int ii=0;ii<vectEdge.size();ii++){
//			Edge* ei = vectEdge[ii];
//			std::vector<Face*> vectFace = ei->getFaces();
//			if (vectFace[0]!=f){
//				stackWaitF.push(vectFace[0]);
//				mesh_.mark(vectFace[0],mark_wait);
//				stackWaitN1.push(NodesOf_f[ii]);
//				stackWaitN2.push(NodesOf_f[(ii+1)%3]);
//			}
//			else if(vectFace.size()>1 && vectFace[1]!=f){
//				stackWaitF.push(vectFace[1]);
//				mesh_.mark(vectFace[1],mark_wait);
//				stackWaitN1.push(NodesOf_f[ii]);
//				stackWaitN2.push(NodesOf_f[(ii+1)%3]);
//			}
//		}
//			int cmpt = 0;
//		while(!stackWaitF.empty()){
//			Face* f = stackWaitF.top();
//			mesh_.mark(f,mark_done);
//			nbFacesDone++;
//			Node* n1 = stackWaitN1.top();
//			Node* n2 = stackWaitN2.top();
//
//			stackWaitF.pop();
//			stackWaitN1.pop();
//			stackWaitN2.pop();
//
//			// check if the orientation of face f is OK
//			std::vector<Node*> NodesOf_f = f->getNodes();
//
//			Node* after_n1=0;
//			if(NodesOf_f[0]==n1)
//				after_n1 = NodesOf_f[1];
//			else if (NodesOf_f[1]==n1)
//				after_n1 = NodesOf_f[2];
//			else
//				after_n1 = NodesOf_f[0];
//
//			if (after_n1==n2)
//			{
//				// Reorientation of the face
//				std::vector<id> lids = f->getNodeIDs();
//				id id_tmp = lids[0];
//				lids[0] = lids[1];
//				lids[1] = id_tmp;
//				f->setNodes(lids);
//				// reorient the edges
//				f->arrangeEdges();
//				cmpt++;
//			}
//			std::vector<Edge*> edgesFace = f->getEdges();
//			NodesOf_f = f->getNodes();
//			for (unsigned int ii=0;ii<edgesFace.size();ii++){
//				Edge* ei = edgesFace[ii];
//				std::vector<Face*> vectFace = ei->getFaces();
//				if (vectFace[0]!=f &&
//					!mesh_.isMarked(vectFace[0],mark_wait))
//				{
//					stackWaitF.push(vectFace[0]);
//					mesh_.mark(vectFace[0],mark_wait);
//					stackWaitN1.push(NodesOf_f[ii]);
//					stackWaitN2.push(NodesOf_f[(ii+1)%3]);
//				}
//				else if(vectFace.size()>1 &&
//						vectFace[1]!=f &&
//						!mesh_.isMarked(vectFace[1],mark_wait))
//				{
//					stackWaitF.push(vectFace[1]);
//					mesh_.mark(vectFace[1],mark_wait);
//					stackWaitN1.push(NodesOf_f[ii]);
//					stackWaitN2.push(NodesOf_f[(ii+1)%3]);
//				}
//			}
//
//		}
//		std::cout<<"Nb of reoriented faces "<<cmpt<<std::endl;
//		mesh_.unmarkAll(mark_wait);
//		mesh_.freeMark(mark_wait);
//
//	}
//
//	mesh_.unmarkAll(mark_done);
//	mesh_.freeMark(mark_done);
//
}
/*----------------------------------------------------------------------------*/
void FacetedGeomManager::OutwardNormal( Face* f)
{
	throw GMDSException("Not yet implemented");
//	TCoord r = 5; // Radius of the neighborhood cylinder
//	int nbOfDrilling = 0; //Number of time the half-line meet the surface
//	int mark_neighbour = mesh_.getNewMark();
//	mesh_.mark(f,mark_neighbour);
//	int mark_drilled = mesh_.getNewMark();
//	mesh_.mark(f,mark_drilled);
//	std::stack<Face*> stackDrilledF;
//
//	// COMPUTATION OF THE COORDINATES OF THE LINE.
//	math::VectorDyn normalRef(3);
//	normalRef.set(0,f->getNormal().get(0));
//	normalRef.set(1,f->getNormal().get(1));
//	normalRef.set(2,f->getNormal().get(2));
//	// Computation of the coordinates of the center of the triangle
//	std::vector<Node*> VNodes = f->getNodes();
//	math::VectorDyn V1(VNodes[0]->getX(),VNodes[0]->getY(),VNodes[0]->getZ());
//	math::VectorDyn V2(VNodes[1]->getX(),VNodes[1]->getY(),VNodes[1]->getZ());
//	math::VectorDyn V3(VNodes[2]->getX(),VNodes[2]->getY(),VNodes[2]->getZ());
//	math::VectorDyn VCenter = (1.0/3.0)*(V1+V2+V3);
//
//	// CREATION OF THE STACK OF ELEMENTS CLOSE TO THE LINE
//	std::stack<Face*> stackNeighbourF;
//	typename Mesh<mask>::nodes_iterator it_nodes = mesh_.nodes_begin();
//	for(;!it_nodes->isDone();it_nodes->next()){
//	  Node* n = it_nodes->currentItem();
//	  math::VectorDyn Xn(n->getX(),n->getY(),n->getZ());
//	  // Now, we compute the distance from node n to the line
//	  TCoord d = ((Xn-VCenter).cross(normalRef)).norm();
//
//	  if(d<r){
//		  std::vector<Face*> VFaces = n->getFaces();
//		  for(unsigned int i=0;i<VFaces.size();i++){
//			  Face* fi = VFaces[i];
//			  if(!mesh_.isMarked(fi,mark_neighbour)){
//				  stackNeighbourF.push(fi);
//				  mesh_.mark(fi,mark_neighbour);
//				  //Return value true if the intersection line/plane
//				  //of the face fi is inside the triangle fi
//				  bool isInsidefi = IsInsideFace(mark_drilled,normalRef,VCenter,fi);
//				  if(isInsidefi){
//					  stackDrilledF.push(fi);
//					  std::cout<<"Drilled face : "<<fi->getID()<<std::endl;
//					  mesh_.mark(fi,mark_drilled);
//					  nbOfDrilling++;
//				  }
//			  }
//		  }
//	  }
//	}
//	std::cout<<"Nb of drillings : "<<nbOfDrilling<<std::endl;
//	if(nbOfDrilling%2==1){
//		// Reorientation of the face
//		std::vector<id> lids = f->getNodeIDs();
//		id id_tmp = lids[0];
//		lids[0] = lids[1];
//		lids[1] = id_tmp;
//		f->setNodes(lids);
//		// reorient the edges
//		f->arrangeEdges();
//	}
//	// Unmark all the faces which was marked
//	while(!stackNeighbourF.empty()){
//		Face* f = stackNeighbourF.top();
//		mesh_.unmark(f,mark_neighbour);
//		stackNeighbourF.pop();
//	}
//
//	mesh_.unmarkAll(mark_neighbour);
//	mesh_.freeMark(mark_neighbour);
//	mesh_.unmarkAll(mark_drilled);
//	mesh_.freeMark(mark_drilled);
}
/*----------------------------------------------------------------------------*/
bool FacetedGeomManager::IsInsideFace(const int& mark_drilled,
									  const math::VectorDyn& normalRef,
									  const math::VectorDyn& VCenter,
									   Face fi)
{
throw GMDSException("Not yer implmented");
return true;
//	// Get the coordinates of the triangle fi
//	std::vector<Edge*> VEdges = fi->getEdges();
//	std::vector<Node*> VNodes = fi->getNodes();
//	math::VectorDyn V1(VNodes[0]->getX(),VNodes[0]->getY(),VNodes[0]->getZ());
//	math::VectorDyn V2(VNodes[1]->getX(),VNodes[1]->getY(),VNodes[1]->getZ());
//	math::VectorDyn V3(VNodes[2]->getX(),VNodes[2]->getY(),VNodes[2]->getZ());
//
//	//PARAMETERIZATION OF THE PLANE OF fi
//	// using the Cartesian equation :a*x + b*y + c*z + d = 0.
//	// a,b and c are directly given by the normal
//	math::VectorDyn normal_fi(3);
//	normal_fi.set(0,fi->getNormal().get(0));
//	normal_fi.set(1,fi->getNormal().get(1));
//	normal_fi.set(2,fi->getNormal().get(2));
//	TCoord a = normal_fi.get(0);
//	TCoord b = normal_fi.get(1);
//	TCoord c = normal_fi.get(2);
//	// Compute the coefficient d
//	TCoord d = -(normal_fi.dot(V1));
//
//	// Check the case the line is parallel to the plane-there is no intersection
//	if (normalRef.dot(normal_fi) == 0.0)
//		return false;
//
//	// Coordinates of the point I, intersection between plane/Line
//	// solving a*xI + b*yI + c*zI + d = 0
//	// with xI=xC+k*xU; yI=yC+k*yU; zI=zC+k*zU
//	TCoord k = -(d + a*VCenter.get(0) + b*VCenter.get(1) +c*VCenter.get(2))
//				/(a*normalRef.get(0) + b*normalRef.get(1) + c*normalRef.get(2));
//	math::VectorDyn I(VCenter.get(0)+k*normalRef.get(0),
//				VCenter.get(1)+k*normalRef.get(1),
//				VCenter.get(2)+k*normalRef.get(2));
//
//	// Case I is on the wrong side of the half-line
//	if (k<0.0)
//		return false;
//
//	// COMPUTATION OF THE BARYCENTRIC COORDINATES
//	math::VectorDyn Coefs(2);
//	Coefs = GetBarycentricCoefs(V1,V2,V3,I,normal_fi);
//
//	TCoord alpha = Coefs.get(0);
//	TCoord beta  = Coefs.get(1);
//	TCoord gamma = 1.0-alpha-beta;
//
//	if (alpha>0.0 && beta>0.0 && gamma>0.0)
//		return true;
//	else if(alpha==0.0 && beta!=0.0){
//		std::vector<Face*> AdjFaces = VEdges[1]->getFaces();
//		if(AdjFaces.size()>0
//		   && !mesh_.isMarked(AdjFaces[0],mark_drilled)
//		   && !mesh_.isMarked(AdjFaces[1],mark_drilled))
//			return true;
//		else if(AdjFaces.size()==1)
//			return true;
//		else
//		return false;
//	}
//	else if(beta==0.0 && alpha!=0.0){
//		std::vector<Face*> AdjFaces = VEdges[2]->getFaces();
//		if(AdjFaces.size()>0
//		   && !mesh_.isMarked(AdjFaces[0],mark_drilled)
//		   && !mesh_.isMarked(AdjFaces[1],mark_drilled))
//			return true;
//		else if(AdjFaces.size()==1)
//			return true;
//		else
//		return false;
//	}
//	else if(gamma==0.0 && alpha!=0.0){
//		std::vector<Face*> AdjFaces = VEdges[0]->getFaces();
//		if(AdjFaces.size()>0
//		   && !mesh_.isMarked(AdjFaces[0],mark_drilled)
//		   && !mesh_.isMarked(AdjFaces[1],mark_drilled))
//			return true;
//		else if(AdjFaces.size()==1)
//			return true;
//		else
//		return false;
//	}
//	else if(alpha==1.0 && beta==0.0){
//		std::vector<Face*> AdjFaces = VNodes[0]->getFaces();
//		for(unsigned int ii=0;ii<AdjFaces.size();ii++){
//			if(mesh_.isMarked(AdjFaces[ii],mark_drilled)){
//				return false;
//			}
//		}
//	}
//	else if(beta==1.0 && alpha==0.0){
//		std::vector<Face*> AdjFaces = VNodes[1]->getFaces();
//		for(unsigned int ii=0;ii<AdjFaces.size();ii++){
//			if(mesh_.isMarked(AdjFaces[ii],mark_drilled)){
//				return false;
//			}
//		}
//	}
//	else if(gamma==1.0 && alpha==0.0){
//		std::vector<Face*> AdjFaces = VNodes[2]->getFaces();
//		for(unsigned int ii=0;ii<AdjFaces.size();ii++){
//			if(mesh_.isMarked(AdjFaces[ii],mark_drilled)){
//				return false;
//			}
//		}
//	}
//	else
//		return false;
//
//	return false;
}
/*----------------------------------------------------------------------------*/
math::VectorDyn FacetedGeomManager::GetBarycentricCoefs(
		math::VectorDyn& V1, math::VectorDyn& V2,
		math::VectorDyn& V3, math::VectorDyn& I,
		math::VectorDyn& normal){
throw GMDSException("Not yet implemented");

//	/* F. Ledoux. La methode a ete changee pour tester tous les plans de projection
//	 * et ne pas se limiter a un plan dans lequel la matrice T aurait pu ne pas
//	 * etre inversible
//	 */
//	// For that, we have to solve the system : T*coefs = RHS
//	// get the longuest component of the normal
//	TCoord compomax = max(fabs(normal[0]),fabs(normal[1]));
//	compomax = max(compomax,fabs(normal[2]));
//	math::MatrixDyn T(2,2); 	    // Initialization of the matrix T
//	TCoord T_det = 0.0;
//	//if (compomax==fabs(normal.get(1))){ // projection sur le plan Oyz
//	T.set(0,0,V1.get(1)-V3.get(1));
//	T.set(0,1,V2.get(1)-V3.get(1));
//	T.set(1,0,V1.get(2)-V3.get(2));
//	T.set(1,1,V2.get(2)-V3.get(2));
//	T_det = T.det();
//	//}
//	if (T_det==0.0){//  && compomax==fabs(normal.get(2))){   // projection sur le plan Ozx
//	    T.set(0,0,V1.get(2)-V3.get(2));
//	    T.set(0,1,V2.get(2)-V3.get(2));
//	    T.set(1,0,V1.get(0)-V3.get(0));
//	    T.set(1,1,V2.get(0)-V3.get(0));
//	    T_det = T.det();
//	}
//	if(T_det==0.0)
//	{  // projection sur le plan Oxy
//	    T.set(0,0,V1.get(0)-V3.get(0));
//	    T.set(0,1,V2.get(0)-V3.get(0));
//	    T.set(1,0,V1.get(1)-V3.get(1));
//	    T.set(1,1,V2.get(1)-V3.get(1));
//	    T_det = T.det();
//	}
//	if(T_det==0.0)
//		throw GMDSException("SINGULAR MATRIX IN THE COMPUTATION OF BARYCENTRIC COEFS");
//
//	math::MatrixDyn Tinv = T.inverse();
//	math::VectorDyn RHS(I.get(0)-V3.get(0),I.get(1)-V3.get(1));
//	math::VectorDyn Coefs = Tinv*RHS;
//	return Coefs;
}
/*----------------------------------------------------------------------------*/
void FacetedGeomManager::updateFromMesh(){
	std::map<TCellID,FacetedPoint* > map_node2point;
	std::map<TCellID,FacetedCurve* > map_edge2curve;

	// vertices
	IGMesh::clouds_iterator itc = mesh_.clouds_begin();
	for(; itc != mesh_.clouds_end(); itc++) {

		IGMesh::cloud current_cloud = *itc;
		std::vector<TCellID>& nodeIDs= current_cloud.cellIDs();

		if(0 == nodeIDs.size()) {
			throw GMDSException("FacetedGeomManager::importVTK a cloud is of size 0.");
		}

		// vertices are clouds with only one node
		if(1 == nodeIDs.size()) {
			FacetedPoint* p = new FacetedPoint (mesh_.get<Node>(nodeIDs[0]), current_cloud.name());
			points_.push_back(p);
			map_node2point[nodeIDs[0]] = p;
		}

	}

	/* in order to build our curves and surfaces, we need to rebuild some edges
	 * and some connectivities */
	IGMeshDoctor doc(&mesh_);
	doc.buildEdgesAndX2E();
	doc.updateUpwardConnectivity();

	// WARNING numerically unstable
//	reorient();

	// curves

	// we begin by creating the surfaces in order to associate
	// faces to surfaces
	Variable<FacetedSurface* >* surfaceAssociation =
			mesh_.newVariable<FacetedSurface* > (GMDS_FACE,"surfaceAssociation");
	{
		// creating geom surfaces from mesh surfaces (faces groups)
		IGMesh::surfaces_iterator its = mesh_.surfaces_begin();
		for(; its != mesh_.surfaces_end(); its++) {

			IGMesh::surface current_surface = *its;
			std::vector<Face> surf_faces= current_surface.cells();
			if(0 == surf_faces.size()) {
				throw GMDSException("FacetedGeomManager::importVTK a surface is of size 0.");
			}

			std::vector<FacetedPoint* > points;
			std::vector<FacetedCurve* > curves;
			FacetedSurface* s =
							new FacetedSurface(points,curves,surf_faces,current_surface.name());
			surfaces_.push_back(s);

			for(unsigned int iFace=0; iFace<surf_faces.size(); iFace++) {
				(*surfaceAssociation)[surf_faces[iFace].getID()] = s;
			}
		} // for(; its != mesh_.surfaces_end(); its++) {

	}

	// now we get create the curves
	itc = mesh_.clouds_begin();
	for(; itc != mesh_.clouds_end(); itc++) {

		IGMesh::cloud current_cloud = *itc;
		std::vector<Node> nodes= current_cloud.cells();
		std::vector<TCellID>& nodeIDs= current_cloud.cellIDs();
		// curves are clouds with more than one node
		if(1 != nodes.size()) {

			std::vector<TCellID> orderedNodeIDs;

			// we need to order the nodes
			// first we find the two end-nodes
			TCellID firstNodeID = NullID;
			TCellID lastNodeID  = NullID;

			for(unsigned int iNode=0; iNode<nodeIDs.size(); iNode++) {
				if(map_node2point.end() != map_node2point.find(nodeIDs[iNode])) {
					if(NullID == firstNodeID) {
						firstNodeID = nodeIDs[iNode];
					}
					else
					{
						lastNodeID = nodeIDs[iNode];
					}
				}
			}
//			std::cout<<"Curve from "<<firstNodeID<<" to "<<lastNodeID<<std::endl;

			if(NullID == firstNodeID) {
				throw GMDSException(""
						"FacetedGeomManager::importVTK point of firstNode of curve was not found.");
			}
			// in this case we have a loop curve
			if(NullID == lastNodeID) {
				lastNodeID = firstNodeID;
			}

			// mark all the nodes of this curve
			Variable<int>* markCurveNodes = mesh_.newVariable<int>(GMDS_NODE,"MarkCurveNodes");

			for(unsigned int iNode=0; iNode<nodeIDs.size(); iNode++) {
				markCurveNodes->set(nodeIDs[iNode],1);
			}

			// identify if this curve is adjacent to only one surface;
			bool hasOneSurface = true;
			// for example in a cylinder
			for(unsigned int iNode=0; iNode<nodes.size(); iNode++) {

				std::vector<Edge> edges = nodes[iNode].get<Edge>();

				for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {
					std::vector<Node> edgeNodes = edges[iEdge].get<Node>();
					if((*markCurveNodes)[edgeNodes[0].getID()]==1 &&
							(*markCurveNodes)[edgeNodes[1].getID()]==1)
					{
						std::vector<Face> edgeFaces = edges[iEdge].get<Face>();
						if(2 != edgeFaces.size()) {
							throw GMDSException("FacetedGeomManager::importVTK "
									"edge must have two faces");
						}

						if((*surfaceAssociation)[edgeFaces[0].getID()] != (*surfaceAssociation)[edgeFaces[1].getID()]) {
							hasOneSurface = false;
							break;
						}
					}
				}
				if(!hasOneSurface) {
					break;
				}
			}

			std::vector<Edge> edges_path;
			TCellID current_edge_id = NullID;

			std::vector<Node> orderedNodes;
			TCellID current_node_id = firstNodeID;
			Node  current_node = mesh_.get<Node>(current_node_id);
			orderedNodes.push_back(current_node);

			// we build the path of edges
			// we use a do...while here because firstNode can be equal
			// to lastNode in case of a loop curve
			do {
				std::vector<Edge> edges = current_node.get<Edge>();

				// select next edge
				bool edgeFound = false;
				for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {

					// must be different from the current one
					if(edges[iEdge].getID() != current_edge_id) {

						// its two nodes must be on the curve
						std::vector<Node> edgeNodes = edges[iEdge].get<Node>();
						if((*markCurveNodes)[edgeNodes[0].getID()]==1 &&
													(*markCurveNodes)[edgeNodes[1].getID()]==1)
						{
							// this must be a curve edge
							std::vector<Face> edgeFaces = edges[iEdge].get<Face>();

							if(2 != edgeFaces.size()) {
								throw GMDSException("FacetedGeomManager::importVTK "
										"edge must have two faces");
							}

							// in case of a one-surface edge, no need for this test
							if(	hasOneSurface ||
									((*surfaceAssociation)[edgeFaces[0].getID()] != (*surfaceAssociation)[edgeFaces[1].getID()])) {

								edgeFound = true;
								Edge current_edge = edges[iEdge];
								edges_path.push_back(current_edge);
								current_edge_id = current_edge.getID();

								if(edgeNodes[0] == current_node) {
									current_node = edgeNodes[1];
								} else {
									current_node = edgeNodes[0];
								}
								current_node_id = current_node.getID();
								orderedNodes.push_back(current_node);
								break;
							}
						}
					}
				} // for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {

				if(!edgeFound) {
					throw GMDSException("FacetedGeomManager::importVTK "
							"next edge not found.");
				}

			} while (current_node_id != lastNodeID);

			mesh_.deleteVariable(GMDS_NODE,"MarkCurveNodes");

			FacetedPoint* p1 = map_node2point[firstNodeID];
			FacetedPoint* p2 = map_node2point[lastNodeID];

			FacetedCurve* c =
							new FacetedCurve(p1,p2,orderedNodes,edges_path, current_cloud.name());

			curves_.push_back(c);
			p1->add(c);
			p2->add(c);

			for(unsigned int iEdge=0; iEdge<edges_path.size(); iEdge++) {
				map_edge2curve[edges_path[iEdge].getID()] = c;
			}

		} // if(1 != nodes.size()) {

	} // curves

	// surfaces
	// they were previously created;
	// now we will add adjacency relations
	for(unsigned int iSurf=0; iSurf<surfaces_.size(); iSurf++) {
		std::vector<Face> surf_faces;
		surfaces_[iSurf]->getMeshFaces(surf_faces);

		if(0 == surf_faces.size()) {
			throw GMDSException("FacetedGeomManager::importVTK a surface is of size 0.");
		}

		// we get all the curves surrounding this surfaces
		// we get the vertices as well
		std::set<FacetedCurve*> curves_set;
		std::vector<FacetedCurve*> curves;
		std::set<FacetedPoint*> points_set;
		std::vector<FacetedPoint*> points;

		for(unsigned int iFace=0; iFace<surf_faces.size(); iFace++) {
			std::vector<Edge> edges = surf_faces[iFace].get<Edge>();

			for(unsigned int iEdge=0; iEdge<edges.size(); iEdge++) {
				if(map_edge2curve.end() != map_edge2curve.find(edges[iEdge].getID())) {

					FacetedCurve* curv = map_edge2curve[edges[iEdge].getID()];
					curves_set.insert(curv);
					if(NULL != curv->getFirstPoint()) {
						points_set.insert(dynamic_cast<gmds::geom::FacetedPoint* > (curv->getFirstPoint()));
					}
					if(NULL != curv->getSecondPoint()) {
						points_set.insert(dynamic_cast<gmds::geom::FacetedPoint* > (curv->getSecondPoint()));
					}
				}
			}
		}

		curves.insert(curves.begin(),curves_set.begin(),curves_set.end());
		points.insert(points.begin(),points_set.begin(),points_set.end());

		surfaces_[iSurf]->replace(points);
		surfaces_[iSurf]->replace(curves);

		//connectivity curves -> surfaces
		for(unsigned int i=0; i<curves.size(); i++) {
			curves[i]->add(surfaces_[iSurf]);
		}

	} // for(unsigned int iSurf=0; iSurf<surfaces_.size(); iSurf++) {

	mesh_.deleteVariable(GMDS_FACE,"surfaceAssociation");

	FacetedVolume* v = new FacetedVolume(surfaces_);
	volumes_.push_back(v);

	// connectivity surface -> volume
	for(unsigned int i=0; i<surfaces_.size(); i++) {
		surfaces_[i]->add(v);
	}

	// we delete every clouds and surfaces of mesh_

	// vertices
	{
		unsigned int nbClouds = mesh_.getNbClouds();
		for(unsigned int iCloud=0; iCloud<nbClouds; iCloud++) {

			IGMesh::clouds_iterator itc = mesh_.clouds_begin();
			IGMesh::cloud cl = *itc;
			mesh_.deleteCloud(cl);
		}
	}

	// surfaces
	{
		unsigned int nbSurfaces = mesh_.getNbSurfaces();
		for(unsigned int iSurf=0; iSurf<nbSurfaces; iSurf++) {

			IGMesh::surfaces_iterator its = mesh_.surfaces_begin();
			IGMesh::surface surf = *its;
			mesh_.deleteSurface(surf);
		}
	}

#ifdef _DEBUG_
	std::cout<<"We imported "<<std::endl;
	std::cout<<"nb vertices "<<this->getNbPoints()<<std::endl;
	std::cout<<"nb curves "<<this->getNbCurves()<<std::endl;
	std::cout<<"nb surfaces "<<this->getNbSurfaces()<<std::endl;
#endif //_DEBUG_

}
/*----------------------------------------------------------------------------*/
void FacetedGeomManager::buildFromMesh( 
				       const bool ASingleSurface)
{
        std::map<TCellID,FacetedPoint* > map_node2point;
        std::map<TCellID,FacetedCurve* > map_edge2curve;

	if(ASingleSurface) {
		std::vector<FacetedPoint* > points;
        	std::vector<FacetedCurve* > curves;

		std::vector<gmds::Face> surf_faces;
		
		gmds::IGMesh::face_iterator itf  = this->mesh_.faces_begin();
        	for(;!itf.isDone();itf.next()) {
                	gmds::Face current_face = itf.value();
			surf_faces.push_back(current_face);
		}				

		std::string surfaceName("Surf0000");

	        FacetedSurface* s =
       			new FacetedSurface(points,curves,surf_faces,surfaceName);
	        surfaces_.push_back(s);

		FacetedVolume* v = new FacetedVolume(surfaces_);
	        volumes_.push_back(v);

		s->add(v);
	} else {
		throw GMDSException("FacetedGeomManager::buildFromMesh : not implemented when ASingleSurface==false");
	} // if(ASingleSurface)	
}
/*----------------------------------------------------------------------------*/
void 
FacetedGeomManager::reorient(FacetedVolume* AVol) 
{
	AVol->reorient();
}
/*----------------------------------------------------------------------------*/
} // namespace geom
/*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
