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
/** \file    SheetOperator.t.h
 *  \author  F. LEDOUX
 *  \date    08/08/2008
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/Algo/SheetOperator.h>
#include <GMDS/CAD/GeomEntity.h>
/*----------------------------------------------------------------------------*/
#include <sstream>
#include <set>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
SheetOperator::SheetOperator(IGMesh& m_mesh)
:m_mesh(m_mesh)
 {}
/*----------------------------------------------------------------------------*/
SheetOperator::~SheetOperator()
{}
/*----------------------------------------------------------------------------*/
std::vector<Face> SheetOperator::
pillow(
	std::vector<Region>& ARegions,
	std::vector<Face>& AFaces,
	std::vector<Region>& ANewRegions,
	bool ADisplacement)
{
	std::vector<Face> return_faces;
	ANewRegions.clear();
	m_newNodesL.clear();

	//-----------------------------------------------------------------------
	//------------------- STEP 1 --------------------------------------------
	//-----------------------------------------------------------------------
	/*
	 * We get nodes of AFaces and mark all the faces of AFaces
	 */
	int layer_mark = m_mesh.getNewMark<Face>();// for faces in AFaces
	int pillow_mark = m_mesh.getNewMark<Region>();// for regions in ARegions

	std::set<TCellID> layerNodes;
	for(unsigned int i=0;i<AFaces.size();i++){
		// mark the faces
		Face current_face = AFaces[i];
		m_mesh.mark(current_face,layer_mark);

		// get the nodes
		std::vector<Node> current_nodes;
		current_face.get<Node>(current_nodes);
		for(unsigned int ni=0;ni<current_nodes.size();ni++){
			layerNodes.insert(current_nodes[ni].getID());
		}
	}

	for(unsigned int i=0;i<ARegions.size();i++){
		// mark the regions
		m_mesh.mark(ARegions[i],pillow_mark);
	}

	/*
	 * Now for each node in layer_nodes, we indicate the classification
	 * of its two substitute nodes. Node ids are used
	 */
	std::map<TCellID, geom::GeomEntity* > N2GeomOnLayer;
	std::map<TCellID, geom::GeomEntity* > N2GeomInVolume;

	std::set<TCellID>::iterator it_nodes = layerNodes.begin();
	for(;it_nodes!=layerNodes.end();it_nodes++){
		Node currentNode = m_mesh.get<Node>(*it_nodes);

		/* we get all the faces adjacent to currentNode and that
		 * belong to the outer boundary of ARegions
		 */
		std::vector<Face> boundaryFaces;
		std::set<TCellID> allFaces;
		std::vector<Region> allRegions;

		currentNode.get<Region>(allRegions);
		for(unsigned int i=0;i<allRegions.size();i++){
			Region ri=allRegions[i];
			std::vector<Face> ri_faces;
			ri.get<Face>(ri_faces);
			for(unsigned int ir=0; ir<ri_faces.size();ir++){
				Face fir = ri_faces[ir];
				std::vector<Node> fir_nodes = fir.get<Node>();
				for(unsigned int in=0;in<fir_nodes.size();in++){
					Node nin = fir_nodes[in];
					if(nin.getID()==currentNode.getID())
						allFaces.insert(fir.getID());
				}

			}
		}
		std::set<TCellID>::iterator itf = allFaces.begin();

		for(;itf!=allFaces.end();itf++){
			Face fi=m_mesh.get<Face>(*itf);
			std::vector<Region> fi_regions;
			fi.get<Region>(fi_regions);
			int nbPillowReg = 0;
			for(unsigned int ir=0; ir<fi_regions.size();ir++){
				Region rir=fi_regions[ir];
				if(m_mesh.isMarked(rir,pillow_mark))
					nbPillowReg++;
			}
			if(nbPillowReg==1)
				boundaryFaces.push_back(fi);
		}

		//Now, boundaryFaces only contains the expected faces
		std::set<geom::GeomEntity* > associatedGeometry;

		geom::GeomEntity* refGeometry = m_mesh.getGeometricClassification(currentNode);

		for(unsigned int i=0;i<boundaryFaces.size();i++){
			Face fi=boundaryFaces[i];
			associatedGeometry.insert(m_mesh.getGeometricClassification(fi));
		}
		if(associatedGeometry.size()==1){
			// all the boundary faces are associated to the same geometric
			// entity (must be refGeometry)
			N2GeomOnLayer [currentNode.getID()]=refGeometry;
			N2GeomInVolume[currentNode.getID()]=0;
		}
		// now the node is classified on a geom. vertex or a geom. curve
		else if (refGeometry->getDim()==0) { // on a vertex
			//we get the faces adjacent to the node and not marked layer_mark
			std::set<geom::GeomEntity* > local_classif;
			bool allMarked=true;
			for(unsigned int i=0;i<boundaryFaces.size();i++){
				Face fi=boundaryFaces[i];
				if(!m_mesh.isMarked(fi,layer_mark)){
					local_classif.insert(m_mesh.getGeometricClassification(fi));
					allMarked=false;
				}

			}
			if(allMarked){
				N2GeomOnLayer [currentNode.getID()]=refGeometry;
				N2GeomInVolume[currentNode.getID()]=0;
			}
			else if(local_classif.size()==1){
				N2GeomOnLayer [currentNode.getID()]=refGeometry;
				N2GeomInVolume[currentNode.getID()]=*(local_classif.begin());
			}
			else if(local_classif.size()==2){
				//local_classif contains all the boundary geometric surfaces that are not
				// expected to be inflated and they are 2 in this case
				std::cout<<"********************"<<std::endl;

				std::set<geom::GeomEntity* >::iterator itGeo = local_classif.begin();
				geom::GeomEntity* g0 = *itGeo; itGeo++;
				geom::GeomEntity* g1 = *itGeo;


				std::vector<Face> f_on_g0, f_on_g1;
				for(unsigned int i=0;i<boundaryFaces.size();i++){
					Face fi=boundaryFaces[i];
					geom::GeomEntity* g_fi= m_mesh.getGeometricClassification(fi);
					if(g_fi==g0)
						f_on_g0.push_back(fi);
					else if(g_fi==g1)
						f_on_g1.push_back(fi);
				}
				if(f_on_g0.size()==1 && f_on_g1.size()==1){
					std::vector<Node> commonNodes = m_mesh.getCommonNodes(f_on_g0[0],f_on_g1[0]);
					if(commonNodes.size()!=2)
						throw GMDSException("PILLOWING: A not handled configuration in 3 Geometric classification (9)");

					N2GeomOnLayer [currentNode.getID()]=refGeometry;

					if(commonNodes[0]==currentNode){
						N2GeomInVolume[currentNode.getID()]=m_mesh.getGeometricClassification(commonNodes[1]);
					}
					else
						N2GeomInVolume[currentNode.getID()]=m_mesh.getGeometricClassification(commonNodes[0]);
				}
				else
					throw GMDSException("PILLOWING: A not handled configuration in 3 Geometric classification (10)");

			}
			else
				throw GMDSException("PILLOWING: A not handled configuration in 3 Geometric classification (8)");

		}
		else if (refGeometry->getDim()==1) { // on a curve
			//std::cout<<"Node "<<currentNode.getID()<<" on curve"<<std::endl;
			//we get the faces adjacent to the node and not marked layer_mark
			std::set<geom::GeomEntity* > local_classif;
			bool allMarked=true;
			for(unsigned int i=0;i<boundaryFaces.size();i++){
				Face fi=boundaryFaces[i];
				if(!m_mesh.isMarked(fi,layer_mark)){
					local_classif.insert(m_mesh.getGeometricClassification(fi));
					allMarked=false;
				}
			}
			if(allMarked){
				//std::cout<<"\t get into the volume"<<std::endl;
				N2GeomOnLayer [currentNode.getID()]=refGeometry;
				N2GeomInVolume[currentNode.getID()]=0;
			}
			else if(local_classif.size()==1){
				//std::cout<<"\t get onto a surface"<<std::endl;
				N2GeomOnLayer [currentNode.getID()]=refGeometry;
				N2GeomInVolume[currentNode.getID()]=*(local_classif.begin());
			}
			else {
//					std::cout<<"Node "<<currentNode.getID()<<std::endl;
//					for(unsigned int i=0;i<boundaryFaces.size();i++){
//						Face *fi=boundaryFaces[i];
//						std::cout<<"\t boundary face["<<i<<"] : "<<fi.getID()<<std::endl;
//					}
				std::set<geom::GeomEntity* >::iterator itGeo = associatedGeometry.begin();
				geom::GeomEntity* g0 = *itGeo; itGeo++;
				geom::GeomEntity* g1 = *itGeo; itGeo++;
				geom::GeomEntity* g2 = *itGeo;
				geom::GeomEntity *g3, *g4;
				if(g0==refGeometry){
					g3=g1;
					g4=g2;
				}
				else if(g1==refGeometry){
					g3=g0;
					g4=g2;
				}
				else if(g2==refGeometry){
					g3=g0;
					g4=g1;
				}
				else if (refGeometry->getDim()==0) {
					/* currentNode is on a geom vertex incident to 3 surfaces
					 * we must get the boundary face among the 3 adj faces
					 */
					if(m_mesh.isMarked(boundaryFaces[0],layer_mark)){
						g3 = m_mesh.getGeometricClassification(boundaryFaces[1]);
						g4 = m_mesh.getGeometricClassification(boundaryFaces[2]);
					}
					else if(m_mesh.isMarked(boundaryFaces[1],layer_mark)){
						g3 = m_mesh.getGeometricClassification(boundaryFaces[0]);
						g4 = m_mesh.getGeometricClassification(boundaryFaces[2]);
					}
					else if(m_mesh.isMarked(boundaryFaces[2],layer_mark)){
						g3 = m_mesh.getGeometricClassification(boundaryFaces[0]);
						g4 = m_mesh.getGeometricClassification(boundaryFaces[1]);
					}
					else
						throw GMDSException("PILLOWING: A not handled configuration in 3 Geometric classification (4)");

				}
				else
					throw GMDSException("PILLOWING: A not handled configuration in 3 Geometric classification ");

				std::vector<Face> f_on_g3, f_on_g4;
				for(unsigned int i=0;i<boundaryFaces.size();i++){
					Face fi=boundaryFaces[i];
					geom::GeomEntity* g_fi= m_mesh.getGeometricClassification(fi);
					if(g_fi==g3)
						f_on_g3.push_back(fi);
					else if(g_fi==g4)
						f_on_g4.push_back(fi);
				}
				if(f_on_g3.size()==1 && f_on_g4.size()==1){
					std::vector<Node> commonNodes = m_mesh.getCommonNodes(f_on_g3[0],f_on_g4[0]);
					if(commonNodes.size()!=2)
						throw GMDSException("PILLOWING: A not handled configuration in 3 Geometric classification (3)");

					N2GeomOnLayer [currentNode.getID()]=refGeometry;

					if(commonNodes[0]==currentNode){
						N2GeomInVolume[currentNode.getID()]=m_mesh.getGeometricClassification(commonNodes[0]);
					}
					else
						N2GeomInVolume[currentNode.getID()]=m_mesh.getGeometricClassification(commonNodes[1]);
				}
				else
					throw GMDSException("PILLOWING: A not handled configuration in 3 Geometric classification (2)");

			}
		}
		else{
			std::cout<<"Node ID "<<currentNode.getID()<<std::endl;
			std::cout<<"All regions: "<<allRegions.size()<<std::endl;
			std::cout<<"All faces: "<<allFaces.size()<<std::endl;
			std::cout<<"nb boundary faces: "<<boundaryFaces.size()<<std::endl;
			std::cout<<"nb associated geometries: "<<associatedGeometry.size()<<std::endl;
			if(refGeometry==0)
				std::cout<<"geometrie NULL"<<std::endl;
			else
				std::cout<<"geometrie dim : "<<refGeometry->getDim()<<std::endl;
			throw GMDSException("PILLOWING: Error in the geometric classification of faces");
		}

	}

	// For each node on the pillowing surface we will associate a vector normal (outward) to the faces
	std::set<TCellID> pillowingFaces;
	for(unsigned int iFace=0; iFace<AFaces.size(); iFace++) {
                pillowingFaces.insert(AFaces[iFace].getID());
        }

	std::map<TCellID,bool> isFaceOutward;
	for(unsigned int iRegion=0; iRegion<ARegions.size(); iRegion++) {
		gmds::Region current_region = ARegions[iRegion];

		std::vector<gmds::Face> faces = current_region.get<gmds::Face>();
		for(unsigned int iFace=0; iFace<faces.size(); iFace++) {
			if(pillowingFaces.find(faces[iFace].getID()) != pillowingFaces.end()) {
				std::vector<gmds::Node> nodes = faces[iFace].get<gmds::Node>();

				isFaceOutward[faces[iFace].getID()] = current_region.isFaceOrientedOutward(nodes);
			}
		}
	
	}
	std::map<TCellID, gmds::math::Vector> N2NormalVector;
	std::map<TCellID, TCoord> N2VectorFactor;

	for(TCellID iFace=0; iFace<AFaces.size(); iFace++) {

		std::vector<gmds::Node> nodes = AFaces[iFace].get<gmds::Node>();

		TCoord minEdgeLength = HUGE_VAL;
		for(unsigned int iNode=0; iNode<nodes.size(); iNode++) {

			TCoord edgeLength = nodes[iNode].getPoint().distance(nodes[(iNode+1)%nodes.size()].getPoint());
			if(edgeLength < minEdgeLength) {
				minEdgeLength = edgeLength;
			}
		}

		for(unsigned int iNode=0; iNode<nodes.size(); iNode++) {

			gmds::math::Vector normal = AFaces[iFace].normal(nodes[iNode]);
                        if(!isFaceOutward[AFaces[iFace].getID()]) {
                                normal = normal * (-1.);
                        }

			if(N2NormalVector.find(nodes[iNode].getID()) == N2NormalVector.end()) {
				N2NormalVector[nodes[iNode].getID()] = normal;
				N2VectorFactor[nodes[iNode].getID()] = minEdgeLength;
			} else {
				N2NormalVector[nodes[iNode].getID()] = N2NormalVector[nodes[iNode].getID()] + normal;
				if(minEdgeLength<N2VectorFactor[nodes[iNode].getID()]) {
					N2VectorFactor[nodes[iNode].getID()] = minEdgeLength;
				}
			}
		}
	}
	for(std::map<TCellID, gmds::math::Vector>::iterator it = N2NormalVector.begin(); it != N2NormalVector.end(); it++) {
		(it->second).normalize();
	} 

	//-----------------------------------------------------------------------
	//------------------- STEP 2 --------------------------------------------
	//-----------------------------------------------------------------------
	// now we create the new nodes, new faces and new regions to insert
	std::map<TCellID,Node> oldNodeToNewNodeR;
	std::map<TCellID,Node> oldNodeToNewNodeL;
	std::map<TCellID,Face> oldFaceToNewFaceR;
	std::map<TCellID,Face> oldFaceToNewFaceL;
	std::map<TCellID,Region> oldFaceToNewRegion;
	std::map<FakeEdge::EdgeID, Face > newLateralFaces;

	//new nodes creation
	for(std::set<TCellID>::iterator it_setLayerNodes=layerNodes.begin();
			it_setLayerNodes!=layerNodes.end();it_setLayerNodes++)
	{
		Node oldNode = m_mesh.get<Node>(*it_setLayerNodes);
		Node newNodeR = m_mesh.newNode(oldNode.getPoint());
		oldNodeToNewNodeR[oldNode.getID()]= newNodeR;

		
		gmds::math::Point point = oldNode.getPoint();
		if(ADisplacement) {
			// this internal node is created with a small displacement from the original node
			point = point + (-N2VectorFactor[oldNode.getID()] * 0.1) * N2NormalVector[oldNode.getID()];
		}
		Node newNodeL = m_mesh.newNode(point);
		m_newNodesL.insert(newNodeL);

		oldNodeToNewNodeL[oldNode.getID()]= newNodeL;

		//geometric classification of new nodes
//			if(N2GeomOnLayer[oldNode]->getDim()==1)
//				std::cout<<"\t NOEUD DROIT SUR UNE COURBE"<<std::endl;
		m_mesh.setGeometricClassification(newNodeL,N2GeomInVolume[oldNode.getID()]);
		m_mesh.setGeometricClassification(newNodeR,N2GeomOnLayer[oldNode.getID()]);

	}

	//new opposite faces, lateral faces and regions creation
	for(std::vector<Face>::iterator it_setLayerFaces=AFaces.begin();
			it_setLayerFaces!=AFaces.end();it_setLayerFaces++)
	{
		Face oldFace = *it_setLayerFaces;
		std::vector<Node> oldNodes = oldFace.get<Node>();

		std::vector<Node> newNodesR;
		std::vector<Node> newNodesL;
		newNodesR.resize(oldNodes.size());
		newNodesL.resize(oldNodes.size());
		for(unsigned int inodes=0;inodes<oldNodes.size();inodes++){
			newNodesR[inodes] = oldNodeToNewNodeR[oldNodes[inodes].getID()];
			newNodesL[inodes] = oldNodeToNewNodeL[oldNodes[inodes].getID()];
		}

		//new faces
		Face newFaceR= m_mesh.newFace(newNodesR);
		return_faces.push_back(newFaceR);
		Face newFaceL= m_mesh.newFace(newNodesL);

		// the right face is classifed as the old one
		m_mesh.setGeometricClassification(newFaceR,m_mesh.getGeometricClassification(oldFace));
		updateGeomClassification(newFaceL);

		if(m_mesh.getGeometricClassification(newFaceL)==0){}
		else if (m_mesh.getGeometricClassification(newFaceL)->getDim()!=2){}
		else
			std::cout<<"==== A NEW LEFT FACE IS ON A SURFACE"<<std::endl;

		oldFaceToNewFaceR[oldFace.getID()]= newFaceR;
		oldFaceToNewFaceL[oldFace.getID()]= newFaceL;

		//we know that we have hexes, so the following code is dedicated
		// to hexes

		// determine whether oldFace is oriented outward from the input region
		std::vector<Region> regions = oldFace.get<Region>();
		Region old_region;
		if(regions.size() == 0) {
			throw GMDSException("Pilowwing problem.");
		} else if (regions.size() == 1) {
			old_region = regions[0];
		} else {
			if(m_mesh.isMarked(regions[0],pillow_mark)) {
				old_region = regions[0];
			} else {
				old_region = regions[1];
			}
		}
		bool isFaceOutward = old_region.isFaceOrientedOutward(oldFace.get<Node>());

		Node n0 = newNodesR[0];
		Node n1 = newNodesR[1];
		Node n2 = newNodesR[2];
		Node n3 = newNodesR[3];
		Node n4 = newNodesL[0];
		Node n5 = newNodesL[1];
		Node n6 = newNodesL[2];
		Node n7 = newNodesL[3];
		Region r;
		if(!isFaceOutward) r = m_mesh.newHex(n0, n1, n2, n3, n4, n5, n6, n7);
		else r = m_mesh.newHex(n4, n5, n6, n7, n0, n1, n2, n3);
		ANewRegions.push_back(r);

		oldFaceToNewRegion[oldFace.getID()]= r;

		//lateral faces
		Face f01, f12, f23, f30;
		FakeEdge fe01(n0.getID(),n1.getID());
		if( newLateralFaces.find(fe01.getID())!=newLateralFaces.end())
			// the face create from n0 and n1 already exist
			f01 = newLateralFaces[fe01.getID()];
		else{
			f01 = m_mesh.newQuad(n0,n1,n5,n4);
			updateGeomClassification(f01);
			newLateralFaces[fe01.getID()] = f01;
		}
		FakeEdge fe12(n1.getID(),n2.getID());
		if(newLateralFaces.find(fe12.getID())!=newLateralFaces.end())
			// the face create from n1 and n2 already exist
			f12 = newLateralFaces[fe12.getID()];
		else{
			f12 = m_mesh.newQuad(n1,n2,n6,n5);
			updateGeomClassification(f12);
			newLateralFaces[fe12.getID()] = f12;
		}
		FakeEdge fe23(n2.getID(),n3.getID());
		if(newLateralFaces.find(fe23.getID())!=newLateralFaces.end())
			// the face create from n2 and n3 already exist
			f23 = newLateralFaces[fe23.getID()];
		else{
			f23 = m_mesh.newQuad(n2,n3,n7,n6);
			updateGeomClassification(f23);
			newLateralFaces[fe23.getID()] = f23;
		}
		FakeEdge fe30(n3.getID(),n0.getID());
		if(newLateralFaces.find(fe30.getID())!=newLateralFaces.end())
			// the face create from n3 and n0 already exist
			f30 = newLateralFaces[fe30.getID()];
		else{
			f30 = m_mesh.newQuad(n3,n0,n4,n7);
			updateGeomClassification(f30);
			newLateralFaces[fe30.getID()] = f30;
		}

		//R2N is done for r
		//N2R
		n0.add<Region>(r);
		n1.add<Region>(r);
		n2.add<Region>(r);
		n3.add<Region>(r);
		n4.add<Region>(r);
		n5.add<Region>(r);
		n6.add<Region>(r);
		n7.add<Region>(r);
		//R2F
		r.add<Face>(newFaceL);
		r.add<Face>(newFaceR);
		r.add<Face>(f01);
		r.add<Face>(f12);
		r.add<Face>(f23);
		r.add<Face>(f30);
		//F2R
		newFaceL.add<Region>(r);
		newFaceR.add<Region>(r);
		f01.add<Region>(r);
		f12.add<Region>(r);
		f23.add<Region>(r);
		f30.add<Region>(r);
	}

	//-----------------------------------------------------------------------
	//------------------- STEP 5 --------------------------------------------
	//-----------------------------------------------------------------------

	//update R2N and F2R for outer faces
	for(unsigned int i=0;i<AFaces.size();i++){
		Face current_face = AFaces[i];
		Face current_faceR= oldFaceToNewFaceR[current_face.getID()];
		Face current_faceL= oldFaceToNewFaceL[current_face.getID()];
		Region new_region = oldFaceToNewRegion[current_face.getID()];


		std::vector<Node  >   current_nodes  = current_face.get<Node>();
		std::vector<Region> current_regions= current_face.get<Region>();
		Region current_region_right, current_region_left;
		if(current_regions.size()==2){
			if(m_mesh.isMarked(current_regions[0],pillow_mark))
			{
				current_region_right = current_regions[1];
				current_region_left = current_regions[0];
			}
			else
			{
				current_region_right = current_regions[0];
				current_region_left = current_regions[1];
			}
			// keep geometry association
			geom::GeomEntity* ge=  m_mesh.getGeometricClassification(current_region_left);
			m_mesh.setGeometricClassification(new_region,ge);
		}
		else if(current_regions.size()==1){
			current_region_right = current_regions[0];
		}
		else
			throw GMDSException("Topological unconsistance");

		//connection NOW
		if(current_region_left.getID()!=NullID){
			//R2F
			current_region_left.replace(current_face,current_faceL);
			//F2R
			current_faceL.add<Region>(current_region_left);
		}
		if(current_region_right.getID()!=NullID){
			//R2F
			current_region_right.replace(current_face,current_faceR);
			//F2R
			current_faceR.add<Region>(current_region_right);
		}
	}

	//-----------------------------------------------------------------------
	//------------------- STEP 6 --------------------------------------------
	//-----------------------------------------------------------------------
	// finally N2R, R2N and F2R
	for(std::set<TCellID>::iterator it_setLayerNodes=layerNodes.begin();
			it_setLayerNodes!=layerNodes.end();it_setLayerNodes++)
	{
		Node oldNode = m_mesh.get<Node>(*it_setLayerNodes);
		Node newNodeL = oldNodeToNewNodeL[oldNode.getID()];
		Node newNodeR = oldNodeToNewNodeR[oldNode.getID()];
		std::vector<Region> incRegions = oldNode.get<Region>();
		for(unsigned int rI=0;rI<incRegions.size();rI++){
			Region r = incRegions[rI];
			std::vector<Face> facesOfR = r.get<Face>();
			if(!m_mesh.isMarked(r,pillow_mark)){
				r.replace(oldNode,newNodeR);
				newNodeR.add<Region>(r);
				for(unsigned int fI=0;fI<facesOfR.size();fI++){
					Face f = facesOfR[fI];
					f.replace(oldNode,newNodeR);
				}
			}
			else{
				r.replace(oldNode,newNodeL);
				newNodeL.add<Region>(r);
				for(unsigned int fI=0;fI<facesOfR.size();fI++){
					Face f = facesOfR[fI];
					f.replace(oldNode,newNodeL);
				}
			}


		}
	}
	//-----------------------------------------------------------------------
	//------------------- STEP 7 --------------------------------------------
	//-----------------------------------------------------------------------

	// Now we can remove the nodes and the faces of the given surface;
	for(std::set<TCellID>::iterator it_setLayerNodes=layerNodes.begin();
			it_setLayerNodes!=layerNodes.end();it_setLayerNodes++)
		m_mesh.deleteNode(*it_setLayerNodes);
	for(unsigned int iF=0;iF<AFaces.size();iF++)
		m_mesh.deleteFace(AFaces[iF]);
//
	/* free mark, but all the marked cell (nodes and faces in AFaces have
	 * been removed.
	 */
	m_mesh.unmarkAll<Face>(layer_mark);
	m_mesh.unmarkAll<Region>(pillow_mark);

	m_mesh.freeMark<Face>(layer_mark);
	m_mesh.freeMark<Region>(pillow_mark);

	return return_faces;

}
//------------------------------------------------------------------------//
void SheetOperator::
updateGeomClassification(Face& AFace)
{
	std::vector<Node> fnodes = AFace.get<Node>();

	geom::GeomEntity* ge_surf=0;
	for(unsigned int ni=0;ni<fnodes.size();ni++){
		Node n = fnodes[ni];
		geom::GeomEntity* ge=  m_mesh.getGeometricClassification(n);
		if(ge==0){
			m_mesh.setGeometricClassification(AFace,0);
			return; // on est dans le volume
		}
		if(ge->getDim()==2)
			ge_surf=ge;
	}

	m_mesh.setGeometricClassification(AFace,ge_surf);
}
/*----------------------------------------------------------------------------*/
std::set<Node>
SheetOperator::getNewInternalNodes()
{
	return m_newNodesL;
}
/*----------------------------------------------------------------------------*/
