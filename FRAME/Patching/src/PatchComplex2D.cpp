/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: F. Ledoux (2015)
 *
 * franck.ledoux@cea.fr
 *
 * The FRAME software is a computer program whose purpose is to provide a set
 * of algorithms to build 2D and 3D meshes using frame field concept. The main
 * focus of these algorithms is quadrilateral and hexahedral meshing.
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
 * PatchComplex2D.cpp
 *
 *  Created on: Sept 11, 2015
 *      Author: F. Ledoux
 */
/*----------------------------------------------------------------------------*/
// GMDS Headers
#include <GMDS/IO/VTKWriter.h>
#include <GMDS/Math/Line.h>
/*----------------------------------------------------------------------------*/
// STL Headers
#include <sstream>
/*----------------------------------------------------------------------------*/
// Patching Headers
#include <PatchComplex2D.h>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
PatchComplex2D::PatchComplex2D(IGMesh* AMesh)
: m_mesh(AMesh)
{
}
/*----------------------------------------------------------------------------*/
PatchComplex2D::~PatchComplex2D()
{
    std::vector<PatchVertex2D*>::iterator it_v = m_vertices.begin();
    for(;it_v!=m_vertices.end();it_v++){
        if(*it_v != NULL)
            delete *it_v;
    }
    std::vector<PatchEdge2D*>::iterator it_e = m_edges.begin();
    for(;it_e!=m_edges.end();it_e++){
        if(*it_e != NULL)
            delete *it_e;
    }
    std::vector<Patch2D*>::iterator it_f = m_patchs.begin();
    for(;it_f!=m_patchs.end();it_f++){
        if(*it_f != NULL)
            delete *it_f;
    }
}
/*----------------------------------------------------------------------------*/
IGMesh* PatchComplex2D::mesh()
{
    return m_mesh;
}
/*----------------------------------------------------------------------------*/
Patch2D* PatchComplex2D::newPatch(std::vector<PatchEdge2D*>& AEdges)
{
    Patch2D* p = new Patch2D(this,AEdges);
    m_patchs.push_back(p);
    return p;
}
/*----------------------------------------------------------------------------*/
PatchEdge2D* PatchComplex2D::newEdge(PatchVertex2D* AV1, PatchVertex2D* AV2)
{
    PatchEdge2D* e = new PatchEdge2D(AV1,AV2);
    m_edges.push_back(e);
    return e;
}
/*----------------------------------------------------------------------------*/
PatchVertex2D* PatchComplex2D::newVertex(const gmds::math::Point& ALocation,
                                         const PatchVertex2D::Type& AType)
{
    PatchVertex2D* v = new PatchVertex2D(ALocation,AType);
    m_vertices.push_back(v);
    return v;
}
/*----------------------------------------------------------------------------*/
void PatchComplex2D::write(const std::string& AFileName)
{
    throw GMDSException("TO IMPLEMENT");
//    static int out = 0;
//    std::stringstream file_name;
//    file_name << AFileName << "_" << out;
//    gmds::IGMesh m(DIM3 | F | N | F2N);
//    
//    std::map<gmds::TCellID, Patch2D >::iterator it;
//    for(it=m_patchs.begin();it!=m_patchs.end();it++){
//        Patch2D current = it->second;
//        
//        //Add the patch in the output mesh
//        std::vector<math::Point> corners = current.corners();
//        Node n1 = m.newNode(corners[0]);
//        Node n2 = m.newNode(corners[1]);
//        Node n3 = m.newNode(corners[2]);
//        Node n4 = m.newNode(corners[3]);
//        m.newQuad(n1,n2,n3,n4);
//    }
//    
//    VTKWriter<IGMesh> writer(m);
//    writer.write(file_name.str(), DIM3 | F | N);
//    out++;
}

///*----------------------------------------------------------------------------*/
//void PatchComplex2D::alignPatch(const gmds::TCellID AID)
//{
//    Face    current_face  = m_mesh->get<Face>(AID);
//    Patch2D current_patch = m_patchs[AID];
//
//    std::vector<Node> current_nodes = current_face.get<Node>();
//    std::set<TCellID> adj_face_ids;
//
//    // We get the IDs of all the faces sharing a node with current_face,
//    // including current_face
//    for(unsigned int i=0; i<current_nodes.size();i++){
//        std::vector<TCellID> local_face_ids = current_nodes[i].getIDs<Face>();
//        adj_face_ids.insert(local_face_ids.begin(), local_face_ids.end());
//    }
//
//    // the face current_face is removed now
//    adj_face_ids.erase(AID);
//
//    //now the patch is really aligned with its adjacent patch
//    std::cout<<"TO ALIGN WITH: "<<adj_face_ids.size()<<std::endl;
//    //========================================================================
//    // We detect edges to move
//    //========================================================================
//    //Point are supposed to be ordered
//    // We traverse the 4 edges of the patch
//
//    std::vector<math::Point>& current_corners = current_patch.corners();
//    for(unsigned int i=0; i<4; i++){
//
//        // [p1,p2] current edge
//        math::Point p1 = current_corners[i];
//        math::Point p2 = current_corners[(i+1)%4];
//        // [p1opp,p2opp] opposite edge
//        math::Point p2opp = current_corners[(i+2)%4];
//        math::Point p1opp = current_corners[(i+3)%4];
//
//        double distSegment = math::min2(p1.distance(p1opp),p2.distance(p2opp));
//        double dmin = distSegment/4;
//        std::cout<<i<<") dmin = "<<dmin<<std::endl;
//        math::Vector v12(p1,p2);
//        v12.normalize();
//        math::Line l12(p1,p2);
//
//        //We traverse all the adjacent patch to perform some corrections
//        std::set<TCellID>::iterator it;
//        for(it=adj_face_ids.begin();it!=adj_face_ids.end();it++){
//            if(m_patchs.find(*it)==m_patchs.end())
//                continue;
//            Patch2D& pi = m_patchs[*it];
//            std::vector<math::Point>& pi_corners = pi.corners();
//
//            //We traverse the four edges of each patch now
//            for(unsigned int k=0; k<4; k++){
//
//                math::Point p3 = pi_corners[k];
//                math::Point p4 = pi_corners[(k+1)%4];
//                math::Vector v34(p3,p4);
//                v34.normalize();
//
//                if(math::near(fabs(v12.dot(v34)),1,0.2)) {
//
//                    //almost parallel edges
//                    math::Segment s34(p3,p4);
//                    //we compute the orthogonal distance
//                    math::Point onL12, onS34;
//                    TCoord dist = l12.distance2D(s34, onL12, onS34);
//                    std::cout<<"   --> almost aligned at "<<dist<<std::endl;
//                    if(dist<dmin && dist!=0.0){
//                        //lines to align
//                        //we look for a medial line between our two segments.
//                        math::Point line_p1, line_p2;
//                        if(v12.dot(v34)>0){
//                            line_p1 = 0.5*(p1 + p3);
//                            line_p2 = 0.5*(p2 + p4);
//                        }
//                        else {
//                            line_p1 = 0.5*(p1 + p4);
//                            line_p2 = 0.5*(p2 + p3);
//
//                        }
//                        math::Line medial(line_p1,line_p2);
//                     //   math::Line medial(p3,p4);
//                        std::cout<<"\t ALIGNEMENT "<<std::endl;
//                        current_patch.setCorner(i, medial.project(p1));
//                        current_patch.setCorner((i+1)%4, medial.project(p2));
////                        std::cout<<p3<<" -- "<<medial.project(p3)<<std::endl;
//                        pi.setCorner(k, medial.project(p3));
//                       pi.setCorner((k+1)%4, medial.project(p4));
//                        m_patchs[*it] = pi;
//                        m_patchs[AID] = current_patch;
//
//                    }
//                }//if(math::near(fabs(v12.dot(v34)),1,0.2))
//
//            }//for(unsigned int k=0; k<4; k++)
//
//        }//for (unsigned int j=0; j<adj_patchs.size(); j++)
//
//    }//for(unsigned int i=0; i<4; i++)
//}


/*----------------------------------------------------------------------------*/
