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
 * Patching2D.cpp
 *
 *  Created on: June 17, 2015
 *      Author: F. Ledoux
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IG.h>
#include <GMDS/IG/IGMeshDoctor.h>
#include <GMDS/IO/MeditReader.h>
#include <GMDS/IO/VTKReader.h>
#include <GMDS/IO/VTKWriter.h>
#include <GMDS/Math/Cross.h>
#include <GMDS/Math/Quaternion.h>
#include <GMDS/Math/Numerics.h>
#include <GMDS/Math/Chart.h>
#include <GMDS/Math/Ray.h>
#include <GMDS/Math/Triangle.h>
#include <GMDS/Math/Line.h>
#include <GMDS/Math/Segment.h>
#include <GMDS/Algo/BoundaryOperator.h>
/*----------------------------------------------------------------------------*/
#include "Patching2D.h"
#include "Patch2D.h"
/*----------------------------------------------------------------------------*/
#include <sstream>
#include <set>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
Patching2D::
Patching2D(IGMesh* AMesh, Variable<math::Cross2D>* AField):
m_mesh(AMesh), m_complex(AMesh),
m_quad_mesh( IGMesh(DIM3 | F | N | F2N | N2F) ),
m_field(AField), m_output_directory_name("")
{
    m_target_size = 0;
}
/*----------------------------------------------------------------------------*/
Patching2D::~Patching2D()
{}
/*----------------------------------------------------------------------------*/
void Patching2D::execute()
{
    //==================================================================
    // Boolean marks initialization
    //==================================================================
    m_mark_faces_with_sing = m_mesh->getNewMark<Face>();
    
    
    //==================================================================
    // GEOMETRY VARIABLE FOR DEBUG
    //==================================================================
    Variable<int>* geom_var = m_mesh->newVariable<int>(GMDS_NODE, "geometry");
    
    
    m_index = m_mesh->newVariable<int>(GMDS_FACE, "index");
    
    IGMesh::node_iterator it_nodes = m_mesh->nodes_begin();
    int nb_on_curves=0;
    for (; !it_nodes.isDone(); it_nodes.next())  {
        Node n = it_nodes.value();
        if (m_mesh->isMarked(n, m_mark_nodes_on_point)){
            m_seed_dims.push_back(0);
            m_seed_ids.push_back(n.getID());
            (*geom_var)[n.getID()] = 2;
        }
        else if (m_mesh->isMarked(n,m_mark_nodes_on_curve)){
            (*geom_var)[n.getID()] = 1;
            nb_on_curves++;
        }
        else
            (*geom_var)[n.getID()] = 0;
    }
    VTKWriter<IGMesh> writer(*m_mesh);
    writer.write("geom_mesh", DIM3 | F | N );
    std::cout << "\t DONE" << std::endl;
    
    
    //========================================================================
    // STEP 1 - Detection of singular triangles and storage
    //========================================================================
    std::cout << "============================================================="
    << std::endl;
    std::cout << "Detection of singular triangles" << std::endl;
    detectSingularTriangles();
    writer.write("singular_triangles", DIM3 | F | N | R);
    std::cout << "\t DONE" << std::endl;
    //========================================================================
    // STEP 3 - Growing expansion of patchs
    //========================================================================
    std::cout << "============================================================="
    << std::endl;
    std::cout << "Growing expansion of patchs" << std::endl;
    growPatchs();
    std::cout << "\t DONE" << std::endl;
    
//    //========================================================================
//    // STEP 3 - Individual patch building
//    //========================================================================
//    std::cout << "============================================================="
//    << std::endl;
//    std::cout << "Patch building for regular triangles" << std::endl;
//    buildRegularPatchs();
//    std::cout << "\t DONE" << std::endl;
//    
//    //========================================================================
//    // STEP 4 - Patch connection - all the patches that overlap each others
//    // are going to be connected
//    //========================================================================
//    std::cout << "============================================================="
//    << std::endl;
//    std::cout << "Overlapping Patch connection building" << std::endl;
//    m_patchs.connectPatches();
//    std::cout << "\t DONE" << std::endl;
//    
//    
//    //========================================================================
//    // STEP 5 - build boundary lines of patches
//    //========================================================================
//    std::cout << "============================================================="
//    << std::endl;
//    std::cout << "Boundary patch building" << std::endl;
//    computePatchBoundaryPoints();
//    std::cout << "\t DONE" << std::endl;
//    //========================================================================
//    // STEP 6 - build the final mesh
//    //========================================================================
//    std::cout << "============================================================="
//    << std::endl;
//    std::cout << "Final Grid building" << std::endl;
//    buildMesh();
//    std::cout << "\t DONE" << std::endl;
//    
    //========================================================================
    // Boolean marks cleaning
    //========================================================================
    m_mesh->unmarkAll<Face>(m_mark_faces_with_sing);
    m_mesh->freeMark<Face>(m_mark_faces_with_sing);
    
}
/*---------------------------------------------------------------------------*/
void Patching2D::growPatchs(){
    
}
/*---------------------------------------------------------------------------*/
void Patching2D::computeTargetSize()
{
    double min_size = 1000000;
    IGMesh::face_iterator it_faces = m_mesh->faces_begin();
    
    for (; !it_faces.isDone(); it_faces.next()) {
        Face current = it_faces.value();
        std::vector<Node> nodes = current.get<Node>();
        
        math::Point p0 = nodes[0].getPoint();
        math::Point p1 = nodes[1].getPoint();
        math::Point p2 = nodes[2].getPoint();
        double d01 = p0.distance(p1);
        double d02 = p0.distance(p2);
        double d12 = p1.distance(p2);
        double current_min = math::min3(d01,d02,d12);
        if(current_min<min_size)
            min_size = current_min;
    }
    
    m_target_size = min_size/8.0;
    std::cout<<"===> Computed target size: "<<m_target_size<<std::endl;
}
/*---------------------------------------------------------------------------*/
void Patching2D::detectSingularTriangles()
{
    IGMesh::face_iterator it_faces = m_mesh->faces_begin();
    m_singular_triangles.clear();
    
    for (; !it_faces.isDone(); it_faces.next()) {
        Face current = it_faces.value();
        std::vector<TCellID> nodeIDs = current.getIDs<Node>();
        int ID1 = nodeIDs[0];
        int ID2 = nodeIDs[1];
        int ID3 = nodeIDs[2];
        
        math::Cross2D c1 = (*m_field)[ID1];
        math::Cross2D c2 = (*m_field)[ID2];
        math::Cross2D c3 = (*m_field)[ID3];
        
        int index = math::Cross2D::index(c1,c2,c3);
        
        (*m_index)[current.getID()] = index;
        
        if (index ==1 || index == -1 ) {
            m_singular_triangles.push_back(current);
            
            m_seed_dims.push_back(2);
            m_seed_ids.push_back(current.getID());
            
            m_mesh->mark(current, m_mark_faces_with_sing);
            //all the faces sharing a node with it are marked too !!!!
            
            std::vector<Node> nodes = current.get<Node>();
            for(unsigned int i_node = 0; i_node<nodes.size(); i_node++) {
                Node ni = nodes[i_node];
                std::vector<Face> ni_faces = ni.get<Face>();
                for(unsigned i_face=0; i_face<ni_faces.size(); i_face++){
                    m_mesh->mark(ni_faces[i_face], m_mark_faces_with_sing);
                }
            }
        }
        
    } //for (; !it_regions.isDone(); it_regions.next())
    
    std::cout << "Nb singular triangles = " << m_singular_triangles.size() << std::endl;
}
/*---------------------------------------------------------------------------*/
void Patching2D::buildRegularPatchs()
{
    double alignment_tol = 0.001;
    const int nb_faces = m_mesh->getNbFaces();
    IGMesh::face_iterator it_faces = m_mesh->faces_begin();
    int compt=0;
    
    for (; !it_faces.isDone(); it_faces.next()) {
        Face current = it_faces.value();
        compt++;
        double percentage = 100.0*compt/nb_faces;
        std::cout<<percentage<<"%"<<std::endl;
        
        //Singular faces are not taken into account now
        if(m_mesh->isMarked(current,m_mark_faces_with_sing))
            continue;
        
        std::vector<Node> nodes = current.get<Node>();
        
        math::Cross2D c[3];
        c[0] = (*m_field)[nodes[0].getID()];
        c[1] = (*m_field)[nodes[1].getID()];
        c[2] = (*m_field)[nodes[2].getID()];
        
        //========================================================================
        //we build a map which provides for each node the vector to consider with
        // an adjacent node
        //========================================================================
        std::map<TCellID, std::map<TCellID, math::Vector> > vec_map;
        for(int i=0;i<3;i++) {
            int j=(i+1)%3;
            int k=(i==0)?2:i-1;
            
            math::Point pi = nodes[i].getPoint();
            math::Point pj = nodes[j].getPoint();
            math::Point pk = nodes[k].getPoint();
            bool align_ij=false;
            //===============================================
            //DIRECTION 1 FROM I TO J
            math::Vector vij(pi,pj);
            vij.normalize();
            std::vector<math::Vector> ci_vectors = c[i].componentVectors();
            std::vector<math::Vector> ci_candidates;
            math::Vector dir;
            for(unsigned int x=0; x<ci_vectors.size(); x++) {
                math::Vector current_vec = ci_vectors[x];
                current_vec.normalize();
                double val = current_vec.dot(vij);
                if(math::near(val-1,0,alignment_tol)){
                    //aligned
                    math::Point witness = pi+current_vec;
                    bool left_w = witness.isStrictlyOnLeft2D(pi,pj);
                    bool left_k = pk.isStrictlyOnLeft2D(pi,pj);
                    if(left_w!=left_k)
                        ci_candidates.push_back(current_vec);
                    else
                        ci_candidates.push_back(vij);
                }
                else if(val>0){
                    ci_candidates.push_back(current_vec);
                }
            }
            if(ci_candidates.size()>2 || ci_candidates.size()<1)
                throw GMDSException("Error in patch cross building");
            else if(ci_candidates.size()==1)
                dir = ci_candidates[0];
            else {
                math::Vector c0 = ci_candidates[0];
                math::Vector c1 = ci_candidates[1];
                if(math::near(c0.dot(vij)-1,0,alignment_tol)){
                    dir=c0;
                    align_ij=true;
                }
                else if(math::near(c1.dot(vij)-1,0,alignment_tol)){
                    dir = c1;
                    align_ij=true;
                }
                else 	{
                    //we need the one on the right side
                    math::Point witness = pi+c0;
                    bool left_w = witness.isStrictlyOnLeft2D(pi,pj);
                    bool left_k = pk.isStrictlyOnLeft2D(pi,pj);
                    if(left_w!=left_k)
                        dir = c0;
                    else
                        dir = c1;
                }
            }
            
            vec_map[nodes[i].getID()][nodes[j].getID()] = dir;
            //===============================================
            //DIRECTION 2 FROM J TO I
            if(align_ij){
                vec_map[nodes[j].getID()][nodes[i].getID()] = dir.opp();
            }
            else {
                math::Vector vji(pj,pi);
                vji.normalize();
                std::vector<math::Vector> cj_vectors = c[j].componentVectors();
                std::vector<math::Vector> cj_candidates;
                for(unsigned int x=0; x<cj_vectors.size(); x++) {
                    math::Vector current_vec = cj_vectors[x];
                    current_vec.normalize();
                    double val = current_vec.dot(vji);
                    if(math::near(val-1,0,alignment_tol)){
                        //aligned
                        math::Point witness = pj+current_vec;
                        bool left_w = witness.isStrictlyOnLeft2D(pj,pi);
                        bool left_k = pk.isStrictlyOnLeft2D(pj,pi);
                        if(left_w!=left_k)
                            cj_candidates.push_back(current_vec);
                        else
                            cj_candidates.push_back(vji);
                    }
                    else if(val>0){
                        cj_candidates.push_back(current_vec);
                    }
                }
                if(cj_candidates.size()>2 || cj_candidates.size()<1)
                    throw GMDSException("Error in patch cross building");
                else if(cj_candidates.size()==1)
                    dir = cj_candidates[0];
                else {
                    math::Vector c0 = cj_candidates[0];
                    math::Vector c1 = cj_candidates[1];
                    if(math::near(c0.dot(vji)-1,0,alignment_tol)){
                        dir=c0;
                    }
                    else if(math::near(c1.dot(vji)-1,0,alignment_tol)){
                        dir = c1;
                    }
                    else
                    {
                        //we need the one on the right side
                        math::Point witness = pi+c0;
                        bool left_w = witness.isStrictlyOnLeft2D(pi,pj);
                        bool left_k = pk.isStrictlyOnLeft2D(pi,pj);
                        if(left_w!=left_k)
                            dir = c0;
                        else
                            dir = c1;
                    }
                }
                
                
                vec_map[nodes[j].getID()][nodes[i].getID()] = dir;
            }
            
        } //for(int i=0;i<3;i++)
        
        //========================================================================
        // Now, we build the patch structure associate to this triangle
        //========================================================================
        std::vector<math::Point> patch_corners;
        math::Point p0 = nodes[0].getPoint();
        math::Point p1 = nodes[1].getPoint();
        math::Point p2 = nodes[2].getPoint();
        double d01 = p0.distance(p1);
        double d02 = p0.distance(p2);
        double d12 = p1.distance(p2);
        double tol = 0.1*math::min3(d01,d02,d12);
        for(int i=0;i<3;i++) {
            
            int j = (i+1)%3;
            int k = (i==0)?2:i-1;
            
            TCellID id_i = nodes[i].getID();
            TCellID id_j = nodes[j].getID();
            TCellID id_k = nodes[k].getID();
            
            math::Point pi = nodes[i].getPoint();
            math::Point pj = nodes[j].getPoint();
            math::Point pk = nodes[k].getPoint();
            
            math::Vector vij = vec_map[id_i][id_j];
            math::Vector vik = vec_map[id_i][id_k];
            //    std::cout<<"In "<<id_i<<" with "<<vij<<" and "<<vik<<" give "<<vij.dot(vik)<<std::endl;
            bool colinear = math::near(fabs(vij.dot(vik)),1.0,alignment_tol);
            if(!colinear){
                //	std::cout<<"Add corner for "<<id_i<<std::endl;
                tryToAdd(pi, patch_corners, tol);
            }
            if(!m_mesh->isMarked(nodes[i], m_mark_nodes_on_curve) ||
               !m_mesh->isMarked(nodes[j], m_mark_nodes_on_curve)) {
                //do we have an intersection point with the next vector
                math::Vector vji = vec_map[id_j][id_i];
                //Warning numeric instabilities can lead to compare the wrong vector
                // We eliminate such degenerate cases here
                math::Vector edge_ij(pi,pj);
                if( math::near(fabs(vij.dot(edge_ij)),0,0.01) && 
                   math::near(fabs(vji.dot(edge_ij)),0,0.01) ) {
                    tryToAdd(pi, patch_corners, tol);
                    tryToAdd(pj, patch_corners, tol);
                }
                else {
                    colinear = math::near(fabs(vij.dot(vji)),1.0,alignment_tol);
                    if(!colinear) {
                        math::Ray ri(pi,vij);
                        math::Ray rj(pj,vji);
                        math::Point p_ij;
                        //	    std::cout<<i<<" - Rays ("<<ri.getPoint()<<", "<<ri.getPoint()+ri.getDir()<<") and ("
                        //		     <<rj.getPoint()<<", "<<ri.getPoint()+rj.getDir()<<") ";
                        if(ri.intersect2D(rj,p_ij)){
                            //	      std::cout<<"intersect"<<std::endl;
                            tryToAdd(p_ij, patch_corners, tol);
                        }
                    }
                }
            }
        }
        if(patch_corners.size()<4){
            std::cout<<"Wrong number of corners: "<<patch_corners.size()<<std::endl;
            throw GMDSException("Impossible to build a patch without 4 corners");
        }
        
        
        //========================================================================
        // We build the patch grid
        //========================================================================
        buildLocalPatch(current, patch_corners);
        
    } //for (; !it_faces.isDone(); it_faces.next())
    
}
/*---------------------------------------------------------------------------*/
  bool Patching2D::tryToAdd(const math::Point& AP, 
			    std::vector<math::Point>& AV,
			    const double& ATol) const
  {
    double tol = ATol;
    bool found = false;
  for(unsigned int i=0;i<AV.size() && !found; i++){
    math::Point pi = AV[i];
    if(pi.distance(AP)<tol)
      found = true;
  }
  if(!found) {
    AV.push_back(AP);
    return true;
  }
  return false;
}
/*---------------------------------------------------------------------------*/
void Patching2D::buildLocalPatch(Face& AFace,
                                 const std::vector<math::Point>& ACorners)
{
    std::cout<<"====== LOCAL PATCH FOR "<<AFace.getID()<<" ========"<<std::endl;
    std::vector<math::Point> corners;
    if(ACorners.size()==4){
        corners = ACorners;
    }
    else {
        //We have more points, useless points are removed;
        for(unsigned int i=0; i<ACorners.size();i++){
            int j=(i+1)%ACorners.size();
            int k=(i==0)?ACorners.size()-1:i-1;
            
            math::Point pi = ACorners[i];
            math::Point pj = ACorners[j];
            math::Point pk = ACorners[k];
            math::Vector vij(pi,pj);
            math::Vector vik(pi,pk);
            vij.normalize();
            vik.normalize();
            if(fabs(vij.dot(vik))<0.4)
                corners.push_back(pi);
        }
    }
    
    if(corners.size()!=4){
        std::cout<<"Face "<<AFace.getID()<<std::endl;
        std::cout<<"Unable to get a 4-sided patch!!! "<<corners.size()<<std::endl;
        for(unsigned int i=0; i<ACorners.size();i++){
            int j=(i+1)%ACorners.size();
            int k=(i==0)?ACorners.size()-1:i-1;
            
            math::Point pi = ACorners[i];
            math::Point pj = ACorners[j];
            math::Point pk = ACorners[k];
            math::Vector vij(pi,pj);
            math::Vector vik(pi,pk);
            vij.normalize();
            vik.normalize();
            std::cout<<"\t"<< pi<<" --> "<<fabs(vij.dot(vik))<<std::endl;
        }
        std::cout<<std::endl;
        exit(0);
    }
    
    
    //Point are supposed to be ordered
    math::Point c0 = corners[0];
    math::Point c1 = corners[1];
    math::Point c2 = corners[2];
    math::Point c3 = corners[3];
    
//    m_patchs.addPatch(AFace.getID(),corners);
//    p.adjustSize(m_target_size);
  
 
    m_complex.write("output_add");
    
    
}

/*---------------------------------------------------------------------------*/
void Patching2D::initMarks(const int AMarkNodePnt,
                           const int AMarkNodeCrv,
                           const int AMarkEdgeCrv)
{
    m_mark_nodes_on_point = AMarkNodePnt;
    m_mark_nodes_on_curve = AMarkNodeCrv;
    m_mark_edges_on_curve = AMarkEdgeCrv;
}

/*----------------------------------------------------------------------------*/
void Patching2D::
writeOutput(const std::string& AFileName)
{
    static int out = 0;
    std::stringstream file_name;
    file_name << m_output_directory_name << "/" << AFileName << "_" << out;
    writeOutputSingle(file_name.str());
    out++;
}
/*----------------------------------------------------------------------------*/
void Patching2D::
writeOutputSingle(const std::string& AFileName)
{
    VTKWriter<IGMesh> writer(m_quad_mesh);
    writer.write(AFileName, DIM3 | F | N);
}
/*----------------------------------------------------------------------------*/
