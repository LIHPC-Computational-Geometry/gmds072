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
 * PatchBuilder.cpp
 *
 *  Created on: sept. 18 2015
 *      Author: Franck Ledoux
 */
/*---------------------------------------------------------------------------*/
// STL Header
#include <iostream>
#include <sstream>
/*---------------------------------------------------------------------------*/
// GMDS Header
#include <GMDS/Algo/BoundaryOperator.h>
#include <GMDS/Algo/DistanceFieldBuilder2D.h>
#include "GMDS/Math/Numerics.h"
#include <GMDS/IO/VTKWriter.h>
/*----------------------------------------------------------------------------*/
// FRAME Header
#include "Tools.h"
#include "PatchBuilder.h"
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
PatchBuilder::PatchBuilder(PatchComplex2D*                AComplex,
                           IGMesh*                        AMesh,
                           Variable<gmds::math::Cross2D>* AField):
m_complex(AComplex),
m_mesh(AMesh),
m_field(AField)
{}
/*---------------------------------------------------------------------------*/
void PatchBuilder::initMarks(const int AMarkNodePnt,
                             const int AMarkNodeCrv,
                             const int AMarkEdgeCrv)
{
    m_mark_nodes_on_point = AMarkNodePnt;
    m_mark_nodes_on_curve = AMarkNodeCrv;
    m_mark_edges_on_curve = AMarkEdgeCrv;
    //================================================================
    // We work on all faces
    //================================================================
    m_mark_face = m_mesh->getNewMark<Face>();
    IGMesh::face_iterator it_f = m_mesh->faces_begin();
    while(!it_f.isDone()){
        m_mesh->mark(it_f.value(),m_mark_face);
        it_f.next();
    }

    m_distance_field = 0;
}

/*----------------------------------------------------------------------------*/
void PatchBuilder::execute()
{
    DistanceFieldBuilder2D distanceBuilder(m_mesh);
    std::cout << "=======================================" << std::endl;
    std::cout << "Computation of a 2D distance field" << std::endl;
    if (!distanceBuilder.isValid())
    {
        std::cout << "Invalid model for distance computation" << std::endl;
        throw GMDSException("Invalid model for distance computation");
    }
    
    //================================================================
    //starting nodes are those on geometrical vertices
    //================================================================
    std::vector<Node> from_nodes, to_nodes;
    IGMesh::node_iterator it = m_mesh->nodes_begin();
    int k=0;
    while(!it.isDone()){
        Node current = it.value();
        if(m_mesh->isMarked(current, m_mark_nodes_on_point) && k<1){
            from_nodes.push_back(current);
            k++;

        }
        else
            to_nodes.push_back(current);
        it.next();
    }
    std::cout<<"From: "<<from_nodes.size()<<std::endl;
    std::cout<<"To:   "<<to_nodes.size()<<std::endl;
    
    m_distance_field =
    distanceBuilder.computeLInfDistanceFromCrossField(from_nodes,
                                                      to_nodes,
                                                      m_field,
                                                      m_mark_face);
    
    
    
    //==================================================================
    // CLEANING - Boolean marks are cleaned
    //==================================================================

    m_mesh->unmarkAll<Face>(m_mark_face);
    m_mesh->freeMark<Face>(m_mark_face);
}
/*----------------------------------------------------------------------------*/
void PatchBuilder::createPatch(Node& AOrigin, TCoord ADistance)
{
    if(m_mesh->isMarked(AOrigin, m_mark_nodes_on_point)){
        createPatchFromAGeometricPoint(AOrigin,ADistance);
    }
    else if (m_mesh->isMarked(AOrigin, m_mark_nodes_on_curve)){
        createPatchFromAGeometricCurve(AOrigin,ADistance);
    }
    else
        createPatchFromAnInnerPoint(AOrigin,ADistance);

}
/*----------------------------------------------------------------------------*/
void PatchBuilder::
createPatchFromAGeometricPoint(Node& AOrigin, TCoord ADistance)
{
    //====================================================================
    // We start from a node, so we are in a regular configuration with 4
    // component vectors
    // (except if we are along a singular triangle)
    //====================================================================
    math::Cross2D c = (*m_field)[AOrigin.getID()];
    std::vector<math::Vector> c_vectors = c.componentVectors();
    
}
/*----------------------------------------------------------------------------*/
void PatchBuilder::
createPatchFromAGeometricCurve(Node& AOrigin, TCoord ADistance)
{
    //====================================================================
    // We start from a node, so we are in a regular configuration with 4
    // component vectors
    // (except if we are along a singular triangle)
    //====================================================================
    math::Cross2D c = (*m_field)[AOrigin.getID()];
    std::vector<math::Vector> c_vectors = c.componentVectors();
    
}
/*----------------------------------------------------------------------------*/
void PatchBuilder::
createPatchFromAnInnerPoint(Node& AOrigin, TCoord ADistance)
{
    //=====================================================================
    // We start from a node, so we are in a regular configuration with 4
    // component vectors
    // (except if we are along a singular triangle)
    //=====================================================================
    math::Cross2D c = (*m_field)[AOrigin.getID()];
    std::vector<math::Vector> c_vectors = c.componentVectors();
    
    //=====================================================================
    // For each vector V, we compute the point at a distance of ADistance
    // from AOrigin in direction V.
    //
    // We can meet the mesh boundary before reaching the expected distance.
    // For each, we store the reached point, and the information about the
    // cell we arrive in.
    //=====================================================================
    std::vector<math::Point>  reached_points  (c_vectors.size());
    std::vector<math::Vector> reached_vectors (c_vectors.size());
    std::vector<int>          reached_cell_dim(c_vectors.size());
    std::vector<TCellID>      reached_cell_id (c_vectors.size());
    std::vector<int>          reached_boundary(c_vectors.size());

    for(unsigned int i=0; i<c_vectors.size(); i++){

        getPoint(AOrigin.getPoint(),
                 c_vectors[i],
                 0,
                 AOrigin.getID(),
                 ADistance,
                 reached_points[i],
                 reached_vectors[i],
                 reached_cell_dim[i],
                 reached_cell_id[i],
                 reached_boundary[i]);
    }
    //=====================================================================
    // For each reached point we compute tow half-line that will enclosed
    // the patch
    //=====================================================================
    for(unsigned int i=0;i<reached_points.size();i++){
        math::Point  pi = reached_points[i];
        math::Vector vi = reached_vectors[i];
        int cell_dim_i = reached_cell_dim[i];
        TCellID cell_id_i = reached_cell_id[i];
        if(cell_dim_i==0){
            Node cell_i;
        }
        else if(cell_dim_i==1){
            
        }
        else{//==2 necessary
            
        }
    }

}

/*----------------------------------------------------------------------------*/
void PatchBuilder::
getPoint(const math::Point&  AFromPnt,
         const math::Vector& AFromDir,
         const int&          AFromCellDim,
         const TCellID&      AFromCellID,
         const double        ADistanceMax,
         math::Point&        AToPnt,
         math::Vector&       AToDir,
         int&                AToCellDim,
         TCellID&            AToCellID,
         int&                AEndOnBnd)
{
    std::cout << "======== Streamline computation ========"
    << std::endl;
    
    Tools tool_object(m_mesh,m_field);
//    ATriangles.clear();
//    APoints.clear();
    
    math::Point  start_pnt = AFromPnt; //starting point
    math::Vector start_dir = AFromDir; //starting direction

    math::Vector  prev_dir = AFromDir; //starting point

    TCellID start_cell_id  = AFromCellID;
    int     start_cell_dim = AFromCellDim;
    
    math::Point  current_pnt = start_pnt;
    math::Vector current_vec = start_dir;
    
    std::cout << "START PNT: " << start_pnt << std::endl;
    math::Point start_dirPnt(start_pnt.X() + start_dir.X(),
                             start_pnt.Y() + start_dir.Y(),
                             start_pnt.Z() + start_dir.Z());
    
    std::cout << "START DIR PNT: " << start_dirPnt << std::endl;
    
    
    bool find_end = false;
    /* indicates that we reach a boundary point or line */
    bool end_on_boundary = false;
    /* indicates that we reach an existing singularity point*/
    //bool end_on_field_singularity = false;
    
    // We check some termination conditions on the boundary.
    if (start_cell_dim==0){
        Node current_node = m_mesh->get<Node>(start_cell_id);
        if(m_mesh->isMarked(current_node, m_mark_nodes_on_point) ||
           m_mesh->isMarked(current_node, m_mark_nodes_on_curve)){
            find_end        = true;
            end_on_boundary = true;
        }
    }
    else { //we have necessarry start_cell_dim=1
        Edge current_edge = m_mesh->get<Edge>(start_cell_id);
        if(m_mesh->isMarked(current_edge, m_mark_edges_on_curve)){
            find_end        = true;
            end_on_boundary = true;
        }
    }
    
    //========================================================================
    // Main loop to create the singularity line
    //========================================================================
    while(!find_end) {
        TCellID next_cell_id  = NullID;
        int     next_cell_dim = -1;
        
        tool_object.findNextCell(start_pnt, start_dir,
                                 start_cell_dim, start_cell_id,
                                 next_cell_dim, next_cell_id);
        

        if (next_cell_dim == -1){
            // The cell defined by (start_cell_dim, start_cell_id) is on the boundary.
            find_end = true;
            end_on_boundary = true;
        }
        else if (next_cell_dim == 1){
            // we are going along an edge.
            // Our simple assumption is to follow this edge until reaching
            // one of its end points and to compute the next direction at
            // this point.
            //
            // As we will arrive into a point, we cannot meet a singularity
            // point
            
            Edge current_edge = m_mesh->get<Edge>(next_cell_id);

//            std::vector<TCellID> adj_faces = current_edge.getIDs<Face>();
//            ATriangles.insert(ATriangles.end(),adj_faces.begin(),adj_faces.end());
            
            std::vector<Node> current_nodes = current_edge.get<Node>();
            //Do we go to the first node of current_edge or to the second?
            math::Vector v0(start_pnt, current_nodes[0].getPoint());
            math::Vector v1(start_pnt, current_nodes[1].getPoint());
            Node next_node;
            if(math::near(v0.norm(),0.0))
                next_node = current_nodes[1];
            else if(math::near(v1.norm(),0.0))
                next_node = current_nodes[0];
            else if(v0.dot(start_dir) > v1.dot(start_dir))
                next_node = current_nodes[0];
            else
                next_node = current_nodes[1];
            
            math::Vector next_dir;
            tool_object.computeOutVectorAtPoint(next_node, start_dir, next_dir);
            
            
            // We assign the new value for the next step
            start_dir = next_dir;
            start_pnt = next_node.getPoint();
            start_cell_dim = 0;
            start_cell_id = next_node.getID();
            find_end = false;
            
  //          APoints.push_back(start_pnt);
        }
        else { //general case, we are in a face
            Face current_face = m_mesh->get<Face>(next_cell_id);
  //          ATriangles.push_back(current_face.getID());
            //==============================================================
            // WE ARE NEVER IN A FACE CONTAINING A SING. POINT, SINCE THE
            // ALGORITHM STARTS FROM SING.
            //
            // CASE 2: SO, WE JUST HAVE TO CROSS THE TRIANGLE (GENERAL CONF)
            //==============================================================
            //Does the current triangle has the same classif
            math::Point  out_pnt;
            math::Vector out_vec;
            TCellID out_cell_id;
            int out_cell_dim;
            tool_object.traverseTriangle(current_face,/* the face we work on*/
                                         start_pnt,/* the point we start from */
                                         start_dir,/* the direction to follow*/
                                         start_cell_dim, /* the dimension of the cell start_pnt is located */
                                         start_cell_id,  /* the id of the cell start_pnt is located on*/
                                         out_pnt, /* the pnt where we go out */
                                         out_vec, /* the direction to follow after*/
                                         out_cell_dim,/* the dim. of the out cell*/
                                         out_cell_id);/* the id of the out cell*/
            
            //we keep the point toPnt to define the line
//            APoints.push_back(out_pnt);
            // we progress to the next point, next vector and so next face too
            prev_dir = start_dir; //we store the prev direction for slot
            //reconnection with balls
            start_pnt = out_pnt;
            start_dir = out_vec;
            start_cell_dim = out_cell_dim;
            start_cell_id = out_cell_id;
            
            
            //post process, we just look we are not arrived onto a geometric boundary
            if (start_cell_dim==0){
                Node current_node = m_mesh->get<Node>(start_cell_id);
                if(m_mesh->isMarked(current_node, m_mark_nodes_on_point) ||
                   m_mesh->isMarked(current_node, m_mark_nodes_on_curve)){
                    find_end        = true;
                    end_on_boundary = true;
                }
            }
            else { //we have necessarry start_cell_dim=1 
                Edge current_edge = m_mesh->get<Edge>(start_cell_id);
                if(m_mesh->isMarked(current_edge, m_mark_edges_on_curve)){
                    find_end        = true;
                    end_on_boundary = true;
                }
            }
        } // else { //general case, we are in a face
    } //while(!find_end)
    
    //==============================================================
    // Update of out parameters
    //==============================================================
    //last followed direction
    AToPnt = start_pnt;
    AToDir = start_dir;
    
    AEndOnBnd = end_on_boundary;
    
    //singularity point data if we found an end point
    AToCellDim = start_cell_dim;
    AToCellID = start_cell_id;
    
}