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
 * IGMeshDoctor.cpp
 *
 *  Created on: 22 mai 2014
 *      Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IGMeshDoctor.h>
#include <GMDS/Utils/CommonTypes.h>
/*----------------------------------------------------------------------------*/
#include <map>
/*----------------------------------------------------------------------------*/
namespace gmds {
    /*----------------------------------------------------------------------------*/
    IGMeshDoctor::IGMeshDoctor(IGMesh* AMesh)
    : m_mesh(AMesh)
    {}
    /*----------------------------------------------------------------------------*/
    IGMeshDoctor::~IGMeshDoctor()
    {}
    
    /*----------------------------------------------------------------------------*/
    void IGMeshDoctor::setMesh(IGMesh* AMesh)
    {
        m_mesh=AMesh;
    }
    /*----------------------------------------------------------------------------*/
    TCoord IGMeshDoctor::isLeft(Node& AN1, Node& AN2, Node& AN3)
    {
        return ( (AN2.X()-AN1.X()) * (AN3.Y()-AN1.Y()) -
                (AN3.X()-AN1.X()) * (AN2.Y()-AN1.Y()) );
    }
    /*----------------------------------------------------------------------------*/
    int IGMeshDoctor::orient2DFaces()
    {
        
        IGMesh::face_iterator it_f = m_mesh->faces_begin();
        TInt nb_reorientation =0;
        
        for(;!it_f.isDone();it_f.next())
        {
            Face f = it_f.value();
            if (orient2DFace(f))
                nb_reorientation++;
        }
        
        return nb_reorientation;
    }
    /*------------------------------------------------------------------------*/
    bool IGMeshDoctor::orient2DFace(Face& AF)
    {
        bool isReoriented = false;
        std::vector<Node> nodes = AF.get<Node>();
        TCoord orientation=0;
        if(AF.getType()==GMDS_TRIANGLE) {
            orientation = isLeft(nodes[0],nodes[1],nodes[2]);
        }
        else {
            //find the rightmost lowest vertex of the polygon
            int index_min=0;
            TCoord x_min = nodes[0].X();
            TCoord y_min = nodes[0].Y();
            for(unsigned int i=0;i<nodes.size();i++)
            {
                if(nodes[i].Y()>y_min)
                    continue;
                if(nodes[i].Y()==y_min) {	// just as low
                    if(nodes[i].X()<x_min)  // and to left
                        continue;
                }
                
                index_min =i;
                x_min = nodes[i].X();
                y_min = nodes[i].Y();
            }
            
            /* Orientation is tested in this vertex */
            if(index_min==0)
                orientation = isLeft(nodes[nodes.size()-1],nodes[0],nodes[1]);
            else if (index_min==nodes.size()-1)
                orientation = isLeft(nodes[index_min-1],nodes[index_min],nodes[0]);
            else
                orientation = isLeft(nodes[index_min-1],nodes[index_min],nodes[index_min+1]);
        }
        if(orientation>0.0) // clockwise or degenerated (=0)
        {
            isReoriented= true;
            std::vector<Node> nodes_inv;
            nodes_inv.resize(nodes.size());
            unsigned int node_size = nodes.size();
            for(unsigned int i=0;i<node_size;i++)
                nodes_inv[i] = nodes[node_size-1-i];
            AF.set<Node>(nodes_inv);
        }
        
        return isReoriented;
    }
    /*----------------------------------------------------------------------------*/
    void IGMeshDoctor::buildFacesAndR2F() const
    {
        if(m_mesh->getModel().has(F) && m_mesh->getModel().has(R2F)) {
            buildFAndR2F();
            //      buildR2F(m_mesh->getModel());
        } else {
            throw GMDSException("IGMeshDoctor::buildFacesAndR2F "
                                "Can not be called when mesh model does not have F && R2F");
        }
    }
    /*----------------------------------------------------------------------------*/
    void IGMeshDoctor::buildEdgesAndX2E() const
    {
        if(m_mesh->getModel().has(E) && (m_mesh->getModel().has(R2E) || m_mesh->getModel().has(F2E))) {
            buildE();
            if(m_mesh->getModel().has(R2E)) {
                buildR2E(m_mesh->getModel());
            }
            if(m_mesh->getModel().has(F2E)) {
                buildF2E(m_mesh->getModel());
            }
        } else {
            throw GMDSException("IGMeshDoctor::buildFacesAndR2E "
                                "Can not be called when mesh model does not have E && (R2E || F2E)");
        }
    }
    /*----------------------------------------------------------------------------*/
    void IGMeshDoctor::updateUpwardConnectivity() const
    {
        MeshModel mod = m_mesh->getModel();
        if(mod.has(N2F) && mod.has(F2N))
        {
            IGMesh::face_iterator it_f = m_mesh->faces_begin();
            for(;!it_f.isDone();it_f.next())
            {
                Face f = it_f.value();
                std::vector<Node> nodes = f.get<Node>();
                
                for(unsigned int i=0; i<nodes.size();i++)
                    if(!nodes[i].has<Face>(f)) {
                        nodes[i].add<Face>(f);
                    }
            }
            
        }
        if(mod.has(E2F) && mod.has(F2E))
        {
            IGMesh::face_iterator it_f = m_mesh->faces_begin();
            for(;!it_f.isDone();it_f.next())
            {
                Face f = it_f.value();
                std::vector<Edge> edges = f.get<Edge>();
                
                for(unsigned int i=0; i<edges.size();i++)
                    if(!edges[i].has<Face>(f)) {
                        edges[i].add<Face>(f);
                    }
            }
        }
        if(mod.has(N2R) && mod.has(R2N))
        {
            IGMesh::region_iterator it_r = m_mesh->regions_begin();
            for(;!it_r.isDone();it_r.next())
            {
                Region r = it_r.value();
                std::vector<Node> nodes = r.get<Node>();
                
                for(unsigned int i=0; i<nodes.size();i++)
                    if(!nodes[i].has<Region>(r)) {
                        nodes[i].add<Region>(r);
                    }
            }
        }
        
        if(mod.has(E2R) && mod.has(R2E))
        {
            IGMesh::region_iterator it_r = m_mesh->regions_begin();
            for(;!it_r.isDone();it_r.next())
            {
                Region r = it_r.value();
                std::vector<Edge> edges = r.get<Edge>();
                
                for(unsigned int i=0; i<edges.size();i++)
                    if(!edges[i].has<Region>(r)) {
                        edges[i].add<Region>(r);
                    }
            }
        }
        if(mod.has(N2E) && mod.has(E2N))
        {
            IGMesh::edge_iterator it_e = m_mesh->edges_begin();
            for(;!it_e.isDone();it_e.next())
            {
                Edge e = it_e.value();
                std::vector<Node> nodes = e.get<Node>();
                
                for(unsigned int i=0; i<nodes.size();i++)
                    if(!nodes[i].has<Edge>(e)) {
                        nodes[i].add<Edge>(e);
                    }
            }
        }
        
        if(mod.has(F2R) && mod.has(R2F))
        {
            IGMesh::region_iterator it_r = m_mesh->regions_begin();
            for(;!it_r.isDone();it_r.next())
            {
                Region r = it_r.value();
                std::vector<Face> faces = r.get<Face>();
                
                for(unsigned int i=0; i<faces.size();i++) {
                    if(!faces[i].has<Region>(r.getID())) {
                        faces[i].add<Region>(r);
                    }
                }
            }
        }
    }
    /*----------------------------------------------------------------------------*/
    void IGMeshDoctor::buildE() const
    {
        MeshModel mod = m_mesh->getModel();
        
        std::map<FakeEdge::EdgeID, TCellID> tmp_edges;
        
        // first register the existing edges
        if(mod.has(E) && mod.has(E2N)) {
            IGMesh::edge_iterator ite = m_mesh->edges_begin();
            for(;!ite.isDone();ite.next()) {
                std::vector<TCellID> e_nodes = ite.value().getIDs<Node>();
                FakeEdge fke(e_nodes[0],e_nodes[1]);
                tmp_edges[fke.getID()] = ite.value().getID();
            }
        } else {
            throw GMDSException("IGMeshDoctor::buildE E or E2N missing.");
        }
        
        if(mod.has(R) && mod.has(R2N)){
            IGMesh::region_iterator it = m_mesh->regions_begin();
            
            for(;!it.isDone();it.next()){
                Region r = it.value();
                ECellType r_type = r.getType();
                std::vector<TCellID> r_nodes = r.getIDs<Node>();
                if(r_type==GMDS_TETRA){
                    //EDGE 1
                    addEdge(r_nodes[0],r_nodes[1],tmp_edges);
                    //EDGE 2
                    addEdge(r_nodes[0],r_nodes[2],tmp_edges);
                    //EDGE 3
                    addEdge(r_nodes[0],r_nodes[3],tmp_edges);
                    //EDGE 4
                    addEdge(r_nodes[1],r_nodes[2],tmp_edges);
                    //EDGE 5
                    addEdge(r_nodes[1],r_nodes[3],tmp_edges);
                    //EDGE 6
                    addEdge(r_nodes[2],r_nodes[3],tmp_edges);
                }
                else if(r_type==GMDS_HEX){
                    addEdge(r_nodes[0],r_nodes[1],tmp_edges);
                    addEdge(r_nodes[1],r_nodes[2],tmp_edges);
                    addEdge(r_nodes[2],r_nodes[3],tmp_edges);
                    addEdge(r_nodes[3],r_nodes[0],tmp_edges);
                    addEdge(r_nodes[4],r_nodes[5],tmp_edges);
                    addEdge(r_nodes[5],r_nodes[6],tmp_edges);
                    addEdge(r_nodes[6],r_nodes[7],tmp_edges);
                    addEdge(r_nodes[7],r_nodes[4],tmp_edges);
                    addEdge(r_nodes[0],r_nodes[4],tmp_edges);
                    addEdge(r_nodes[1],r_nodes[5],tmp_edges);
                    addEdge(r_nodes[2],r_nodes[6],tmp_edges);
                    addEdge(r_nodes[3],r_nodes[7],tmp_edges);
                }
                else if(r_type==GMDS_PYRAMID){
                    addEdge(r_nodes[0],r_nodes[1],tmp_edges);
                    addEdge(r_nodes[1],r_nodes[2],tmp_edges);
                    addEdge(r_nodes[2],r_nodes[3],tmp_edges);
                    addEdge(r_nodes[3],r_nodes[0],tmp_edges);
                    addEdge(r_nodes[0],r_nodes[4],tmp_edges);
                    addEdge(r_nodes[1],r_nodes[4],tmp_edges);
                    addEdge(r_nodes[2],r_nodes[4],tmp_edges);
                    addEdge(r_nodes[3],r_nodes[4],tmp_edges);
                }
                else if(r_type==GMDS_PRISM3){
                    addEdge(r_nodes[0],r_nodes[1],tmp_edges);
                    addEdge(r_nodes[1],r_nodes[2],tmp_edges);
                    addEdge(r_nodes[2],r_nodes[0],tmp_edges);
                    addEdge(r_nodes[3],r_nodes[4],tmp_edges);
                    addEdge(r_nodes[4],r_nodes[5],tmp_edges);
                    addEdge(r_nodes[5],r_nodes[3],tmp_edges);
                    addEdge(r_nodes[0],r_nodes[3],tmp_edges);
                    addEdge(r_nodes[1],r_nodes[4],tmp_edges);
                    addEdge(r_nodes[2],r_nodes[5],tmp_edges);
                }
                else
                    throw GMDSException("IGMeshDoctor::buildE Cell type unknown.");
            }
        }
        else if(mod.has(F) && mod.has(F2N)){
            IGMesh::face_iterator it = m_mesh->faces_begin();
            
            for(;!it.isDone();it.next()){
                Face f = it.value();
                ECellType f_type = f.getType();
                std::vector<TCellID> f_nodes = f.getIDs<Node>();
                if(f_type==GMDS_TRIANGLE){
                    //EDGE 1
                    addEdge(f_nodes[0],f_nodes[1],tmp_edges);
                    //EDGE 2
                    addEdge(f_nodes[1],f_nodes[2],tmp_edges);
                    //EDGE 3
                    addEdge(f_nodes[2],f_nodes[0],tmp_edges);
                }
                else if(f_type==GMDS_QUAD){
                    //EDGE 1
                    addEdge(f_nodes[0],f_nodes[1],tmp_edges);
                    //EDGE 2
                    addEdge(f_nodes[1],f_nodes[2],tmp_edges);
                    //EDGE 3
                    addEdge(f_nodes[2],f_nodes[3],tmp_edges);
                    //EDGE 4
                    addEdge(f_nodes[3],f_nodes[0],tmp_edges);
                }
                else if(f_type==GMDS_POLYGON){
                    unsigned int nb_nodes = f_nodes.size();
                    for(unsigned int i_node=0;i_node<nb_nodes;i_node++){
                        addEdge(f_nodes[i_node],f_nodes[(i_node+1)%nb_nodes],tmp_edges);
                    }
                }
                else
                    throw GMDSException("IGMeshDoctor::buildE Cell type unknown.");
            }
        }
        else
            throw GMDSException("IGMeshDoctor::buildE Not yet implemented");
    }
    /*----------------------------------------------------------------------------*/
    void IGMeshDoctor::addEdge(TCellID AN1, TCellID AN2,
                               std::map<FakeEdge::EdgeID, TCellID>& AFakeEdgeMap) const
    {
        FakeEdge fke(AN1,AN2);
        if(AFakeEdgeMap.find(fke.getID())==AFakeEdgeMap.end())
        {
            //New face to add
            Edge e = m_mesh->newEdge(AN1,AN2);
            AFakeEdgeMap[fke.getID()]=e.getID();
        }
    }
    /*----------------------------------------------------------------------------*/
    TCellID IGMeshDoctor::
    addFace(std::vector<TCellID>& ANodeIDs,
            std::map<FakeFace::FaceID, TCellID>& AFakeFaceMap) const
    {
        TCellID f_id;

        FakeFace fkf(ANodeIDs);
        auto itf =AFakeFaceMap.find(fkf.getID());

        if(itf==AFakeFaceMap.end()){
            //New face to add
            Face f = m_mesh->newFace(ANodeIDs);
            AFakeFaceMap[fkf.getID()]=f.getID();
            f_id=f.getID();
        }
        else{
            f_id= itf->second;
        }
        return f_id;
    }
    /*----------------------------------------------------------------------------*/
    TCellID IGMeshDoctor::
    addFace(Face& AFace,
            std::map<FakeFace::FaceID, TCellID>& AFakeFaceMap) const
    {
        TCellID f_id;
        std::vector<TCellID> f_nodes = AFace.getIDs<Node>();
        
        FakeFace fkf(f_nodes);
        auto itf =AFakeFaceMap.find(fkf.getID());
        if(itf==AFakeFaceMap.end())
        {
            //New face to add
            AFakeFaceMap[fkf.getID()]=AFace.getID();
            f_id = AFace.getID();
        }
        else{
            f_id = itf->second;
        }
        return f_id;
    }
    /*----------------------------------------------------------------------------*/
    void IGMeshDoctor::buildFAndR2F() const
    {
        MeshModel mod = m_mesh->getModel();
        
        std::map<FakeFace::FaceID, TCellID> tmp_faces;
        //==============================================================
        // First we put existing faces into our tmp_faces container
        //==============================================================
        if(mod.has(F) && mod.has(F2N)){
            IGMesh::face_iterator it = m_mesh->faces_begin();
            
            for(;!it.isDone();it.next()){
                Face f = it.value();
                std::vector<TCellID> f_nodes = f.getIDs<Node>();
                addFace(f,tmp_faces);
            }
        }
        //==============================================================
        // Second we look for missing faces
        //==============================================================
        if(mod.has(R) && mod.has(R2N)){
            IGMesh::region_iterator it = m_mesh->regions_begin();
            
            for(;!it.isDone();it.next()){
                Region r = it.value();
                ECellType r_type = r.getType();
                std::vector<TCellID> r_nodes = r.getIDs<Node>();
                std::vector<TCellID> r_faces;
                if(r_type==GMDS_TETRA){
                    r_faces.resize(4);
                    std::vector<TCellID> f_nodes;
                    f_nodes.resize(3);
                    //FACE 1
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[2];
                    f_nodes[2] = r_nodes[1];
                    r_faces[0] = addFace(f_nodes,tmp_faces);
                    
                    //FACE 2
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[1];
                    f_nodes[2] = r_nodes[3];
                    r_faces[1] = addFace(f_nodes,tmp_faces);
                    
                    //FACE 3
                    f_nodes[0] = r_nodes[1];
                    f_nodes[1] = r_nodes[2];
                    f_nodes[2] = r_nodes[3];
                    r_faces[2] = addFace(f_nodes,tmp_faces);
                    
                    //FACE 4
                    f_nodes[0] = r_nodes[2];
                    f_nodes[1] = r_nodes[0];
                    f_nodes[2] = r_nodes[3];
                    r_faces[3] = addFace(f_nodes,tmp_faces);
                }
                else if(r_type==GMDS_HEX){
                    std::vector<TCellID> f_nodes;
                    f_nodes.resize(4);
                    r_faces.resize(6);
                    //FACE 1
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[3];
                    f_nodes[2] = r_nodes[2];
                    f_nodes[3] = r_nodes[1];
                    r_faces[0] = addFace(f_nodes,tmp_faces);
                    //FACE 2
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[1];
                    f_nodes[2] = r_nodes[5];
                    f_nodes[3] = r_nodes[4];
                    r_faces[1] = addFace(f_nodes,tmp_faces);
                    //FACE 3
                    f_nodes[0] = r_nodes[1];
                    f_nodes[1] = r_nodes[2];
                    f_nodes[2] = r_nodes[6];
                    f_nodes[3] = r_nodes[5];
                    r_faces[2] = addFace(f_nodes,tmp_faces);
                    //FACE 4
                    f_nodes[0] = r_nodes[2];
                    f_nodes[1] = r_nodes[3];
                    f_nodes[2] = r_nodes[7];
                    f_nodes[3] = r_nodes[6];
                    r_faces[3] = addFace(f_nodes,tmp_faces);
                    //FACE 5
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[4];
                    f_nodes[2] = r_nodes[7];
                    f_nodes[3] = r_nodes[3];
                    r_faces[4] = addFace(f_nodes,tmp_faces);
                    //FACE 6
                    f_nodes[0] = r_nodes[4];
                    f_nodes[1] = r_nodes[5];
                    f_nodes[2] = r_nodes[6];
                    f_nodes[3] = r_nodes[7];
                    r_faces[5] = addFace(f_nodes,tmp_faces);
                }
                else if(r_type==GMDS_PRISM3){
                    std::vector<TCellID> f_nodes;
                    f_nodes.resize(4);
                    r_faces.resize(5);
                    //FACE 1
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[1];
                    f_nodes[2] = r_nodes[4];
                    f_nodes[3] = r_nodes[3];
                    r_faces[0] = addFace(f_nodes,tmp_faces);
                    //FACE 2
                    f_nodes[0] = r_nodes[1];
                    f_nodes[1] = r_nodes[2];
                    f_nodes[2] = r_nodes[5];
                    f_nodes[3] = r_nodes[4];
                    r_faces[1] = addFace(f_nodes,tmp_faces);
                    //FACE 3
                    f_nodes[0] = r_nodes[2];
                    f_nodes[1] = r_nodes[0];
                    f_nodes[2] = r_nodes[3];
                    f_nodes[3] = r_nodes[5];
                    r_faces[2] = addFace(f_nodes,tmp_faces);
                    //FACE 4
                    f_nodes.resize(3);
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[2];
                    f_nodes[2] = r_nodes[1];
                    r_faces[3] = addFace(f_nodes,tmp_faces);
                    //FACE 5
                    f_nodes[0] = r_nodes[3];
                    f_nodes[1] = r_nodes[4];
                    f_nodes[2] = r_nodes[5];
                    r_faces[4] = addFace(f_nodes,tmp_faces);
                }
                else if(r_type==GMDS_PYRAMID){
                    std::vector<TCellID> f_nodes;
                    f_nodes.resize(4);
                    r_faces.resize(5);
                    //FACE 1
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[3];
                    f_nodes[2] = r_nodes[2];
                    f_nodes[3] = r_nodes[1];
                    r_faces[0] = addFace(f_nodes,tmp_faces);
                    //FACE 2
                    f_nodes.resize(3);
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[1];
                    f_nodes[2] = r_nodes[4];
                    r_faces[1] = addFace(f_nodes,tmp_faces);
                    //FACE 3
                    f_nodes[0] = r_nodes[1];
                    f_nodes[1] = r_nodes[2];
                    f_nodes[2] = r_nodes[4];
                    r_faces[2] = addFace(f_nodes,tmp_faces);
                    //FACE 4
                    f_nodes[0] = r_nodes[2];
                    f_nodes[1] = r_nodes[3];
                    f_nodes[2] = r_nodes[4];
                    r_faces[3] = addFace(f_nodes,tmp_faces);
                    //FACE 5
                    f_nodes[0] = r_nodes[3];
                    f_nodes[1] = r_nodes[0];
                    f_nodes[2] = r_nodes[4];
                    r_faces[4] = addFace(f_nodes,tmp_faces);
                }
                else
                    throw GMDSException("IGMeshDoctor::buildF Not yet implemented");
                r.set<Face>(r_faces);
            }//for(;!it.isDone();it.next())
        }//if(mod.has(R) && mod.has(R2N))
        else
            throw GMDSException("IGMeshDoctor::buildF Not yet implemented");
    }    /*----------------------------------------------------------------------------*/
    void IGMeshDoctor::buildF() const
    {
        MeshModel mod = m_mesh->getModel();
        
        std::map<FakeFace::FaceID, TCellID> tmp_faces;
        //==============================================================
        // First we put existing faces into our tmp_faces container
        //==============================================================
        if(mod.has(F) && mod.has(F2N)){
            IGMesh::face_iterator it = m_mesh->faces_begin();
            
            for(;!it.isDone();it.next()){
                Face f = it.value();
                std::vector<TCellID> f_nodes = f.getIDs<Node>();
                addFace(f,tmp_faces);
            }
        }
        //==============================================================
        // Second we look for missing faces
        //==============================================================
        if(mod.has(R) && mod.has(R2N)){
            IGMesh::region_iterator it = m_mesh->regions_begin();
            
            for(;!it.isDone();it.next()){
                Region r = it.value();
                ECellType r_type = r.getType();
                std::vector<TCellID> r_nodes = r.getIDs<Node>();
                if(r_type==GMDS_TETRA){
                    std::vector<TCellID> f_nodes;
                    f_nodes.resize(3);
                    //FACE 1
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[2];
                    f_nodes[2] = r_nodes[1];
                    addFace(f_nodes,tmp_faces);
                    
                    //FACE 2
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[1];
                    f_nodes[2] = r_nodes[3];
                    addFace(f_nodes,tmp_faces);
                    
                    //FACE 3
                    f_nodes[0] = r_nodes[1];
                    f_nodes[1] = r_nodes[2];
                    f_nodes[2] = r_nodes[3];
                    addFace(f_nodes,tmp_faces);
                    
                    //FACE 4
                    f_nodes[0] = r_nodes[2];
                    f_nodes[1] = r_nodes[0];
                    f_nodes[2] = r_nodes[3];
                    addFace(f_nodes,tmp_faces);
                }
                else if(r_type==GMDS_HEX){
                    std::vector<TCellID> f_nodes;
                    f_nodes.resize(4);
                    
                    //FACE 1
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[3];
                    f_nodes[2] = r_nodes[2];
                    f_nodes[3] = r_nodes[1];
                    addFace(f_nodes,tmp_faces);
                    //FACE 2
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[1];
                    f_nodes[2] = r_nodes[5];
                    f_nodes[3] = r_nodes[4];
                    addFace(f_nodes,tmp_faces);
                    //FACE 3
                    f_nodes[0] = r_nodes[1];
                    f_nodes[1] = r_nodes[2];
                    f_nodes[2] = r_nodes[6];
                    f_nodes[3] = r_nodes[5];
                    addFace(f_nodes,tmp_faces);
                    //FACE 4
                    f_nodes[0] = r_nodes[2];
                    f_nodes[1] = r_nodes[3];
                    f_nodes[2] = r_nodes[7];
                    f_nodes[3] = r_nodes[6];
                    addFace(f_nodes,tmp_faces);
                    //FACE 5
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[4];
                    f_nodes[2] = r_nodes[7];
                    f_nodes[3] = r_nodes[3];
                    addFace(f_nodes,tmp_faces);
                    //FACE 6
                    f_nodes[0] = r_nodes[4];
                    f_nodes[1] = r_nodes[5];
                    f_nodes[2] = r_nodes[6];
                    f_nodes[3] = r_nodes[7];
                    addFace(f_nodes,tmp_faces);
                }
                else if(r_type==GMDS_PRISM3){
                    std::vector<TCellID> f_nodes;
                    f_nodes.resize(4);
                    
                    //FACE 1
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[1];
                    f_nodes[2] = r_nodes[4];
                    f_nodes[3] = r_nodes[3];
                    addFace(f_nodes,tmp_faces);
                    //FACE 2
                    f_nodes[0] = r_nodes[1];
                    f_nodes[1] = r_nodes[2];
                    f_nodes[2] = r_nodes[5];
                    f_nodes[3] = r_nodes[4];
                    addFace(f_nodes,tmp_faces);
                    //FACE 3
                    f_nodes[0] = r_nodes[2];
                    f_nodes[1] = r_nodes[0];
                    f_nodes[2] = r_nodes[3];
                    f_nodes[3] = r_nodes[5];
                    addFace(f_nodes,tmp_faces);
                    //FACE 4
                    f_nodes.resize(3);
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[2];
                    f_nodes[2] = r_nodes[1];
                    addFace(f_nodes,tmp_faces);
                    //FACE 5
                    f_nodes[0] = r_nodes[3];
                    f_nodes[1] = r_nodes[4];
                    f_nodes[2] = r_nodes[5];
                    addFace(f_nodes,tmp_faces);
                }
                else if(r_type==GMDS_PYRAMID){
                    std::vector<TCellID> f_nodes;
                    f_nodes.resize(4);
                    //FACE 1
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[3];
                    f_nodes[2] = r_nodes[2];
                    f_nodes[3] = r_nodes[1];
                    addFace(f_nodes,tmp_faces);
                    //FACE 2
                    f_nodes.resize(3);
                    f_nodes[0] = r_nodes[0];
                    f_nodes[1] = r_nodes[1];
                    f_nodes[2] = r_nodes[4];
                    addFace(f_nodes,tmp_faces);
                    //FACE 3
                    f_nodes[0] = r_nodes[1];
                    f_nodes[1] = r_nodes[2];
                    f_nodes[2] = r_nodes[4];
                    addFace(f_nodes,tmp_faces);
                    //FACE 4
                    f_nodes[0] = r_nodes[2];
                    f_nodes[1] = r_nodes[3];
                    f_nodes[2] = r_nodes[4];
                    addFace(f_nodes,tmp_faces);
                    //FACE 5
                    f_nodes[0] = r_nodes[3];
                    f_nodes[1] = r_nodes[0];
                    f_nodes[2] = r_nodes[4];
                    addFace(f_nodes,tmp_faces);
                }
                else
                    throw GMDSException("IGMeshDoctor::buildF Not yet implemented");
            }//for(;!it.isDone();it.next())
        }//if(mod.has(R) && mod.has(R2N))
        else
            throw GMDSException("IGMeshDoctor::buildF Not yet implemented");
    }
    /*----------------------------------------------------------------------------*/
    void IGMeshDoctor::buildN2N(const MeshModel &ARefModel) const
    {
        throw GMDSException("IGMeshDoctor::buildN2N Not yet implemented");
    }
    /*----------------------------------------------------------------------------*/
    void IGMeshDoctor::buildN2E(const MeshModel &ARefModel) const
    {
        if(ARefModel.has(E2N)){
            IGMesh::edge_iterator it = m_mesh->edges_begin();
            for(;!it.isDone();it.next()){
                Edge e = it.value();
                std::vector<Node> e_nodes;
                e.get<Node>(e_nodes);
                for(unsigned int i=0;i<e_nodes.size();i++)
                {
                    Node n_i = e_nodes[i];
                    n_i.add<Edge>(e);
                }
            }
        }
        else
            throw GMDSException("IGMeshDoctor::buildN2E Not yet implemented");}
    /*----------------------------------------------------------------------------*/
    void IGMeshDoctor::buildN2F(const MeshModel &ARefModel) const
    {
        if(ARefModel.has(F2N)){
            IGMesh::face_iterator it = m_mesh->faces_begin();
            for(;!it.isDone();it.next()){
                Face f = it.value();
                std::vector<Node> f_nodes;
                f.get<Node>(f_nodes);
                for(unsigned int i=0;i<f_nodes.size();i++)
                {
                    Node n_i = f_nodes[i];
                    n_i.add<Face>(f);
                }
            }
        }
        else
            throw GMDSException("IGMeshDoctor::buildN2F Not yet implemented");
    }
    /*----------------------------------------------------------------------------*/
    void IGMeshDoctor::buildN2R(const MeshModel &ARefModel) const
    {
        if(ARefModel.has(R2N)){
            IGMesh::region_iterator it = m_mesh->regions_begin();
            for(;!it.isDone();it.next()){
                Region r = it.value();
                std::vector<Node> r_nodes;
                r.get<Node>(r_nodes);
                for(unsigned int i=0;i<r_nodes.size();i++)
                {
                    Node n_i = r_nodes[i];
                    n_i.add<Region>(r);
                }
            }
        }
        else
            throw GMDSException("IGMeshDoctor::buildN2R Not yet implemented");
    }
    /*----------------------------------------------------------------------------*/
    void IGMeshDoctor::buildE2N(const MeshModel &ARefModel) const
    {
        if(ARefModel.has(F) && ARefModel.has(F2N)){
            IGMesh::face_iterator it = m_mesh->faces_begin();
            std::map<FakeEdge::EdgeID,TCellID> fakeEdgeMap;
            for(;!it.isDone();it.next()){
                Face f = it.value();
                std::cout<<"FACE "<<f.getID()<<std::endl;
                std::vector<TCellID> f_node_ids;
                f.getIDs<Node>(f_node_ids);
                for(unsigned int i=0;i<f_node_ids.size();i++)
                {
                    TCellID  id1 = f_node_ids[i];
                    TCellID  id2 = f_node_ids[(i+1)%(f_node_ids.size())];
                    addEdge(id1,id2,fakeEdgeMap);
                }
            }
        }
        else
            throw GMDSException("It is impossible to retrieve the E2N adjacency in a generic manner");
    }
    /*----------------------------------------------------------------------------*/
    void IGMeshDoctor::buildE2E(const MeshModel &ARefModel) const
    {
        if(ARefModel.has(E2N) && ARefModel.has(N2E)){
            IGMesh::edge_iterator it = m_mesh->edges_begin();
            for(;!it.isDone();it.next()){
                Edge e = it.value();
                TCellID e_id = e.getID();
                std::vector<Node> e_nodes;
                e.get<Node>(e_nodes);
                for(unsigned int i=0;i<e_nodes.size();i++)
                {
                    Node n_i = e_nodes[i];
                    std::vector<Edge> adj_edges;
                    n_i.get<Edge>(adj_edges);
                    for(unsigned int j=0;j<adj_edges.size();j++){
                        Edge adj_edge = adj_edges[j];
                        if(adj_edge.getID()!=e_id){
                            e.add<Edge>(adj_edge);
                        }
                    }
                }
            }
        }
        else if (ARefModel.has(E2N)){
            IGMesh::edge_iterator it = m_mesh->edges_begin();
            //we build a tmp relations from nodes to edges
            std::map<TCellID,std::vector<TCellID> > tmp_N2E;
            for(;!it.isDone();it.next()){
                Edge e = it.value();
                TCellID e_id = e.getID();
                std::vector<TCellID> node_ids;
                e.getIDs<Node>(node_ids);
                for(unsigned int i=0;i<node_ids.size();i++)
                    tmp_N2E[node_ids[i]].push_back(e_id);
            }
            //we traverse again the edges to update E2E
            for(it = m_mesh->edges_begin();!it.isDone();it.next()){
                Edge e = it.value();
                TCellID e_id = e.getID();
                std::vector<TCellID> node_ids;
                e.getIDs<Node>(node_ids);
                for(unsigned int i=0;i<node_ids.size();i++)
                    tmp_N2E[node_ids[i]].push_back(e_id);
                if(node_ids.size()!=2)
                    throw GMDSException("Error: An edge does not have two end nodes in the E2E building process");
                TCellID n1 = node_ids[0];
                TCellID n2 = node_ids[1];
                
                std::vector<TCellID> edges1 = tmp_N2E[n1];
                std::vector<TCellID> edges2 = tmp_N2E[n2];
                std::vector<TCellID> common_edges = getCommonBut(edges1,edges2,e_id);
                
                for(unsigned int i_common = 0; i_common<common_edges.size(); i_common++){
                    e.add<Edge>(common_edges[i_common]);
                }
            }
        }
        else
            throw GMDSException("Impossible to retrieve the E2E adjacency with this ARefModelel");
    }
    /*----------------------------------------------------------------------------*/
    void IGMeshDoctor::buildE2F(const MeshModel &ARefModel) const
    {
        if(ARefModel.has(F2E)){
            IGMesh::face_iterator it = m_mesh->faces_begin();
            for(;!it.isDone();it.next()){
                Face f = it.value();
                std::vector<Edge> f_edges;
                f.get<Edge>(f_edges);
                for(unsigned int i=0;i<f_edges.size();i++)
                {
                    Edge e_i = f_edges[i];
                    e_i.add<Face>(f);
                }
            }
        }
        else
            throw GMDSException("IGMeshDoctor::buildE2F Not yet implemented");
    }
    /*----------------------------------------------------------------------------*/
    void IGMeshDoctor::buildE2R(const MeshModel &ARefModel) const
    {
        if(ARefModel.has(R2E)){
            IGMesh::region_iterator it = m_mesh->regions_begin();
            for(;!it.isDone();it.next()){
                Region r = it.value();
                std::vector<Edge> r_edges;
                r.get<Edge>(r_edges);
                for(unsigned int i=0;i<r_edges.size();i++)
                {
                    Edge e_i = r_edges[i];
                    e_i.add<Region>(r);
                }
            }
        }
        else
            throw GMDSException("IGMeshDoctor::buildE2R Not yet implemented");
    }
    /*----------------------------------------------------------------------------*/
    void IGMeshDoctor::buildF2N(const MeshModel &ARefModel) const
    {
        throw GMDSException("IGMeshDoctor::buildF2N Not yet implemented");
    }
    /*----------------------------------------------------------------------------*/
    void IGMeshDoctor::buildF2E(const MeshModel &ARefModel) const
    {
        //we just have F2N and E2N
        std::map<TCellID,std::vector<TCellID> > tmp_N2E;
        IGMesh::edge_iterator it_e = m_mesh->edges_begin();
        for(;!it_e.isDone();it_e.next()){
            Edge e= it_e.value();
            TCellID e_id = e.getID();
            std::vector<TCellID> e_node_ids;
            e.getIDs<Node>(e_node_ids);
            for(unsigned int i=0;i<e_node_ids.size();i++)
            {
                tmp_N2E[e_node_ids[i]].push_back(e_id);
            }
        }
        
        IGMesh::face_iterator it_f = m_mesh->faces_begin();
        for(;!it_f.isDone();it_f.next()){
            Face f = it_f.value();
            std::vector<TCellID> f_node_ids;
            f.getIDs<Node>(f_node_ids);
            std::vector<TCellID> multiset_edge;
            for(unsigned int i=0;i<f_node_ids.size();i++)
            {
                std::vector<TCellID> current_edges= tmp_N2E[f_node_ids[i]];
                multiset_edge.insert(multiset_edge.end(),current_edges.begin(),current_edges.end());
            }
            std::vector<TCellID> result = keepFilter(multiset_edge,2);
            f.set<Edge>(result);
        }
    }
    /*----------------------------------------------------------------------------*/
    void IGMeshDoctor::buildF2F(const MeshModel &ARefModel) const
    {
        if(ARefModel.has(F2E) && ARefModel.has(E2F)){
            IGMesh::face_iterator it = m_mesh->faces_begin();
            for(;!it.isDone();it.next()){
                Face f = it.value();
                TCellID f_id = f.getID();
                std::vector<Edge> f_edges;
                f.get<Edge>(f_edges);
                for(unsigned int i=0;i<f_edges.size();i++)
                {
                    Edge e_i = f_edges[i];
                    std::vector<Face> adj_faces;
                    e_i.get<Face>(adj_faces);
                    for(unsigned int j=0;j<adj_faces.size();j++){
                        Face adj_face = adj_faces[j];
                        if(adj_face.getID()!=f_id){
                            f.add<Face>(adj_face);
                        }
                    }
                }
            }
        }
        else if (ARefModel.has(F2N) && ARefModel.has(N2F)){
            IGMesh::face_iterator it = m_mesh->faces_begin();
            for(it=m_mesh->faces_begin();!it.isDone();it.next())
            {
                Face    f    = it.value();
                TCellID f_id = f.getID();
                std::vector<Node> f_nodes;
                f.get<Node>(f_nodes);
                //now we get the set of faces adjacent to the edges of f
                for(unsigned int i=0;i<f_nodes.size();i++)
                {
                    Node n1 = f_nodes[i];
                    Node n2 = f_nodes[(i+1)%f_nodes.size()];
                    std::vector<TCellID> faces1 = n1.getIDs<Face>();
                    std::vector<TCellID> faces2 = n2.getIDs<Face>();
                    std::vector<TCellID> common_faces = getCommonBut(faces1,faces2,f_id);
                    
                    for(unsigned int i_common = 0; i_common<common_faces.size(); i_common++){
                        f.add<Face>(common_faces[i_common]);
                    }
                }
            }
        }
        else if (ARefModel.has(F2E)){
            IGMesh::face_iterator it = m_mesh->faces_begin();
            // map an edge id to the incident face ids.
            std::map<TCellID,std::vector<TCellID> > tmpE2F;
            for(;!it.isDone();it.next()){
                Face    f    = it.value();
                TCellID f_id = f.getID();
                std::vector<TCellID> f_edges;
                f.getIDs<Edge>(f_edges);
                for(unsigned int i=0;i<f_edges.size();i++)
                {
                    tmpE2F[f_edges[i]].push_back(f_id);
                }
            }
            //Now we connect each face to adjacent ones using tmpE2F
            for(it=m_mesh->faces_begin();!it.isDone();it.next())
            {
                Face    f    = it.value();
                TCellID f_id = f.getID();
                std::vector<TCellID> f_edges;
                f.getIDs<Edge>(f_edges);
                //now we get the set of faces adjacent to the edges of f
                for(unsigned int i=0;i<f_edges.size();i++)
                {
                    TCellID edge_i= f_edges[i];
                    std::vector<TCellID> faces_i = tmpE2F[edge_i];
                    for(unsigned int j = 0; j<faces_i.size(); j++){
                        if(faces_i[j]!=f_id)
                            f.add<Face>(faces_i[j]);
                    }
                }
            }
            
        }
        else if (ARefModel.has(F2N)){
            IGMesh::face_iterator it = m_mesh->faces_begin();
            // map a node id to the incident face ids.
            std::map<TCellID,std::vector<TCellID> > tmpN2F;
            for(;!it.isDone();it.next()){
                Face    f    = it.value();
                TCellID f_id = f.getID();
                std::vector<TCellID> f_nodes;
                f.getIDs<Node>(f_nodes);
                for(unsigned int i=0;i<f_nodes.size();i++)
                {
                    tmpN2F[f_nodes[i]].push_back(f_id);
                }
            }
            //Now we connect each face to adjacent ones using tmpN2F
            for(it=m_mesh->faces_begin();!it.isDone();it.next())
            {
                Face    f    = it.value();
                TCellID f_id = f.getID();
                std::vector<TCellID> f_nodes;
                f.getIDs<Node>(f_nodes);
                //now we get the set of faces adjacent to the edges of f
                for(unsigned int i=0;i<f_nodes.size();i++)
                {
                    TCellID node_id1 = f_nodes[i];
                    TCellID node_id2 = f_nodes[(i+1)%f_nodes.size()];
                    std::vector<TCellID> faces1 = tmpN2F[node_id1];
                    std::vector<TCellID> faces2 = tmpN2F[node_id2];
                    std::vector<TCellID> common_faces = getCommonBut(faces1,faces2,f_id);
                    
                    for(unsigned int i_common = 0; i_common<common_faces.size(); i_common++){
                        f.add<Face>(common_faces[i_common]);
                    }
                    
                }
            }
            
        }
        else
            throw GMDSException("IGMeshDoctor::buildF2F Not yet implemented");
    }
    /*----------------------------------------------------------------------------*/
    void IGMeshDoctor::buildF2R(const MeshModel &ARefModel) const
    {
        if(ARefModel.has(R2F)){
            IGMesh::region_iterator it = m_mesh->regions_begin();
            for(;!it.isDone();it.next()){
                Region r = it.value();
                std::vector<Face> r_faces;
                r.get<Face>(r_faces);
                for(unsigned int i=0;i<r_faces.size();i++)
                {
                    Face f_i = r_faces[i];
                    f_i.add<Region>(r);
                }
            }
        }
        else
            throw GMDSException("IGMeshDoctor::buildF2R Not yet implemented");
    }
    /*----------------------------------------------------------------------------*/
    void IGMeshDoctor::buildR2N(const MeshModel &ARefModel) const
    {
        throw GMDSException("IGMeshDoctor::buildR2N Not yet implemented");
    }
    /*----------------------------------------------------------------------------*/
    void IGMeshDoctor::buildR2E(const MeshModel &ARefModel) const
    {
        //	if(ARefModel.has(R2F) && ARefModel.has(F2E)){
        //		throw GMDSException("IGMeshDoctor::buildR2E Not yet implemented");
        //	}
        //	else if(ARefModel.has(R2F) && ARefModel.has(E2F)){
        //		throw GMDSException("IGMeshDoctor::buildR2E Not yet implemented");
        //	}
        //	else if(ARefModel.has(F2R) && ARefModel.has(F2E)){
        //		throw GMDSException("IGMeshDoctor::buildR2E Not yet implemented");
        //	}
        //	else if(ARefModel.has(F2R) && ARefModel.has(E2F)){
        //		throw GMDSException("IGMeshDoctor::buildR2E Not yet implemented");
        //	}
        //	else if (ARefModel.has(F2R)){
        //		throw GMDSException("IGMeshDoctor::buildR2E Not yet implemented");
        //	}
        //	else if(ARefModel.has(R2N) && ARefModel.has(E2N)){
        //		throw GMDSException("IGMeshDoctor::buildR2E Not yet implemented");
        //	}
        //	else if(ARefModel.has(R2N) && ARefModel.has(N2E)){
        //		throw GMDSException("IGMeshDoctor::buildR2E Not yet implemented");
        //	}
        //	else if(ARefModel.has(N2R) && ARefModel.has(E2N)){
        //		throw GMDSException("IGMeshDoctor::buildR2E Not yet implemented");
        //	}
        //	else if(ARefModel.has(N2R) && ARefModel.has(N2E)){
        //		throw GMDSException("IGMeshDoctor::buildR2E Not yet implemented");
        //	}
        //	else
        {
            //we just have R2N and E2N
            std::map<TCellID,std::vector<TCellID> > tmp_N2E;
            IGMesh::edge_iterator it_e = m_mesh->edges_begin();
            for(;!it_e.isDone();it_e.next()){
                Edge e= it_e.value();
                TCellID e_id = e.getID();
                std::vector<TCellID> e_node_ids;
                e.getIDs<Node>(e_node_ids);
                for(unsigned int i=0;i<e_node_ids.size();i++)
                {
                    tmp_N2E[e_node_ids[i]].push_back(e_id);
                }
            }
            IGMesh::region_iterator it_r = m_mesh->regions_begin();
            for(;!it_r.isDone();it_r.next()){
                Region r = it_r.value();
                std::vector<TCellID> r_node_ids;
                r.getIDs<Node>(r_node_ids);
                std::vector<TCellID> multiset_edges;
                for(unsigned int i=0;i<r_node_ids.size();i++)
                {
                    std::vector<TCellID> current_edges = tmp_N2E[r_node_ids[i]];
                    multiset_edges.insert(multiset_edges.end(),current_edges.begin(),current_edges.end());
                }
                std::vector<TCellID> r_face_ids = keepFilter(multiset_edges,2);
                r.set<Edge>(r_face_ids);
            }
            
            
        }
        
    }
    /*----------------------------------------------------------------------------*/
    void IGMeshDoctor::buildR2F(const MeshModel &ARefModel) const
    {
        //faces have been built before calling this method, thus F and F2N exist.
        
        //	if(ARefModel.has(R2E) && ARefModel.has(F2E)){
        //		throw GMDSException("IGMeshDoctor::buildR2F Not yet implemented");
        //	}
        //	else if(ARefModel.has(R2E) && ARefModel.has(E2F)){
        //		throw GMDSException("IGMeshDoctor::buildR2F Not yet implemented");
        //	}
        //	else if(ARefModel.has(E2R) && ARefModel.has(F2E)){
        //		throw GMDSException("IGMeshDoctor::buildR2F Not yet implemented");
        //	}
        //	else if(ARefModel.has(E2R) && ARefModel.has(E2F)){
        //		throw GMDSException("IGMeshDoctor::buildR2F Not yet implemented");
        //	}
        //	else if(ARefModel.has(R2E) && ARefModel.has(F2E)){
        //		throw GMDSException("IGMeshDoctor::buildR2F Not yet implemented");
        //	}
        //	else if(ARefModel.has(R2N) && ARefModel.has(E2F)){
        //		throw GMDSException("IGMeshDoctor::buildR2F Not yet implemented");
        //	}
        //	else if(ARefModel.has(N2R) && ARefModel.has(F2E)){
        //		throw GMDSException("IGMeshDoctor::buildR2F Not yet implemented");
        //	}
        //	else if(ARefModel.has(N2R) && ARefModel.has(E2F)){
        //		throw GMDSException("IGMeshDoctor::buildR2F Not yet implemented");
        //	}
        //	else if (ARefModel.has(F2R)){
        //		throw GMDSException("Not yet implemented");
        //	}
        //else
        {
            //we just have R2N and F2N
            std::map<TCellID,std::vector<TCellID> > tmp_N2F;
            IGMesh::face_iterator it_f = m_mesh->faces_begin();
            for(;!it_f.isDone();it_f.next()){
                Face f= it_f.value();
                TCellID f_id = f.getID();
                std::vector<TCellID> f_node_ids;
                f.getIDs<Node>(f_node_ids);
                for(unsigned int i=0;i<f_node_ids.size();i++)
                {
                    tmp_N2F[f_node_ids[i]].push_back(f_id);
                }
            }
            IGMesh::region_iterator it_r = m_mesh->regions_begin();
            for(;!it_r.isDone();it_r.next()){
                Region r = it_r.value();
                std::vector<TCellID> r_node_ids;
                r.getIDs<Node>(r_node_ids);
                std::vector<TCellID> multiset_faces;
                for(unsigned int i=0;i<r_node_ids.size();i++)
                {
                    std::vector<TCellID> current_faces = tmp_N2F[r_node_ids[i]];
                    multiset_faces.insert(multiset_faces.end(),current_faces.begin(),current_faces.end());
                }
                std::vector<TCellID> r_face_ids = keepFilter(multiset_faces,3);
                r.set<Face>(r_face_ids);
            }
            
        }
    }
    /*----------------------------------------------------------------------------*/
    void IGMeshDoctor::buildR2R(const MeshModel &ARefModel) const
    {
        if(ARefModel.has(R2F) && ARefModel.has(F2R)){
            throw GMDSException("IGMeshDoctor::buildR2R Not yet implemented");
        }
        else if(ARefModel.has(R2N) && ARefModel.has(N2R)){
            throw GMDSException("IGMeshDoctor::buildR2R Not yet implemented");
        }
        else if(ARefModel.has(R2E) && ARefModel.has(E2R)){
            throw GMDSException("IGMeshDoctor::buildR2R Not yet implemented");
        }
        else if (ARefModel.has(F2R)){
            throw GMDSException("IGMeshDoctor::buildR2R Not yet implemented");
        }
        else
            throw GMDSException("IGMeshDoctor::buildR2R Not yet implemented");
    }
    /*----------------------------------------------------------------------------*/
} // namespace gmds
/*----------------------------------------------------------------------------*/
