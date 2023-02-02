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
/** \file    Log.cpp
 *  \author  F. LEDOUX
 *  \date    09/17/2009
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/Utils/OrientedGraph.h>
/*----------------------------------------------------------------------------*/
namespace gmds{
    /*----------------------------------------------------------------------------*/
    GraphNode::GraphNode(const TCellID ANb)
    :m_id(ANb)
    {}
    /*----------------------------------------------------------------------------*/
    bool GraphNode::addOutEdge(GraphEdge* AEdge)
    {
        if(AEdge->tail()!=this)
            return false;
        
        m_out_edges.push_back(AEdge);
        return true;
    }
    /*----------------------------------------------------------------------------*/
    bool GraphNode::addInEdge(GraphEdge* AEdge)
    {
        if(AEdge->head()!=this)
            return false;
        
        m_in_edges.push_back(AEdge);
        return  true;
    }
    /*----------------------------------------------------------------------------*/
    TCellID GraphNode::id()
    {
        return m_id;
    }
    
    /*----------------------------------------------------------------------------*/
    std::vector<GraphEdge*>& GraphNode::outEdges()
    {
        return m_out_edges;
    }
    /*----------------------------------------------------------------------------*/
    std::vector<GraphEdge*>& GraphNode::inEdges()
    {
        return m_in_edges;
    }
    /*----------------------------------------------------------------------------*/
    GraphEdge::GraphEdge(const TCellID AID,  GraphNode* ATail,  GraphNode* AHead)
    :m_id(AID), m_tail(ATail),m_head(AHead){
    }
    /*----------------------------------------------------------------------------*/
    TCellID GraphEdge::id()
    {
        return m_id;
    }
    /*----------------------------------------------------------------------------*/
    std::vector<GraphEdge*> GraphEdge::getEdgesStartingFrom(GraphNode* ANode){
        std::vector<GraphEdge*> edges;
        
        if(ANode!=m_head && ANode!=m_tail){
            return edges;
        }
        
        std::vector<GraphEdge*> all_edges = ANode->outEdges();
        for(int e=0; e<all_edges.size(); e++){
            if(all_edges[e]->id()!=this->id()){
//                std::cout<<"Sharing edge: "<<e->tail()->id()<<", "
//                <<e->head()->id()<<std::endl;
                edges.push_back(all_edges[e]);
            }
        }
        
        return edges;
    }
    
    /*----------------------------------------------------------------------------*/
    GraphNode* GraphEdge::tail()
    {
        return m_tail;
    }
    /*----------------------------------------------------------------------------*/
    GraphNode* GraphEdge::head()
    {
        return m_head;
    }
    /*----------------------------------------------------------------------------*/
    OrientedGraph::OrientedGraph(const int ANbNodes)
    {
        m_nodes.resize(ANbNodes);
        for(auto i=0; i<ANbNodes;i++){
            m_nodes[i].m_id=i;
            m_map_nodes[i]=&(m_nodes[i]);
        }
    }
    /*----------------------------------------------------------------------------*/
    OrientedGraph::OrientedGraph(const std::set<TCellID>& AIDS)
    {
        m_nodes.resize(AIDS.size());
        auto i=0;
        for(std::set<TCellID>::iterator id=AIDS.begin(); id!=AIDS.end(); id++){
            m_nodes[i].m_id=*id;
            m_map_nodes[*id]=&(m_nodes[i]);
            i++;
        }
    }
    /*----------------------------------------------------------------------------*/
    OrientedGraph::~OrientedGraph()
    {}
    /*----------------------------------------------------------------------------*/
    bool OrientedGraph::addEdge(const TCellID AID,
                                const TCellID ATail,
                                const TCellID AHead)
    {
        m_edges.push_back(GraphEdge(AID,
                                    m_map_nodes[ATail],
                                    m_map_nodes[AHead]));
        return true;
    }
    
    /*----------------------------------------------------------------------------*/
    void OrientedGraph::updateNodes()
    {
        
        for(int i=0; i<m_edges.size(); i++){
            GraphEdge* ei = &m_edges[i];
            ei->tail()->addOutEdge(ei);
            ei->head()->addInEdge(ei);
            
        }
    }
    /*----------------------------------------------------------------------------*/
    GraphNode* OrientedGraph::node(const TCellID AIndex)
    {
        if(m_map_nodes.find(AIndex)==m_map_nodes.end())
            throw GMDSException("Not found");
        
        return m_map_nodes[AIndex];
        
    }
    /*----------------------------------------------------------------------------*/
    GraphEdge* OrientedGraph::edge(const TCellID AIndex)
    {
        return &(m_edges[AIndex]);
    }
    /*----------------------------------------------------------------------------*/
    int OrientedGraph::nbNodes() const
    {
        return m_nodes.size();
    }
    /*----------------------------------------------------------------------------*/
    int OrientedGraph::nbEdges() const
    {
        return m_edges.size();
    }
}
/*----------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& AStr,  gmds::OrientedGraph& AG)
{
    AStr<<"Oriented Graph ("<<AG.nbNodes()<<", "<<AG.nbEdges()<<")"<<std::endl;
    std::cout<<" Nodes: ";
    for(int i=0; i<AG.m_nodes.size(); i++)
        std::cout<<AG.m_nodes[i].id()<<" ";
    std::cout<<"\n Edges: ";
    for(int i=0; i<AG.m_edges.size(); i++)
        std::cout<<"["<<AG.m_edges[i].head()->id()<<", "<<AG.m_edges[i].tail()->id()<<"] ";
    std::cout<<std::endl;
    return AStr;
}

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
