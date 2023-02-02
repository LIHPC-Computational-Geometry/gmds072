/*----------------------------------------------------------------------------*/
/* Copyright: CEA
 * contributors: R. Viertel (2016)
 *
 * rvierte@sandia.gov
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
/*---------------------------------------------------------------------------*/
#include <THexer.h>

using namespace fhedo;

void THexer::execute()
{
  initMarks();

  buildNodesAndEdges();

  buildHexFacesOnTetFaces();

  buildHexFacesInsideTetRegions();

  buildHexRegions();

  buildNodeMap();

  cleanMarks();
}

void THexer::buildNodesAndEdges()
{
  IGMesh::region_iterator r_it = m_Tmesh->regions_begin();

  for(; !r_it.isDone();r_it.next())
  {
    Region tet = r_it.value();

    //build hex node at center of tet
    Node hex_center_of_tet = m_Hmesh->newNode(tet.center());
    m_tets_to_hex_nodes[tet.getID()] = hex_center_of_tet;

    std::vector<Face> tet_faces = tet.get<Face>();

    for(Face tet_f: tet_faces)
    {
      //get hex node corresponding to center of tet face
      if(m_Tmesh->isMarked(tet_f,m_face_on_hmesh))
      {
        Node hex_center_of_tet_face = m_tet_faces_to_hex_nodes[tet_f.getID()];
        Edge new_edge = m_Hmesh->newEdge(hex_center_of_tet,hex_center_of_tet_face);
        hex_center_of_tet.add<Edge>(new_edge);
        hex_center_of_tet_face.add<Edge>(new_edge);
      }
      else
      {
        Node hex_center_of_tet_face = m_Hmesh->newNode(tet_f.center());
        m_tet_faces_to_hex_nodes[tet_f.getID()] = hex_center_of_tet_face;
        m_Tmesh->mark(tet_f,m_face_on_hmesh);
        Edge new_edge = m_Hmesh->newEdge(hex_center_of_tet,hex_center_of_tet_face);
        hex_center_of_tet.add<Edge>(new_edge);
        hex_center_of_tet_face.add<Edge>(new_edge);


        std::vector<Edge> tet_face_edges = tet_f.get<Edge>();

        for(Edge tet_e: tet_face_edges)
        {
          //get hex node corresponding to center of tet edge
          if(m_Tmesh->isMarked(tet_e,m_edge_on_hmesh))
          {
            Node hex_center_of_tet_edge = m_tet_edges_to_hex_nodes[tet_e.getID()];
            Edge new_edge = m_Hmesh->newEdge(hex_center_of_tet_face,hex_center_of_tet_edge);
            hex_center_of_tet_face.add<Edge>(new_edge);
            hex_center_of_tet_edge.add<Edge>(new_edge);
          }
          else
          {
            Node hex_center_of_tet_edge = m_Hmesh->newNode(tet_e.center());
            m_tet_edges_to_hex_nodes[tet_e.getID()] = hex_center_of_tet_edge;
            m_Tmesh->mark(tet_e,m_edge_on_hmesh);
            Edge new_edge = m_Hmesh->newEdge(hex_center_of_tet_face,hex_center_of_tet_edge);
            hex_center_of_tet_face.add<Edge>(new_edge);
            hex_center_of_tet_edge.add<Edge>(new_edge);

            std::vector<Node> tet_edge_nodes = tet_e.get<Node>();

            for(Node tet_n: tet_edge_nodes)
            {
              if(m_Tmesh->isMarked(tet_n,m_node_on_hmesh))
              {
                Node hex_corner_tet_node = m_tet_nodes_to_hex_nodes[tet_n.getID()];
                Edge new_edge = m_Hmesh->newEdge(hex_center_of_tet_edge,hex_corner_tet_node);
                hex_center_of_tet_edge.add<Edge>(new_edge);
                hex_corner_tet_node.add<Edge>(new_edge);
              }
              else
              {
                Node hex_corner_tet_node = m_Hmesh->newNode(tet_n.getPoint());
                m_tet_nodes_to_hex_nodes[tet_n.getID()] = hex_corner_tet_node;
                m_Tmesh->mark(tet_n,m_node_on_hmesh);
                Edge new_edge = m_Hmesh->newEdge(hex_center_of_tet_edge,hex_corner_tet_node);
                hex_center_of_tet_edge.add<Edge>(new_edge);
                hex_corner_tet_node.add<Edge>(new_edge);
              }
            }
          }
        }
      }
    }
  }
}

void THexer::buildHexFacesOnTetFaces()
{
  IGMesh::face_iterator f_it = m_Tmesh->faces_begin();
  for(; !f_it.isDone();f_it.next())
  {
    Face f = f_it.value();

    std::vector<Edge> edges = f.get<Edge>();
    std::vector<Node> nodes = f.get<Node>();

    for(Node n: nodes)
    {
      std::vector<Edge> adj_edges;
      //get edges on face adjacent to node n
      for(Edge e: edges)
      {
        std::vector<Node> edge_nodes = e.get<Node>();
        if(edge_nodes[0].getID() == n.getID() || edge_nodes[1].getID() == n.getID())
        {
          adj_edges.push_back(e);
        }
      }

      Node node0 = m_tet_nodes_to_hex_nodes[n.getID()];
      Node node1 = m_tet_edges_to_hex_nodes[adj_edges[0].getID()];
      Node node2 = m_tet_faces_to_hex_nodes[f.getID()];
      Node node3 = m_tet_edges_to_hex_nodes[adj_edges[1].getID()];

      Face new_quad = m_Hmesh->newQuad(node0,node1,node2,node3);
      node0.add<Face>(new_quad);
      node1.add<Face>(new_quad);
      node2.add<Face>(new_quad);
      node3.add<Face>(new_quad);

      //connect edges to faces
      std::vector<Edge> node1_edges = node1.get<Edge>();

      for(Edge e: node1_edges)
      {
        std::vector<Node> edge_nodes = e.get<Node>();
        for(Node en: edge_nodes)
        {
          if(en.getID() == node0.getID() || en.getID() == node3.getID())
          {
            new_quad.add<Edge>(e);
            e.add<Face>(new_quad);
          }
        }
      }

      std::vector<Edge> node3_edges = node3.get<Edge>();

      for(Edge e: node3_edges)
      {
        std::vector<Node> edge_nodes = e.get<Node>();
        for(Node en: edge_nodes)
        {
          if(en.getID() == node0.getID() || en.getID() == node2.getID())
          {
            new_quad.add<Edge>(e);
            e.add<Face>(new_quad);
          }
        }
      }
    }
  }
}

void THexer::buildHexFacesInsideTetRegions()
{
  IGMesh::region_iterator r_it = m_Tmesh->regions_begin();

  for(;!r_it.isDone();r_it.next())
  {
    Region tet = r_it.value();

    std::vector<Edge> tet_edges = tet.get<Edge>();
    std::vector<Face> tet_faces = tet.get<Face>();

    for(Edge tet_e: tet_edges)
    {
      std::vector<Face> adj_faces;
      //get faces on tet adjacent to edge tet_e
      for(Face tet_f: tet_faces)
      {
        std::vector<Edge> face_edges = tet_f.get<Edge>();
        if( face_edges[0].getID() == tet_e.getID() ||
            face_edges[1].getID() == tet_e.getID() ||
            face_edges[2].getID() == tet_e.getID()    )
        {
          adj_faces.push_back(tet_f);
        }
      }
      Node node0 = m_tets_to_hex_nodes[tet.getID()];
      Node node1 = m_tet_faces_to_hex_nodes[adj_faces[0].getID()];
      Node node2 = m_tet_edges_to_hex_nodes[tet_e.getID()];
      Node node3 = m_tet_faces_to_hex_nodes[adj_faces[1].getID()];

      Face new_quad = m_Hmesh->newQuad(node0,node1,node2,node3);
      node0.add<Face>(new_quad);
      node1.add<Face>(new_quad);
      node2.add<Face>(new_quad);
      node3.add<Face>(new_quad);

      //connect edges to faces
      std::vector<Edge> node1_edges = node1.get<Edge>();

      for(Edge e: node1_edges)
      {
        std::vector<Node> edge_nodes = e.get<Node>();
        for(Node en: edge_nodes)
        {
          if(en.getID() == node0.getID() || en.getID() == node3.getID())
          {
            new_quad.add<Edge>(e);
            e.add<Face>(new_quad);
          }
        }
      }

      std::vector<Edge> node3_edges = node3.get<Edge>();

      for(Edge e: node3_edges)
      {
        std::vector<Node> edge_nodes = e.get<Node>();
        for(Node en: edge_nodes)
        {
          if(en.getID() == node0.getID() || en.getID() == node2.getID())
          {
            new_quad.add<Edge>(e);
            e.add<Face>(new_quad);
          }
        }
      }
    }
  }
}

void THexer::buildHexRegions()
{
  IGMesh::region_iterator r_it = m_Tmesh->regions_begin();

  for(;!r_it.isDone();r_it.next())
  {
    Region tet = r_it.value();

    std::vector<Node> nodes = tet.get<Node>();
    std::vector<Face> faces = tet.get<Face>();

    //build hex on each node
    for(Node n: nodes)
    {
      std::vector<Face> adj_faces;
      for(Face f: faces)
      {
        std::vector<Node> face_nodes = f.get<Node>();
        if( face_nodes[0].getID() == n.getID() ||
            face_nodes[1].getID() == n.getID() ||
            face_nodes[2].getID() == n.getID()    )
        {
          adj_faces.push_back(f);
        }
      }

      Face base = adj_faces[0];
      Edge edge01;
      Edge edge12;
      Edge edge20;

      std::vector<Edge> base_edges = base.get<Edge>();
      for(Edge base_e: base_edges)
      {
        std::vector<Face> edge_faces = base_e.get<Face>();

        for(Face edge_f: edge_faces)
        {
          if(edge_f.getID() == adj_faces[1].getID())
          {
            edge01 = base_e;
          }

          if(edge_f.getID() == adj_faces[2].getID())
          {
            edge20 = base_e;
          }
        }
      }

      std::vector<Edge> face1_edges = adj_faces[1].get<Edge>();
      for(Edge face1_e: face1_edges)
      {
        std::vector<Face> edge_faces = face1_e.get<Face>();
        for(Face edge_f: edge_faces)
        {
          if(edge_f.getID() == adj_faces[2].getID())
          {
            edge12 = face1_e;
            break;
          }
        }
      }

      Node node0 = m_tet_faces_to_hex_nodes[base.getID()];
      Node node1 = m_tet_edges_to_hex_nodes[edge01.getID()];
      Node node2 = m_tet_nodes_to_hex_nodes[n.getID()];
      Node node3 = m_tet_edges_to_hex_nodes[edge20.getID()];

      Node node4 = m_tets_to_hex_nodes[tet.getID()];
      Node node5 = m_tet_faces_to_hex_nodes[adj_faces[1].getID()];
      Node node6 = m_tet_edges_to_hex_nodes[edge12.getID()];
      Node node7 = m_tet_faces_to_hex_nodes[adj_faces[2].getID()];

      Region new_hex = m_Hmesh->newHex(node0,node1,node2,node3,
                                      node4,node5,node6,node7);
      node0.add<Region>(new_hex);
      node1.add<Region>(new_hex);
      node2.add<Region>(new_hex);
      node3.add<Region>(new_hex);
      node4.add<Region>(new_hex);
      node5.add<Region>(new_hex);
      node6.add<Region>(new_hex);
      node7.add<Region>(new_hex);


      //connect faces and hexes
      std::vector<Face> node2_faces = node2.get<Face>();

      for(Face f: node2_faces)
      {
        std::vector<Node> face_nodes = f.get<Node>();
        for(Node fn: face_nodes)
        {
          if(fn.getID() == node0.getID() || fn.getID() == node5.getID() ||
             fn.getID() == node7.getID()                                   )
          {
            new_hex.add<Face>(f);
            f.add<Region>(new_hex);
            break;
          }
        }
      }

      std::vector<Face> node4_faces = node4.get<Face>();

      for(Face f: node4_faces)
      {
        std::vector<Node> face_nodes = f.get<Node>();
        for(Node fn: face_nodes)
        {
          if(fn.getID() == node1.getID() || fn.getID() == node3.getID() ||
             fn.getID() == node6.getID()                                   )
          {
            new_hex.add<Face>(f);
            f.add<Region>(new_hex);
            break;
          }
        }
      }

      //connect edges and hexes
      std::vector<Edge> node2_edges = node2.get<Edge>();
      for(Edge e: node2_edges)
      {
        std::vector<Node> edge_nodes = e.get<Node>();
        for(Node n: edge_nodes)
        {
          if(n.getID() == node1.getID() || n.getID() == node3.getID() ||
             n.getID() == node6.getID()                                  )
          {
            new_hex.add<Edge>(e);
            e.add<Region>(new_hex);
          }
        }
      }

      std::vector<Edge> node0_edges = node0.get<Edge>();
      for(Edge e: node0_edges)
      {
        std::vector<Node> edge_nodes = e.get<Node>();
        for(Node n: edge_nodes)
        {
          if(n.getID() == node1.getID() || n.getID() == node3.getID() ||
             n.getID() == node4.getID()                                  )
          {
            new_hex.add<Edge>(e);
            e.add<Region>(new_hex);
          }
        }
      }

      std::vector<Edge> node5_edges = node5.get<Edge>();
      for(Edge e: node5_edges)
      {
        std::vector<Node> edge_nodes = e.get<Node>();
        for(Node n: edge_nodes)
        {
          if(n.getID() == node1.getID() || n.getID() == node6.getID() ||
             n.getID() == node4.getID()                                  )
          {
            new_hex.add<Edge>(e);
            e.add<Region>(new_hex);
          }
        }
      }

      std::vector<Edge> node7_edges = node7.get<Edge>();
      for(Edge e: node7_edges)
      {
        std::vector<Node> edge_nodes = e.get<Node>();
        for(Node n: edge_nodes)
        {
          if(n.getID() == node3.getID() || n.getID() == node6.getID() ||
             n.getID() == node4.getID()                                  )
          {
            new_hex.add<Edge>(e);
            e.add<Region>(new_hex);
          }
        }
      }
    }
  }
}

void THexer::buildNodeMap()
{
  Variable<TCellID>* H_node_map;
  Variable<TCellID>* T_node_map;
  H_node_map = m_Hmesh->newVariable<TCellID>(GMDS_NODE, "node_map");
  T_node_map = m_Tmesh->newVariable<TCellID>(GMDS_NODE, "node_map");

  IGMesh::node_iterator hex_n_it = m_Hmesh->nodes_begin();

  for(;!hex_n_it.isDone();hex_n_it.next())
  {
    Node hex_n = hex_n_it.value();
    (*H_node_map)[hex_n.getID()] = -99;
  }


  IGMesh::node_iterator n_it = m_Tmesh->nodes_begin();

  for(;!n_it.isDone();n_it.next())
  {
    Node n = n_it.value();
    Node hex_node = m_tet_nodes_to_hex_nodes[n.getID()];
    (*H_node_map)[hex_node.getID()] = n.getID();
    (*T_node_map)[n.getID()] = hex_node.getID();
  }
}

void THexer::initMarks()
{
  m_node_on_hmesh = m_Tmesh->getNewMark<Node>();
  m_edge_on_hmesh = m_Tmesh->getNewMark<Edge>();
  m_face_on_hmesh = m_Tmesh->getNewMark<Face>();
  m_edge_on_hex   = m_Hmesh->getNewMark<Edge>();
}

void THexer::cleanMarks()
{
  m_Tmesh->unmarkAll<Node>(m_node_on_hmesh);
  m_Tmesh->unmarkAll<Edge>(m_edge_on_hmesh);
  m_Tmesh->unmarkAll<Face>(m_face_on_hmesh);
  m_Hmesh->unmarkAll<Edge>(m_edge_on_hex);

  m_Tmesh->freeMark<Node>(m_node_on_hmesh);
  m_Tmesh->freeMark<Edge>(m_edge_on_hmesh);
  m_Tmesh->freeMark<Face>(m_face_on_hmesh);
  m_Hmesh->freeMark<Edge>(m_edge_on_hex);
}
