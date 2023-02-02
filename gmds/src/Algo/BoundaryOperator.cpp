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
/** \file    BoundaryOperator.t.h
 *  \author  F. LEDOUX
 *  \date    08/08/2008
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/Algo/BoundaryOperator.h>
#include <GMDS/CAD/GeomEntity.h>
/*----------------------------------------------------------------------------*/
#include <sstream>
#include <set>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
BoundaryOperator::BoundaryOperator(IGMesh* AMesh)
:m_mesh(AMesh)
{}
/*----------------------------------------------------------------------------*/
BoundaryOperator::~BoundaryOperator()
{}
/*----------------------------------------------------------------------------*/
bool BoundaryOperator::isValid() const
{
  MeshModel model = m_mesh->getModel();
  if(model.has(R)) {
    if (!model.has(F2R))
      return false;
    if (!model.has(F2E))
      return false;
    if (!model.has(E2F))
      return false;

    return true;
  }
  else if(model.has(F)) {
    if (!model.has(F2E))
      return false;
    if (!model.has(E2F))
      return false;

    return true;
  }

  return false;
}
/*----------------------------------------------------------------------------*/
void BoundaryOperator::getBoundaryNodes(std::vector<TCellID>& ANodeIDs)
{
  std::set<TCellID> elts;
  MeshModel model = m_mesh->getModel();
  if(model.has(R)) {
    throw GMDSException("Not yet implemented in 3D");
  }
  else { 
    IGMesh::edge_iterator it = m_mesh->edges_begin();
    
    for (; !it.isDone(); it.next())  {
      Edge e = it.value();
      std::vector<TCellID> faces_e = e.getIDs<Face>();
      if(faces_e.size()==1){
	std::vector<TCellID> nodes_e = e.getIDs<Node>();
	elts.insert(nodes_e.begin(),nodes_e.end());
      }
    }

  } 
  ANodeIDs.clear();
  ANodeIDs.insert(ANodeIDs.end(),elts.begin(),elts.end());
}
/*----------------------------------------------------------------------------*/
void BoundaryOperator::
markCellOnGeometry(const int AMarkFOnSurf,
                   const int AMarkEOnSurf,
                   const int AMarkNOnSurf,
                   const int AMarkEOnCurve,
                   const int AMarkNOnCurve,
                   const int AMarkNOnPnt,
                   const int AMarkAN)
{   
  MeshModel model = m_mesh->getModel();
  if(model.has(R)) {
    markCellsOnSurfaces(AMarkFOnSurf, AMarkEOnSurf, AMarkNOnSurf); 
    markCellsOnCurves(AMarkFOnSurf, AMarkEOnSurf, AMarkEOnCurve, AMarkNOnCurve);
    markNodesOnPoint(AMarkEOnCurve, AMarkNOnCurve,AMarkNOnPnt);
    markAloneNodes(AMarkAN);
    colorFaces(AMarkFOnSurf,AMarkEOnCurve);
  }
  else {  
   
    //all the faces, edges and nodes are on the surface, so we invert the mesh 
    //mask
    m_mesh->negateMaskMark<Node>(AMarkNOnSurf);
    m_mesh->negateMaskMark<Edge>(AMarkEOnSurf);
    m_mesh->negateMaskMark<Face>(AMarkFOnSurf);
  
    markCellsOnCurves(AMarkEOnCurve, AMarkNOnCurve);
      
    markNodesOnPoint(AMarkEOnCurve, AMarkNOnCurve,AMarkNOnPnt);
      
    markAloneNodes(AMarkAN);
    
    colorFaces(AMarkFOnSurf,AMarkEOnCurve); 
  }
}

/*----------------------------------------------------------------------------*/
void BoundaryOperator::markAloneNodes(const int AMarkAlone) 
{
    IGMesh::node_iterator it = m_mesh->nodes_begin();

  int cpt = 0;
  for (; !it.isDone(); it.next())
    {
      Node n = it.value();
      if (n.getNbEdges() == 0) {
	m_mesh->mark(n, AMarkAlone);
	cpt++;
      }

    }
  std::cout << "Nodes alones : " << cpt << std::endl;
}
/*----------------------------------------------------------------------------*/
void BoundaryOperator::
markCellsOnSurfaces(const int AMarkBF, const int AMarkBE, const int AMarkBN)
{
    IGMesh::face_iterator it = m_mesh->faces_begin();
    
    int cpt1 = 0, cpt2 = 0;
    
    for (; !it.isDone(); it.next())
    {
        Face f = it.value();
        cpt1++;
        if (f.get<Region>().size() == 1)
        {
            m_mesh->mark(f, AMarkBF);
            std::vector<Edge> f_edges = f.get<Edge>();
            for (unsigned int i = 0; i < f_edges.size(); i++)
            {
                m_mesh->mark(f_edges[i], AMarkBE);
            }
            std::vector<Node> f_nodes = f.get<Node>();
            for (unsigned int i = 0; i < f_nodes.size(); i++)
            {
                m_mesh->mark(f_nodes[i], AMarkBN);
            }
            cpt2++;
        }
    }
    std::cout << cpt2 << " marked faces on " << cpt1 << std::endl;
}
/*----------------------------------------------------------------------------*/
void BoundaryOperator::
markCellsOnCurves(const int AMarkBF,    //mark for faces on surfaces //IN
                  const int AMarkBE,            //mark for edges on surfaces //IN
                  const int AMarkCE,    //mark for edges on curves //OUT
                  const int AMarkCN)            //mark for nodes on curves //OUT
{
    IGMesh::edge_iterator it = m_mesh->edges_begin();
    
    int cpt1 = 0, cpt2 = 0;
    for (; !it.isDone(); it.next())  {
        Edge e = it.value();
        
        if (m_mesh->isMarked(e, AMarkBE)) {
            cpt1++;
            //on va calculer les normales des faces adjacentes et au bord
            //et leur prod scalaire
            std::vector<Face> adj_faces = e.get<Face>();
            if(adj_faces.size()==1) {
                //2D boundary edge
                m_mesh->mark(e, AMarkCE);
                std::vector<Node> e_nodes = e.get<Node>();
                m_mesh->mark(e_nodes[0], AMarkCN);
                m_mesh->mark(e_nodes[1], AMarkCN);
            } // if(adj_faces.size()==1) {
            else{
                std::vector<Face> boundary_adj_faces;
                
                for (unsigned int i = 0; i < adj_faces.size(); i++) {
                    Face current_face = adj_faces[i];
                    if (m_mesh->isMarked(current_face, AMarkBF))
                        boundary_adj_faces.push_back(current_face);
                }
                
                if (boundary_adj_faces.size() != 2){
                    throw GMDSException("a boundary edge should be adjacent to 2 boundary faces!!!");
                }
                
                Face f0 = boundary_adj_faces[0];
                Face f1 = boundary_adj_faces[1];
                
                //LA SUITE DES VECTEURS A CONSTRUIRE A BASE DE POINTS
                math::Vector n0 = getOutputNormalOfABoundaryFace(f0);
                math::Vector n1 = getOutputNormalOfABoundaryFace(f1);
                
                
                double dotProduct = n0.dot(n1);
                double angle_dev = 45*math::Constants::PI/180;
                if (dotProduct < (sqrt(2.0) / 2.0)){
                    cpt2++;
                    m_mesh->mark(e, AMarkCE);
                    std::vector<Node> e_nodes = e.get<Node>();
                    m_mesh->mark(e_nodes[0], AMarkCN);
                    m_mesh->mark(e_nodes[1], AMarkCN);
                }
            }
        }//else{
    } //for (; !it.isDone(); it.next())  {
    std::cout << cpt2 << " marked edges on " << cpt1 <<" boundary edges"<< std::endl;
}
/*----------------------------------------------------------------------------*/
void BoundaryOperator::
markCellsOnCurves(const int AMarkCE, //mark for edges on curves //OUT
                  const int AMarkCN) //mark for nodes on curves //OUT
{
  IGMesh::edge_iterator it = m_mesh->edges_begin();

  int cpt1 = 0, cpt2 = 0;
  for (; !it.isDone(); it.next())  {
    Edge e = it.value();
    cpt1++; 
    std::vector<Face> adj_faces = e.get<Face>();
    if(adj_faces.size()==1) {
      //2D boundary edge
      m_mesh->mark(e, AMarkCE);
      std::vector<Node> e_nodes = e.get<Node>();
      m_mesh->mark(e_nodes[0], AMarkCN);
      m_mesh->mark(e_nodes[1], AMarkCN);
      cpt2++;
    }
  } //for (; !it.isDone(); it.next())  {
  std::cout << cpt2 << " marked edges on " << cpt1 <<" boundary edges"<< std::endl;
}
/*----------------------------------------------------------------------------*/
void BoundaryOperator::markNodesOnPoint(const int AMarkCE,// edge on curve IN
                                        const int AMarkCN,// node on curve IN
                                        const int AMarkPN)// node on vertex OUT
{
    int cpt1=0, cpt2=0;
    
    for (IGMesh::node_iterator it = m_mesh->nodes_begin();
         !it.isDone(); it.next()) {
        Node n = it.value();
        cpt1++;
        if (m_mesh->isMarked(n, AMarkCN)){
            //We have a node on curve
            std::vector<Edge> adj_edges = n.get<Edge>();
            int cpt_tmp = 0;
            for (int ei=0; ei<adj_edges.size(); ei++){
                if (m_mesh->isMarked(adj_edges[ei], AMarkCE))
                    cpt_tmp++;
            }
            if (cpt_tmp > 2) {
                m_mesh->mark(n, AMarkPN);
                cpt2++;
            }
            else if (cpt_tmp == 2){
                //check if we have a brutal normal change
                
                //First we get all the nodes connected to n by a
                // boundary edge
                std::vector<Node> connected_nodes;
                
                for (int ei=0; ei<adj_edges.size(); ei++){
                    if (m_mesh->isMarked(adj_edges[ei], AMarkCE)){
                        std::vector<Node> edge_nodes = adj_edges[ei].get<Node>();
                        for (int nj=0; nj<edge_nodes.size(); nj++){
                            if (edge_nodes[nj] != n)
                                connected_nodes.push_back(edge_nodes[nj]);
                        }
                    }
                }
                
                Node n0 = connected_nodes[0];
                Node n1 = connected_nodes[1];
                
                math::Vector v0(n.getPoint(), n0.getPoint());
                math::Vector v1(n.getPoint(), n1.getPoint());
                v0.normalize();
                v1.normalize();
                double dotProduct = v0.dot(v1);
                // if we have a brutal normal change
                if (dotProduct > (-sqrt(2.0) / 2.0)){
                    m_mesh->mark(n, AMarkPN);
                    cpt2++;
                }

            }
        }
    }
    std::cout << cpt2 << " marked nodes among " << cpt1 << std::endl;
}
/*----------------------------------------------------------------------------*/
void BoundaryOperator::
colorFaces(const int AMarkFOnSurf, const int AMarkEOnCurv)
{
    
    m_var_color_surf = m_mesh->newVariable<int>(GMDS_FACE, "BND_SURFACE_COLOR");
    
    int color = 1; //Default value is 0
    int markDone = m_mesh->getNewMark<Face>();
    IGMesh::face_iterator it = m_mesh->faces_begin();
    for (; !it.isDone(); it.next() )
    {
        Face f = it.value();
        //on ne considere que les faces au bord
        //qui n'ont pas encore ete traitees
        if (m_mesh->isMarked(f, AMarkFOnSurf) && !m_mesh->isMarked(f, markDone))
        {
            //new surface
            color++; // so new color
            m_mesh->mark(f, markDone);
            (*m_var_color_surf)[f.getID()] = color;
            
            //on se propage et on marque
            std::vector<Face> next;
            next.push_back(f);
            
            while (!next.empty()){
                Face current = next.back();
                next.pop_back();
                //recuperation des faces voisines non traitees et appartenant a la surface
                std::vector<Edge> current_edges = current.get<Edge>();
                
                for (unsigned int ie = 0; ie < current_edges.size(); ie++)
                {
                    Edge ei = current_edges[ie];
                    
                    if (!m_mesh->isMarked(ei, AMarkEOnCurv))//si ce n'est pas une arete au bord
                    {
                        
                        std::vector<Face> f_edges = ei.get<Face>();
                        
                        for (unsigned int ifa = 0; ifa < f_edges.size(); ifa++){
                            Face fi = f_edges[ifa];
                            if (m_mesh->isMarked(fi, AMarkFOnSurf) &&
                                !m_mesh->isMarked(fi, markDone)){
                                m_mesh->mark(fi, markDone);
                                (*m_var_color_surf)[fi.getID()] = color;
                                next.push_back(fi);
                                
                            }
                        }
                    }
                }
            }
        }
    }
    m_mesh->unmarkAll<Face>(markDone);
    m_mesh->freeMark<Face>(markDone);
}
/*----------------------------------------------------------------------------*/
void BoundaryOperator::
colorEdges(const int AMarkEOnCurv, const int AMarkNOnPnt)
{
    
    m_var_color_curve = m_mesh->newVariable<int>(GMDS_EDGE, "BND_CURVE_COLOR");
    
    int color = 1; //Default value is 0
    int markDone = m_mesh->getNewMark<Edge>();
    IGMesh::edge_iterator it = m_mesh->edges_begin();
    std::vector<Edge> done_edges;
    for (; !it.isDone(); it.next() )
    {
        Edge e = it.value();
        // We only go throug edges classigied on curves and that have not been
        // yet handled
        if ( m_mesh->isMarked(e, AMarkEOnCurv) &&
            !m_mesh->isMarked(e, markDone)) {
            //new curve
            color++; // so new color
            m_mesh->mark(e, markDone);
            (*m_var_color_curve)[e.getID()] = color;
            
            //propagation to curve edges sharing a point with e
            std::vector<Edge> next;
            next.push_back(e);
            
            while (!next.empty()){
                Edge current = next.back();
                next.pop_back();
                //We get the ajacent edges that are on a curve but not yet done
                std::vector<Node> current_nodes = current.get<Node>();
                
                for (int ni=0; ni<current_nodes.size(); ni++) {
                    
                    if (!m_mesh->isMarked(current_nodes[ni], AMarkNOnPnt)){
                        //If it is not a node classified on a point, we can found
                        // a next edge on this curve
                        
                        std::vector<Edge> n_edges = current_nodes[ni].get<Edge>();
                        for (int ej=0; ej<n_edges.size(); ej++){
                            if (m_mesh->isMarked(n_edges[ej], AMarkEOnCurv) &&
                                !m_mesh->isMarked(n_edges[ej], markDone)){
                                m_mesh->mark(n_edges[ej], markDone);
                                (*m_var_color_curve)[n_edges[ej].getID()] = color;
                                next.push_back(n_edges[ej]);
                                
                            }
                        }
                    }
                }
            }
        }
    }
    m_mesh->unmarkAll<Edge>(markDone);
    m_mesh->freeMark<Edge>(markDone);
}
/*----------------------------------------------------------------------------*/
void BoundaryOperator::
colorNodes(const int AMarkNOnPnt)
{
    Variable<int>* v_color=0;
    try{
    v_color = m_mesh->newVariable<int>(GMDS_NODE, "BND_VERTEX_COLOR");
    }
    catch(GMDSException& e){
        v_color = m_mesh->getVariable<int>(GMDS_NODE, "BND_VERTEX_COLOR");
        
    }
        
    int color = 0; //Default value is 0
    IGMesh::node_iterator it = m_mesh->nodes_begin();
    for (; !it.isDone(); it.next() )
    {
        Node n = it.value();
        // We only go throug nodes classigied on point
        
        if (m_mesh->isMarked(n, AMarkNOnPnt)) {
            //new point
            color++; // so new color
            (*v_color)[n.getID()] = color;
        }
        
    }
}
/*----------------------------------------------------------------------------*/
math::Vector BoundaryOperator::getOutputNormal(Face& AFace, Region& ARegion)
{
  std::vector<Node> region_nodes = ARegion.get<Node>();
  std::vector<Node> face_nodes = AFace.get<Node>();

  if (region_nodes.size() != 4)
    throw GMDSException("SingularityGraphBuilder::getOutputNormal can only be used on tetrahedral regions");
  if (face_nodes.size() != 3)
    throw GMDSException("SingularityGraphBuilder::getOutputNormal can only be used on triangular faces");

  //we go through all the nodes of ARegion to find the one that do not belong
  //to AFAce
  for (unsigned int i = 0; i<region_nodes.size(); i++){
    Node n = region_nodes[i];
    if (n != face_nodes[0] && n != face_nodes[1] && n != face_nodes[2])
      {
	//n is the node opposite to the face AFace
	Node n0 = face_nodes[0];
	Node n1 = face_nodes[1];
	Node n2 = face_nodes[2];
	math::Vector normal_to_face = AFace.normal();
	math::Vector in_vector(n0.getPoint(), n.getPoint());
	if (normal_to_face.dot(in_vector)>0.0)
	  {
	    return math::Vector(
				-normal_to_face.X(),
				-normal_to_face.Y(),
				-normal_to_face.Z());
	  }
	else
	  {
	    return normal_to_face;
	  }
      }
  }

  throw GMDSException("BoundaryOperator::getOutputNormal face not found");
}
/*----------------------------------------------------------------------------*/
math::Vector BoundaryOperator::
getOutputNormalOfABoundaryFace(Face& AFace)
{
    std::vector<Region> adj_regions = AFace.get<Region>();
    if (adj_regions.size() != 1)
        throw GMDSException("A boundary face must be adjacent to 1 region!!!");
    
    return getOutputNormal(AFace, adj_regions[0]);
}
/*----------------------------------------------------------------------------*/
math::Vector BoundaryOperator::
getOutputNormalOfABoundaryNode(Node& ANode)
{
    std::vector<Face> adj_faces = ANode.get<Face>();
    std::vector<math::Vector> weighted_normals;
    for(unsigned int i=0; i<adj_faces.size(); i++){
        Face fi = adj_faces[i];
        std::vector<Region> fi_regions = fi.get<Region>();
        if (fi_regions.size() == 1){
            // we have a boundary face
            math::Vector vi = getOutputNormal(fi, fi_regions[0]);
            TCoord ai = fi.area();
            weighted_normals.push_back(ai*vi);
            
        }
    }
    // now we compute the normal vector at ANode
    math::Vector n = weighted_normals[0];
    for(unsigned int i=1; i <weighted_normals.size();i++){
        n = n+weighted_normals[i];
    }
    n.normalize();
    return n;
}
/*----------------------------------------------------------------------------*/
