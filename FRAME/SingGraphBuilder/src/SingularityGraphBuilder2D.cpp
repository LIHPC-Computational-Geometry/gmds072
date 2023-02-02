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
 * SingularityGraphBuilder2D.cpp
 *
 *  Created on: 13 juil. 2014
 *      Author: bibi
 */
/*----------------------------------------------------------------------------*/
#include <GMDS/IG/IG.h>
#include <GMDS/IG/IGMeshDoctor.h>
#include <GMDS/IO/MeditReader.h>
#include <GMDS/IO/VTKReader.h>
#include <GMDS/IO/VTKWriter.h>
#include <GMDS/Math/Cross.h>
#include <GMDS/Math/Quaternion.h>
#include <GMDS/Math/Chart.h>
#include <GMDS/Math/Ray.h>
#include <GMDS/Math/Triangle.h>
#include <GMDS/Math/Segment.h>
#include <GMDS/Algo/BoundaryOperator.h>
#include <GMDS/Math/Numerics.h>
#include "SingularityGraphBuilder2D.h"
/*----------------------------------------------------------------------------*/
#include <sstream>
#include <set>
/*----------------------------------------------------------------------------*/
using namespace gmds;
/*----------------------------------------------------------------------------*/
SingularityGraphBuilder2D::
SingularityGraphBuilder2D(IGMesh* AMesh,
                          Variable<math::Cross2D>* AField,
                          const double ATolerance,
                          const bool ABuildGeomSing)
:m_mesh(AMesh), m_field(AField),
m_tool(AMesh,AField),
m_output_directory_name(""),
m_graph(AMesh)
{
    m_build_geometric_singularities= ABuildGeomSing;
    
    if(ATolerance<0.01 || ATolerance>0.1)
        throw GMDSException("SingularityGraphBuilder2D: Tolerance must be taken in [0.01,0.1]");
    IGMesh::node_iterator it_nodes = m_mesh->nodes_begin();
    double x_min, y_min, z_min, x_max, y_max, z_max;
    math::Point current_pnt = it_nodes.value().getPoint();
    x_min = current_pnt.X();
    x_max = current_pnt.X();
    y_min = current_pnt.Y();
    y_max = current_pnt.Y();
    z_min = current_pnt.Z();
    z_max = current_pnt.Z();
    it_nodes.next();
    while(!it_nodes.isDone()){
        math::Point current_pnt = it_nodes.value().getPoint();
        
        if(current_pnt.X()<x_min)
            x_min = current_pnt.X();
        else  if(current_pnt.X()>x_max)
            x_max = current_pnt.X();
        
        if(current_pnt.Y()<y_min)
            y_min = current_pnt.Y();
        else  if(current_pnt.Y()>y_max)
            y_max = current_pnt.Y();
        
        if(current_pnt.Z()<z_min)
            z_min = current_pnt.Z();
        else  if(current_pnt.Z()>z_max)
            z_max = current_pnt.Z();
        
        it_nodes.next();
    }
    math::Point p_min(x_min, y_min, z_min);
    math::Point p_max(x_max, y_max, z_max);
    m_mesh_radius = p_min.distance(p_max);
    
    m_confusing_distance = ATolerance*m_mesh_radius;
}
/*----------------------------------------------------------------------------*/
void SingularityGraphBuilder2D::execute()
{
    //==================================================================
    // Boolean marks initialization
    //==================================================================
    m_mark_faces_with_sing_point = m_mesh->getNewMark<Face>();
    m_mark_faces_with_sing_line  = m_mesh->getNewMark<Face>();
    
    //==================================================================
    // GEOMETRY VARIABLE FOR DEBUG
    //==================================================================
    Variable<int>* geom_var = m_mesh->newVariable<int>(GMDS_NODE, "geometry");
    
    m_index = m_mesh->newVariable<int>(GMDS_FACE, "index");
    
    IGMesh::node_iterator it_nodes = m_mesh->nodes_begin();
    for (; !it_nodes.isDone(); it_nodes.next())
    {
        Node n = it_nodes.value();
        if (m_mesh->isMarked(n, m_mark_nodes_on_point))
            (*geom_var)[n.getID()] = 2;
        else if (m_mesh->isMarked(n,m_mark_nodes_on_curve))
            (*geom_var)[n.getID()] = 1;
        else
            (*geom_var)[n.getID()] = 0;
    }
    VTKWriter<IGMesh> writer(*m_mesh);
    writer.write("geom_mesh", DIM3 | F | N | R);
    std::cout << "\t DONE" << std::endl;
    
    
    //========================================================================
    // STEP 1 - Detection of singular tetrahedra and storage
    //========================================================================
    std::cout << "============================================================="
    << std::endl;
    std::cout << "Detection of singular triangles" << std::endl;
    detectSingularTriangles();
    writer.write("index_mesh", DIM3 | F | N | R);
    std::cout << "\t DONE" << std::endl;
    
    //========================================================================
    // STEP 2 - Creation of singularity points and slots
    //========================================================================
    // For 3-valent singular points ...
    std::list<gmds::Face>::iterator sing_it;
    
    for(sing_it=m_singularities_3.begin();sing_it!=m_singularities_3.end();
        sing_it++) {
        Face current_face = *sing_it;
        createSingPointAndSlots(current_face);
    }
    // ... and 5-valent singular points
    for(sing_it=m_singularities_5.begin();sing_it!=m_singularities_5.end();
        sing_it++){
        Face current_face = *sing_it;
        createSingPointAndSlots(current_face);
    }
    //========================================================================
    // STEP 3 - Geometrical features are added in the singularity graph
    //========================================================================
    std::cout << "Addition of geometrical features" << std::endl;
    addGeometryToSingularityGraph(m_build_geometric_singularities);
    std::cout << "\t DONE" << std::endl;
    writeOutputSingle("with_geom");
    
    writeSingularityPointsAndSlots();
    
    //========================================================================
    // STEP 4 - Initialization of confusing balls around singularity points
    //========================================================================
    std::vector<SingularityPoint*> singularity_points = m_graph.getPoints();
    for (unsigned int i = 0; i < singularity_points.size(); i++) {
        SingularityPoint* pi = singularity_points[i];
        initConfusingBalls(pi);
    }
    //variable for debug purpose
    //DEBUG - DEBUG - DEBUG - DEBUG - DEBUG
    Variable<int>* ball_var =
    m_mesh->newVariable<int>(GMDS_FACE, "sing_ball");
    
    IGMesh::face_iterator it_faces = m_mesh->faces_begin();
    for (; !it_faces.isDone(); it_faces.next()) {
        Face f = it_faces.value();
        SingularityPoint* sing = m_faces_to_singularity_on_surf[f.getID()];
        if (sing==0)
            (*ball_var)[f.getID()] = 0;
        else
            (*ball_var)[f.getID()] = sing->index();
    }
    VTKWriter<IGMesh> writerB(*m_mesh);
    writerB.write("confusing_balls", DIM3 | F | N);
    //DEBUG - DEBUG - DEBUG - DEBUG - DEBUG
    //========================================================================
    // STEP 5 - Singularity line building
    //========================================================================
    createSingularityLines();
    //========================================================================
    // STEP 6 - Detect singularity lines intersection and split them
    //========================================================================
    detectLineIntersections();
    //========================================================================
    // STEP 7 - Build surface patchs
    //========================================================================
    m_graph.buildSurfacePatchs();
    writeOutput("patchs");
    std::vector<SingularityPatch* >  patchs = m_graph.getSurfacePatchs();
    std::cout<<"Nb Built patchs: "<<patchs.size()<<std::endl;
    //========================================================================
    // Boolean marks cleaning
    //========================================================================
    m_mesh->unmarkAll<Face>(m_mark_faces_with_sing_point);
    m_mesh->unmarkAll<Face>(m_mark_faces_with_sing_line);
    m_mesh->freeMark<Face>(m_mark_faces_with_sing_point);
    m_mesh->freeMark<Face>(m_mark_faces_with_sing_line);
    
}
/*---------------------------------------------------------------------------*/
void SingularityGraphBuilder2D::detectLineIntersections()
{
    std::vector<math::Point> added_points;
    IGMesh::face_iterator it_faces = m_mesh->faces_begin();
    std::vector<Face> candidate_faces;
    std::cout<<"======= SURFACE LINE INTERSECTIONS ================"<<std::endl;
    //========================================================================
    // We detect interesting faces before splitting any curves. After, it will
    // too late to detecte them.
    //=======================================================================
    for (; !it_faces.isDone(); it_faces.next()) {
        Face f = it_faces.value();
        
        std::vector<SurfaceSingularityLine*> current_lines = getSingularityLinesIn(f);
        if (current_lines.size()>1){
            candidate_faces.push_back(f);
        }
    } //for (; !it_faces.isDone(); it_faces.next())
    
    for(unsigned int i=0;i<candidate_faces.size();i++){
        Face f = candidate_faces[i];
        
        std::vector<SurfaceSingularityLine*> current_lines = getSingularityLinesIn(f);
        if(current_lines.size()==2)
            createLineIntersection(current_lines[0],current_lines[1],f, added_points);
        else
            createLineIntersection(current_lines,f, added_points);
    }//for(unsigned int i=0;i<candidate_faces.size();i++)
    
    std::cout<<"NB NEW POINTS: "<<added_points.size()<<std::endl;
}
/*---------------------------------------------------------------------------*/
void SingularityGraphBuilder2D::
createLineIntersection(SurfaceSingularityLine *ALine1,
                       SurfaceSingularityLine *ALine2,
                       Face& AFace,
                       std::vector<gmds::math::Point>& AAddedPoints)
{
    // We look for the geometrical intersection point of the curve lines
    SurfaceSingularityLine *l0 = ALine1;
    SurfaceSingularityLine *l1 = ALine2;
    math::Point p;
    
    math::Point face_center = AFace.center();
    std::vector<Node> face_nodes = AFace.get<Node>();
    double face_radius =0;
    int nb_face_nodes = face_nodes.size();
    for(unsigned int i=0; i < nb_face_nodes; i++) {
        math::Point pi = face_nodes[i].getPoint();
        math::Point pj = face_nodes[(i+1)%nb_face_nodes].getPoint();
        double distance_ij = pi.distance(pj);
        if(distance_ij>face_radius)
            face_radius = distance_ij;
    }
    
    if(l0->getIntersectionPoint(l1,face_center,face_radius,p)) {
        bool already_added = false;
        for(unsigned int i=0;!already_added && i<AAddedPoints.size();i++){
            if(math::near(p.distance(AAddedPoints[i]),0.0))
                already_added = true;
        }
        if(!already_added){
            AAddedPoints.push_back(p);
            //==============================================================
            // Creation of the singularity point
            //==============================================================*
            SurfaceSingularityPoint* new_pnt = m_graph.newSurfacePoint();
            new_pnt->setLocation(p);
            new_pnt->addMeshFace(AFace);
            //==================================================================
            // Creation of a new sing point
            // And splitting of intersected sing. lines
            //==================================================================
            m_graph.splitSurfaceLine(new_pnt,l0);
            m_graph.splitSurfaceLine(new_pnt,l1);
            writeOutput("boundary_line");
        }
    }//for(unsigned int i=0;i<candidate_faces.size();i++)
    
}
/*---------------------------------------------------------------------------*/
void SingularityGraphBuilder2D::
createLineIntersection(std::vector<SurfaceSingularityLine *>& ALines,
                       Face& AFace,
                       std::vector<gmds::math::Point>& AAddedPoints)
{
    // We look for the geometrical intersection point of the curve lines
    math::Point p;
    
    math::Point face_center = AFace.center();
    std::vector<Node> face_nodes = AFace.get<Node>();
    double face_radius =0;
    int nb_face_nodes = face_nodes.size();
    for(unsigned int i=0; i < nb_face_nodes; i++) {
        math::Point pi = face_nodes[i].getPoint();
        math::Point pj = face_nodes[(i+1)%nb_face_nodes].getPoint();
        double distance_ij = pi.distance(pj);
        if(distance_ij>face_radius)
            face_radius = distance_ij;
    }
    // WARNING THIS DOUBLE LOOP DOES NOT ENSURE TO GET ALL THE POSSIBLE
    // INTERSECTIONS!!!! THIS SHOULD BE IMPROVED.
    for(unsigned int i=0;i<ALines.size()-1; i++){
        
        SurfaceSingularityLine *li = ALines[i];
        
        for(unsigned int j=i+1; j<ALines.size(); j++){
            
            SurfaceSingularityLine *lj = ALines[j];
            
            if(li->getIntersectionPoint(lj,face_center,face_radius,p)) {
                bool already_added = false;
                for(unsigned int i=0;!already_added && i<AAddedPoints.size();i++){
                    if(math::near(p.distance(AAddedPoints[i]),0.0))
                        already_added = true;
                }
                if(!already_added){
                    AAddedPoints.push_back(p);
                    //==============================================================
                    // Creation of the singularity point
                    //==============================================================*
                    SurfaceSingularityPoint* new_pnt = m_graph.newSurfacePoint();
                    new_pnt->setLocation(p);
                    new_pnt->addMeshFace(AFace);
                    //==================================================================
                    // Creation of a new sing point
                    // And splitting of intersected sing. lines
                    //==================================================================
                    m_graph.splitSurfaceLine(new_pnt,li);
                    m_graph.splitSurfaceLine(new_pnt,lj);  
                    writeOutput("boundary_line");
                }
            }
            
        }
    }
    
    
}
/*---------------------------------------------------------------------------*/
void SingularityGraphBuilder2D::detectSingularTriangles()
{
    IGMesh::face_iterator it_faces = m_mesh->faces_begin();
    m_singularities_3.clear();
    m_singularities_5.clear();
    
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
        
        if (index ==1) {
            m_singularities_5.push_back(current);
            m_mesh->mark(current, m_mark_faces_with_sing_point);
        }
        else if (index == -1) {
            m_singularities_3.push_back(current);
            m_mesh->mark(current, m_mark_faces_with_sing_point);
            
        }
        
    } //for (; !it_regions.isDone(); it_regions.next())
    
    std::cout << "Nb 3-singular simplices: " << m_singularities_3.size() << std::endl;
    std::cout << "Nb 5-singular simplices: " << m_singularities_5.size() << std::endl;
}
/*---------------------------------------------------------------------------*/
void SingularityGraphBuilder2D::initMarks(const int AMarkNodePnt,
                                          const int AMarkNodeCrv, 
                                          const int AMarkEdgeCrv)
{
    m_mark_nodes_on_point = AMarkNodePnt;
  m_mark_nodes_on_curve = AMarkNodeCrv;
  m_mark_edges_on_curve = AMarkEdgeCrv;
}
/*----------------------------------------------------------------------------*/
void SingularityGraphBuilder2D::
initConfusingBalls(SingularityPoint* APnt)
{
  math::Point sing_location = APnt->getLocation();
  //========================================================================
  // Faces
  //========================================================================
  int nb_faces = 0;
  IGMesh::face_iterator it_faces = m_mesh->faces_begin();
  for (; !it_faces.isDone(); it_faces.next()) {
      Face current_face = it_faces.value();
      math::Point center = current_face.center();
      if(center.distance(sing_location)< m_confusing_distance) { 
	m_faces_to_singularity_on_surf[current_face.getID()] = APnt;
	nb_faces++;
	}
  }

  //========================================================================
  // Edges
  //========================================================================
  IGMesh::edge_iterator it_edges = m_mesh->edges_begin();
  for (; !it_edges.isDone(); it_edges.next()) {
      Edge current_edge = it_edges.value();
      math::Point center = current_edge.center();
      if(center.distance(sing_location)< m_confusing_distance) { 
	m_edges_to_singularity_on_surf[current_edge.getID()] = APnt;
	}
  }
  //========================================================================
  // Nodes
  //========================================================================
  IGMesh::node_iterator it_nodes = m_mesh->nodes_begin();
  for (; !it_nodes.isDone(); it_nodes.next()) {
      Node current_node = it_nodes.value();
      math::Point center = current_node.getPoint();
      if(center.distance(sing_location)< m_confusing_distance) { 
	m_nodes_to_singularity_on_surf[current_node.getID()] = APnt;
	}
  }
}
/*----------------------------------------------------------------------------*/
void SingularityGraphBuilder2D::computeSingPointInfo(Face& AFace,
						     math::Point& APosSing)
{
    std::cout << "== Singularity point info for face "<<AFace.getID()<<" =="<<std::endl;

  std::vector<Node> current_nodes = AFace.get<Node>();
  Node node1 = current_nodes[0];
  Node node2 = current_nodes[1];
  Node node3 = current_nodes[2];

  math::Vector v1 =(*m_field)[node1.getID()].referenceVector();
  math::Vector v2 =(*m_field)[node2.getID()].referenceVector();
  math::Vector v3 =(*m_field)[node3.getID()].referenceVector();
 
  //We solve a 2x2 Ax=b system with x = (alpha,beta)
  double alpha=0, beta=0, gamma=0;
  math::Vector A = v1-v3;
  math::Vector B = v2-v3;
  math::Vector C = v3.opp();
  
  double dA = A[0]*B[1]-A[1]*C[0];
  if(dA!=0){
    double Dx = C[0]*B[1]-C[1]*B[0];
    double Dy = A[0]*C[1]-A[1]*C[0];
    alpha = Dx/dA;
    beta  = Dy/dA;
    gamma = 1-alpha-beta;
  }
  else{
    throw GMDSException("Null Determinant in the computation of a singularity point location");
  }
  if(true){//gamma<0){ 
    alpha = 0.333;
    beta  = 0.333;
    gamma = 1-alpha-beta;
  }
  std::cout<<"(a,b,g) = "<<alpha<<", "<<beta<<", "<<gamma<<std::endl;
  //(alpha, beta, gamma) are the barycentric coordinates where the cross field vanishes
  //We can then compute the singularity point location
  double x = alpha * node1.X() + beta * node2.X() + gamma * node3.X();
  double y = alpha * node1.Y() + beta * node2.Y() + gamma * node3.Y();

  APosSing.setXYZ(x,y,0);

}
/*----------------------------------------------------------------------------*/
void SingularityGraphBuilder2D::createSingPointAndSlots(Face& AFace) 
{
  //=========================================================================
  // WE BEGIN BY COMPUTING THE LOCATION OF THE SINGULARITY POINT IN AFace
  //=========================================================================
  math::Point s;
  computeSingPointInfo(AFace, s);  
  double xs = s.X();
  double ys = s.Y();

  //=========================================================================
  // NOW, WE COMPUTE THE INTERSECTIONS POINT ALONG EACH EDGE
  //========================================================================= 
  //for each detected slot, we store its location, its direction, its out cell
  // id, and the dimension of this out cell.

  std::vector<math::Point>  slot_points; 
  std::vector<double>       slot_param; 
  std::vector<math::Vector> slot_dirs; 
  std::vector<int>          slot_cell_dim;
  std::vector<TCellID>      slot_cell_id;

  std::vector<Edge> edges = AFace.get<Edge>();
  for(int i=0;i<3;i++){ //we walk along each edge of AFace
    Edge ei = edges[i];
    std::vector<Node> nodes_ei = ei.get<Node>();
    Node ni = nodes_ei[0];
    Node nj = nodes_ei[1];
    std::vector<math::Vector> vectors_i = (*m_field)[ni.getID()].componentVectors();
    math::Cross2D cross_j = (*m_field)[nj.getID()];

    math::Point pi = ni.getPoint();
    math::Point pj = nj.getPoint();

    double xi = pi.X();
    double yi = pi.Y();
    double xj = pj.X();
    double yj = pj.Y();

    double xij = xi - xj;
    double yij = yi - yj;

    double xjs = xj - xs;
    double yjs = yj - ys;

    for(int k=0;k<2;k++){
      math::Vector cik =  vectors_i[k];
      math::Vector cjk =  cross_j.closestComponentVector(cik);

      double x_cik = cik.X();
      double y_cik = cik.Y();
      double x_cjk = cjk.X();
      double y_cjk = cjk.Y();

      double x_cijk = x_cik - x_cjk;
      double y_cijk = y_cik - y_cjk;

      double a = (xij* y_cijk) - (yij   * x_cijk);
      double b = (y_cjk * xij) + (xjs   * y_cijk) - (x_cjk * yij) - (yjs * x_cijk);
      double c = (y_cjk * xjs) - (x_cjk * yjs);

      std::vector<double> solutions;
      math::solve2ndDegreePolynomial(a,b,c,solutions);
      std::cout<<"2nd degree eq: "<<solutions.size()<<std::endl;

      for(unsigned int i_sol=0; i_sol<solutions.size();i_sol++){

	  double alpha = solutions[i_sol];
	  if(alpha>1 || alpha<0.0)
	    continue;

	  std::cout<<"\t "<<alpha<<std::endl;
	  math::Point p = alpha*pi + (1-alpha)*pj;
	  slot_param.push_back(alpha);
	  slot_points.push_back(p);
	  slot_dirs.push_back(math::Vector(s,p));
	  if(alpha==0){ //we go out from a node
	    slot_cell_dim.push_back(0);
	    slot_cell_id.push_back(ni.getID());
	  } 
	  else if(alpha==0){ //we go out from a node
	    slot_cell_dim.push_back(0);
	    slot_cell_id.push_back(nj.getID());
	  }
	  else { // general case, the edge
	    slot_cell_dim.push_back(1);
	    slot_cell_id.push_back(ei.getID());
	  }
      }

    }//for(int k=0;k<2;k++)

  }// for(int i=0;i<3;i++)
  std::cout<<"Nb slots = "<< slot_points.size()<<std::endl;

  // SINGULARITY POINT CREATION 
  SurfaceSingularityPoint*  singularity = m_graph.newSurfacePoint();
  singularity->setLocation(s); 
  singularity->addMeshFace(AFace);
  for(unsigned int i=0;i<slot_points.size(); i++){
    singularity->newSlot(slot_points[i],
			 slot_dirs[i],
			 slot_cell_id[i],
			 slot_cell_dim[i], 
			 true, 
			 0);
  }

}
/*----------------------------------------------------------------------------*/
void SingularityGraphBuilder2D::
addGeometryToSingularityGraph(const bool ABuildGeomSlots)
{
  //Now we add singularty points and lines for the corner and edges of the geometry
  IGMesh::node_iterator it_node = m_mesh->nodes_begin();
  //============================================================================
  //	GEOMETRIC  SINGULARITY POINTS
  //============================================================================
  for (; !it_node.isDone(); it_node.next()) {
      Node current_node = it_node.value();
      if (m_mesh->isMarked(current_node,  m_mark_nodes_on_point))
	{
	  VertexSingularityPoint* sing_point = m_graph.newVertexPoint();
	  sing_point->setXYZ(current_node.X(), current_node.Y(), current_node.Z());
	  sing_point->addMeshNode(current_node);
	}
    }

  //============================================================================
  //	GEOMETRIC  SINGULARITY LINES
  //============================================================================
  //Now we have all the corner of the geom as singularity
  //We will create separatrices based on the geometric edges

  int mark_geom_edges = m_mesh->getNewMark<Node>();
  std::vector<SingularityLine*> added_geom_lines;

  it_node = m_mesh->nodes_begin();
  for (; !it_node.isDone(); it_node.next()){
    Node current_node = it_node.value();
    if ((!m_mesh->isMarked(current_node, m_mark_nodes_on_point)) &&
	(m_mesh->isMarked(current_node, m_mark_nodes_on_curve)) &&
	(!m_mesh->isMarked(current_node, mark_geom_edges)))  {
      //new singularity line to create here
      std::vector<int> listOfNodesInSingLeft;
      std::vector<int> listOfNodesInSingRight;
      m_mesh->mark(current_node, mark_geom_edges);
      Node NodeLeft;
      Node NodeRight;
      int gotFirstOne = 0;
      for (unsigned int i = 0; i < current_node.get<Edge>().size(); i++) {
	if (m_mesh->isMarked(current_node.get<Edge>()[i], m_mark_edges_on_curve)) {
	  for (unsigned int j = 0; j < 2; j++){
	    if (current_node.get<Edge>()[i].get<Node>()[j].getID() != current_node.getID()){
	      if (!gotFirstOne){
		gotFirstOne = 1;
		NodeLeft = current_node.get<Edge>()[i].get<Node>()[j];
	      }
	      else{
		NodeRight = current_node.get<Edge>()[i].get<Node>()[j];
	      }
	    }
	  }
	}
      }
      Node Ncurr = current_node;
 //     int numberOfTheFirstSing = 0;
 //     int numberOfTheSecondSing = 0;
      //Here we have the two nodes left and right that continue the line
      //		std::cout << "entree dans premiere boucle while" << std::endl;
      //------------------------------------------------------------------------------------
      //la premiere condition de la boucle while empeche la fermeture des cercles
      while ((NodeLeft.getID() != NodeRight.getID()) &&
	     (!m_mesh->isMarked(NodeLeft, m_mark_nodes_on_point))) {
	//	std::cout << "first boucle: on est dans le node " << NodeLeft.getID() << std::endl;
	m_mesh->mark(NodeLeft, mark_geom_edges);
	listOfNodesInSingLeft.push_back(NodeLeft.getID());
	//	if (m_mesh->isMarked(NodeLeft, m_mark_nodes_on_point)){
	//		std::cout << "this is a corner" << std::endl;
	//	}
	/*	if (m_mesh->isMarked(NodeLeft, m_mark_nodes_on_curve)){
		std::cout << "this is an edge" << std::endl;
		}
		std::cout << "on a un nombre de edges voisines de "
		<< NodeLeft.get<Edge>().size() << std::endl;*/
	Node NodeNext = NodeLeft;
	for (unsigned int i = 0; i < NodeLeft.get<Edge>().size(); i++){
	  if (m_mesh->isMarked(NodeLeft.get<Edge>()[i], m_mark_edges_on_curve)){
	    //std::cout << "start of edge" << std::endl;
	    for (unsigned int j = 0; j < NodeLeft.get<Edge>()[i].get<Node>().size(); j++){


	      if (NodeLeft.get<Edge>()[i].get<Node>()[j].getID() != NodeLeft.getID()){
		if (NodeLeft.get<Edge>()[i].get<Node>()[j].getID() != Ncurr.getID()){
		  NodeNext = NodeLeft.get<Edge>()[i].get<Node>()[j];
		}
	      }
	    }
	    //		std::cout << "out of edge" << std::endl;
	  }
	}
	Ncurr = NodeLeft;
	NodeLeft = NodeNext;
      }
      //	std::cout << "sortie de premiere boucle while" << std::endl;
      //------------------------------------------------------------------------------------
      m_mesh->mark(NodeLeft, mark_geom_edges);
      listOfNodesInSingLeft.push_back(NodeLeft.getID());
      if (NodeLeft.getID() != NodeRight.getID()){
	//we need to get on the right too
	Node Ncurr = current_node;
	//	std::cout << "entree dans seconde boucle while" << std::endl;
	while (!m_mesh->isMarked(NodeRight, m_mark_nodes_on_point)) {
	  //	std::cout << "seconde boucle: on est dans le node " << NodeRight.getID() << std::endl;
	  m_mesh->mark(NodeRight, mark_geom_edges);
	  listOfNodesInSingRight.push_back(NodeRight.getID());
	  //if (m_mesh->isMarked(NodeRight, m_mark_nodes_on_point))
	  //{
	  //	std::cout << "this is a corner" << std::endl;
	  //}
	  Node NodeNext = NodeRight;
	  for (unsigned int i = 0; i < NodeRight.get<Edge>().size(); i++){
	    if (m_mesh->isMarked(NodeRight.get<Edge>()[i], m_mark_edges_on_curve)){
	      for (unsigned int j = 0; j < 2; j++){
		//std::cout << "on test le node " << NodeRight.get<Edge>()[i].get<Node>()[j].getID() << std::endl;
		if (NodeRight.get<Edge>()[i].get<Node>()[j].getID() != NodeRight.getID()){
		  if (NodeRight.get<Edge>()[i].get<Node>()[j].getID() != Ncurr.getID()){
		    NodeNext = NodeRight.get<Edge>()[i].get<Node>()[j];
		  }
		}
	      }
	    }
	  }
	  Ncurr = NodeRight;
	  NodeRight = NodeNext;
	}
	//	std::cout << "sortie de seconde boucle while" << std::endl;
	m_mesh->mark(NodeRight, mark_geom_edges);
	listOfNodesInSingRight.push_back(NodeRight.getID());
      }
      else{ //cycle
	listOfNodesInSingRight.push_back(NodeLeft.getID());
      }
      //Here we have the list of points to insert in the separatrix
      CurveSingularityLine* new_line = m_graph.newCurveLine();
      //we give the line id
      new_line->setNumber(m_graph.getNbLines());

      std::vector<Node> curve_nodes;
      for (unsigned int i = 0; i < listOfNodesInSingLeft.size(); i++) {
	Node NAtThisPoint = m_mesh->get<Node>(listOfNodesInSingLeft[listOfNodesInSingLeft.size() - 1 - i]);
	new_line->addDiscretizationPoint(NAtThisPoint.getPoint());
	curve_nodes.push_back(NAtThisPoint);
      }

      new_line->addDiscretizationPoint(current_node.getPoint());
      curve_nodes.push_back(current_node);

      for (unsigned int i = 0; i < listOfNodesInSingRight.size(); i++)
	{
	  Node NAtThisPoint = m_mesh->get<Node>(listOfNodesInSingRight[i]);
	  new_line->addDiscretizationPoint(NAtThisPoint.getPoint());
	  curve_nodes.push_back(NAtThisPoint);
	}

      //now we add the edges from the nodes in the ordered way
      std::vector<Edge> curve_edges;
      for (unsigned int i = 0; i < curve_nodes.size() - 1; i++)
	{
	  Node current = curve_nodes[i];
	  Node next = curve_nodes[i + 1];
	  std::vector<Edge> current_edges = current.get<Edge>();

	  bool found_edge = false;
	  for (unsigned int j = 0; j < current_edges.size() && !found_edge; j++)
	    {
	      Edge ej = current_edges[j];
	      if (m_mesh->isMarked(ej, m_mark_edges_on_curve)){
		std::vector<Node> ej_nodes = ej.get<Node>();
		if (ej_nodes[0].getID() == next.getID() ||
		    ej_nodes[1].getID() == next.getID())
		  {
		    curve_edges.push_back(ej);
		    found_edge = true;
		  }
	      }
	    }
	}
      //TODO Attention aux courbes cycliques, a priori non traitees
      new_line->setMeshEdges(curve_edges);
      added_geom_lines.push_back(new_line);
    }
  }

  // Geometrical curves made of only one mesh edges are missing. We add them now.
  IGMesh::edge_iterator itEdgeGeo = m_mesh->edges_begin();

  for (; !itEdgeGeo.isDone(); itEdgeGeo.next()){
    Edge currentEdge = itEdgeGeo.value();
    // on ne traite que les aretes sur une courbe geometrique
    if (m_mesh->isMarked(currentEdge, m_mark_edges_on_curve)){
      std::vector<Node> currentNodes = currentEdge.get<Node>();
      Node n1 = currentNodes[0];
      Node n2 = currentNodes[1];
      // on regarde si les deux noeuds correspondent ˆ des sommets geometriques
      if (m_mesh->isMarked(n1, m_mark_nodes_on_point) && m_mesh->isMarked(n2, m_mark_nodes_on_point))
	{ // On cree donc la separatrice reliant les singularites associŽes ˆ n1 et n2

	  SingularityLine* new_line = m_graph.newCurveLine();
	  new_line->setNumber(m_graph.getNbLines());

	  new_line->addDiscretizationPoint(n1.getPoint());
	  new_line->addDiscretizationPoint(n2.getPoint());
	  added_geom_lines.push_back(new_line);

	}
    }
  }
  //============================================================================
  // GEOM SINGULARITY POINTS AND GEOM SINGULARITY LINES MUST BE CONNECTED
  //============================================================================
  std::vector<VertexSingularityPoint*> geom_points = m_graph.getVertexPoints();
  for (unsigned int i = 0; i < geom_points.size(); i++)
    {
      VertexSingularityPoint* pi = geom_points[i];
      gmds::math::Point pi_point = pi->getLocation();
      for (unsigned int j = 0; j < added_geom_lines.size(); j++)
	{
	  SingularityLine* lj = added_geom_lines[j];
	  std::vector<SingularityPoint*>& lj_points = lj->getEndPoints();
	  if (lj_points.size() == 2)
	    continue; //the line is connected to its both end points

	  std::vector<gmds::math::Point >&
	    lj_discretization = lj->getDiscretizationPoints();

	  gmds::math::Point  lj_end_loc[2];
	  gmds::math::Vector lj_end_dir[2];
	  lj_end_loc[0] = lj_discretization[0];
	  lj_end_loc[1] = lj_discretization[lj_discretization.size() - 1];
	  lj_end_dir[0] = gmds::math::Vector(lj_end_loc[0], lj_discretization[1]);
	  lj_end_dir[1] = gmds::math::Vector(lj_end_loc[1], 
					     lj_discretization[lj_discretization.size() - 2]);
	  for (int k = 0; k < 2; k++){
	    gmds::math::Point pk = lj_end_loc[k];
	    gmds::math::Vector v(pk, pi_point);
	    if (v.norm2() < 1e-5){
	      //pi and lj must be connected
	      pi->newSlot(lj_end_loc[k], lj_end_dir[k], 
			  pi->getMeshNode().getID() /*starting cell id*/,
			  0 /*starting cell dim*/,
			  true /*on surface*/, lj);
	      lj->addSingularityPoint(pi);
	    }
	  }
	}
    }
  //============================================================================
  // CYCLE GEOMETRIC CURVE WILL BE NOT CONNECTED TO GEOM POINT. WE CREATE SUCH 
  // POINT IN AN ARTIFICAL WAY
  //============================================================================

  std::vector<CurveSingularityLine*> geom_curves = m_graph.getCurveLines();
  for (unsigned int i = 0; i < geom_curves.size(); i++){
    CurveSingularityLine* ci =  geom_curves[i];
    if (ci->getEndPoints().empty()){ 
      gmds::math::Point loc = ci->getDiscretizationPoints()[0];
      CurveSingularityPoint*	p= m_graph.newCurvePoint();
      p->setLocation(loc);

      std::vector<math::Point > line_disc =  ci->getDiscretizationPoints();
      math::Point p_begin = line_disc[0];
      math::Point p_end   = line_disc[line_disc.size()-1];

      math::Vector slot_dir1 = 	math::Vector(p_begin, line_disc[1]);
      math::Vector slot_dir2 = 	math::Vector(p_end, line_disc[line_disc.size()-2]);
      
      SingularityPoint::Slot* s = p->newSlot(p->getLocation(), // slot location
					     slot_dir1, // slot direction
					     0,    // No linked cell (id)
					     0,    // No linked cell (dim)
					     true, // Always on surface (Maybe false in the future)
					     ci,// Connected line
					     slot_dir1.opp());// Line direction is the slot direction
      s->isFreeze = true;
      s =p->newSlot(p->getLocation(), // slot location
		    slot_dir2, // slot direction
		    0,    // No linked cell (id)
		    0,    // No linked cell (dim)
		    true, // Always on surface (Maybe false in the future)
		    ci,// Connected line
		    slot_dir2.opp());// Line direction is the slot direction
      s->isFreeze = true;
      ci->addSingularityPoint(p);
      ci->addSingularityPoint(p);
    }
  }
  m_mesh->unmarkAll<Node>(mark_geom_edges);
  m_mesh->freeMark<Node>(mark_geom_edges);


  //================================================================================
  //	FINALLY, WE TRAVERSE GEOM SING POINT TO DEFINE THEIR SLOTS AT NON-CONVEX AREAS
  //================================================================================
  if( ABuildGeomSlots) {
    std::vector<VertexSingularityPoint* >  vertex_points = m_graph.getVertexPoints();
    for(unsigned int i=0;i<vertex_points.size();i++) {
      VertexSingularityPoint* current_point = vertex_points[i];
      Node current_node = current_point->getMeshNode();
      math::Point in_pnt = current_node.getPoint();
      std::cout<<"============================"<<std::endl
	       <<"Current node "<<current_node.getID()<<" at "
	       <<current_node.getPoint()<<std::endl;

      //Volume singularity lines must be added for non-convex surface points
      std::vector<Edge> current_edges = current_node.get<Edge>();
      std::vector<Face> current_faces = current_node.get<Face>();
      for (unsigned int j = 0; j < current_edges.size(); j++){
	Edge ej = current_edges[j];

	//only edges classified on geometrical curves are taken into account
	if (!m_mesh->isMarked(ej, m_mark_edges_on_curve))
	  continue;

	std::vector<Node> ej_nodes = ej.get<Node>();
	Node other_node;
	if (ej_nodes[0] == current_node)
	  other_node = ej_nodes[1];
	else
	  other_node = ej_nodes[0];
	std::cout<<"edge "<<j<<" from "<< other_node.getID()<<" to "
		  <<current_node.getID()<<std::endl;
       	math::Vector vec(other_node.getPoint(), current_node.getPoint());
	//	math::Vector vec(current_node.getPoint(), other_node.getPoint());
	vec.normalize();
	
	// Now for each face adjacent to current_node, we look for one intersected
	// by vec
	bool found = false; //indique si on a trouve une intersection
	for (unsigned int k = 0; k < current_faces.size() && !found; k++){
	  Face fk = current_faces[k];

	  Node node_from = current_node;  
	  //=====================================================================
	  // We look for the opposite edge  
	  //=====================================================================
	  Edge opposite_edge;
	  std::vector<Edge> current_edges = fk.get<Edge>();
	  for (unsigned int i = 0; i < current_edges.size(); i++) {
	    Edge ei = current_edges[i];
	    std::vector<Node> ei_nodes = ei.get<Node>();
	    if (node_from != ei_nodes[0] && node_from != ei_nodes[1])
	      opposite_edge = ei;
	  }


	  //=====================================================================
	  // We compute an out pnt and vector
	  //=====================================================================
	  math::Point  out_pnt;
	  math::Vector out_vec;
	  int out_dim =0;
	  TCellID out_id = NullID;
	  math::Vector in_vec = vec;
	  bool found_out = false;
	  std::vector<Node> other_nodes = opposite_edge.get<Node>();
	  //===========================================================
	  // Go through the first opposite node if it not on a curve ?
	  //===========================================================
	  if(!m_mesh->isMarked(other_nodes[0],m_mark_nodes_on_curve)) {
	    math::Point opp_node_loc1 = other_nodes[0].getPoint();

	    math::Vector v_opp1(in_pnt, opp_node_loc1);
	    v_opp1.normalize();
	      
	    if (math::near(v_opp1.dot(in_vec) - 1,0)) {
	      found_out = true;
	      out_pnt = opp_node_loc1;    
	      m_tool.computeOutVectorAtPoint(other_nodes[0], in_vec, out_vec);
	      out_dim = 0;
	      out_id  = other_nodes[0].getID(); 
	    }
	  }
	  if(!found_out && 
	     !m_mesh->isMarked(other_nodes[1],m_mark_nodes_on_curve)) {
	    math::Point opp_node_loc2 = other_nodes[1].getPoint();

	    math::Vector v_opp2(in_pnt, opp_node_loc2);
	    v_opp2.normalize();
	      
	    if (math::near(v_opp2.dot(in_vec) - 1,0)) {
	      found_out = true;
	      out_pnt = opp_node_loc2;    
	      m_tool.computeOutVectorAtPoint(other_nodes[1], in_vec, out_vec);
	      out_dim = 0;
	      out_id  = other_nodes[1].getID(); 
	    }
	  }

	  //================================================
	  // Go through the opposite edge
	  // And not through one of its end points due to 
	  // previous tests.
	  //================================================ 
	  if(!found_out) {
	    found_out = m_tool.computeOutVectorFromRayAndEdge(opposite_edge,
						       in_pnt,
						       in_vec,
						       out_pnt,		
						       out_vec);
	    if(found_out){ 
	      out_dim = 1;
	      out_id  = opposite_edge.getID(); 
	    }
	  }

	  if(found_out){
	    std::cout<<"   -> out in ("<<out_dim<<", "<<out_id<<")"<<std::endl;
	    current_point->newSlot(out_pnt, out_vec, out_id, out_dim, false);
	  }


	}// for (unsigned int k = 0; k < current_faces.size() && !found; k++)

      }//for (unsigned int j = 0; j < current_edges.size(); j++)


    }// for(unsigned int i=0;i<vertex_points.size();i++)

  }// if( ABuildGeomSlots) 
  
}
/*----------------------------------------------------------------------------*/
void SingularityGraphBuilder2D::createSingularityLines()
{
  //Now we need to create separatrices on the faces
  //we will start by launching them from the slot of the singularities
  //for each one we have the pos and dir to launch the sep into, as well as the first triangle
  std::cout << "== Singularity lines extraction ==" << std::endl;

  //========================================================================
  // At the beginning, only the singularity point of the cross field are
  // known. Some other points may be created on the fly, when we interset
  // boundary, singularity lines, etc.
  //========================================================================
  std::vector<SingularityPoint*> singularity_points = m_graph.getPoints();

  for (unsigned int i = 0; i < singularity_points.size(); i++) {
    SingularityPoint* pi = singularity_points[i];
    std::vector<SingularityPoint::Slot*> pi_slots = pi->getSlots();

    for (unsigned int j = 0; j < pi_slots.size(); j++) {
      if (!pi_slots[j]->isLaunched)
	m_free_slots.push_back(pi_slots[j]);
    }
    
  }
  //========================================================================
  // Creation of singularity lines from cross field singular points
  //========================================================================
  while(!m_free_slots.empty()){
    SingularityPoint::Slot* current_slot = m_free_slots.front();
    m_free_slots.pop_front();
    if(current_slot->isLaunched!=true) {
      std::cout<<"CURRENT NON LAUNCHED SLOT -- loc: ";
      std::cout<<current_slot->location<<", dir: "<<std::endl;
      std::cout<<current_slot->direction<<std::endl;
	      
      current_slot->isLaunched = true;
      computeSingularityLine(current_slot->from_point, current_slot);    
      writeOutput("boundary_line");
    }
  } // while(!m_free_slots.empy())
}
/*----------------------------------------------------------------------------*/
void SingularityGraphBuilder2D::
computeSingularityLine(SingularityPoint*       AFromPoint,
		       SingularityPoint::Slot* AFromSlot)
{
  std::cout << "======== Computation of a single singularity line ========" 
	    << std::endl;

  //========================================================================
  // Data initialization for line building
  //========================================================================
  bool end_on_bnd = false;
  bool end_on_free_slot = false;
  bool must_create_pnt = false;
  SingularityPoint       *to_sing_pnt  = 0;
  SingularityPoint::Slot *to_slot = 0;
  TCellID to_cell_id;
  int to_cell_dim;
  math::Point  to_pnt; 
  math::Vector to_dir; 
  std::vector<math::Point> line_discretization;
  std::vector<TCellID>     line_triangles;

  //========================================================================
  // Stream line computation
  //========================================================================
  computeStreamLine(AFromPoint,
		    AFromSlot,
		    to_sing_pnt,
		    to_slot,
		    to_pnt,
		    to_dir,
		    line_discretization,
		    line_triangles,
		    to_cell_dim,
		    to_cell_id,
		    end_on_bnd,
		    end_on_free_slot,
		    must_create_pnt);

  //========================================================================
  // Stream line analysis and creation
  //========================================================================

  //line creation
  SurfaceSingularityLine* surf_line = m_graph.newSurfaceLine();
  int sepNumberTmp = m_graph.getNbLines();
  surf_line->setNumber(sepNumberTmp);

  //connect line to initial singularity point
  SingularityPoint* from_sing_pnt = AFromSlot->from_point;
  surf_line->addSingularityPoint(from_sing_pnt);
  surf_line->addDiscretizationPoint(from_sing_pnt->getLocation());
  //and initial singularity point to the line
  AFromSlot->line = surf_line;
  AFromSlot->line_direction =  AFromSlot->direction;
  AFromSlot->isLaunched = true;

  //Insertion of line points
  for(unsigned int i=0; i<line_discretization.size(); i++) {
    surf_line->addDiscretizationPoint(line_discretization[i]);
  }

  for(unsigned int i=0; i<line_triangles.size(); i++) {
    surf_line->addTraversedFace(line_triangles[i]);
  }

  //========================================================================
  // Termination of the streamline
  //========================================================================
  if(end_on_bnd) {  
    
    std::cout<<"END ON BOUNDARY"<<std::endl;
    //======================================================================
    // CASE 1 - We finish on the boundary. A geometric point must be created
    // The cell defined by (start_cell_dim, start_cell_id) is on the boundary.
    // We have to create a geometric singularity point so.

    SingularityPoint* geom_pnt;
    SingularityPoint::Slot* incoming_slot;
    createGeometricSingularityPoint( to_pnt,      // the point we are
				     to_dir,      // the direction we come from
				     to_cell_dim, // the dim. of the cell 
				     to_cell_id,  // the id of the cell
				     geom_pnt,       // the created point
				     incoming_slot); // and the slot 

   

    surf_line->addSingularityPoint(geom_pnt);
    surf_line->addDiscretizationPoint(geom_pnt->getLocation());
    //    geom_pnt->connectLine(surf_line, to_dir); 

  } // if(to_cell_dim==0)
  else { 
    if(!end_on_free_slot){  
      //======================================================================
      // CASE 3 - We finish on a field singularity where a slot is not free
      /// The slot must be free before connection. 
      m_graph.removeLine(to_slot->line);
    }
    //======================================================================
    // CASE 2 - We finish on a field singularity where a slot is free NOW
    std::cout<<"END ON A FREE SLOT"<<std::endl;

    to_slot->isLaunched = true;
    to_slot->line = surf_line;
    to_slot->line_direction = to_dir;

    // We need to compute the discretization point, which are in the 
    // triangles of the confusing ball.
    surf_line->addSingularityPoint(to_sing_pnt);
    surf_line->addDiscretizationPoint(to_sing_pnt->getLocation());

    // We have the two singularities and the corresponding slots.
    // In this process, we interpolate the 
     backtrackSingularityLine(surf_line, // the line we modify 
     			     to_sing_pnt, // the point we start from
     			     to_slot,// the slot we start from
     			     AFromPoint,// the point we go to
     			     AFromSlot);  // the slot we go to
    

  } //else { 

}
/*----------------------------------------------------------------------------*/
void SingularityGraphBuilder2D::
computeStreamLine(SingularityPoint*         AFromPnt,
		  SingularityPoint::Slot*   AFromSlot,
		  SingularityPoint*&        AToSingPnt,
		  SingularityPoint::Slot*&  AToSlot,
		  math::Point&              AToPnt,
		  math::Vector&             AToDir,
		  std::vector<math::Point>& APoints, 
		  std::vector<TCellID>&     ATriangles,
		  int&                      AToCellDim,
		  TCellID&                  AToCellID,
		  bool&                     AEndOnBnd,
		  bool&                     AToSlotIsFree,
		  bool&                     APntToCreate)
{
  std::cout << "======== Streamline computation ========" 
	    << std::endl;

  ATriangles.clear();
  APoints.clear();

  math::Point  start_pnt = AFromSlot->location ; //starting point
  math::Vector start_dir = AFromSlot->direction; //starting direction
  math::Vector prev_dir  = AFromSlot->direction; //prev direction used in the
						 //termination of extrapolation
						 //process (when we get into a
						 //confusing ball)
  
  TCellID start_cell_id  = AFromSlot->starting_cell_id ;
  int     start_cell_dim = AFromSlot->starting_cell_dim;


  math::Point  current_pnt = start_pnt;
  math::Vector current_vec = start_dir;

  std::cout << "START PNT: " << start_pnt << std::endl;
  math::Point start_dirPnt(start_pnt.X() + start_dir.X(),
			   start_pnt.Y() + start_dir.Y(),
			   start_pnt.Z() + start_dir.Z());

  std::cout << "START DIR PNT: " << start_dirPnt << std::endl;
  //========================================================================
  // Singularity point we will be connecting to. It can be another 
  // singularity point of the cross field or a geometric singularity point
  // created on the fly.
  //========================================================================
  SingularityPoint*       found_pnt = 0;
  SingularityPoint::Slot* found_slot =0;

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

      m_tool.findNextCell(start_pnt, start_dir,
                          start_cell_dim, start_cell_id,
                          next_cell_dim, next_cell_id);

    /*    std::cout<<"NEXT CELL ("
	     <<next_cell_dim<<", "
	     <<next_cell_id<<")"<<std::endl;
    */
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
      //
      // WARNING, we do not check the confusing ball!!!!!

      Edge current_edge = m_mesh->get<Edge>(next_cell_id);
      std::vector<TCellID> adj_faces = current_edge.getIDs<Face>();
      ATriangles.insert(ATriangles.end(),adj_faces.begin(),adj_faces.end());
      
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
      m_tool.computeOutVectorAtPoint(next_node, start_dir, next_dir);

   
      // We assign the new value for the next step
      start_dir = next_dir;
      start_pnt = next_node.getPoint();
      start_cell_dim = 0;
      start_cell_id = next_node.getID();
      find_end = false;

      APoints.push_back(start_pnt);
    }
    else { //general case, we are in a face
      Face current_face = m_mesh->get<Face>(next_cell_id);
      ATriangles.push_back(current_face.getID());
      //==============================================================
      // CASE 1: DO WE ARE IN A FACE CONTAINING A SING. POINT?
      //==============================================================
      bool intersect_sing_point = false;

      SingularityPoint* next_sing_point = m_faces_to_singularity_on_surf[current_face.getID()];
      bool must_try_to_connect = false;
      if (next_sing_point != NULL && //face in a confusing ball ...
	  next_sing_point != AFromPnt) {// ... of another singularity point
	  
	must_try_to_connect = true;
      }
      else if (next_sing_point != NULL && //face in the confusing ball ...
	       next_sing_point == AFromPnt) {//... of the incoming singularity point
	  
	//Warning: completly empiric, we just try to detect cyclic lines
	if (APoints.size() >= 100)
	  must_try_to_connect = true;
      }
	
      if(must_try_to_connect) {
	std::cout << "\n --> triangle "  << current_face.getID() 
		  << " with sing. point" <<next_sing_point->getLocation()<< std::endl;
	std::cout << "START PNT: " << start_pnt << std::endl;
	math::Point start_dirPnt(start_pnt.X() + start_dir.X(),
				 start_pnt.Y() + start_dir.Y(),
				 start_pnt.Z() + start_dir.Z());

	std::cout << "START DIR PNT: " << start_dirPnt << std::endl;
	//Now, we look for a compatible slot
	std::vector<SingularityPoint::Slot*>& cur_slots = next_sing_point->getSlots();

	//==============================================================
	// We look for a free slot and we connect the line to it
	//==============================================================
	bool found_compatible_slot = false; 
	bool found_free_slot = false;
	found_pnt = next_sing_point;
	double slot_epsilon = 0.9;
	//we normalize the vector we arrive with
	math::Vector current_vec = start_dir;
	current_vec.normalize();
	while (!found_compatible_slot && slot_epsilon > 0.4) {
	  double best_deviation = -2;
	  double best_slot_id = 0;

	  std::cout<<"Current vec: "<<current_vec<<" and epsilon: "
		   <<slot_epsilon<<" and nb slots: "<<cur_slots.size()<<std::endl;
	  for (unsigned int i_slot = 0; i_slot < cur_slots.size(); i_slot++) {

	    SingularityPoint::Slot* current_slot = cur_slots[i_slot];
	    if(current_slot->isFreeze)
	      continue;
	    math::Vector            slot_opp_dir = current_slot->direction.opp();
	    slot_opp_dir.normalize();
	    std::cout<<"Slot "<<i_slot<<": "<<slot_opp_dir<<std::endl;
	    double slot_deviation =slot_opp_dir.dot(current_vec);
	    
	    if (slot_deviation > slot_epsilon && slot_deviation > best_deviation) {
	      best_deviation = slot_deviation;
	      best_slot_id=i_slot;
	    } //if (slot_deviation < slot_epsilon) {

	    
	  } //for (unsigned int i_slot = 0; !found_free_slot && i_slot < ....
	  
	  if(best_deviation!=-2 &&  !cur_slots.empty()){
	
	    SingularityPoint::Slot* best_slot = cur_slots[best_slot_id];
	    math::Vector            slot_opp_dir = best_slot->direction.opp();
	    std::cout<<"Best slot "<<best_slot_id<<": "<<slot_opp_dir<<std::endl;
	    if ( best_slot->isLaunched) {
	      //slot already assigned with a prvious line (and direction)
	      math::Vector prev_line_dir = best_slot->line_direction; 
	      double prev_deviation =slot_opp_dir.dot(prev_line_dir);
	      if(best_deviation>prev_deviation) {
		//the new alignment is better than the previous one 
		found_free_slot = false;	
		found_compatible_slot = true;
		found_slot = best_slot;
		intersect_sing_point = true;
		std::cout<<" WE CHOOSE A NEW CONNECTION"<<std::endl;
	      }
	      else {
		// We keep the previous association
		found_compatible_slot = false;
		std::cout<<" WE KEEP THE PREV CONNECTION"<<std::endl;
	      }
		  
	    } //if ( current_slot>isLaunched) {
	    else {  //WE HAVE A FREE SLOT
	      // We keep the previous association
	      found_free_slot = true;	
	      found_compatible_slot = true;
	      found_slot = best_slot;
	      intersect_sing_point = true;
	      std::cout<<" WE FOUND A FREE SLOT"<<std::endl;
	    }
	    
	    // HAVE WE FIND THE END OF THE LINE??
	    if(found_compatible_slot) {
	      // COMPATIBLE AND FREE SLOT 
	      
	      find_end = true;
	    }
	  } 
	  slot_epsilon -=0.1;
	    	    
	}//while (!found_compatible_slot && slot_epsilon < -0.4) 

      }//	if(must_try_to_connect) {

      
      //==============================================================
      // CASE 2: SO, WE JUST HAVE TO CROSS THE TRIANGLE (GENERAL CONF)
      //==============================================================
      //Does the current triangle has the same classif
      if ( !intersect_sing_point) {

	math::Point  out_pnt;
	math::Vector out_vec;
	TCellID out_cell_id;
	int out_cell_dim;
	m_tool.traverseTriangle(current_face,   /* the face we work on*/
			 start_pnt,      /* the geometric point we start from */
			 start_dir,      /* the geometric direction to follow*/
			 start_cell_dim, /* the dimension of the cell start_pnt is located */ 
			 start_cell_id,  /* the id of the cell start_pnt is located on*/ 
			 out_pnt,        /* the geometric point where we go out */
			 out_vec,        /* the geometric direction to follow after*/
			 out_cell_dim,   /* the dimension of the out cell (0 or 1) */
			 out_cell_id);   /* the id of the out cell*/ 

	//we keep the point toPnt to define the line
	APoints.push_back(out_pnt);
	// we progress to the next point, next vector and so next face too
	prev_dir = start_dir; //we store the prev direction for slot
			      //reconnection with balls
	start_pnt = out_pnt;
	start_dir = out_vec;
	start_cell_dim = out_cell_dim;
	start_cell_id = out_cell_id;

      } //if (!intersect_line && !intersect_sing_point) {

      if (intersect_sing_point)
	find_end = true;
  

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
  AToDir = start_dir;
  AToPnt = start_pnt;


  AEndOnBnd = end_on_boundary;
  
  //the end point must be created if it has not been found
  APntToCreate = (found_pnt==0);

  //singularity point data if we found an end point
  AToSingPnt = found_pnt;
  AToSlot = found_slot;
  AToCellDim = start_cell_dim;
  AToCellID = start_cell_id;
 
}
/*----------------------------------------------------------------------------*/
void SingularityGraphBuilder2D::
createGeometricSingularityPoint(const math::Point&      AInPnt,
				const math::Vector&     AInVec,
				const int               ACellDim,
				const TCellID           ACellID,
				SingularityPoint*&       APnt,
				SingularityPoint::Slot*& ASlot)
{
  if(ACellDim==0){

    Node n = m_mesh->get<Node>(ACellID);

    if(m_mesh->isMarked(n,m_mark_nodes_on_point)){
      //We arrive onto a geometric point !!!!
      // It means that a singularity geometric point already exist for it
      // We look for it
      std::vector<VertexSingularityPoint* >  
	geom_sing_points = m_graph.getVertexPoints();

      bool found_pnt = false;
      for(unsigned int i=0; i<geom_sing_points.size() && !found_pnt; i++){
	VertexSingularityPoint* current_pnt = geom_sing_points[i];
	Node current_node = current_pnt->getMeshNode();
	if(current_node.getID()==ACellID){
	  found_pnt = true;
	}
      }

      if(!found_pnt){
	throw GMDSException("createGeometricSingularityPoint: Error, no geometric point to be connected to!");
      }
    } // if(m_mesh->isMarked(n,m_mark_nodes_on_point))
    else{
      // We are on a geometrical curve
      std::cout<<"Split a boundary line at node "<<n.getID()<<std::endl;
      m_graph.splitCurveLine(AInPnt, AInVec, n, APnt, ASlot);
      
    }
  } // if(ACellDim==0)
  else {
    Edge e = m_mesh->get<Edge>(ACellID);
    std::cout<<"Split a boundary line at edge "
	     <<e.get<Node>()[0].getID()<<", "
	     <<e.get<Node>()[1].getID()<<std::endl;
    m_graph.splitCurveLine(AInPnt, AInVec, e, APnt, ASlot);

  }
}
/*----------------------------------------------------------------------------*/
void SingularityGraphBuilder2D::
backtrackSingularityLine(SurfaceSingularityLine* ALine,
			 SingularityPoint*       AFromPnt,
			 SingularityPoint::Slot* AFromSlot,
			 SingularityPoint*       AToPnt,
			 SingularityPoint::Slot* AToSlot)
{
  //==============================================================
  // We keep in mind the length of the line 
  //==============================================================
  double line_length = ALine->length(); 
  double length_tolerance = 1.2*line_length;
  
  SingularityPoint*        arrival_sing;
  SingularityPoint::Slot*  arrival_slot;
  math::Point              arrival_pnt;
  math::Vector             arrival_dir;
  int                      arrival_cell_dim;
  TCellID                  arrival_cell_id;
  std::vector<math::Point> line_pnts;
  std::vector<TCellID>     new_traversed_triangles;
  bool arrival_on_bnd;
  bool arrival_on_free_slot;
  bool arrival_pnt_to_create;
  
  computeStreamLine(AFromPnt,
		    AFromSlot,
		    arrival_sing,
		    arrival_slot,
		    arrival_pnt,
		    arrival_dir,
		    line_pnts,
		    new_traversed_triangles,
		    arrival_cell_dim,
		    arrival_cell_id,
		    arrival_on_bnd,
		    arrival_on_free_slot,
		    arrival_pnt_to_create);

  if(arrival_sing!=AToPnt)
    throw GMDSException("Backtraking issue: We don't reach the departure point");



  std::vector<gmds::math::Point> old_pnts = ALine->getDiscretizationPoints(); 
  //first and last points must not be taken into account, they are on the
  //singularity pnts.
  //
  // WARNING
  // old_pnts and line_pnts are traversed in opposite directions. Moreover,
  // they don't terminate at the same location. One goes from the singularity
  // point 1 to the boundary of the  confusing ball 2 and the other does the 
  // opposite.
  std::vector<gmds::math::Point> new_pnts;

  new_pnts.push_back(old_pnts[0]);
  math::Point last_point = old_pnts[old_pnts.size()-1];

   
  for(unsigned int i_old = 1; i_old<old_pnts.size()-1;i_old++){

    math::Point current_pnt = old_pnts[i_old];

    math::Segment seg(line_pnts[0],line_pnts[1]);
    math::Point proj_pnt = seg.project(current_pnt);
    double proj_dist = current_pnt.distance(proj_pnt);
    
    for(unsigned int i_new = 2; i_new<line_pnts.size();i_new++){
      math::Segment seg_i(line_pnts[i_new-1],line_pnts[i_new]);
      math::Point proj_pnt_i = seg_i.project(current_pnt);
      double proj_dist_i = current_pnt.distance(proj_pnt_i);
      if(proj_dist_i<proj_dist){
	proj_dist = proj_dist_i;
	proj_pnt = proj_pnt_i;
      }
    }
    // We have our new point
    new_pnts.push_back(0.5*current_pnt+0.5*proj_pnt);
  }

  new_pnts.push_back(last_point);
  //  double line_nb_discretization_pnts = line_points.size();
  ALine->setDiscretizationPoints(new_pnts);

  //==============================================================
  // Now we recompute the intersected faces. It is quite complicated
  // and could be improved. But up to now, it works.
  //==============================================================

  // We keep in mind the faces traversed the first time and during the
  // backtracking process. Indeed, we don't know which faces will be really
  // traversed at the end.
  std::vector<TCellID> all_traversed_faces = ALine->getTraversedFaces();
  all_traversed_faces.insert(all_traversed_faces.end(),
			     new_traversed_triangles.begin(),
			     new_traversed_triangles.end());
  std::set<TCellID> candidates;
  for(unsigned int i=0; i<all_traversed_faces.size();i++) {
    Face f=m_mesh->get<Face>(all_traversed_faces[i]);
    std::vector<Node> f_nodes = f.get<Node>();
    for(unsigned int j=0;j<f_nodes.size();j++){
      Node nj = f_nodes[j];
      std::vector<TCellID> nj_faces = nj.getIDs<Face>();
      for(unsigned int k=0;k<nj_faces.size();k++){
	candidates.insert(nj_faces[k]);
      }
    }
  }
  
  std::set<TCellID> set_of_traversed_faces;
  
  for(std::set<TCellID>::iterator it = candidates.begin(); 
      it!=candidates.end(); it++) {
    Face f=m_mesh->get<Face>(*it); 

    std::vector<Node> f_nodes = f.get<Node>();
    math::Triangle t(f_nodes[0].getPoint(),
		     f_nodes[1].getPoint(),
		     f_nodes[2].getPoint());
    bool found_pnt=false;
    for(unsigned int j=0;j<new_pnts.size() && !found_pnt; j++){
      math::Point pj = new_pnts[j];
      if(t.isIn(pj)){
	found_pnt = true;
	set_of_traversed_faces.insert(f.getID());
      }
      if(j!=0){ 
	math::Point pk = new_pnts[j-1];
	if(t.intersect(math::Segment(pj,pk))){
	  found_pnt = true;
	  set_of_traversed_faces.insert(f.getID());
	}
	
      }
    }
  }

  std::vector<TCellID> final_traversed_faces;
  final_traversed_faces.insert(final_traversed_faces.end(),
			       set_of_traversed_faces.begin(),
			       set_of_traversed_faces.end());
  ALine->setTraversedFaces(final_traversed_faces);
}
/*----------------------------------------------------------------------------*/
std::vector<SurfaceSingularityLine*> 
SingularityGraphBuilder2D::getSingularityLinesIn(const Face& AFace)
{
  std::vector<SurfaceSingularityLine* > found_lines;

  if(m_faces_to_singularity_on_surf[AFace.getID()]!=0) 
    return found_lines;
  std::vector<SurfaceSingularityLine* > lines = m_graph.getSurfaceLines();

  TCellID face_id = AFace.getID();
  for(unsigned int i=0;i<lines.size();i++){
    SurfaceSingularityLine* line_i = lines[i];
    if(line_i->isTraversed(face_id)){
      found_lines.push_back(line_i);
    }
  }

  return found_lines;
}

/*----------------------------------------------------------------------------*/
void SingularityGraphBuilder2D::
writeOutput(const std::string& AFileName)
{
  static int out = 0;
  std::stringstream file_name;
  file_name << m_output_directory_name << "/" << AFileName << "_" << out;
  writeOutputSingle(file_name.str());
  out++;
}
/*----------------------------------------------------------------------------*/
void SingularityGraphBuilder2D::
writeOutputSingle(const std::string& AFileName)
{
  m_graph.createVTKOutputFile(AFileName);
  //std::stringstream file_name;
  //file_name << AFileName<<"_FROM_MESH";
  VTKWriter<IGMesh> writer(*m_mesh);
  //writer.write(file_name.str(), DIM3 | R | F | N);
}
/*-----------------------------------------------------------------*/
void SingularityGraphBuilder2D::
writeSingularityPointsAndSlots() 
{

  gmds::IGMesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
  std::vector<SingularityPoint*> points = m_graph.getPoints();
  std::vector<SingularityPoint*>::iterator it = points.begin();
  for(;it!=points.end();it++){
    SingularityPoint* sing = *it;
    math::Point loc = sing->getLocation(); 
    gmds::Node center = m.newNode(loc.X(), loc.Y(), loc.Z());
    std::vector< SingularityPoint::Slot*>& slots = sing->getSlots();
   
    for (unsigned int i = 0; i <slots.size(); i++){
      SingularityPoint::Slot* si = slots[i];
      
      math::Point si_loc = si->location;
      math::Vector si_dir = si->direction;
      si_dir.normalize();
      gmds::Node si_departure = m.newNode(si_loc.X(), si_loc.Y(), si_loc.Z());
      gmds::Node si_end = m.newNode(si_loc.X() + si_dir.X(), 
				    si_loc.Y() + si_dir.Y(), 
				    si_loc.Z() + si_dir.Z() );
      
      m.newTriangle(center, si_departure, si_departure);
      m.newTriangle(si_end, si_departure, si_departure);
    }
  }
  
  gmds::VTKWriter<IGMesh> w(m);
  //<gmds::DIM3 | gmds::F | gmds::N | gmds::F2N> w(m);
  w.write("points_and_slots", gmds::F | gmds::N);

}

/*----------------------------------------------------------------------------*/
