/*---------------------------------------------------------------------------*/
/*
 * SingularityGraph.cpp
 *
 *  Created on: 13 juil. 2014
 *      Author: bibi
 */
/*---------------------------------------------------------------------------*/
#include "SingularityGraph.h"
/*---------------------------------------------------------------------------*/
#include <GMDS/IG/IGMesh.h>
#include <GMDS/IO/VTKWriter.h>
#include <GMDS/Utils/Exception.h>
#include <GMDS/Math/Segment.h>
#include <GMDS/Math/Triangle.h>
#include <GMDS/Math/Numerics.h>
/*---------------------------------------------------------------------------*/
using namespace gmds;
/*---------------------------------------------------------------------------*/
SingularityGraph::SingularityGraph(IGMesh *AMesh)
 :m_mesh(AMesh)
{}
/*---------------------------------------------------------------------------*/
SingularityGraph::~SingularityGraph()
{
  for (unsigned int i = 0; i < m_points.size(); i++)
    delete m_points[i];
  for (unsigned int i = 0; i < m_lines.size(); i++)
    delete m_lines[i];
}
/*---------------------------------------------------------------------------*/
std::vector<SingularityPoint*>&
SingularityGraph::getPoints(){
  return m_points;
}
/*---------------------------------------------------------------------------*/
std::vector<VolumeSingularityLine*>
SingularityGraph::getVolumeLines(){
  std::vector<VolumeSingularityLine*> l;
  for (unsigned int i = 0; i < m_lines.size(); i++)
    {
      if (m_lines[i]->getType() == SingularityLine::VOLUME)
	l.push_back(dynamic_cast<VolumeSingularityLine*>(m_lines[i]));
    }
  return l;
}
/*---------------------------------------------------------------------------*/
std::vector<SurfaceSingularityLine*>
SingularityGraph::getSurfaceLines(){
  std::vector<SurfaceSingularityLine*> l;
  for (unsigned int i = 0; i < m_lines.size(); i++)
    {
      if (m_lines[i]->getType() == SingularityLine::SURFACE)
	l.push_back(dynamic_cast<SurfaceSingularityLine*>(m_lines[i]));
    }
  return l;
}
/*---------------------------------------------------------------------------*/
std::vector<CurveSingularityLine*>
SingularityGraph::getCurveLines(){
  std::vector<CurveSingularityLine*> l;
  for (unsigned int i = 0; i < m_lines.size(); i++)
    {
      if (m_lines[i]->getType() == SingularityLine::CURVE)
	l.push_back(dynamic_cast<CurveSingularityLine*>(m_lines[i]));
    }
  return l;
}

/*---------------------------------------------------------------------------*/
std::vector<VolumeSingularityPoint*>
SingularityGraph::getVolumePoints(){
  std::vector<VolumeSingularityPoint*> vol_points;
  for (unsigned int i = 0; i < m_points.size(); i++)
    {
      if (m_points[i]->getGeomType() == SingularityPoint::VOLUME)
	{
	  VolumeSingularityPoint* vol_pnt = dynamic_cast<VolumeSingularityPoint*>(m_points[i]);
	  vol_points.push_back(vol_pnt);
	}
	
    }
  return vol_points;
}
/*---------------------------------------------------------------------------*/
std::vector<VertexSingularityPoint*>
SingularityGraph::getVertexPoints()
{
  std::vector<VertexSingularityPoint*> points;
  for (unsigned int i = 0; i < m_points.size(); i++)
    {
      if (m_points[i]->getGeomType() == SingularityPoint::VERTEX)
	points.push_back(dynamic_cast<VertexSingularityPoint*>(m_points[i]));
    }
  return points;
}
/*---------------------------------------------------------------------------*/
std::vector<CurveSingularityPoint*>
SingularityGraph::getCurvePoints()
{
  std::vector<CurveSingularityPoint*> points;
  for (unsigned int i = 0; i < m_points.size(); i++)
    {
      if (m_points[i]->getGeomType() == SingularityPoint::CURVE)
	points.push_back(dynamic_cast<CurveSingularityPoint*>(m_points[i]));
    }
  return points;
}
/*---------------------------------------------------------------------------*/
std::vector<SurfaceSingularityPoint*>
SingularityGraph::getSurfacePoints()
{
  std::vector<SurfaceSingularityPoint*> points;
  for (unsigned int i = 0; i < m_points.size(); i++)
    {
      if (m_points[i]->getGeomType() == SingularityPoint::SURFACE)
	points.push_back(dynamic_cast<SurfaceSingularityPoint*>(m_points[i]));
    }
  return points;
}
/*---------------------------------------------------------------------------*/
std::vector<SingularityPatch*>
SingularityGraph::getSurfacePatchs()
{
  return  m_patchs;
}
/*---------------------------------------------------------------------------*/
std::vector<SingularityLine*>&
SingularityGraph::getLines(){
  return m_lines;
}
/*---------------------------------------------------------------------------*/
SurfaceSingularityLine* SingularityGraph::newSurfaceLine(){
  SurfaceSingularityLine* l = new SurfaceSingularityLine();
  m_lines.push_back(l);
  l->setNumber(getNbLines());
  return l;
}
/*---------------------------------------------------------------------------*/
SingularityPatch* SingularityGraph::newSurfacePatch() {
  SingularityPatch* p = new SingularityPatch();
  m_patchs.push_back(p);
  return p;
}
/*---------------------------------------------------------------------------*/
SurfaceSingularityPoint* 
SingularityGraph::newSurfacePoint()
{
  SurfaceSingularityPoint* p = new SurfaceSingularityPoint(m_mesh);
  m_points.push_back(p);
  return p;
}
/*---------------------------------------------------------------------------*/
VolumeSingularityLine* SingularityGraph::newVolumeLine(){
  VolumeSingularityLine* l = new VolumeSingularityLine();
  m_lines.push_back(l);
  l->setNumber(getNbLines());
  return l;
}
/*---------------------------------------------------------------------------*/
VolumeSingularityPoint* 
SingularityGraph::newVolumePoint()
{
  VolumeSingularityPoint* p = new VolumeSingularityPoint(m_mesh);
  m_points.push_back(p);

  return p;
}
/*---------------------------------------------------------------------------*/
CurveSingularityPoint* 
SingularityGraph::newCurvePoint(){
  CurveSingularityPoint* p = new CurveSingularityPoint(m_mesh);
  m_points.push_back(p);
  return p;
}
/*---------------------------------------------------------------------------*/
CurveSingularityLine* SingularityGraph::newCurveLine(){
  CurveSingularityLine* l = new CurveSingularityLine();
  m_lines.push_back(l);
  l->setNumber(getNbLines());
  return l;
}
/*---------------------------------------------------------------------------*/
VertexSingularityPoint* 
SingularityGraph::newVertexPoint(){
  VertexSingularityPoint* p = new VertexSingularityPoint(m_mesh);
  m_points.push_back(p);
  return p;
}
/*---------------------------------------------------------------------------*/
int SingularityGraph::getNbPoints() const
{
  return m_points.size();
}
/*---------------------------------------------------------------------------*/
int SingularityGraph::getNbLines() const
{
  return m_lines.size();
}
/*---------------------------------------------------------------------------*/
int SingularityGraph::getNbVolumePoints() const
{
  int nbPnts = 0;
  std::cout << "nb points = " << m_points.size();
  for (unsigned int i = 0; i < m_points.size(); i++)
    {
      SingularityPoint* p = m_points[i];
      if (p->getGeomType() == SingularityPoint::VOLUME)
	{
	  nbPnts++;
	}
    }
  return nbPnts;
}
/*---------------------------------------------------------------------------*/
int SingularityGraph::getNbSurfacePoints() const
{
  int nbPnts = 0;
  for (unsigned int i = 0; i < m_points.size(); i++)
    {
      SingularityPoint* p = m_points[i];
      if (p->getGeomType() == SingularityPoint::SURFACE)
	nbPnts++;
    }
  return nbPnts;
}
/*---------------------------------------------------------------------------*/
int SingularityGraph::getNbCurvePoints() const
{
  int nbPnts = 0;
  for (unsigned int i = 0; i < m_points.size(); i++)
    {
      SingularityPoint* p = m_points[i];
      if (p->getGeomType() == SingularityPoint::CURVE)
	nbPnts++;
    }
  return nbPnts;
}
/*---------------------------------------------------------------------------*/
int SingularityGraph::getNbVertexPoints() const
{
  int nbPnts = 0;
  for (unsigned int i = 0; i < m_points.size(); i++)
    {
      SingularityPoint* p = m_points[i];
      if (p->getGeomType() == SingularityPoint::VERTEX)
	nbPnts++;
    }
  return nbPnts;
}
/*---------------------------------------------------------------------------*/
int SingularityGraph::getNbVolumeLines() const
{
  int nb = 0;
  for (unsigned int i = 0; i < m_lines.size(); i++)
    {
      SingularityLine* l = m_lines[i];
      if (l->getType() == SingularityLine::VOLUME)
	nb++;
    }
  return nb;
}
/*---------------------------------------------------------------------------*/
int SingularityGraph::getNbSurfaceLines() const
{
  int nb = 0;
  for (unsigned int i = 0; i < m_lines.size(); i++)
    {
      SingularityLine* l = m_lines[i];
      if (l->getType() == SingularityLine::SURFACE)
	nb++;
    }
  return nb;
}
/*---------------------------------------------------------------------------*/
int SingularityGraph::getNbCurveLines() const
{
  int nb = 0;
  for (unsigned int i = 0; i < m_lines.size(); i++)
    {
      SingularityLine* l = m_lines[i];
      if (l->getType() == SingularityLine::CURVE)
	nb++;
    }
  return nb;
}
/*---------------------------------------------------------------------------*/
void  SingularityGraph::
splitCurveLine(const gmds::math::Point&  APnt, 
	       const gmds::math::Vector& AVec, 
	       const gmds::Edge&         AEdge,
	       SingularityPoint*&        ASing,
	       SingularityPoint::Slot*&  ASlot)
{
  std::cout<<"splitting in point "<<APnt<<std::endl;
  //==============================================================
  // STEP 1 - We find the line that must be splitted
  //==============================================================
  std::vector<CurveSingularityLine*> lines = getCurveLines();

  bool line_found = false;
  CurveSingularityLine* split_curve = 0;
  int curve_index = -1;

  for (unsigned int i = 0;/* !line_found &&*/ i < lines.size(); i++)  {
    CurveSingularityLine* l = lines[i];
    l->healOrientation();
    std::vector<gmds::Edge> l_edges = l->getMeshEdges();
    for (unsigned int j = 0; !line_found && j < l_edges.size(); j++)
      {
	if (l_edges[j].getID() == AEdge.getID())
	  {
	    line_found = true;

	    curve_index = j;
	    split_curve = l;
	  }
      }
  }
  
  if (!line_found)
    throw gmds::GMDSException("ERROR: No curve line found");

  //==============================================================
  // STEP 2 - CREATION OF THE SINGULARITY POINT
  //==============================================================

  ASing = newCurvePoint();
  ASing->setLocation(APnt);
  ASing->addMeshEdge(AEdge);

  //==============================================================
  // STEP 3 - CONNECTION OF LINES TO POINT
  //==============================================================
  SingularityPoint* last_end_point = split_curve->removeSingularityPoint();
  split_curve->addSingularityPoint(ASing);

  //new line from last_end_point to ASing (inverted direction)
  CurveSingularityLine* new_half_line = newCurveLine();
  new_half_line->addSingularityPoint(last_end_point);
  new_half_line->addSingularityPoint(ASing);

  std::cout<<"NEW HALF FROM "<<last_end_point->getLocation()
	   <<" TO "<<ASing->getLocation()<<std::endl;

  //OPEN SLOTS OF ASing
  //no associated line for the second 

    ASing->newSlot(APnt, AVec.opp(),      // slot geometrical data
		 AEdge.getID(),1, //cell owner info (id+dim) 
		 true);           //it is on the surface
  

  //curve mesh association and discretization

  std::vector<Edge> first_part, second_part;
  std::vector<Edge> init_edges = split_curve->getMeshEdges();
  for (unsigned int j = 0; j < init_edges.size(); j++)  {
    if (j <= curve_index)
      first_part.push_back(init_edges[j]);
  }

  for (unsigned int j = init_edges.size()-1; j >0; j--) {
    if (j >= curve_index)
      second_part.push_back(init_edges[j]);
  }
  //The edge were new curves meet is shared by the 2 new curves
  

  split_curve->setMeshEdges(first_part);
  new_half_line->setMeshEdges(second_part);

  //now the discretization is built on the edges traversal
  std::vector<gmds::math::Point > new_discretization;
  //=========================================================
  //1st curve, i.e. split_curve
  //=========================================================
  SingularityPoint* first_end = split_curve->getEndPoints()[0];
  gmds::math::Point prev_pnt  = first_end->getLocation();
  std::cout<<"FIRST CURVE from "<<prev_pnt
	   <<" with edge "
	   <<first_part[0].get<Node>()[0]<<", "
	   <<first_part[0].get<Node>()[1]<<std::endl;
  bool first_end_on_vertex = (first_end->getGeomType() == SingularityPoint::VERTEX);
  new_discretization.push_back(prev_pnt);
  if (first_end_on_vertex){
    Edge current_edge = first_part[0];
    Edge next_edge = first_part[1];

    std::vector<Node> current_nodes = current_edge.get<Node>();
    std::vector<Node> next_nodes = next_edge.get<Node>();

    Node common_node;
    for (unsigned int i = 0; i < current_nodes.size(); i++){
      for (unsigned int j = 0; j < next_nodes.size(); j++){
	if (current_nodes[i] == next_nodes[j])
	  common_node = current_nodes[i];
      }
    }
   
    prev_pnt = common_node.getPoint();
    new_discretization.push_back(prev_pnt);
    
  }
  else { // initial sing. point on a mesh edge 
    gmds::Edge current_edge = first_part[1];
    std::vector<Node> current_nodes = current_edge.get<Node>();

    gmds::math::Point p0 = current_nodes[0].getPoint();
    gmds::math::Point p1 = current_nodes[1].getPoint();

    double d0 = prev_pnt.distance(p0);
    double d1 = prev_pnt.distance(p1);

    if (d0 < d1){
      prev_pnt = p1;
    }
    else{
      prev_pnt = p0;
    }
    new_discretization.push_back(prev_pnt);
  }
  for (unsigned int i = 2; i < first_part.size()-1; i++){
    Edge current_edge = first_part[i];
    std::vector<Node> current_nodes = current_edge.get<Node>();
		
    gmds::math::Point p0 = current_nodes[0].getPoint();
    gmds::math::Point p1 = current_nodes[1].getPoint();
		
    double d0 = prev_pnt.distance(p0);
    double d1 = prev_pnt.distance(p1);

    if (d0 < d1){
      new_discretization.push_back(p1);
      prev_pnt = p1;
    }
    else{
      new_discretization.push_back(p0);
      prev_pnt = p0;
    }

  }
  std::cout<<"here"<<std::endl;
  //we have reached the end of the first line
  //we add the last point
  new_discretization.push_back(split_curve->getEndPoints()[1]->getLocation());

  split_curve->setDiscretizationPoints(new_discretization);
  split_curve->healOrientation();
 
  std::cout<<"here 2"<<std::endl;

  //=========================================================
  //2nd curve, i.e. new_half_line 
  //=========================================================

  //prev_pnt is the AEdge point that we must not use to discretize the second
  //line
  new_discretization.clear();
  first_end = new_half_line->getEndPoints()[0];
  prev_pnt = first_end->getLocation();
  std::cout<<"We are in "<<prev_pnt<<std::endl;
  first_end_on_vertex = (first_end->getGeomType() == SingularityPoint::VERTEX);
  new_discretization.push_back(prev_pnt);

  if (first_end_on_vertex){
    std::cout<<"Start from vertex"<<std::endl;
    Edge current_edge = second_part[0];
    Edge next_edge = second_part[1];

    std::vector<Node> current_nodes = current_edge.get<Node>();
    std::vector<Node> next_nodes = next_edge.get<Node>();

    Node common_node;
    for (unsigned int i = 0; i < current_nodes.size(); i++){
      for (unsigned int j = 0; j < next_nodes.size(); j++){
	if (current_nodes[i] == next_nodes[j])
	  common_node = current_nodes[i];
      }
    }
    prev_pnt = common_node.getPoint();
    new_discretization.push_back(prev_pnt);

  }
  else { // initial sing. point on a mesh edge  
    std::cout<<"No from an edge???"<<std::endl;

    Edge current_edge = second_part[1];
    std::vector<Node> current_nodes = current_edge.get<Node>();

    gmds::math::Point p0 = current_nodes[0].getPoint();
    gmds::math::Point p1 = current_nodes[1].getPoint();

    double d0 = prev_pnt.distance(p0);
    double d1 = prev_pnt.distance(p1);

    if (d0 < d1){
      prev_pnt = p1;
    }
    else{
      prev_pnt = p0;
    }
    new_discretization.push_back(prev_pnt);
  }
  for (unsigned int i = 2; i < second_part.size() - 1; i++){
    Edge current_edge = second_part[i];
    std::vector<Node> current_nodes = current_edge.get<Node>();

    gmds::math::Point p0 = current_nodes[0].getPoint();
    gmds::math::Point p1 = current_nodes[1].getPoint();

    double d0 = prev_pnt.distance(p0);
    double d1 = prev_pnt.distance(p1);

    if (d0 < d1){
      new_discretization.push_back(p1);
      prev_pnt = p1;
    }
    else{
      new_discretization.push_back(p0);
      prev_pnt = p0;
    }

  }  std::cout<<"here 3"<<std::endl;

  //we have reached the end of the second line
  //we add the last point
  new_discretization.push_back(new_half_line->getEndPoints()[1]->getLocation());
  new_half_line->setDiscretizationPoints(new_discretization);
  new_half_line->healOrientation();  

  //INVERSE CONNECTIVITY
  std::cout<<"####### SPLIT OF CURVE LINES FOR AN EDGE ############"<<std::endl;
  ASing->addLine(split_curve);
  ASing->addLine(new_half_line);

}
/*---------------------------------------------------------------------------*/
void SingularityGraph::
splitCurveLine(const gmds::math::Point&  APnt, 
	       const gmds::math::Vector& AVec, 
	       const gmds::Node&         ANode,
	       SingularityPoint*&        ASing,
	       SingularityPoint::Slot*&  ASlot)
{
  
  //==============================================================
  // STEP 1 - We find the line that must be splitted
  //==============================================================
  std::vector<CurveSingularityLine*> lines = getCurveLines();

  bool line_found = false;
  CurveSingularityLine* split_curve = 0;
  int curve_index = -1;
  int node_index_in_edge = -1;

  for (unsigned int i = 0; !line_found && i < lines.size(); i++)  {

    CurveSingularityLine* l = lines[i];
    std::vector<gmds::Edge> l_edges = l->getMeshEdges();

    for (unsigned int j = 0; !line_found && j < l_edges.size(); j++) {

      std::vector<gmds::Node> nodes_j = l_edges[j].get<Node>();
      if (nodes_j[0].getID() == ANode.getID()) {
	line_found = true;
	node_index_in_edge = 0;
	curve_index = j;
	split_curve = l;
      }
      else if (nodes_j[1].getID() == ANode.getID()) {
	line_found = true;
	node_index_in_edge = 1;
	curve_index = j;
	split_curve = l;
      }

    } //for (unsigned int j = 0; !line_found && j < l_edges.size(); j++)

  }//for (unsigned int i = 0; !line_found && i < lines.size(); i++)

  std::cout << "Curve found? " << line_found << std::endl;
  if (!line_found)
    throw gmds::GMDSException("ERROR: No curve line found");

  //==============================================================
  // STEP 2 - CREATION OF THE SINGULARITY POINT
  //==============================================================
  ASing = newCurvePoint();
  ASing->setLocation(APnt);
  ASing->addMeshNode(ANode);

  //==============================================================
  // STEP 3 - CONNECTION OF LINES TO POINT
  //==============================================================
  //old line keeps its first part
  SingularityPoint* last_end_point =
    split_curve->removeSingularityPoint();

  split_curve->addSingularityPoint(ASing);
  
  //new line from last_end_point to ASing (inverted direction)
  CurveSingularityLine* new_half_line = newCurveLine();
  new_half_line->addSingularityPoint(last_end_point);
  new_half_line->addSingularityPoint(ASing);

  //OPEN SLOTS OF ASing
  //no associated line for the second 

  ASing->newSlot(APnt, AVec,       // slot geometrical informaiton
		 ANode.getID(), 0, //cell owner info (id+dim) 
		 true);            //it is on the surface
 
  //curve mesh association and discretization

  std::vector<Edge> first_part, second_part;
  std::vector<Edge> init_edges = split_curve->getMeshEdges();
  for (unsigned int j = 0; j < init_edges.size(); j++) {
      if (j <= curve_index)
	first_part.push_back(init_edges[j]);
    }
  if( node_index_in_edge==0)//the last edge must be removed
    first_part.pop_back();

  for (unsigned int j = init_edges.size()-1; j >0; j--) {
      if (j >= curve_index)
	second_part.push_back(init_edges[j]);
    }
  if( node_index_in_edge==1)//the last edge must be removed
    second_part.pop_back();

  split_curve->setMeshEdges(first_part);

  new_half_line->setMeshEdges(second_part);

  //now the discretization build on the edges traversal
  std::vector<gmds::math::Point > new_discretization;
  //=========================================================
  //1st curve, i.e. split_curve
  SingularityPoint* first_end = split_curve->getEndPoints()[0];
  gmds::math::Point prev_pnt = first_end->getLocation();

  bool first_end_on_vertex = (first_end->getGeomType() == SingularityPoint::VERTEX);
  new_discretization.push_back(prev_pnt);

  if (!first_end_on_vertex){
    Edge current_edge = first_part[0];
    Edge next_edge = first_part[1];

    std::vector<Node> current_nodes = current_edge.get<Node>();
    std::vector<Node> next_nodes = next_edge.get<Node>();

    Node common_node;
    for (unsigned int i = 0; i < current_nodes.size(); i++){
      for (unsigned int j = 0; j < next_nodes.size(); j++){
	if (current_nodes[i] == next_nodes[j])
	  common_node = current_nodes[i];
      }
    }
    prev_pnt = common_node.getPoint();
    new_discretization.push_back(prev_pnt);
		
  }
  else { // initial sing. point on a mesh edge 
    gmds::Edge current_edge = first_part[0];
    std::vector<Node> current_nodes = current_edge.get<Node>();

    gmds::math::Point p0 = current_nodes[0].getPoint();
    gmds::math::Point p1 = current_nodes[1].getPoint();

    double d0 = prev_pnt.distance(p0);
    double d1 = prev_pnt.distance(p1);

    if (d0 < d1){
      prev_pnt = p1;
    }
    else{
      prev_pnt = p0;
    }
    new_discretization.push_back(prev_pnt);
  }
  for (unsigned int i = 1; i < first_part.size()-1; i++){
    Edge current_edge = first_part[i];
    std::vector<Node> current_nodes = current_edge.get<Node>();
		
    gmds::math::Point p0 = current_nodes[0].getPoint();
    gmds::math::Point p1 = current_nodes[1].getPoint();
		
    double d0 = prev_pnt.distance(p0);
    double d1 = prev_pnt.distance(p1);

    if (d0 < d1){
      new_discretization.push_back(p1);
      prev_pnt = p1;
    }
    else{
      new_discretization.push_back(p0);
      prev_pnt = p0;
    }

  }
  //we have reached the end of the firs line
  //we add the last point
  new_discretization.push_back(split_curve->getEndPoints()[1]->getLocation());

  split_curve->setDiscretizationPoints(new_discretization);

  //prev_pnt is the AEdge point that we must not use to discretize the second line
  new_discretization.clear();

  //=========================================================
  //2nd curve, i.e. new_half_line
  first_end = new_half_line->getEndPoints()[0];
  prev_pnt = first_end->getLocation();

  first_end_on_vertex = (first_end->getGeomType() == SingularityPoint::VERTEX);
  new_discretization.push_back(prev_pnt);

  if (!first_end_on_vertex){
    Edge current_edge = second_part[0];
    Edge next_edge = second_part[1];

    std::vector<Node> current_nodes = current_edge.get<Node>();
    std::vector<Node> next_nodes = next_edge.get<Node>();

    Node common_node;
    for (unsigned int i = 0; i < current_nodes.size(); i++){
      for (unsigned int j = 0; j < next_nodes.size(); j++){
	if (current_nodes[i] == next_nodes[j])
	  common_node = current_nodes[i];
      }
    }
    prev_pnt = common_node.getPoint();
    new_discretization.push_back(prev_pnt);

  }
  else { // initial sing. point on a mesh edge 
    Edge current_edge = second_part[0];
    std::vector<Node> current_nodes = current_edge.get<Node>();

    gmds::math::Point p0 = current_nodes[0].getPoint();
    gmds::math::Point p1 = current_nodes[1].getPoint();

    double d0 = prev_pnt.distance(p0);
    double d1 = prev_pnt.distance(p1);

    if (d0 < d1){
      prev_pnt = p1;
    }
    else{
      prev_pnt = p0;
    }
    new_discretization.push_back(prev_pnt);
  }
  for (unsigned int i = 1; i < second_part.size() - 1; i++){
    Edge current_edge = second_part[i];
    std::vector<Node> current_nodes = current_edge.get<Node>();

    gmds::math::Point p0 = current_nodes[0].getPoint();
    gmds::math::Point p1 = current_nodes[1].getPoint();

    double d0 = prev_pnt.distance(p0);
    double d1 = prev_pnt.distance(p1);

    if (d0 < d1){
      new_discretization.push_back(p1);
      prev_pnt = p1;
    }
    else{
      new_discretization.push_back(p0);
      prev_pnt = p0;
    }

  }
  //we have reached the end of the second line
  //we add the last point
  new_discretization.push_back(new_half_line->getEndPoints()[1]->getLocation());

  new_half_line->setDiscretizationPoints(new_discretization);

  //INVERSE CONNECTIVITY
  std::cout<<"####### SPLIT OF CURVE LINES FOR A NODE ############"<<std::endl;
  ASing->addLine(split_curve);
  ASing->addLine(new_half_line);
}
/*---------------------------------------------------------------------------*/
void SingularityGraph::
splitSurfaceLine(SurfaceSingularityPoint*& APnt,
		 SurfaceSingularityLine*& ALine) 
{
  gmds::math::Point pnt_loc = APnt->getLocation();
  //==============================================================
  // We found the interesected segment
  //==============================================================
  std::vector<gmds::math::Point>& pnts = ALine->getDiscretizationPoints();
  int  index_pnt = 0;
  bool found_pnt = false;

  for (unsigned int i = 1; !found_pnt && i < pnts.size()-2; i++) {
    gmds::math::Segment sij(pnts[i], pnts[i+1]);
    if(sij.isIn(pnt_loc)){
      found_pnt = true;
      index_pnt = i;
    }
  }

  if(!found_pnt)
    throw gmds::GMDSException("Error during surface line splitting");

  
  //==============================================================
  //old line keeps its first part
  SingularityPoint* snd_end_point = ALine->removeSingularityPoint();
  ALine->addSingularityPoint(APnt);
  
  //==============================================================
  //new line from last_end_point to ASing (inverted direction)
  SurfaceSingularityLine* new_half_line = newSurfaceLine();
  new_half_line->addSingularityPoint(snd_end_point);
  new_half_line->addSingularityPoint(APnt);

 
  //==============================================================
  //discretization are computed
  //==============================================================

  std::vector<gmds::math::Point> first_part, second_part;
  first_part.push_back(ALine->getEndPoints()[0]->getLocation());
  for (unsigned int i = 0; i < index_pnt; i++) {
    first_part.push_back(pnts[i]);
  }

  first_part.push_back(ALine->getEndPoints()[1]->getLocation());
  second_part.push_back(new_half_line->getEndPoints()[0]->getLocation());
  for (unsigned int i = pnts.size()-1; i > index_pnt; i--) {
    second_part.push_back(pnts[i]);
  }
  second_part.push_back(new_half_line->getEndPoints()[1]->getLocation());

  ALine->setDiscretizationPoints(first_part);
  new_half_line->setDiscretizationPoints(second_part);

  ALine->healOrientation();
  new_half_line->healOrientation();

  //==============================================================
  // We recompute traversed faces
  //==============================================================
  std::vector<gmds::TCellID> prev_faces = ALine->getTraversedFaces();
  std::vector<gmds::TCellID> new_faces_1, new_faces_2;

  for(unsigned int i=0; i<prev_faces.size();i++){
    Face f = m_mesh->get<Face>(prev_faces[i]);
    std::vector<Node> f_nodes = f.get<Node>();
    math::Triangle t(f_nodes[0].getPoint(),
		     f_nodes[1].getPoint(),
		     f_nodes[2].getPoint());
    bool found = false;
    for(unsigned int j=1; !found && j<first_part.size();j++){
         math::Point pj = first_part[j];
	 if(t.isIn(pj)){
	   found = true;
	   new_faces_1.push_back(f.getID());
	 }
	 if(j!=0){ 
	   math::Point pk = first_part[j-1];
	   if(t.intersect(math::Segment(pj,pk))){
	     found_pnt = true;
	     new_faces_1.push_back(f.getID());
	   }
	 }
    }
    found = false;
    for(unsigned int j=1; !found && j<second_part.size();j++){
         math::Point pj = second_part[j];
	 if(t.isIn(pj)){
	   found = true;
	   new_faces_2.push_back(f.getID());
	 }
	 if(j!=0){ 
	   math::Point pk = second_part[j-1];
	   if(t.intersect(math::Segment(pj,pk))){
	     found_pnt = true;
	     new_faces_2.push_back(f.getID());
	   }
	 }
    }
  }

  ALine->setTraversedFaces(new_faces_1);
  new_half_line->setTraversedFaces(new_faces_2);
  
  //==============================================================
  //Connectivity from points to curves
  std::cout<<"####### SPLIT OF SURF LINES ############"<<std::endl;
  APnt->addLine(ALine);
  APnt->addLine(new_half_line);
}
/*---------------------------------------------------------------------------*/
void SingularityGraph::removeLine(SingularityLine* ALine)
{
  if(ALine==0)
    return;

  bool found_line = false;
  for (unsigned int i = 0; i < m_lines.size() && !found_line; i++)
    {
      SingularityLine* current_line = m_lines[i];
      if(current_line==ALine) {
	found_line = true;
	m_lines[i] = m_lines.back();
	m_lines.pop_back();
	delete ALine;
      }
    }
}
/*---------------------------------------------------------------------------*/
void SingularityGraph::buildSurfacePatchs()
{  
  //========================================================================
  //We start from a set of lines and points that must be connected
  // in an acceptable configuration. Wrong configurations will not
  // be detected leading to a wrong patch decomposition.  
  //========================================================================

  for(unsigned int i=0;i<m_points.size(); i++) {
    SingularityPoint *pi = m_points[i];
    pi->clearSlots();
  }

  for(unsigned int i=0;i<m_lines.size(); i++) {
    SingularityLine *li = m_lines[i];
    li->removeAllSingularityPoints();
  }

  //========================================================================
  // TOPOLOGY CLEANING
  //========================================================================

  // For each point, we remove the connected lines
  for(unsigned int i=0;i<m_points.size(); i++) {
    SingularityPoint *pi = m_points[i];
    pi->clearSlots();
  }

  //For each line, we remove connected points.
  for(unsigned int i=0;i<m_lines.size(); i++) {
    SingularityLine *li = m_lines[i];
    li->removeAllSingularityPoints();
  }

  //========================================================================
  // TOPOLOGY REBUILDING
  //========================================================================
  // Now we reconnnect to have a clean configuration. Each curve must have 
  // two end points.

  for(unsigned int i=0;i<m_lines.size(); i++) {
    SingularityLine *li = m_lines[i];
    std::vector<math::Point > li_points = li->getDiscretizationPoints();
    math::Point p1 = li_points[0];
    math::Point p2 = li_points[li_points.size()-1];
    double nb_connections=0;
    for(unsigned int j=0; j<m_points.size(); j++) {
      SingularityPoint *pj = m_points[j];
      if(math::near(p1.distance(pj->getLocation()),0.0) ||
	 math::near(p2.distance(pj->getLocation()),0.0) ) {
	pj->addLine(li);
	li->addSingularityPoint(pj);
	nb_connections++;
      }
    }
    if(nb_connections!=2){
      std::cout<<"nb connections: "<<nb_connections<<std::endl;
      throw GMDSException("Error during topoloy rebuilding: a line is not adjacent to 2 points");
    }
      
  }

  for(unsigned int i=0;i<m_points.size(); i++) {
    SingularityPoint *pi = m_points[i];
    std::cout<<"Point "<<pi->getLocation()<<":";
    for(unsigned int j=0;j<pi->getSlots().size();j++) {
      SingularityPoint::Slot* sj = pi->getSlots()[j];
      std::cout<<" "<<sj->line->getNumber();
    }
    std::cout<<std::endl;
  }

   
  // We store an integer for each curve to indicate how much time it has
  // been traversed. It must be 1 (on the boundary) or 2 (inside). To detect
  // if a curve is on the boundary, we check if it classified onto a curve.
  std::map<SingularityLine*, bool> direct, reverse;
  for(unsigned int i=0;i<m_lines.size(); i++) {
    direct [m_lines[i]] = false;
    reverse[m_lines[i]] = false;
  }

  //========================================================================
  // LOOP to build patches from lines
  //========================================================================
  for(unsigned int i=0;i<m_lines.size(); i++) {
    SingularityLine* li = m_lines[i];
   
    if(!direct[li]){ // first time at least one patch to build
      //We never work with curve line as first
      if(li->getType()==SingularityLine::CURVE)
	continue;
      //=====================================================
      // First traversal
      //=====================================================
      direct[li] = true;
      std::cout<<"Patch with lines ";
      // A patch is created
      SingularityPatch* patch = newSurfacePatch();

      std::vector<SingularityPoint*> points = li->getEndPoints();

      SingularityPoint* first_point = points[0];

      patch->addPoint(first_point);
      patch->addLine(li);
      std::cout<<" from "<<li->getNumber()<<" ";
      SingularityPoint* current_point = points[1];
      SingularityLine* current_line = current_point->nextLine(li);

      while(current_line->getNumber()!=li->getNumber()){
	patch->addPoint(current_point);
	patch->addLine(current_line);

	std::cout<<current_line->getNumber()<<" ";
	points = current_line->getEndPoints();
	if(current_point == points[0]){
	  current_point = points[1];
	  direct[current_line]=true;
	}
	else{
	  current_point = points[0];
	  reverse[current_line]=true;
	}

	current_line = current_point->nextLine(current_line);

      }
      
      //we decrement the index to traverse it again
      i--;
    }
    else if(!reverse[li]) { 
      //=====================================================
      // Second traversal
      //=====================================================
      //second time, a patch can be built if the curve is not on the boundary.
      // Second time we go from the 2nd end point towards the 1st one A patch is
      // created
      if(li->getType()==SingularityLine::CURVE)
	continue;

      std::cout<<"Patch with lines ";
      SingularityPatch* patch = newSurfacePatch();

      std::vector<SingularityPoint*> points = li->getEndPoints();

      SingularityPoint* first_point = points[1];

      patch->addPoint(first_point);
      patch->addLine(li);

      std::cout<<" from "<<li->getNumber()<<" ";

      SingularityPoint* current_point = points[0];
      SingularityLine* current_line = current_point->nextLine(li);
   
      while(current_line->getNumber()!=li->getNumber()) {

	patch->addPoint(current_point);
	patch->addLine(current_line);

	std::cout<<current_line->getNumber()<<" ";

	points = current_line->getEndPoints();
	if(current_point == points[0]){
	  current_point = points[1];
	  direct[current_line]=true;
	}
	else{
	  current_point = points[0];
	  reverse[current_line]=true;
	}

	current_line = current_point->nextLine(current_line);

      }
      
    }
    std::cout<<std::endl;
  }
}
/*---------------------------------------------------------------------------*/
void SingularityGraph::
createVTKOutputFile(const std::string& AFileName) const
{
  gmds::IGMesh m(gmds::MeshModel(gmds::DIM3 | gmds::F | gmds::N | gmds::F2N));
  if(m_patchs.empty()) {

    gmds::Variable<int>* var=0;
    try {
      var = m.getVariable<int>(GMDS_FACE, "lines");
    }
    catch(gmds::GMDSException&){
      var = m.newVariable<int>(GMDS_FACE, "lines");
    }

    for (unsigned int i = 0; i < m_lines.size(); i++)
      {
	SingularityLine* current_line = m_lines[i];
	std::vector<gmds::math::Point >&
	  points = current_line->getDiscretizationPoints();
	for (unsigned int j = 0; j < points.size() - 1; j++)
	  {
	    gmds::Node n1 = m.newNode(points[j].X(), points[j].Y(), points[j].Z());
	    gmds::Node n2 = m.newNode(points[j + 1].X(), points[j + 1].Y(), points[j + 1].Z());
	    gmds::Face f  =  m.newTriangle(n1, n1, n2);
	    (*var)[f.getID()] = current_line->getNumber();
	  }

      }
  }
  else {

    for (unsigned int i = 0; i < m_patchs.size(); i++){

      SingularityPatch* current_patch = m_patchs[i];

      std::vector<SingularityPoint* > current_points;
      current_patch->getPoints(current_points);

      std::vector<Node> current_nodes;
      for (unsigned int j = 0; j < current_points.size(); j++) {
	gmds::math::Point p = current_points[j]->getLocation();
	current_nodes.push_back(m.newNode(p.X(), p.Y(), p.Z()));
      }
      m.newFace(current_nodes);
    }

  }
  std::cout
    << "We write a mesh with " << m.getNbNodes() << " nodes and "
    << m.getNbFaces() << " faces in file "<<AFileName<< std::endl;
  
  gmds::VTKWriter<IGMesh> w(m);
  //<gmds::DIM3 | gmds::F | gmds::N | gmds::F2N> w(m);
  w.write(AFileName, gmds::F | gmds::N);
  std::cout << "\t writing done in " << AFileName << std::endl;


}
/*---------------------------------------------------------------------------*/
