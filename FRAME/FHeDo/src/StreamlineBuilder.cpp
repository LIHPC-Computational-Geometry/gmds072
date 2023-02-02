#include <StreamlineBuilder.h>
#include <GMDS/Algo/BoundaryOperator.h>
#include <GMDS/IO/VTKWriter.h>
#include <GMDS/Math/Ray.h>
#include <GMDS/Math/Triangle.h>
#include <float.h>
#include <iostream>
#include <cmath>

void StreamlineBuilder::execute()
{
  getCrossField();

  initMeshMarks();

  storeBoundaryData();

  findConcavities();

  // sendSheetsFromConcavities

  for(Node n: m_surface_nodes)
  {
    if( n.getID() % 50 == 0 )
    {
      math::Vector3d dir = m_bnd_normals[n.getID()];
      dir = -1*dir;
      math::Vector3d z = {0,1,0};
      // if( dir.dot(z) < .9 && dir.dot(z) > -.9 )
        computeStreamlineAtNode(n,dir);
    }
  }

  writeSolution();

  cleanMeshMarks();
}

void StreamlineBuilder::getCrossField()
{
    m_chart_field = m_mesh.newVariable<math::Chart>(GMDS_NODE, "chart");
    m_harmonic_field = m_mesh.getVariable<math::SHarmonicL4>(GMDS_NODE,"SHL4");

    IGMesh::node_iterator node_it = m_mesh.nodes_begin();

    for(; !node_it.isDone(); node_it.next())
    {
      Node n = node_it.value();
      (*m_chart_field)[n.getID()] = (*m_harmonic_field)[n.getID()].chart();
    }
}

void StreamlineBuilder::computeStreamlineAtNode(Node& current_node, math::Vector3d& current_dir)
{
  //===========================================
  // set up initial conditions for extend function
  //===========================================
  bool region_set = false;

  math::Point current_point = current_node.getPoint();
  Face current_face;
  Region current_region;

  std::vector<Face> adj_faces;
  current_node.get<Face>(adj_faces);

  for(Face f: adj_faces)
  {
    if ( m_mesh.isMarked(f, m_bm.mark_face_on_surf) )
    {
      current_face = f; // it doesn't matter which face we use when starting at a node
      break;
    }
  }

  std::vector<Region> adj_regions;
  current_node.get<Region>(adj_regions);

  double length = 1;

  for(Region r: adj_regions)
  {
    math::Vector3d b(current_point,r.center());
    length = min(length,b.norm());
  }

  math::Chart frame = (*m_chart_field)[current_node.getID()];
  math::Vector3d v = closestComponentVector(current_dir,frame);
  v.normalize();
  v = v*length;

  math::Point tip(current_node.X() + v.X(), current_node.Y() + v.Y(), current_node.Z() + v.Z());

  current_dir = v;

  for(Region r: adj_regions)
  {
    if( isPointinTetrahedron(tip,r) )
    {
      current_region = r;
      region_set = true;
    }
  }

  if(!region_set)
  {
    printf("region not set!\n");
    return;
  }

  if( m_mesh.isMarked(current_region, m_singular_region) )
  {
    printf("streamline ran into singular region\n");
    return;
  }

  //======================================
  // extend until we reach the boundary
  //======================================
  math::Point next_point;
  bool next_point_set = false;
  math::Vector3d next_dir;
  Face next_face;
  Region next_region;

  while( extend(current_point, current_dir, current_face, current_region, next_point, next_dir, next_face, next_region) )
  {
    // if( m_mesh.isMarked(current_region, m_singular_region) )
    // {
    //   printf("streamline ran into singular region\n");
    //   return;
    // }

    drawStreamline(current_point, next_point);


    // printf("current_face: %d next_face: %d\n", current_face.getID(),next_face.getID());
    // printf("current_region: %d next_region: %d\n", current_region.getID(),next_region.getID());

    current_point = next_point;
    next_point_set = true;
    current_dir = next_dir;
    current_face = next_face;
    current_region = next_region;
  }
  // printf("current_face: %d next_face: %d\n", current_face.getID(),next_face.getID());

  if(next_point_set)
    drawStreamline(current_point, next_point);
}

bool StreamlineBuilder::extend(math::Point& current_point, math::Vector3d& current_dir, Face& current_face, Region& current_region,
                                                          math::Point& next_point, math::Vector3d& next_dir, Face& next_face, Region& next_region)
{
  bool next_point_set = false;
  bool next_dir_set = false;
  bool next_face_set = false;
  bool next_region_set = false;

  math::Chart current_chart;

  current_chart = interpolateFrameField(current_point,current_region);

  math::Vector3d current_ref3d = closestComponentVector(current_dir,current_chart);
  math::Vector current_ref(current_ref3d.X(),current_ref3d.Y(), current_ref3d.Z());

  math::Ray ray0(current_point,current_ref);

  std::vector<Face> adj_faces;
  current_region.get<Face> (adj_faces);

  math::Point temp_point;
  bool temp_set = false;

  //======================================================================
  // loop through faces of region until we find one that intersects the ray in the current direction to find temp point
  //======================================================================
  for(Face f: adj_faces)
  {
    if ( f.getID() != current_face.getID() )
    {
      std::vector<Node> nodes;
      f.get<Node>(nodes);

      math::Point intersection_point;

      math::Triangle triangle(nodes[0].getPoint(), nodes[1].getPoint(), nodes[2].getPoint() );

      bool intersected;

      try
      {
        intersected = ray0.intersect3D(triangle, intersection_point);
      }
      catch (GMDSException& e)
      {
        intersected = false;
      }

      if( intersected )
      {
        temp_point = intersection_point;
        temp_set = true;
        break;
      }
    }
  }

  if(!temp_set)
  {
    Region other_region;
    for(Face f: adj_faces)
    {
      if( getOtherRegionOnFace(f, current_region, other_region) )
      {
        std::vector<Face> other_region_faces;
        other_region.get<Face>(other_region_faces);

        for(Face f: other_region_faces)
        {
          if ( f.getID() != current_face.getID() )
          {
            std::vector<Node> nodes;
            f.get<Node>(nodes);

            math::Point intersection_point;

            math::Triangle triangle(nodes[0].getPoint(), nodes[1].getPoint(), nodes[2].getPoint() );

            bool intersected;

            try
            {
              intersected = ray0.intersect3D(triangle, intersection_point);
            }
            catch (GMDSException& e)
            {
              intersected = false;
            }

            if( intersected )
            {
              temp_point = intersection_point;
              temp_set = true;
              break;
            }
          }
        }
      }

      if(temp_set)
        break;
    }
  }

  if(!temp_set)
  {
    printf("temp not set on line 238\n");
    return false;
  }

  //==============================================
  // find the average of the reference directions and extend in that direction
  //==============================================
  math::Chart temp_frame;
  temp_frame = interpolateFrameField(temp_point, current_region);

  math::Vector3d temp_ref3d = closestComponentVector(current_dir,temp_frame);
  math::Vector temp_ref(temp_ref3d.X(),temp_ref3d.Y(),temp_ref3d.Z());

  math::Vector next_dir_for_ray = (current_ref + temp_ref)/2;
  next_dir = (current_ref3d + temp_ref3d)/2;
  next_dir_set = true;

  next_dir_for_ray.normalize();
  math::Point ray_point(current_point.X() + .0001*next_dir_for_ray.X(), current_point.Y() + .0001*next_dir_for_ray.Y(), current_point.Z() + .0001*next_dir_for_ray.Z());

  math::Ray ray1(ray_point,next_dir_for_ray);

  //=======================================================================
  //  loop through facess unti we find the face that intersects the extension direction. Update info for next iteration
  //=======================================================================
  for(Face f: adj_faces)
  {
    if ( f.getID() != current_face.getID() )
    {
      std::vector<Node> nodes;
      f.get<Node>(nodes);

      math::Point intersection_point;

      math::Triangle triangle(nodes[0].getPoint(), nodes[1].getPoint(), nodes[2].getPoint() );

      bool intersected;

      try
      {
        intersected = ray1.intersect3D(triangle, intersection_point);
      }
      catch (GMDSException& e)
      {
        intersected = false;
      }

      if( intersected )
      {
        next_point = intersection_point;
        next_point_set = true;
        next_face = f;
        next_face_set = true;

        if( getOtherRegionOnFace(f, current_region, next_region) )
        {
          next_region_set = true;
          break;
        }
        else
        {
          // printf("no region on opposite face line 287\n");
          return false;
        }

      }
    }
  }

  if(!next_point_set)
  {
    Region other_region;
    for(Face f: adj_faces)
    {
      if( getOtherRegionOnFace(f, current_region, other_region) )
      {
        std::vector<Face> other_region_faces;
        other_region.get<Face>(other_region_faces);

        for(Face f: other_region_faces)
        {
          if ( f.getID() != current_face.getID() )
          {
            std::vector<Node> nodes;
            f.get<Node>(nodes);

            math::Point intersection_point;

            math::Triangle triangle(nodes[0].getPoint(), nodes[1].getPoint(), nodes[2].getPoint() );

            bool intersected;

            try
            {
              intersected = ray1.intersect3D(triangle, intersection_point);
            }
            catch (GMDSException& e)
            {
              intersected = false;
            }

            if( intersected )
            {
              next_point = intersection_point;
              next_point_set = true;
              next_face = f;
              next_face_set = true;

              if( getOtherRegionOnFace(f, current_region, next_region) )
              {
                next_region_set = true;
                break;
              }
              else
              {
                // printf("no region on opposite of face line 330\n");
                return false;
              }

            }
          }
        }
      }

      if(next_point_set)
        break;
    }
  }

  if(!next_point_set || !next_dir_set || !next_face_set || !next_region_set)
  {
    printf("something didn't get set from within extend function\n");
    return false;
    // throw GMDSException("something didn't get set from within extend function");
  }

  return true;
}

void StreamlineBuilder::initMeshMarks()
{
    m_mark_node_on_streamline  = m_mesh.getNewMark<Node>();

    m_streamline_region = m_mesh.getNewMark<Region>();

    m_mark_concavity = m_mesh.getNewMark<Node>();
}

void StreamlineBuilder::cleanMeshMarks()
{
    // ==================================================================
    // clean marks
    m_mesh.unmarkAll<Node>(m_mark_node_on_streamline);

    m_mesh.unmarkAll<Region>(m_streamline_region);

    m_mesh.unmarkAll<Node>(m_mark_concavity);


    //==================================================================
    //free marks
    m_mesh.freeMark<Node>(m_mark_node_on_streamline);

    m_mesh.freeMark<Region>(m_streamline_region);

    m_mesh.freeMark<Node>(m_mark_concavity);

}

void StreamlineBuilder::storeBoundaryData()
{
  //============================================
  // store nodes
  //============================================
    IGMesh::node_iterator node_it = m_mesh.nodes_begin();

    for(; !node_it.isDone(); node_it.next())
    {
      Node n = node_it.value();
      if( m_mesh.isMarked(n, m_bm.mark_node_on_surf))
        m_surface_nodes.push_back(n);
    }
    std::cout << "number of surface nodes: " << m_surface_nodes.size() << std::endl;

    for(Node n: m_surface_nodes)
    {
      if( m_mesh.isMarked(n, m_bm.mark_node_on_curv))
        m_curve_nodes.push_back(n);
    }
    std::cout << "number of curve nodes: " << m_curve_nodes.size() << std::endl;

    for(Node n: m_curve_nodes)
    {
      if( m_mesh.isMarked(n, m_bm.mark_node_on_pnt) )
        m_point_nodes.push_back(n);
    }
    std::cout << "number of point nodes: " << m_point_nodes.size() << std::endl;


  //==================================================
  // store edges
  //==================================================
    IGMesh::edge_iterator edge_it = m_mesh.edges_begin();

    for(; !edge_it.isDone(); edge_it.next())
    {
      Edge e = edge_it.value();
      if( m_mesh.isMarked(e, m_bm.mark_edge_on_surf))
        m_surface_edges.push_back(e);
    }
    std::cout << "number of surface edges: " << m_surface_edges.size() << std::endl;

    for(Edge e: m_surface_edges)
    {
      if( m_mesh.isMarked(e, m_bm.mark_edge_on_curv))
        m_curve_edges.push_back(e);
    }
    std::cout << "number of curve edges: " << m_curve_edges.size() << std::endl;

  //==================================================
  // store faces
  //==================================================
    IGMesh::face_iterator face_it = m_mesh.faces_begin();

    for(; !face_it.isDone(); face_it.next())
    {
      Face f = face_it.value();
      if( m_mesh.isMarked(f, m_bm.mark_face_on_surf))
        m_surface_faces.push_back(f);
    }
    std::cout << "number of surface faces: " << m_surface_faces.size() << std::endl;

    BoundaryOperator boundaryOp(&m_mesh);

    for(IGMesh::node_iterator itn = m_mesh.nodes_begin(); !itn.isDone();
        itn.next()){
        Node n = itn.value();
        if(m_mesh.isMarked(n,m_bm.mark_node_on_surf)){
            math::Vector3d nv= boundaryOp.getOutputNormalOfABoundaryNode(n);
            m_bnd_normals[n.getID()]=math::Vector3d(nv.X(),nv.Y(),nv.Z());
        }
    }
}

void StreamlineBuilder::writeSolution()
{
    static int nb_file = 0;

    VTKWriter<IGMesh> writer(m_streamline_mesh);

    writer.write(m_global_params.output_dir+"/streamlines", DIM3 | R | N);

    nb_file++;
}

void StreamlineBuilder::drawStreamline(math::Point& point1, math::Point& point2)
{
  Node node1 = m_streamline_mesh.newNode(point1);
  Node node2 = m_streamline_mesh.newNode(point1.X() + .001, point1.Y() + .001  , point1.Z() + .001);
  Node node3 = m_streamline_mesh.newNode(point2.X() + .001, point2.Y() + .001, point2.Z() + .001);
  Node node4 = m_streamline_mesh.newNode(point2);

  Region tet = m_streamline_mesh.newTet(node1,node2,node3,node4);
}

math::Vector3d StreamlineBuilder::closestComponentVector(math::Vector3d& dir, math::Chart& frame)
{
  math::Vector3d closest_vec;
  double max  = 0;
  for(int i = 0;i < 3; i++)
  {
    math::Vector3d vec = frame[i];
    double product = vec.dot(dir);
    product = abs(product);
    if(product > max)
    {
      max = product;
      closest_vec = vec;
    }
  }
  if( closest_vec.dot(dir) < 0 )
    closest_vec = -1*closest_vec;

  return closest_vec;
}

int StreamlineBuilder::findSingularRegions()
{
  int num_sing_regions = 0;
  IGMesh::region_iterator it = m_mesh.regions_begin();
  for(; !it.isDone(); it.next() )
  {
    Region current_region = it.value();
    std::vector<TCellID> node_ids = current_region.getIDs<Node>();
    math::Quaternion q[4];
    for(int i_n=0; i_n<4; i_n++){
        math::Chart ci= (*m_chart_field)[node_ids[i_n]];
        q[i_n]= math::Quaternion(ci);
    }

    int sing_type = math::Quaternion::testSingularity(q[0],q[1],q[2],q[3]);
    if (sing_type != 0)
    {
      m_mesh.mark(current_region, m_singular_region);
      num_sing_regions++;
    }
  }
  return num_sing_regions;
}

bool StreamlineBuilder::isPointinTetrahedron(math::Point& point,Region& r)
{
  std::vector<Node> adj_nodes;
  r.get<Node> (adj_nodes);

  std::vector<math::Point> tet_points;
  std::vector<TCoord> barcycentric_coordinates;

  for(Node n: adj_nodes)
  {
    tet_points.push_back(n.getPoint());
  }

  point.computeBarycentric(tet_points,point,barcycentric_coordinates);

  for(TCoord t: barcycentric_coordinates)
  {
    if(t < 0)
      return false;
  }

  return true;
}

math::Chart StreamlineBuilder::interpolateFrameField(math::Point& point, Region& r)
{
  std::vector<Node> adj_nodes;
  r.get<Node> (adj_nodes);

  std::vector<math::Point> tet_points;
  std::vector<TCoord> barcycentric_coordinates;

  for(Node n: adj_nodes)
  {
    tet_points.push_back(n.getPoint());
  }

  point.computeBarycentric(tet_points,point,barcycentric_coordinates);

  math::SHarmonicL4 frame;

  TCoord d1 = barcycentric_coordinates[0];
  TCoord d2 = barcycentric_coordinates[1];
  TCoord d3 = barcycentric_coordinates[2];
  TCoord d4 = barcycentric_coordinates[3];

  TCellID id1 = adj_nodes[0].getID();
  TCellID id2 = adj_nodes[1].getID();
  TCellID id3 = adj_nodes[2].getID();
  TCellID id4 = adj_nodes[3].getID();

  frame = d1*(*m_harmonic_field)[id1] + d2*(*m_harmonic_field)[id2] + d3*(*m_harmonic_field)[id3] + d4*(*m_harmonic_field)[id4];

  return frame.chart();
}

bool StreamlineBuilder::getOtherRegionOnFace(Face& face, Region& this_region, Region& other_region)
{
  std::vector<Region> adj_regions;
  face.get<Region>(adj_regions);

  if(adj_regions.size() == 1)
  {
    return false;
  }
  else if(adj_regions.size() == 2)
  {
    for(Region r: adj_regions)
      if( r.getID() != this_region.getID() )
      {
        other_region = r;
      }
    return true;
  }
  else
    throw GMDSException("face has invalid number of regions");
}

void StreamlineBuilder::streamlineAtNodeID(int id)
{
  Node n = m_mesh.get<Node>(id);
  math::Vector3d dir = m_bnd_normals[id];
  dir = -1*dir;

  computeStreamlineAtNode(n,dir);
}

void StreamlineBuilder::findConcavities()
{
  for(Edge e: m_curve_edges)
  {
    math::Point e_center = e.center();

    //get faces of each curve edge
    std::vector<Face> adj_faces;
    e.get<Face>(adj_faces);

    std::vector<Face> surface_faces;

    double totalAngle = 0;

    for(int i=0; i < adj_faces.size() - 1; i++)
    {
      Face current_face = adj_faces[i];
      Face next_face = adj_faces[i+1];
      math::Point current_center = current_face.center();
      math::Point next_center = next_face.center();
      math::Vector3d current_vector(e_center, current_center);
      math::Vector3d next_vector(e_center, next_center);

      current_vector.normalize();
      next_vector.normalize();
      totalAngle += current_vector.angle(next_vector);
    }

    for(Face f: adj_faces)
    {
      if( m_mesh.isMarked(f, m_bm.mark_face_on_surf) )
      {
        surface_faces.push_back(f);
      }
    }

    if( totalAngle > 180 )
    {
      std::vector<Node> adj_nodes;
      e.get<Node> (adj_nodes);

      m_mesh.mark(adj_nodes[0], m_mark_concavity);
      m_mesh.mark(adj_nodes[1], m_mark_concavity);
      m_concave_nodes.insert(adj_nodes[0]);
      m_concave_nodes.insert(adj_nodes[1]);

      math::Chart chart0 = (*m_chart_field)[adj_nodes[0].getID()];
      math::Chart chart1 = (*m_chart_field)[adj_nodes[1].getID()];

      math::Vector3d dir(e_center, surface_faces[0].center());
      math::Vector3d other_dir(e_center, surface_faces[1].center());

      m_concavity_map_1[adj_nodes[0]] = -dir;
      m_concavity_map_2[adj_nodes[0]] = -other_dir;
      // m_concavity_map_1[adj_nodes[1]] = dir;
      // m_concavity_map_2[adj_nodes[1]] = other_dir;
    }
  }

}

void StreamlineBuilder::sendSheetsFromConcavities()
{
  for(Node n: m_concave_nodes)
  {
    math::Vector3d dir;

    dir = m_concavity_map_1[n];
    computeStreamlineAtNode(n,dir);

    dir = m_concavity_map_2[n];
    computeStreamlineAtNode(n,dir);
  }
}

void StreamlineBuilder::writeFrameField()
{
  IGMesh::node_iterator it = m_mesh.nodes_begin();
  for(; !it.isDone(); it.next() )
  {
    Node n = it.value();
    math::SHarmonicL4 frame = (*m_harmonic_field)[n.getID()];

    ofstream out;
    out.open("FrameField.ffd");

    out << n.getID() << " " << frame << std::endl;
  }
}
