#ifndef STREAMLINE_BUILDER_H_
#define STREAMLINE_BUILDER_H_

#include <GMDS/IG/IGMesh.h>
#include <GMDS/Math/Chart.h>
#include <GMDS/Math/AxisAngleRotation.h>
#include <GMDS/Math/SHarmonicL4.h>
#include <GMDS/Utils/Log.h>
#include <GMDS/Utils/Variable.h>
#include <GMDS/Math/VectorND.h>
#include <Params.h>
#include <map>
#include <set>

using namespace gmds;

class StreamlineBuilder
{
  public:
    StreamlineBuilder(IGMesh& AMesh,const fhedo::ParamsGlobal& AGParam, IGMesh& blank, fhedo::ParamsMark& AMarks)
    :m_mesh(AMesh),m_global_params(AGParam),m_streamline_mesh(blank),m_bm(AMarks)
    {
    };

    ~StreamlineBuilder(){};

    void getCrossField();

    void execute();


  private:
    void computeStreamlineAtNode(Node& current_node, math::Vector3d& current_dir);

    void initMeshMarks();

    void cleanMeshMarks();

    void storeBoundaryData();

    void writeSolution();

    void drawStreamline(math::Point& point1, math::Point& point2);

    math::Vector3d closestComponentVector(math::Vector3d& dir, math::Chart& frame);

    int findSingularRegions();

    bool extend(math::Point& current_point, math::Vector3d& current_dir, Face& current_face, Region& current_region,
                                                          math::Point& next_point, math::Vector3d& next_dir, Face& next_face, Region& next_region);

    bool isPointinTetrahedron(math::Point& point,Region& r);

    math::Chart interpolateFrameField(math::Point& point, Region& r);

    bool getOtherRegionOnFace(Face& f, Region& current_region, Region& other_region);

    void streamlineAtNodeID(int id);

    void findConcavities();

    void sendSheetsFromConcavities();

    void writeFrameField();

    IGMesh& m_mesh;

    IGMesh& m_streamline_mesh;

    fhedo::ParamsGlobal m_global_params;

    Variable<math::Chart>* m_chart_field;

    Variable<math::AxisAngleRotation>* m_rot_field;

    Variable<math::SHarmonicL4>* m_harmonic_field;

    fhedo::ParamsMark m_bm;

    std::map<gmds::TCellID, gmds::math::Vector3d> m_bnd_normals;

    std::vector<Node> m_surface_nodes;

    std::vector<Node> m_curve_nodes;

    std::vector<Node> m_point_nodes;

    std::vector<Edge> m_surface_edges;

    std::vector<Edge> m_curve_edges;

    std::vector<Face> m_surface_faces;

    std::set<Node> m_concave_nodes;

    std::map<Node,math::Vector3d> m_concavity_map_1;

    std::map<Node,math::Vector3d> m_concavity_map_2;


    int m_mark_node_on_streamline;

    int m_streamline_region;

    int m_singular_region;

    int m_mark_concavity;

};

#endif //STREAMLINE_BUILDER_H_
