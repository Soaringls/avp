#pragma once

#include <utility>
#include <vector>

#include "Eigen/Dense"
#include "opencv2/core/core.hpp"

namespace qcraft {
namespace localization {
namespace visual {
typedef Eigen::Vector2d Point2d;
typedef Eigen::Vector2f Point2f;
typedef Eigen::Vector3d Point3d;
typedef Eigen::Vector3f Point3f;
typedef Eigen::Vector4d Point4d;
typedef Eigen::Vector4f Point4f;
typedef Eigen::Quaternionf Quaternionf;

struct PerceptionLine {
  // type: Divider-0 Curb-1
  int type;
  // line type: Solid-0 Dashed-1 Virtual-3
  int line_type;
  // position type: LLL-3 LL-2 L-1 R-1 RR-2 RRR-3
  int p_type = 0;
  int id;
  Point2f start_point;
  Point2f end_point;
  float coffes[4];
  std::vector<Point2f> points_uv;
  std::vector<float> points_scores;
  std::vector<Point2f> points_unuv;
  int index_vaild = 0;
  float length_vaild = 0.0;
  float distance = 0.0;
  bool invailed = true;
  std::vector<Point2f> points_vcs;
  std::vector<Point2f> points_vcs_sampling;
};

struct PerceptionObject {
  int type;
  int id;
  float score{0.0};
  std::vector<Point2f> points_uv;
  std::vector<Point2f> points_vcs;
};

struct PerceptionData {
  int64_t timestamp;
  int camera_id;
  std::vector<PerceptionLine> dividers;
  std::vector<PerceptionLine> road_edges;
  std::vector<PerceptionObject> divider_kps;
  std::vector<PerceptionObject> marking_kps;
  std::vector<PerceptionObject> markings;
  std::vector<PerceptionObject> poles;
  std::vector<PerceptionObject> signs;
};

struct ImgData {
  int64_t timestamp;
  int camera_id;
  cv::Mat img;
};

struct MapLine {
  enum Type : uint16_t {
    CURB = 0,
    DIVIDER = 1,
  };
  int line_type;
  // 0 curb 1 divider
  int type;
  int64_t id;
  int64_t section_id;
  std::vector<Point2d> line;
  std::vector<Point3d> line3d;
  std::vector<Point2f> line_vcs;  // actually local-enu
};

struct MapObject {
  int type;
  int64_t id;
  int64_t section_id;
  std::vector<Point2d> kps;
  std::vector<Point2f> kps_vcs;
  int confidence = 0;
};

struct MapData {
  enum Type : uint16_t {
    // motor vechile type
    GENERAL = 0,
    BUS_ONLY = 1,
    RAMP = 2,
    PARKING = 3,
    EMERGENCY = 4,
    MIXED_WITH_CYCLIST = 5,
    // invalid for vehicle
    BICYCLE_ONLY = 6,
    WALKING_STREET = 7,
  };
  Type lane_type = Type::GENERAL;
  bool is_in_intersection = false;
  float dist2intersection = 0.0;
  int loc_lane_num_l = 0;
  int loc_lane_num_r = 0;
  int64_t section_id = 0;
  int64_t lane_id = 0;
  std::vector<MapLine> dividers;
  std::vector<MapLine> road_edges;
  std::vector<MapLine> level_switch_zones;
  std::vector<MapObject> signs;
  std::vector<MapObject> poles;
  std::vector<MapObject> divider_kps;
  std::vector<MapObject> marking_kps;
  std::vector<MapObject> markings;
  std::vector<MapObject> traffic_lights;
};

struct LidarMatrixData {
  LidarMatrixData() = default;
  void clear() {
    num_scans = 0;
    scan_points.clear();
  }
  struct ScanData {
    int scan_id;
    std::vector<Point4d> points;
  };
  int num_scans;
  std::vector<ScanData> scan_points;
};

struct LidarPerceptionData {
  struct Pole {
    bool vaild{false};
    std::vector<Point4d> pts;
    PerceptionObject object;
  };
  std::vector<Point4d> raw_pts;
  int64_t e_timestamp;
  int64_t m_timestamp;
  int64_t s_timestamp;
  std::vector<Pole> poles;
  Point2f delta_xy;
  float delta_yaw;
};

struct GnssData {
  enum GnssType {
    INVALID = 0,
    PROPAGATED = 1,
    SINGLE = 2,
    PSRDIFF = 3,
    PPP = 4,
    RTK_FLOAT = 5,
    RTK_INTEGER = 6
  };
  int64_t timestamp;
  double lla[3];
  double lla_sigma[3];
  double rpy[3];
  double rpy_sigma[3];
  float v[3];
  int num_star;
  float hdop;
  GnssType type;
  bool vaild{false};
};

enum LocCode {
  Loc_Initing,
  Loc_Normal,
  Loc_LowPecersion,
  Loc_Invialed,
};

// algorithm error code
enum PerceptionErr {
  Perception_Err_NoBoundary,
  Perception_Err_NoLane,
  Perception_Err_NoPole,
  Perception_Err_NoSign,
};

enum MapErr {
  Map_Err_NoMap,
  Map_Err_OutDate,
};

enum MatchErr {
  Match_Err_Boundary,
  Match_Err_Lane,
  Match_Err_Pole,
  Match_Err_Sign,
  Match_Err_Marking,
};

enum ResetErr {
  Reset_Err_Dismatch_Too_Long,
  Reset_Err_Large_YawDiff,
  Reset_Err_Match_Avrdist_Too_Large,
};

enum MatchType : uint16_t {
  DIVIDER_PAIR = 0,
  ROAD_EDGE_PAIR = 1,
  POLE_PAIR = 2,
  OBJECT_PAIR = 3,
  SHAPE_PAIR = 4,
};

enum SemanticLaneIdx {
  R1 = -1,
  R2 = -2,
  R3 = -3,
  R4 = -4,
  L1 = 1,
  L2 = 2,
  L3 = 3,
  L4 = 4,
};

struct Node {
  Point2f pt_enu;
  Point2f pt_vcs;
};
struct ConstraintedPair {
  MatchType type;
  float offset_distance;
  bool valid = false;
};

struct PointConstraintedPair : ConstraintedPair {
  Node src;
  Node target;
};

struct PointLineConstraintedPair : ConstraintedPair {
  Node pt;
  std::vector<Node> line;
};

struct ShapeConstraintedPair : ConstraintedPair {
  std::vector<Node> src_shape;
  std::vector<Node> target_shape;
};

struct Pose {
  Point3f p;
  Quaternionf q;
};

struct LocateData {
  int64_t timestamp;
  Pose pose;
  double lla[3];
  float rpy[3];
  LocCode status{LocCode::Loc_Invialed};
  // perception map macth
  uint16_t error[4];
  int level_id;
  Point2f xy_offset;
  float yaw_offset;
  float avr_dist;
};

struct MapLaneLine {
  mapping::LaneBoundaryProto::Type type;
  mapping::ElementId id;
  // 2d points in vehicle coordinate system
  std::vector<Vec2d> pts_vcs;

  MapLaneLine(mapping::LaneBoundaryProto::Type type_, mapping::ElementId id_,
              std::vector<Vec2d> pts_vcs_)
      : type(type_), id(id_), pts_vcs(std::move(pts_vcs_)) {}
};

}  // namespace visual
}  // namespace localization
}  // namespace qcraft
