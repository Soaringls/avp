#ifndef AVP_FRAME_H_
#define AVP_FRAME_H_
#include "data_types.h"
namespace avp {
class Frame {
public:
  double ts;
  int64_t frame_id;
  OdometryData odom;
  ImuData imu;
  WheelData wheel;
  LidarData lidar;
  VisionData vision;
};
} // namespace avp
#endif
