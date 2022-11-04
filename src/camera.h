#ifndef AVP_CAMERA_H_
#define AVP_CAMERA_H_
#include <vector>

#include "Eigen/Core"

namespace avp {
class Camera {
 public:
  typedef Eigen::Matrix4f Mat4f;
  typedef Eigen::Matrix3f Mat3f;
  typedef Eigen::Matrix2f Mat2f;
  typedef Eigen::Vector4f Point4f;
  typedef Eigen::Vector3f Point3f;
  typedef Eigen::Vector2f Point2f;
  Camera() {}
  ~Camera() {}
  static Camera::Point3f Body2Camera(const Camera::Mat4f& T_c_b,
                                     const Camera::Point3f& p_b);
  static Camera::Point3f Camera2Body(const Camera::Mat4f& T_b_c,
                                     const Camera::Point3f& p_c);
  static Camera::Point3f Pixel2Camera(const Camera::Mat3f& K,
                                      const Camera::Point2f& p_p,
                                      const float depth = 1.0);
  static Camera::Point3f Pixel2Body(const Camera::Mat4f& T_c_b,
                                    const Camera::Mat3f& K,
                                    const Camera::Point2f& p_p,
                                    const float depth = 1.0);
  static Camera::Point2f Body2Pixel(const Camera::Mat4f& T_c_b,
                                    const Camera::Mat3f& K,
                                    const Camera::Point3f& p_b);
  static Camera::Point2f Camera2Pixel(const Camera::Mat3f& K,
                                      const Camera::Point3f& p_c);
  static Camera::Point2f Undistort(const float* distceoffs,
                                   const Camera::Mat3f& K,
                                   const Camera::Point2f& p);
  static void Undistort(const float* distceoffs, const Camera::Mat3f& K,
                        const std::vector<Camera::Point2f>& raw_ps,
                        std::vector<Camera::Point2f>* un_ps);
};
}  // namespace avp

#endif
