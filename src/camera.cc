#include "camera.h"

#include <stdlib.h>

#include <cmath>

// #define USE_OPENCV

#ifdef USE_OPENCV
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/core/core.hpp"
#endif

namespace avp {
#define EPSILON (1e-6)
Camera::Point3f Camera::Body2Camera(const Camera::Mat4f& T_c_b,
                                    const Camera::Point3f& p_b) {
  Camera::Point4f _p_b{p_b[0], p_b[1], p_b[2], 1.0};
  return (T_c_b * _p_b).head(3);
}

Camera::Point3f Camera::Camera2Body(const Camera::Mat4f& T_b_c,
                                    const Camera::Point3f& p_c) {
  Camera::Point4f _p_c{p_c[0], p_c[1], p_c[2], 1.0};
  return (T_b_c * _p_c).head(3);
}

Camera::Point3f Camera::Pixel2Camera(const Camera::Mat3f& K,
                                     const Camera::Point2f& p_p,
                                     const float depth) {
  const float fx = K(0, 0);
  const float fy = K(1, 1);
  const float cx = K(0, 2);
  const float cy = K(1, 2);
  if (std::abs(fx) < EPSILON || std::abs(fy) < EPSILON) abort();
  return Camera::Point3f((p_p[0] - cx) * depth / fx, (p_p[1] - cy) * depth / fy,
                         depth);
}

Camera::Point2f Camera::Camera2Pixel(const Camera::Mat3f& K,
                                     const Camera::Point3f& p_c) {
  Camera::Point3f _p_c{p_c[0] / p_c[2], p_c[1] / p_c[2], 1.0};
  return (K * _p_c).head(2);
}

Camera::Point2f Camera::Undistort(const float* distceoffs,
                                  const Camera::Mat3f& K,
                                  const Camera::Point2f& p) {
  const float k1 = distceoffs[0];
  const float k2 = distceoffs[1];
  const float p1 = distceoffs[2];
  const float p2 = distceoffs[3];
  const float k3 = distceoffs[4];

  const float fx = K(0, 0);
  const float fy = K(1, 1);
  const float cx = K(0, 2);
  const float cy = K(1, 2);

  if (std::abs(fx) < EPSILON || std::abs(fy) < EPSILON) abort();

  float m = 0.0, n = 0.0, a = 0.0, b = 0.0, r = 0.0;
  m = (p[0] - cx) / fx;
  n = (p[1] - cy) / fy;
  r = std::sqrt(std::pow(m, 2) + std::pow(n, 2));
  a = m * (1 + k1 * std::pow(r, 2) + k2 * std::pow(r, 4) +
           k3 * std::pow(r, 6)) +
      2 * p1 * m * n + p2 * (std::pow(r, 2) + 2 * std::pow(m, 2));
  b = n * (1 + k1 * std::pow(r, 2) + k2 * std::pow(r, 4) +
           k3 * std::pow(r, 6)) +
      p1 * (std::pow(r, 2) + 2 * std::pow(n, 2)) + 2 * p2 * m * n;
  return Camera::Point2f(fx * a + cx, fy * b + cy);
}

void Camera::Undistort(const float* distceoffs, const Camera::Mat3f& K,
                       const std::vector<Camera::Point2f>& raw_ps,
                       std::vector<Camera::Point2f>* un_ps) {
#ifndef USE_OPENCV
  const float k1 = distceoffs[0];
  const float k2 = distceoffs[1];
  const float p1 = distceoffs[2];
  const float p2 = distceoffs[3];
  const float k3 = distceoffs[4];

  const float fx = K(0, 0);
  const float fy = K(1, 1);
  const float cx = K(0, 2);
  const float cy = K(1, 2);
  for (const auto& raw_p : raw_ps) {
    if (std::abs(fx) < EPSILON || std::abs(fy) < EPSILON) abort();

    float m = 0.0, n = 0.0, a = 0.0, b = 0.0, r = 0.0;
    m = (raw_p[0] - cx) / fx;
    n = (raw_p[1] - cy) / fy;
    r = std::sqrt(std::pow(m, 2) + std::pow(n, 2));
    a = m * (1 + k1 * std::pow(r, 2) + k2 * std::pow(r, 4) +
             k3 * std::pow(r, 6)) +
        2 * p1 * m * n + p2 * (std::pow(r, 2) + 2 * std::pow(m, 2));
    b = n * (1 + k1 * std::pow(r, 2) + k2 * std::pow(r, 4) +
             k3 * std::pow(r, 6)) +
        p1 * (std::pow(r, 2) + 2 * std::pow(n, 2)) + 2 * p2 * m * n;
    un_ps->emplace_back(Camera::Point2f(fx * a + cx, fy * b + cy));
  }
#else
  const cv::Mat _K = (cv::Mat_<float>(3, 3) << K(0, 0), 0.0, K(0, 2), 0.0,
                      K(1, 1), K(1, 2), 0.0, 0.0, 1.0);
  std::vector<float> v_dist(distceoffs, distceoffs + 8);
  std::vector<cv::Point2f> raw_pts, un_pts;
  for (const auto& pt : raw_ps) {
    raw_pts.emplace_back(cv::Point2f(pt.x(), pt.y()));
  }

  cv::undistortPoints(raw_pts, un_pts, _K, v_dist, cv::Mat(), _K);

  for (const auto& pt : un_pts) {
    un_ps->emplace_back(Camera::Point2f(pt.x, pt.y));
  }
#endif
}

Camera::Point3f Camera::Pixel2Body(const Camera::Mat4f& T_c_b,
                                   const Camera::Mat3f& K,
                                   const Camera::Point2f& p_p,
                                   const float depth) {
  Eigen::Matrix<float, 3, 4> t_tmp = K * T_c_b.block<3, 4>(0, 0);
  Camera::Mat3f t = t_tmp.block<3, 3>(0, 0);
  t.col(2) = t_tmp.col(3);
  Camera::Point3f pw = t.inverse() * Camera::Point3f(p_p[0], p_p[1], 1.0);
  return Camera::Point3f(pw[0] / pw[2], pw[1] / pw[2], 0.0);
}

Camera::Point2f Camera::Body2Pixel(const Camera::Mat4f& T_c_b,
                                   const Camera::Mat3f& K,
                                   const Camera::Point3f& p_b) {
  return Camera2Pixel(K, Body2Camera(T_c_b, p_b));
}

}  // namespace avp
