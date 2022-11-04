#ifndef ONBOARD_ROTATION_H_
#define ONBOARD_ROTATION_H_

#include "data_types.h"
#include "fast_math.h"

namespace avp {
// 2d
inline void RotatePoint(const float* pt, const float theta, float* p,
                        float* j_p_theta = nullptr) {
  const float cos = fast_math::Cos(theta);
  const float sin = fast_math::Sin(theta);
  p[0] = cos * pt[0] - sin * pt[1];
  p[1] = sin * pt[0] + cos * pt[1];
  if (j_p_theta) {
    j_p_theta[0] = -sin * pt[0] - cos * pt[1];
    j_p_theta[1] = cos * pt[0] - sin * pt[1];
  }
}

inline Point2f RotatePoint(const Point2f& pt, const float theta,
                           float* j_p_theta = nullptr) {
  Point2f point_tf;
  RotatePoint(reinterpret_cast<const float*>(pt.data()), theta,
              reinterpret_cast<float*>(point_tf.data()), j_p_theta);

  return point_tf;
}

// yaw(z), pitch(y), roll(x)
struct EulerAngles {
  float yaw;
  float pitch;
  float roll;
  EulerAngles(float _yaw, float _pitch, float _roll) {
    yaw = _yaw;
    pitch = _pitch;
    roll = _roll;
  }
  EulerAngles() { yaw = pitch = roll = 0.0; }
};

struct Quaternion {
  float w;
  float x;
  float y;
  float z;
  Quaternion(float _x, float _y, float _z, float _w) {
    x = _x;
    y = _y;
    z = _z;
    w = _w;
  }
  Quaternion() {
    x = y = z = 0.0;
    w = 1.0;
  }
};
Quaternion ToQuaternion(float yaw, float pitch, float roll, float* j_q_e = 0);
Quaternion ToQuaternion(const EulerAngles& angles, float* j_q_e = 0);
EulerAngles ToEulerAngles(const Quaternion& q);

// quaternion use matrix to rotate
void dfQuaternionToRMatrix(const float* q, float* r, float* j_q_r = 0);
void dfRMatrixRotatePoint(const float* r, const float* p1, float* p,
                          float* j_p_r = 0);
void dfQuaternionRotatePointUseRMatrix(const float* q, const float* p1,
                                       float* p, float* j_p_r = 0);

// quaternion wxyz
void dfQuaternionRotatePoint(const float* q, const float* p1, float* p,
                             float* j_p_q = 0);
void dfUnitQuaternionRotatePoint(const float* q, const float* p1, float* p,
                                 float* j_p_uq = 0);
void dfQuaternionInvRotatePoint(const float* q, const float* p1, float* p2,
                                float* j_p2_q = 0, float* r_inv = 0);
void dfUnitQuaternionInvRotatePoint(const float* q, const float* p1, float* p,
                                    float* j_p_uq = 0, float* r_inv = 0);
void dfQuaternionMultiply(const float* q1, const float* q2, float* q12,
                          float* j_q12_q1 = 0, float* j_q12_q2 = 0);
void dfQuaternionConjugate(const float* q, float* qc, float* j_qc_q = 0);
}  // namespace avp
#endif
