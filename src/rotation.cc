#include "rotation.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmath>

#include "Eigen/Core"

namespace avp {

Quaternion ToQuaternion(float yaw, float pitch, float roll,
                        float* j_q_e) {  // yaw(z), pitch(y), roll(x)
  float cy = cos(yaw * 0.5);
  float sy = sin(yaw * 0.5);
  float cp = cos(pitch * 0.5);
  float sp = sin(pitch * 0.5);
  float cr = cos(roll * 0.5);
  float sr = sin(roll * 0.5);
  Quaternion q;
  q.w = cr * cp * cy + sr * sp * sy;
  q.x = sr * cp * cy - cr * sp * sy;
  q.y = cr * sp * cy + sr * cp * sy;
  q.z = cr * cp * sy - sr * sp * cy;
  if (j_q_e) {
#define _J4x3Get(ptr, row, col) ((ptr)[row * 3 + col])  // NOLINT
    // j_qw_e
    _J4x3Get(j_q_e, 0, 0) = 0.5 * (cr * cp * (-sy) + sr * sp * cy);  // NOLINT
    _J4x3Get(j_q_e, 0, 1) = 0.5 * (cr * (-sp) * cy + sr * cp * sy);  // NOLINT
    _J4x3Get(j_q_e, 0, 2) = 0.5 * ((-sr) * cp * cy + cr * sp * sy);  // NOLINT
    // j_qx_e
    _J4x3Get(j_q_e, 1, 0) = 0.5 * (sr * cp * (-sy) - cr * sp * cy);  // NOLINT
    _J4x3Get(j_q_e, 1, 1) = 0.5 * (sr * (-sp) * cy - cr * cp * sy);  // NOLINT
    _J4x3Get(j_q_e, 1, 2) = 0.5 * (cr * cp * cy - (-sr) * sp * sy);  // NOLINT
    // j_qy_e
    _J4x3Get(j_q_e, 2, 0) = 0.5 * (cr * sp * (-sy) + sr * cp * cy);  // NOLINT
    _J4x3Get(j_q_e, 2, 1) = 0.5 * (cr * cp * cy + sr * (-sp) * sy);  // NOLINT
    _J4x3Get(j_q_e, 2, 2) = 0.5 * ((-sr) * sp * cy + cr * cp * sy);  // NOLINT
    // j_qy_e
    _J4x3Get(j_q_e, 3, 0) = 0.5 * (cr * cp * cy - sr * sp * (-sy));  // NOLINT
    _J4x3Get(j_q_e, 3, 1) = 0.5 * (cr * (-sp) * sy - sr * cp * cy);  // NOLINT
    _J4x3Get(j_q_e, 3, 2) = 0.5 * ((-sr) * cp * sy - cr * sp * cy);  // NOLINT
#undef _J4x3Get                                                      // NOLINT
  }
  return q;
}

Quaternion ToQuaternion(const EulerAngles& angles, float* j_q_e) {
  return ToQuaternion(angles.yaw, angles.pitch, angles.roll, j_q_e);
}

EulerAngles ToEulerAngles(const Quaternion& q) {
  EulerAngles angles;

  // roll (x-axis rotation)
  float sinr_cosp = 2 * (q.w * q.x + q.y * q.z);
  float cosr_cosp = 1 - 2 * (q.x * q.x + q.y * q.y);
  angles.roll = std::atan2(sinr_cosp, cosr_cosp);

  // pitch (y-axis rotation)
  float sinp = 2 * (q.w * q.y - q.z * q.x);
  if (std::abs(sinp) >= 1)
    angles.pitch = std::copysign(M_PI / 2, sinp);
  else
    angles.pitch = std::asin(sinp);

  // yaw (z-axis rotation)
  float siny_cosp = 2 * (q.w * q.z + q.x * q.y);
  float cosy_cosp = 1 - 2 * (q.y * q.y + q.z * q.z);
  angles.yaw = std::atan2(siny_cosp, cosy_cosp);
  return angles;
}

void dfQuaternionToRMatrix(const float* q, float* r, float* j_r_q) {
  const float w = q[0];
  const float x = q[1];
  const float y = q[2];
  const float z = q[3];
  // const float ww = w * w;
  const float wx = w * x;
  const float wy = w * x;
  const float wz = w * x;
  const float xx = -x * x;
  const float xy = x * y;
  const float xz = x * z;
  const float yy = -y * y;
  const float yz = y * z;
  const float zz = -z * z;
  const float scale = 1.f / std::sqrt(w * w + x * x + y * y + z * z);

#define _R(row, col) r[(row)*3 + (col)]  // NOLINT
  _R(0, 0) = 1 - 2 * (yy + zz);          // NOLINT
  _R(0, 1) = 2 * (xy - wz);              // NOLINT
  _R(0, 2) = 2 * (xz + wy);              // NOLINT

  _R(1, 0) = 2 * (xy + wz);      // NOLINT
  _R(1, 1) = 1 - 2 * (xx + zz);  // NOLINT
  _R(1, 2) = 2 * (yz - wx);      // NOLINT

  _R(2, 0) = 2 * (xz - wy);      // NOLINT
  _R(2, 1) = 2 * (yz + wx);      // NOLINT
  _R(2, 2) = 1 - 2 * (xx + yy);  // NOLINT

  if (j_r_q) {                                          // NOLINT
#define _J4x9Get(ptr, row, col) ((ptr)[row * 4 + col])  // NOLINT
    const float _2ss = 2 * scale * scale;               // NOLINT// NOLINT
    const float _2ws = 2 * w * scale;                   // NOLINT// NOLINT
    const float _2xs = 2 * x * scale;                   // NOLINT// NOLINT
    const float _2ys = 2 * y * scale;                   // NOLINT// NOLINT
    const float _2zs = 2 * z * scale;                   // NOLINT// NOLINT

    const float _2ssw = _2ss * w;  // NOLINT// NOLINT
    const float _2ssx = _2ss * x;  // NOLINT// NOLINT
    const float _2ssy = _2ss * y;  // NOLINT// NOLINT
    const float _2ssz = _2ss * z;  // NOLINT// NOLINT
    // j_r(0,0)_q// NOLINT
    _J4x9Get(j_r_q, 0, 0) = _2ws - _R(0, 0) * _2ssw;   // NOLINT
    _J4x9Get(j_r_q, 0, 1) = _2xs - _R(0, 0) * _2ssx;   // NOLINT
    _J4x9Get(j_r_q, 0, 2) = -_2ys - _R(0, 0) * _2ssy;  // NOLINT
    _J4x9Get(j_r_q, 0, 3) = -_2zs - _R(0, 0) * _2ssz;  // NOLINT
    // j_r(0,1)_q// NOLINT
    _J4x9Get(j_r_q, 1, 0) = -_2zs - _R(0, 1) * _2ssw;  // NOLINT
    _J4x9Get(j_r_q, 1, 1) = _2ys - _R(0, 1) * _2ssx;   // NOLINT
    _J4x9Get(j_r_q, 1, 2) = _2xs - _R(0, 1) * _2ssy;   // NOLINT
    _J4x9Get(j_r_q, 1, 3) = -_2ws - _R(0, 1) * _2ssz;  // NOLINT
    // j_r(0,2)_q// NOLINT
    _J4x9Get(j_r_q, 2, 0) = _2ys - _R(0, 2) * _2ssw;  // NOLINT
    _J4x9Get(j_r_q, 2, 1) = _2zs - _R(0, 2) * _2ssx;  // NOLINT
    _J4x9Get(j_r_q, 2, 2) = _2ws - _R(0, 2) * _2ssy;  // NOLINT
    _J4x9Get(j_r_q, 2, 3) = _2xs - _R(0, 2) * _2ssz;  // NOLINT
    // j_r(1,0)_q// NOLINT
    _J4x9Get(j_r_q, 3, 0) = _2zs - _R(1, 0) * _2ssw;  // NOLINT
    _J4x9Get(j_r_q, 3, 1) = _2ys - _R(1, 0) * _2ssx;  // NOLINT
    _J4x9Get(j_r_q, 3, 2) = _2xs - _R(1, 0) * _2ssy;  // NOLINT
    _J4x9Get(j_r_q, 3, 3) = _2ws - _R(1, 0) * _2ssz;  // NOLINT
    // j_r(1,1)_q// NOLINT
    _J4x9Get(j_r_q, 4, 0) = _2ws - _R(1, 1) * _2ssw;   // NOLINT
    _J4x9Get(j_r_q, 4, 1) = -_2xs - _R(1, 1) * _2ssx;  // NOLINT
    _J4x9Get(j_r_q, 4, 2) = _2ys - _R(1, 1) * _2ssy;   // NOLINT
    _J4x9Get(j_r_q, 4, 3) = -_2zs - _R(1, 1) * _2ssz;  // NOLINT
    // j_r(1,2)_q// NOLINT
    _J4x9Get(j_r_q, 5, 0) = -_2xs - _R(1, 2) * _2ssw;  // NOLINT
    _J4x9Get(j_r_q, 5, 1) = -_2ws - _R(1, 2) * _2ssx;  // NOLINT
    _J4x9Get(j_r_q, 5, 2) = _2zs - _R(1, 2) * _2ssy;   // NOLINT
    _J4x9Get(j_r_q, 5, 3) = _2ys - _R(1, 2) * _2ssz;   // NOLINT
    // j_r(2,0)_q// NOLINT
    _J4x9Get(j_r_q, 6, 0) = -_2ys - _R(2, 0) * _2ssw;  // NOLINT
    _J4x9Get(j_r_q, 6, 1) = _2zs - _R(2, 0) * _2ssx;   // NOLINT
    _J4x9Get(j_r_q, 6, 2) = -_2ws - _R(2, 0) * _2ssy;  // NOLINT
    _J4x9Get(j_r_q, 6, 3) = _2xs - _R(2, 0) * _2ssz;   // NOLINT
    // j_r(2,1)_q// NOLINT
    _J4x9Get(j_r_q, 7, 0) = _2xs - _R(2, 1) * _2ssw;  // NOLINT
    _J4x9Get(j_r_q, 7, 1) = _2ws - _R(2, 1) * _2ssx;  // NOLINT
    _J4x9Get(j_r_q, 7, 2) = _2zs - _R(2, 1) * _2ssy;  // NOLINT
    _J4x9Get(j_r_q, 7, 3) = _2ys - _R(2, 1) * _2ssz;  // NOLINT
    // j_r(2,2)_q// NOLINT
    _J4x9Get(j_r_q, 8, 0) = _2ws - _R(2, 2) * _2ssw;   // NOLINT
    _J4x9Get(j_r_q, 8, 1) = -_2xs - _R(2, 2) * _2ssx;  // NOLINT
    _J4x9Get(j_r_q, 8, 2) = -_2ys - _R(2, 2) * _2ssy;  // NOLINT
    _J4x9Get(j_r_q, 8, 3) = _2zs - _R(2, 2) * _2ssz;   // NOLINT

#undef _J4x9Get  // NOLINT
  }              // NOLINT
#undef _R        // NOLINT
  for (int i = 0; i < 9; i++) r[i] *= scale;
}

void dfRMatrixRotatePoint(const float* r, const float* p1, float* p,
                          float* j_p_r) {
#define _R(i, j) r[(i)*3 + (j)]  // NOLINT
  p[0] = _R(0, 0) * p1[0] + _R(0, 1) * p1[1] + _R(0, 2) * p1[2];
  p[1] = _R(1, 0) * p1[0] + _R(1, 1) * p1[1] + _R(1, 2) * p1[2];
  p[2] = _R(2, 0) * p1[0] + _R(2, 1) * p1[1] + _R(2, 2) * p1[2];
#undef _R
  if (j_p_r) {
#define _J_P_R(i, j) j_p_r[(i)*9 + (j)]  // NOLINT
    memset(j_p_r, 0.0, sizeof(float) * 3 * 9);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) _J_P_R(i, i * 3 + j) = p1[j];
#undef _J_P_R
  }
}

void dfQuaternionRotatePointUseRMatrix(const float* q, const float* p1,
                                       float* p, float* j_p_q) {
  float R[9];  // NOLINT
  float j_r_q[9 * 4];
  float j_p_r[3 * 9];
  dfQuaternionToRMatrix(q, R, (j_p_q) ? j_r_q : 0);
  dfRMatrixRotatePoint(R, p1, p, (j_p_q) ? j_p_r : 0);

  if (j_p_q) {
    Eigen::Map<const Eigen::Matrix<float, 9, 4, Eigen::RowMajor>> m_j_r_q(
        j_r_q);
    Eigen::Map<const Eigen::Matrix<float, 3, 9, Eigen::RowMajor>> m_j_p_r(
        j_p_r);
    Eigen::Map<Eigen::Matrix<float, 3, 4, Eigen::RowMajor>> m_j_p_q(j_p_q);
    m_j_p_q = m_j_p_r * m_j_r_q;
  }
}

// quaternion wxyz
void dfUnitQuaternionRotatePoint(const float* q, const float* p1, float* p,
                                 float* j_p_q) {
  const float w = q[0];
  const float x = q[1];
  const float y = q[2];
  const float z = q[3];
  // const float ww = w * w;
  const float wx = w * x;
  const float wy = w * y;
  const float wz = w * z;
  const float xx = -x * x;
  const float xy = x * y;
  const float xz = x * z;
  const float yy = -y * y;
  const float yz = y * z;
  const float zz = -z * z;

  const float a = p1[0], b = p1[1], c = p1[2];
  p[0] = 2 * ((yy + zz) * a + (xy - wz) * b + (xz + wy) * c) + a;
  p[1] = 2 * ((xy + wz) * a + (xx + zz) * b + (yz - wx) * c) + b;
  p[2] = 2 * ((xz - wy) * a + (yz + wx) * b + (xx + yy) * c) + c;

  if (j_p_q) {
#define _J3x4Get(ptr, row, col) ((ptr)[row * 4 + col])         // NOLINT
    _J3x4Get(j_p_q, 0, 0) = 2 * (-z * b + y * c);              // NOLINT
    _J3x4Get(j_p_q, 0, 1) = 2 * (y * b + z * c);               // NOLINT
    _J3x4Get(j_p_q, 0, 2) = 2 * (-2 * y * a + x * b + w * c);  // NOLINT
    _J3x4Get(j_p_q, 0, 3) = 2 * (-2 * z * a - w * b + x * c);  // NOLINT

    _J3x4Get(j_p_q, 1, 0) = 2 * (z * a - x * c);              // NOLINT
    _J3x4Get(j_p_q, 1, 1) = 2 * (y * a - 2 * x * b - w * c);  // NOLINT
    _J3x4Get(j_p_q, 1, 2) = 2 * (x * a + z * c);              // NOLINT
    _J3x4Get(j_p_q, 1, 3) = 2 * (w * a - 2 * z * b + y * c);  // NOLINT

    _J3x4Get(j_p_q, 2, 0) = 2 * (-y * a + x * b);              // NOLINT
    _J3x4Get(j_p_q, 2, 1) = 2 * (z * a + w * b - 2 * x * c);   // NOLINT
    _J3x4Get(j_p_q, 2, 2) = 2 * (-w * a + z * b - 2 * y * c);  // NOLINT
    _J3x4Get(j_p_q, 2, 3) = 2 * (x * a + y * b);               // NOLINT
#undef _J3x4Get                                                // NOLINT
  }
}

void dfQuaternionRotatePoint(const float* q, const float* p1, float* p,
                             float* j_p_q) {
  const float w = q[0];
  const float x = q[1];
  const float y = q[2];
  const float z = q[3];
  const float scale = 1.f / std::sqrt(w * w + x * x + y * y + z * z);
  if (std::isinf(scale)) {
    printf("scale is inf\n");
    abort();
  }
  float unit_q[4] = {scale * w, scale * x, scale * y, scale * z};
  float j_p_uq[12];
  dfUnitQuaternionRotatePoint(unit_q, p1, p, j_p_q ? j_p_uq : 0);

  if (j_p_q) {
    Eigen::Matrix4f j_uq_q(Eigen::Matrix4f::Zero());
    float scale3 = scale * scale * scale;
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) j_uq_q(i, j) = -scale3 * q[i] * q[j];
    for (int i = 0; i < 4; i++) j_uq_q(i, i) += scale;

    Eigen::Map<Eigen::Matrix<float, 3, 4, Eigen::RowMajor>> m_j_p_q(j_p_q);
    Eigen::Map<Eigen::Matrix<float, 3, 4, Eigen::RowMajor>> m_j_p_uq(j_p_uq);
    m_j_p_q = m_j_p_uq * j_uq_q;
  }
}

void dfQuaternionInvRotatePoint(const float* q, const float* p1, float* p,
                                float* j_p_q, float* r_inv) {
  const float w = q[0];
  const float x = q[1];
  const float y = q[2];
  const float z = q[3];
  const float scale = 1.0 / std::sqrt(w * w + x * x + y * y + z * z);
  if (std::isinf(scale)) {
    printf("scale is inf\n");
    abort();
  }
  float unit_q[4] = {scale * w, scale * x, scale * y, scale * z};
  float j_p_uq[9];
  dfUnitQuaternionInvRotatePoint(unit_q, p1, p, j_p_q ? j_p_uq : 0, r_inv);

  if (j_p_q) {
    Eigen::Matrix4f j_uq_q(Eigen::Matrix4f::Zero());
    float scale3 = scale * scale * scale;
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) j_uq_q(i, j) = -scale3 * q[i] * q[j];
    for (int i = 0; i < 4; i++) j_uq_q(i, i) += scale;

    Eigen::Map<Eigen::Matrix<float, 3, 4, Eigen::RowMajor>> m_j_p_q(j_p_q);
    Eigen::Map<const Eigen::Matrix<float, 3, 4, Eigen::RowMajor>> m_j_p_uq(
        j_p_uq);
    m_j_p_q = m_j_p_uq * j_uq_q;
  }
}

void dfUnitQuaternionInvRotatePoint(const float* q, const float* p1, float* p,
                                    float* j_p_q, float* r_inv) {
  const float w = q[0];
  const float x = -q[1];
  const float y = -q[2];
  const float z = -q[3];
  const float ww = w * w;
  const float wx = w * x;
  const float wy = w * y;
  const float wz = w * z;
  const float xx = -x * x;
  const float xy = x * y;
  const float xz = x * z;
  const float yy = -y * y;
  const float yz = y * z;
  const float zz = -z * z;

  if (r_inv) {
#define _R(row, col) r_inv[(row)*3 + (col)]  // NOLINT
    _R(0, 0) = ww - xx + yy + zz;            // NOLINT
    _R(0, 1) = 2 * (xy + wz);                // NOLINT
    _R(0, 2) = 2 * (xz - wy);                // NOLINT

    _R(1, 0) = 2 * (xy - wz);      // NOLINT
    _R(1, 1) = ww + xx - yy + zz;  // NOLINT
    _R(1, 2) = 2 * (yz + wx);      // NOLINT

    _R(2, 0) = 2 * (xz + wy);      // NOLINT
    _R(2, 1) = 2 * (yz - wx);      // NOLINT
    _R(2, 2) = ww + xx + yy - zz;  // NOLINT
#undef _R                          // NOLINT
  }

  const float a = p1[0], b = p1[1], c = p1[2];
  p[0] = 2 * ((yy + zz) * a + (xy + wz) * b + (xz - wy) * c) + a;
  p[1] = 2 * ((xy - wz) * a + (xx + zz) * b + (yz + wx) * c) + b;
  p[2] = 2 * ((xz + wy) * a + (yz - wx) * b + (xx + yy) * c) + c;

  if (j_p_q) {
#define _J3x4Get(ptr, row, col) ((ptr)[row * 4 + col])         // NOLINT
    _J3x4Get(j_p_q, 0, 0) = 2 * (z * b - y * c);               // NOLINT
    _J3x4Get(j_p_q, 0, 1) = 2 * (y * b + z * c);               // NOLINT
    _J3x4Get(j_p_q, 0, 2) = 2 * (-2 * y * a + x * b - w * c);  // NOLINT
    _J3x4Get(j_p_q, 0, 3) = 2 * (-2 * z * a + w * b + x * c);  // NOLINT

    _J3x4Get(j_p_q, 1, 0) = 2 * (-z * a + x * c);              // NOLINT
    _J3x4Get(j_p_q, 1, 1) = 2 * (y * a - 2 * x * b + w * c);   // NOLINT
    _J3x4Get(j_p_q, 1, 2) = 2 * (x * a + z * c);               // NOLINT
    _J3x4Get(j_p_q, 1, 3) = 2 * (-w * a - 2 * z * b + y * c);  // NOLINT

    _J3x4Get(j_p_q, 2, 0) = 2 * (y * a - x * b);              // NOLINT
    _J3x4Get(j_p_q, 2, 1) = 2 * (z * a - w * b - 2 * x * c);  // NOLINT
    _J3x4Get(j_p_q, 2, 2) = 2 * (w * a + z * b - 2 * y * c);  // NOLINT
    _J3x4Get(j_p_q, 2, 3) = 2 * (z * a + y * b);              // NOLINT
#undef _J3x4Get                                               // NOLINT
  }
}

void dfQuaternionMultiply(const float* q1, const float* q2, float* q12,
                          float* j_q12_q1, float* j_q12_q2) {
  const float w1 = q1[0], x1 = q1[1], y1 = q1[2], z1 = q1[3];
  const float w2 = q2[0], x2 = q2[1], y2 = q2[2], z2 = q2[3];

  q12[0] = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2;
  q12[1] = x1 * w2 + w1 * x2 - z1 * y2 + y1 * z2;
  q12[2] = y1 * w2 + z1 * x2 + w1 * y2 - x1 * z2;
  q12[3] = z1 * w2 - y1 * x2 + x1 * y2 + w1 * z2;

#define _J4x4Get(ptr, row, col) ((ptr)[row * 4 + col])  // NOLINT
  if (j_q12_q1) {                                       // NOLINT
    _J4x4Get(j_q12_q1, 0, 0) = w2;                      // NOLINT
    _J4x4Get(j_q12_q1, 0, 1) = -x2;                     // NOLINT
    _J4x4Get(j_q12_q1, 0, 2) = -y2;                     // NOLINT
    _J4x4Get(j_q12_q1, 0, 3) = -z2;                     // NOLINT
                                                        // NOLINT
    _J4x4Get(j_q12_q1, 1, 0) = x2;                      // NOLINT
    _J4x4Get(j_q12_q1, 1, 1) = w2;                      // NOLINT
    _J4x4Get(j_q12_q1, 1, 2) = z2;                      // NOLINT
    _J4x4Get(j_q12_q1, 1, 3) = -y2;                     // NOLINT
                                                        // NOLINT
    _J4x4Get(j_q12_q1, 2, 0) = y2;                      // NOLINT
    _J4x4Get(j_q12_q1, 2, 1) = -z2;                     // NOLINT
    _J4x4Get(j_q12_q1, 2, 2) = w2;                      // NOLINT
    _J4x4Get(j_q12_q1, 2, 3) = x2;                      // NOLINT
                                                        // NOLINT
    _J4x4Get(j_q12_q1, 3, 0) = z2;                      // NOLINT
    _J4x4Get(j_q12_q1, 3, 1) = y2;                      // NOLINT
    _J4x4Get(j_q12_q1, 3, 2) = -x2;                     // NOLINT
    _J4x4Get(j_q12_q1, 3, 3) = w2;                      // NOLINT
  }                                                     // NOLINT
                                                        // NOLINT
  if (j_q12_q2) {                                       // NOLINT
    _J4x4Get(j_q12_q2, 0, 0) = w1;                      // NOLINT
    _J4x4Get(j_q12_q2, 0, 1) = -x1;                     // NOLINT
    _J4x4Get(j_q12_q2, 0, 2) = -y1;                     // NOLINT
    _J4x4Get(j_q12_q2, 0, 3) = -z1;                     // NOLINT
                                                        // NOLINT
    _J4x4Get(j_q12_q2, 1, 0) = x1;                      // NOLINT
    _J4x4Get(j_q12_q2, 1, 1) = w1;                      // NOLINT
    _J4x4Get(j_q12_q2, 1, 2) = -z1;                     // NOLINT
    _J4x4Get(j_q12_q2, 1, 3) = y1;                      // NOLINT
                                                        // NOLINT
    _J4x4Get(j_q12_q2, 2, 0) = y1;                      // NOLINT
    _J4x4Get(j_q12_q2, 2, 1) = z1;                      // NOLINT
    _J4x4Get(j_q12_q2, 2, 2) = w1;                      // NOLINT
    _J4x4Get(j_q12_q2, 2, 3) = -x1;                     // NOLINT
                                                        // NOLINT
    _J4x4Get(j_q12_q2, 3, 0) = z1;                      // NOLINT
    _J4x4Get(j_q12_q2, 3, 1) = -y1;                     // NOLINT
    _J4x4Get(j_q12_q2, 3, 2) = x1;                      // NOLINT
    _J4x4Get(j_q12_q2, 3, 3) = w1;                      // NOLINT
  }                                                     // NOLINT
#undef _J4x4Get                                         // NOLINT
}

void dfQuaternionConjugate(const float* q, float* qc, float* j_qc_q) {
  qc[0] = q[0];
  qc[1] = -q[1];
  qc[2] = -q[2];
  qc[3] = -q[3];

  if (j_qc_q) {
    memset(j_qc_q, 0.0, sizeof(float) * 16);
    j_qc_q[0] = 1.0;
    j_qc_q[5] = -1.0;
    j_qc_q[10] = -1.0;
    j_qc_q[15] = -1.0;
  }
}

}  // namespace avp
