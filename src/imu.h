#ifndef AVP_IMU_H_
#define AVP_IMU_H_

namespace avp {
class IMU {
public:
  IMU() = default;
  ~IMU() = default;
  static void PreInteraciton();
  static void Interaciton();
};
} // namespace avp
#endif
