#ifndef AVP_ODOMETRY_H_
#define AVP_ODOMETRY_H_
#include "frame.h"

namespace avp {
class Odometry {
public:
    Odometry();
    ~Odometry();
    void process(Frame* frame);
private:
    class Impl;
    Impl* p_impl;
};
} // namespace avp
#endif