#ifndef AVP_SOLVER_H_
#define AVP_SOLVER_H_

#include <cmath>

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/SparseQR"

namespace avp {
template <typename _Scalar>
class DenseFunctor {
 public:
  typedef _Scalar Scalar;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> InputType;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> ValueType;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> JacobianType;
  typedef Eigen::ColPivHouseholderQR<JacobianType> QRSolver;

  DenseFunctor(int intputs, int values) : m_inputs(intputs), m_values(values) {}
  ~DenseFunctor() {}
  virtual int df(const InputType& x, JacobianType* fjac) = 0;
  virtual int fdf(const InputType& x, ValueType* err, JacobianType* fjac) = 0;
  virtual int operator()(const InputType& x, ValueType* err) = 0;
  virtual float error(const ValueType& err) = 0;
  int m_inputs, m_values, m_value_dim;
};

template <typename _Scalar>
class SparseFunctor {
 public:
  typedef _Scalar Scalar;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> InputType;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> ValueType;
  typedef Eigen::SparseMatrix<Scalar> JacobianType;
  typedef Eigen::SparseQR<JacobianType, Eigen::AMDOrdering<int>> QRSolver;

  SparseFunctor(int intputs, int values)
      : m_inputs(intputs), m_values(values) {}
  ~SparseFunctor() {}
  virtual int df(const InputType& x, JacobianType* fjac) = 0;
  virtual int fdf(const InputType& x, ValueType* err, JacobianType* fjac) = 0;
  virtual int operator()(const InputType& x, ValueType* err) = 0;
  virtual float error(const ValueType& err) = 0;
  int m_inputs, m_values, m_value_dim;
};

template <typename CostFunction>
class GNSolver {
 public:
  typedef typename CostFunction::InputType InputType;
  typedef typename CostFunction::ValueType ValueType;
  typedef typename CostFunction::JacobianType JacobianType;

  GNSolver(CostFunction cf, int max_iter) : cf_(cf), max_iter_(max_iter) {}
  ~GNSolver() {}
  bool Minimize(InputType* x) {
    bool ret = false;
    float init_err = 0.0;
    ValueType e(cf_.m_values, 1);
    JacobianType j(cf_.m_values, cf_.m_inputs);
    for (size_t iter = 0; iter < max_iter_; iter++) {
      cf_.fdf(*x, &e, &j);
      if (iter == 0) {
        init_err = cf_.error(e);
      }
      JacobianType jtj = j.transpose() * j;
      qr_.compute(jtj);
      InputType delta = qr_.solve(j.transpose() * e);
      (*x) -= qr_.solve(j.transpose() * e);
      if (Eigen::ComputationInfo::Success != qr_.info()) {
        return ret;
      }
    }
    cf_(*x, &e);
    float iter_err = cf_.error(e);
    ret = (iter_err < init_err);
    return ret;
  }

 private:
  CostFunction cf_;
  const int max_iter_;
  typename CostFunction::QRSolver qr_;
};

}  // namespace avp
#endif
