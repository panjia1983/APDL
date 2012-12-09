#ifndef CSPACE_LEARNING_EVAL_H
#define CSPACE_LEARNING_EVAL_H

#include "hash/metric.h"

namespace cspace_learning
{

// evaluate a database item
namespace eval
{

template<typename Metric>
struct PointEval : public std::unary_function<typename Metric::first_argument_type, float>
{
  typedef typename Metric::first_argument_type Value;

  PointEval(const Metric& metric) : metric_(metric)
  {
  }

  void reset(Value ref)
  {
    ref_ = ref;
  }

  float operator () (Value elem)
  {
    return metric_(ref_, elem);
  }

  unsigned int dim() const
  {
    return metric_.dim();
  }

protected:
  Value ref_; // reference/query point
  Metric metric_;
};


struct PointToLineEval : public std::unary_function<const float*, const float>
{
  typedef const float* Value;
  PointToLineEval(unsigned int dim)
  {
    dim_ = dim;
  }

  void reset(Value v, Value a)
  {
    v_ = v;
    a_ = a;
  }

  float operator () (Value elem)
  {
    float dist = 0, tmp = 0;
    for(unsigned int d = 0; d < dim_; ++d)
    {
      dist += sqr(elem[d] - a_[d]);
      tmp += (elem[d] - a_[d]) * v_[d];
    }

    return dist - tmp * tmp;
  }

  unsigned int dim() const
  {
    return dim_;
  }

protected:
  Value v_;
  Value a_;
  unsigned int dim_;
};

}


}

#endif
