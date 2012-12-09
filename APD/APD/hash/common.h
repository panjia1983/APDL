#ifndef CSPACE_LEARNING_COMMON_H
#define CSPACE_LEARNING_COMMON_H

#include <boost/random.hpp>
#include <map>

namespace cspace_learning
{

typedef boost::mt19937 DefaultRng;

// Some of the frequently used distributions.

/// Gaussian distribution.
typedef boost::normal_distribution<float> Gaussian;
/// Cauchy distribution.
typedef boost::cauchy_distribution<float> Cauchy;
/// Uniform distribution.
typedef boost::uniform_real<float> Uniform;
/// Uniform distribution with int values.
typedef boost::uniform_int<int> UniformInt;
/// Uniform distribution with unsigned values.
typedef boost::uniform_int<unsigned> UniformUnsigned;


template <typename T> T sqr (const T &x) { return x * x; }


typedef std::map<unsigned int, float> SparseVector;

template<typename T>
class dotprod : public std::binary_function<const T*, const T*, float>
{
  unsigned int dim_;
public:
  dotprod(unsigned int dim) : dim_(dim) {}

  unsigned int dim() const { return dim_; }

  float operator () (const T* elem1, const T* elem2) const
  {
    float r = 0.0;
    for(unsigned int i = 0; i < dim_; ++i)
      r += elem1[i] * elem2[i];
    return r;
  }

  float operator () (const SparseVector& v1, const SparseVector& v2) const
  {
    float r = 0.0;
    typename SparseVector::const_iterator pos;
    for(typename SparseVector::const_iterator it = v1.begin(); it != v1.end(); ++it)
    {
      unsigned int id = it->first;
      if((pos = v2.find(id)) != v2.end())
      {
        r += it->second * pos->second;
      }
    }
    return r;
  }

  float operator () (const SparseVector& v1, const T* v2) const
  {
    float r = 0.0;
    for(typename SparseVector::const_iterator it = v1.begin(); it != v1.end(); ++it)
    {
      unsigned int id = it->first;
      r += it->second * v2[id];
    }
    return r;
  }
};


template<>
class dotprod<float> : public std::binary_function<const float*, const float*, float>
{
  unsigned int dim_;
public:
  dotprod(unsigned int dim) : dim_(dim) {}

  unsigned int dim() const { return dim_; }

  float operator () (const float* elem1, const float* elem2) const
  {
    unsigned int d = dim_ & ~unsigned(7);
    const float *aa = elem1, *end_a = aa + d;
    const float *bb = elem2, *end_b = bb + d;
#ifdef __GNUC__
    __builtin_prefetch(aa, 0, 3);
    __builtin_prefetch(bb, 0, 0);
#endif
    float r = 0.0;
    float r0, r1, r2, r3, r4, r5, r6, r7;

    const float* a = end_a, *b = end_b;

    r0 = r1 = r2 = r3 = r4 = r5 = r6 = r7 = 0.0;

    switch(dim_ & 7)
    {
      case 7: r6 = a[6] * b[6];
      case 6: r5 = a[5] * b[5];
      case 5: r4 = a[4] * b[4];
      case 4: r3 = a[3] * b[3];
      case 3: r2 = a[2] * b[2];
      case 2: r1 = a[1] * b[1];
      case 1: r0 = a[0] * b[0];
    }

    a = aa; b = bb;

    for(; a < end_a; a += 8, b += 8)
    {
#ifdef __GNUC__
        __builtin_prefetch(a + 32, 0, 3);
        __builtin_prefetch(b + 32, 0, 0);
#endif
        r += r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7;
        r0 = a[0] * b[0];
        r1 = a[1] * b[1];
        r2 = a[2] * b[2];
        r3 = a[3] * b[3];
        r4 = a[4] * b[4];
        r5 = a[5] * b[5];
        r6 = a[6] * b[6];
        r7 = a[7] * b[7];
    }

    r += r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7;

    return r;
  }

  float operator () (const SparseVector& v1, const SparseVector& v2) const
  {
    float r = 0.0;
    SparseVector::const_iterator pos;
    for(SparseVector::const_iterator it = v1.begin(); it != v1.end(); ++it)
    {
      unsigned int id = it->first;
      if((pos = v2.find(id)) != v2.end())
      {
        r += it->second * pos->second;
      }
    }
    return r;
  }

  float operator () (const SparseVector& v1, const float* v2) const
  {
    float r = 0.0;
    for(SparseVector::const_iterator it = v1.begin(); it != v1.end(); ++it)
    {
      unsigned int id = it->first;
      r += it->second * v2[id];
    }
    return r;
  }
};


struct Line
{
  Line(const float* v_, const float* a_, unsigned int dim_)
  {
    v = v_;
    a = a_;
    dim = dim_;
  }

  Line() {}

  const float* v;
  const float* a;
  unsigned int dim;
};


}

#endif
