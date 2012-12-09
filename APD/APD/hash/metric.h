#ifndef CSPACE_LEARNING_METRIC_H
#define CSPACE_LEARNING_METRIC_H

#include <functional>
#include <hash/common.h>

namespace cspace_learning
{

namespace metric
{

template<typename T>
class l1 : public std::binary_function<const T*, const T*, float>
{
  unsigned int dim_;
public:
  l1(unsigned int dim) : dim_(dim) {}

  unsigned int dim() const { return dim_; }

  float operator () (const T* elem1, const T* elem2) const
  {
    float r = 0.0;
    for(unsigned int i = 0; i < dim_; ++i)
      r += fabs(elem1[i] - elem2[i]);
    return r;
  }
};

template<typename T>
class l2 : public std::binary_function<const T*, const T*, float>
{
  unsigned int dim_;
public:
  l2(unsigned int dim) : dim_(dim) {}

  unsigned int dim() const { return dim_; }

  float operator () (const T* elem1, const T* elem2) const
  {
    float r = 0.0;
    for(unsigned int i = 0; i < dim_; ++i)
      r += sqr(elem1[i] - elem2[i]);
    return sqrt(r);
  }
};

template<>
class l2<float> : public std::binary_function<const float*, const float*, float>
{
  unsigned int dim_;
public:
  l2(unsigned int dim) : dim_(dim) {}

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
      case 7: r6 = sqr(a[6] - b[6]);
      case 6: r5 = sqr(a[5] - b[5]);
      case 5: r4 = sqr(a[4] - b[4]);
      case 4: r3 = sqr(a[3] - b[3]);
      case 3: r2 = sqr(a[2] - b[2]);
      case 2: r1 = sqr(a[1] - b[1]);
      case 1: r0 = sqr(a[0] - b[0]);
    }

    a = aa; b = bb;

    for(; a < end_a; a += 8, b += 8)
    {
#ifdef __GNUC__
        __builtin_prefetch(a + 32, 0, 3);
        __builtin_prefetch(b + 32, 0, 0);
#endif
        r += r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7;
        r0 = sqr(a[0] - b[0]);
        r1 = sqr(a[1] - b[1]);
        r2 = sqr(a[2] - b[2]);
        r3 = sqr(a[3] - b[3]);
        r4 = sqr(a[4] - b[4]);
        r5 = sqr(a[5] - b[5]);
        r6 = sqr(a[6] - b[6]);
        r7 = sqr(a[7] - b[7]);
    }

    r += r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7;

    return sqrt(r);
  }
};

template<typename T>
class l2sqr : public std::binary_function<const T*, const T*, float>
{
  unsigned int dim_;
public:
  l2sqr(unsigned int dim) : dim_(dim) {}

  unsigned int dim() const { return dim_; }

  float operator () (const T* elem1, const T* elem2) const
  {
    float r = 0.0;
    for(unsigned int i = 0; i < dim_; ++i)
      r += sqr(elem1[i] - elem2[i]);
    return r;
  }
};

template<>
class l2sqr<float> : public std::binary_function<const float*, const float*, float>
{
  unsigned int dim_;
public:
  l2sqr(unsigned int dim) : dim_(dim) {}

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
      case 7: r6 = sqr(a[6] - b[6]);
      case 6: r5 = sqr(a[5] - b[5]);
      case 5: r4 = sqr(a[4] - b[4]);
      case 4: r3 = sqr(a[3] - b[3]);
      case 3: r2 = sqr(a[2] - b[2]);
      case 2: r1 = sqr(a[1] - b[1]);
      case 1: r0 = sqr(a[0] - b[0]);
    }

    a = aa; b = bb;

    for(; a < end_a; a += 8, b += 8)
    {
#ifdef __GNUC__
        __builtin_prefetch(a + 32, 0, 3);
        __builtin_prefetch(b + 32, 0, 0);
#endif
        r += r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7;
        r0 = sqr(a[0] - b[0]);
        r1 = sqr(a[1] - b[1]);
        r2 = sqr(a[2] - b[2]);
        r3 = sqr(a[3] - b[3]);
        r4 = sqr(a[4] - b[4]);
        r5 = sqr(a[5] - b[5]);
        r6 = sqr(a[6] - b[6]);
        r7 = sqr(a[7] - b[7]);
    }

    r += r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7;

    return r;
  }
};

template<typename T>
class max : public std::binary_function<const T*, const T*, float>
{
  unsigned dim_;
public:
  max(unsigned int dim) : dim_(dim) {}

  unsigned int dim() const { return dim_; }

  float operator () (const T* elem1, const T* elem2) const
  {
    double r = std::numeric_limits<T>::min();
    for(unsigned int i = 0; i < dim_; ++i)
    {
      r += std::max(r, fabs(elem1[i] - elem2[i]));
    }
    return (float)sqrtf(r);
  }
};

struct basic_hamming
{
  static unsigned int char_bit_cnt[];

  template<typename B>
  unsigned __hamming(B a, B b)
  {
    B c = a ^ b;
    unsigned char *p = reinterpret_cast<unsigned char *>(&c);
    unsigned int r = 0;
    for(unsigned int i = 0; i < sizeof(B); i++)
    {
      r += char_bit_cnt[*p++];
    }
    return r;
  }

  unsigned int __hamming(unsigned char c1, unsigned char c2) const
  {
      return char_bit_cnt[c1 ^ c2];
  }
};

template<typename T>
struct hamming : public std::binary_function<const T*, const T*, float>, public basic_hamming
{
  unsigned int dim_;

public:
  hamming(unsigned int dim) : dim_(dim) {}

  unsigned int dim() const { return dim_; }

  float operator () (const T* elem1, const T* elem2) const
  {
    unsigned int r = 0;
    for(unsigned int i = 0; i < dim_; ++i)
    {
      r += __hamming(elem1[i], elem2[i]);
    }
    return r;
  }
};


}

}

#endif
