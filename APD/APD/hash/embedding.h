#ifndef CSPACE_LEARNING_EMBEDDING_H
#define CSPACE_LEARNING_EMBEDDING_H

#include "hash/common.h"
#include <assert.h>

namespace cspace_learning
{



class DotProductConservativeRandomProjector
{
  unsigned int dim_;
  std::vector<float> cache_;
  std::vector<float> samples_;
public:
  DotProductConservativeRandomProjector(unsigned int dim)
  {
    dim_ = dim;
    cache_.resize(dim_);
  }

  ~DotProductConservativeRandomProjector()
  {
  }

  template<typename RNG>
  void init(unsigned int N, RNG& rng)
  {
    samples_.resize(N);
    boost::variate_generator<RNG &, Uniform> generator(rng, Uniform(0, 1));
    for(unsigned int i = 0; i < N; ++i)
      samples_[i] = generator();
  }

  SparseVector operator () (const float* in)
  {
    SparseVector v;
    cache_[0] = in[0] * in[0];
    for(unsigned int i = 1; i < dim_; ++i)
    {
      cache_[i] = in[i] * in[i] + cache_[i - 1];
    }

    float sum = cache_[dim_ - 1];

    for(unsigned int i = 0; i < dim_; ++i)
    {
      cache_[i] /= sum;
    }
    cache_[dim_ - 1] = 1;

    for(unsigned int i = 0; i < samples_.size(); ++i)
    {
      unsigned int pos = std::lower_bound(cache_.begin(), cache_.end(), samples_[i]) - cache_.begin();

      if(v.find(pos) != v.end())
      {
        v[pos] += sum / in[pos];
      }
      else
      {
        v[pos] = sum / in[pos];
      }
    }

    for(SparseVector::iterator it = v.begin(); it != v.end(); ++it)
    {
      it->second /= samples_.size();
    }

    return v;
  }
};

class DistanceConservativeRandomProjector
{
  unsigned int dim_in_;
  unsigned int dim_out_;
  float* Pmatrix_;
  float* cache_;

public:
  DistanceConservativeRandomProjector(unsigned int dim_in, unsigned int dim_out)
  {
    dim_in_ = dim_in;
    dim_out_ = dim_out;
    Pmatrix_ = new float[dim_in_ * dim_out_];
    cache_ = new float[dim_out_];
  }

  ~DistanceConservativeRandomProjector()
  {
    delete [] Pmatrix_;
    delete [] cache_;
  }

  void reset(unsigned int dim_in, unsigned int dim_out)
  {
    dim_in_ = dim_in;
    dim_out_ = dim_out;
    delete [] Pmatrix_;
    delete [] cache_;
    Pmatrix_ = new float[dim_in_ * dim_out_];
    cache_ = new float[dim_out_];
  }

  template<typename RNG>
  void init(RNG& rng)
  {
    boost::variate_generator<RNG&, Gaussian> generator(rng, Gaussian());

    float ratio = 1.0 / sqrt((float)dim_out_);
    for(unsigned int i = 0; i < dim_out_; ++i)
    {
      for(unsigned int j = 0; j < dim_in_; ++j)
      {
        Pmatrix_[i * dim_in_ + j] = ratio * generator();
      }
    }
  }

  float* operator () (const float* in)
  {
    dotprod<float> dot(dim_in_);
    for(unsigned int i = 0; i < dim_out_; ++i)
      cache_[i] = dot(&Pmatrix_[i * dim_in_], in);
    return cache_;
  }
};

class DistanceConservativeSparseRandomProjector
{
  unsigned int dim_in_;
  unsigned int dim_out_;
  float* Pmatrix_;
  float* cache_;
  float s_;

public:
  DistanceConservativeSparseRandomProjector(unsigned int dim_in, unsigned int dim_out, float s)
  {
    dim_in_ = dim_in;
    dim_out_ = dim_out;
    Pmatrix_ = new float[dim_in_ * dim_out_];
    cache_ = new float[dim_out_];
    s_ = s;
  }

  ~DistanceConservativeSparseRandomProjector()
  {
    delete [] Pmatrix_;
    delete [] cache_;
  }

  void reset(unsigned int dim_in, unsigned int dim_out)
  {
    dim_in_ = dim_in;
    dim_out_ = dim_out;
    delete [] Pmatrix_;
    delete [] cache_;
    Pmatrix_ = new float[dim_in_ * dim_out_];
    cache_ = new float[dim_out_];
  }

  float s() const { return s_; }
  float& s() { return s_; }

  template<typename RNG>
  void init(RNG& rng)
  {
    boost::variate_generator<RNG&, Uniform> generator(rng, Uniform(0, 1));

    float ratio = 1.0 / sqrt((float)dim_out_);
    float inv_s = 1.0 / s_;
    float sqrts = sqrt(s_);
    for(unsigned int i = 0; i < dim_out_; ++i)
    {
      for(unsigned int j = 0; j < dim_in_; ++j)
      {
        float v = generator();
        if(v < 0.5 * inv_s)
          Pmatrix_[i * dim_in_ + j] = ratio * sqrts;
        else if(v > 1 - 0.5 * inv_s)
          Pmatrix_[i * dim_in_ + j] = -ratio * sqrts;
        else
          Pmatrix_[i * dim_in_ + j] = 0;
      }
    }
  }

  float* operator () (const float* in)
  {
    for(unsigned int i = 0; i < dim_out_; ++i)
    {
      float tmp = 0;
      for(unsigned int j = 0; j < dim_in_; ++j)
      {
        if(Pmatrix_[i * dim_in_ + j] != 0)
          tmp += in[j] * Pmatrix_[i * dim_in_ + j];
      }
      cache_[i] = tmp;
    }
    return cache_;
  }
};


// MUST set t_ before using it.
class Embedding
{
  unsigned int dim_;
  unsigned int span_;
  float t_;
  float* cache_;
public:
  Embedding(unsigned int dim) : dim_(dim), span_(dim + 2)
  {
    cache_ = new float[span_ * span_];
  }

  ~Embedding()
  {
    delete [] cache_;
  }

  void reset(unsigned int dim)
  {
    dim_ = dim;
    span_ = dim + 2;
    delete [] cache_;
    cache_ = new float[span_ * span_];
  }

  float& t()
  {
    return t_;
  }

  float t() const
  {
    return t_;
  }

  float* embedded_datapoint(const float* in);

  float* embedded_line(const float* v, const float* a);
};

class ReducedEmbedding
{
  unsigned int dim_;
  unsigned int span_;
  float t_;
  float* cache_;
public:
  ReducedEmbedding(unsigned int dim) : dim_(dim), span_(dim + 2)
  {
    cache_ = new float[span_ * (span_ + 1) / 2];
  }

  ~ReducedEmbedding()
  {
    delete [] cache_;
  }

  void reset(unsigned int dim)
  {
    dim_ = dim;
    span_ = dim + 2;
    delete [] cache_;
    cache_ = new float[span_ * (span_ + 1) / 2];
  }

  float& t()
  {
    return t_;
  }

  float t() const
  {
    return t_;
  }

  float* embedded_datapoint(const float* in);

  float* embedded_line(const float* v, const float* a);
};



}

#endif
