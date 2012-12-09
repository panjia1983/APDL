#include "hash/embedding.h"

namespace cspace_learning
{
float* Embedding::embedded_datapoint(const float* in)
{
  for(unsigned int i = 0; i < dim_; ++i)
  {
    for(unsigned int j = 0; j < dim_; ++j)
    {
      cache_[i * span_ + j] = in[i] * in[j];
    }
  }

  for(unsigned int i = 0; i < dim_; ++i)
  {
    cache_[i * span_ + dim_] = in[i];
    cache_[i * span_ + dim_ + 1] = in[i] * t_;
  }

  for(unsigned int i = 0; i < dim_; ++i)
  {
    cache_[dim_ * span_ + i] = in[i];
    cache_[(dim_ + 1) * span_ + i] = in[i] * t_;
  }

  cache_[dim_ * span_ + dim_] = 1;
  cache_[dim_ * span_ + dim_ + 1] = t_;
  cache_[(dim_ + 1) * span_ + dim_] = t_;
  cache_[(dim_ + 1) * span_ + dim_ + 1] = t_ * t_;

  return cache_;
}

float* Embedding::embedded_line(const float* v, const float* a)
{
  for(unsigned int i = 0; i < dim_; ++i)
  {
    for(unsigned int j = 0; j < dim_; ++j)
    {
      if(i == j)
        cache_[i * span_ + j] = - 1 + v[i] * v[j];
      else
        cache_[i * span_ + j] = v[i] * v[j];
    }
  }

  float vdota = 0, adota = 0;
  for(unsigned int i = 0; i < dim_; ++i)
  {
    vdota += v[i] * a[i];
    adota += a[i] * a[i];
  }

  for(unsigned int i = 0; i < dim_; ++i)
  {
    float tmp = a[i] - vdota * v[i];
    cache_[i * span_ + dim_] = tmp;
    cache_[dim_ * span_ + i] = tmp;
  }

  cache_[dim_ * span_ + dim_] = vdota * vdota - adota;

  for(unsigned int i = 0; i < dim_ + 1; ++i)
  {
    cache_[i * span_ + dim_ + 1] = 0;
    cache_[(dim_ + 1) * span_ + i] = 0;
  }

  cache_[(dim_ + 1) * span_ + dim_ + 1] = 0;

  return cache_;
}

float* ReducedEmbedding::embedded_datapoint(const float* in)
{
  unsigned int offset = 0;
  float inv_sqrt2 = 1.0 / sqrt(2.0);
  for(unsigned int i = 0; i < span_; ++i)
  {
    for(unsigned int j = i; j < span_; ++j)
    {
      if(i == j)
      {
        if(i < dim_)
          cache_[offset] = in[i] * in[j] * inv_sqrt2;
        else if(i == dim_)
          cache_[offset] = inv_sqrt2;
        else
          cache_[offset] = t_ * t_ * inv_sqrt2;
      }
      else
      {
        if(j < dim_)
          cache_[offset] = in[i] * in[j];
        else if(j == dim_)
          cache_[offset] = in[i];
        else if(j == dim_ + 1)
        {
          if(i < dim_)
            cache_[offset] = in[i] * t_;
          else
            cache_[offset] = t_;
        }
      }
      offset++;
    }
  }

  return cache_;
}

float* ReducedEmbedding::embedded_line(const float* v, const float* a)
{
  unsigned int offset = 0;
  float inv_sqrt2 = 1.0 / sqrt(2.0);

  float vdota = 0, adota = 0;
  for(unsigned int i = 0; i < dim_; ++i)
  {
    vdota += v[i] * a[i];
    adota += a[i] * a[i];
  }

  for(unsigned int i = 0; i < span_; ++i)
  {
    for(unsigned int j = i; j < span_; ++j)
    {
      if(i == j)
      {
        if(i < dim_)
          cache_[offset] = (- 1 + v[i] * v[j]) * inv_sqrt2;
        else if(i == dim_)
          cache_[offset] = (vdota * vdota - adota) * inv_sqrt2;
        else
          cache_[offset] = 0;
      }
      else
      {
        if(j < dim_)
          cache_[offset] = v[i] * v[j];
        else if(j == dim_)
          cache_[offset] = a[i] - vdota * v[i];
        else
          cache_[offset] = 0;
      }
      offset++;
    }
  }

  return cache_;
}

}
