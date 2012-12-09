#ifndef CSPACE_LEARNING_DATASET_H
#define CSPACE_LEARNING_DATASET_H

#include <vector>
#include <boost/dynamic_bitset.hpp>
#include <assert.h>

namespace cspace_learning
{

class Dataset
{
public:

  typedef unsigned int Key;
  typedef const float* Value;

  struct DataItem
  {
    DataItem()
    {
      tag = 0;
    }

    std::vector<float> q;
    int tag;
  };


  void free()
  {
    data_.clear();
  }

  Dataset(unsigned int dim) { dim_ = dim; }

  void insert(float* q, unsigned int dim, int tag = 0)
  {
    if(dim != dim_)
    {
      std::cerr << "dimension not match" << std::endl;
      return;
    }

    DataItem item;
    item.q.resize(dim_);
    for(unsigned int i = 0; i < dim; ++i)
    {
      item.q[i] = q[i];
    }
    item.tag = tag;

    data_.push_back(item);
  }

  const float* operator [] (unsigned int i) const
  {
    assert(i < data_.size());
    return &(data_[i].q[0]);
  }

  float* operator [] (unsigned int i)
  {
    assert(i < data_.size());
    return &(data_[i].q[0]);
  }

  int tag(unsigned int i) const
  {
    assert(i < data_.size());
    return data_[i].tag;
  }

  int& tag(unsigned int i)
  {
    assert(i < data_.size());
    return data_[i].tag;
  }

  unsigned int dim() const
  {
    return dim_;
  }

  unsigned int size() const
  {
    return data_.size();
  }

  class Accessor
  {
    const Dataset& dataset_;
    boost::dynamic_bitset<> flags_;

  public:
    typedef unsigned int Key;
    typedef const float* Value;

    Accessor(const Dataset& dataset) : dataset_(dataset), flags_(dataset.size()) {}

    void reset()
    {
      flags_.resize(dataset_.size());
      flags_.reset();
    }

    bool mark(Key key)
    {
      if(flags_[key]) return false;
      flags_.set(key);
      return true;
    }

    const float* operator () (Key key)
    {
      return dataset_[key];
    }

    const Dataset& dataset() const
    {
      return dataset_;
    }
  };

protected:
  std::vector<DataItem> data_;
  unsigned int dim_;
};

}

#endif
