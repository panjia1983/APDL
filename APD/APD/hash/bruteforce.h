#ifndef CSPACE_LEARNING_BRUTEFORCE_H
#define CSPACE_LEARNING_BRUTEFORCE_H

#include "hash/topk.h"
#include "hash/dataset.h"

namespace cspace_learning
{

template<typename Accessor, typename Eval>
class Bruteforce
{
public:

  typedef typename Accessor::Key Key;
  typedef typename Accessor::Value Value;

  Bruteforce(const Accessor& accessor, Eval* eval, unsigned int k, float r = std::numeric_limits<float>::max()) : accessor_(accessor), eval_(eval), knn_k_(k), radius_(r)
  {
  }

  void reset()
  {
    accessor_.reset();
    topk_.reset(knn_k_, radius_);
  }

  const Topk<Key>& topk() const
  {
    return topk_;
  }

  Topk<Key>& topk()
  {
    return topk_;
  }

  void knn()
  {
    if(!eval_)
    {
      std::cerr << "Need to set eval first!" << std::endl;
      return;
    }

    for(unsigned int i = 0; i < accessor_.dataset().size(); ++i)
      topk_ << typename Topk<Key>::Element(i, (*eval_)(accessor_(i)));
  }

protected:
  Accessor accessor_;
  Eval* eval_;
  unsigned int knn_k_;
  float radius_;
  Topk<Key> topk_;
};

}

#endif
