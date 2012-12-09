#ifndef CSPACE_LEARNING_TOPK_H
#define CSPACE_LEARNING_TOPK_H

#include <limits>
#include <vector>
#include <algorithm>
#include <hash/common.h>
#include <hash/metric.h>

namespace cspace_learning
{

/**
  * The entry stored in the top-k data structure.  The class Topk is implemented
  * as a heap of TopkEntry.
  */
template<typename Key>
struct TopkEntry
{
  Key key;
  float dist;
  bool match(const TopkEntry& other) const
  {
    return key == other.key;
  }

  bool match(Key e) const
  {
    return key == e;
  }

  TopkEntry(Key e, float d) : key(e), dist(d) {}
  TopkEntry() : dist(std::numeric_limits<float>::max()) {}

  void reset()
  {
    dist = std::numeric_limits<float>::max();
    key = -1;
  }

  friend bool operator < (const TopkEntry& e1, const TopkEntry& e2)
  {
    return e1.dist < e2.dist;
  }
};

template<typename Key>
class Topk : public std::vector<TopkEntry<Key> >
{
  unsigned int knn_k_;
  float radius_;
  float threshold_;

public:
  typedef TopkEntry<Key> Element;
  typedef typename std::vector<TopkEntry<Key> > Base;

  Topk() {}
  Topk(unsigned int k) { reset(k); }
  ~Topk() {}

  unsigned int nTopk() const
  {
    unsigned int n = 0;
    for(typename Base::const_iterator i = this->begin(); i != this->end(); ++i)
    {
      if(i->key != -1)
        n++;
      else
        break;
    }

    return n;
  }

  void reset(unsigned int k, float r = std::numeric_limits<float>::max())
  {
    knn_k_ = k;
    radius_ = threshold_ = r;
    this->resize(k);

    for(typename Base::iterator i = this->begin(); i != this->end(); ++i)
      i->reset();
  }

  void reset(unsigned int k, Key key, float r = std::numeric_limits<float>::max())
  {
    knn_k_ = k;
    radius_ = threshold_ = r;
    this->resize(k);

    for(typename Base::iterator i = this->begin(); i != this->end(); ++i)
    {
      i->reset();
      i->key = key;
    }
  }

  void reset(float r)
  {
    knn_k_ = 0;
    radius_ = threshold_ = r;
    this->clear();
  }

  float threshold() const
  {
    return threshold_;
  }

  Topk& operator << (const Element& t)
  {
    if(t.dist >= threshold_) return *this;

    // R-NN
    if(knn_k_ == 0)
    {
      this->push_back(t);
      return *this;
    }

    // K-NN
    unsigned int i = this->size() - 1;
    unsigned int j;

    for(;;)
    {
      if(i == 0) break;
      j = i - 1;
      if(this->at(j).match(t)) return *this;
      if(this->at(j) < t) break;
      i = j;
    }

    // i is the place to insert
    j = this->size() - 1;
    for(;;)
    {
      if(j == i) break;
      this->at(j) = this->at(j-1);
      --j;
    }
    this->at(i) = t;
    threshold_ = this->back().dist;
    return *this;
  }

  float recall(const Topk<Key>& topk) const
  {
    unsigned int matched = 0;
    for(typename Base::const_iterator it1 = this->begin(); it1 != this->end(); ++it1)
    {
      for(typename Base::const_iterator it2 = topk.begin(); it2 != topk.end(); ++it2)
      {
        if(it1->match(*it2))
        {
          matched++;
          break;
        }
      }
    }
    return float(matched + 1) / (float)(this->size() + 1);
  }

  float error(const Topk<Key>& topk) const
  {

    float ret = 0.0;
    typename Base::const_iterator it1 = this->begin();
    typename Base::const_iterator it2 = topk.begin();

    unsigned int cnt = 0;
    for(; (it1 != this->end()) && (it2 != topk.end()); ++it1, ++it2, ++cnt)
    {
      ret += (it1->dist / it2->dist);
    }

    ret /= cnt;

    return ret;
  }


  unsigned int getKNN_K() const
  {
    return knn_k_;
  }
};

template<typename Accessor, typename Query, typename Eval>
class TopkScanner
{
public:
  typedef typename Accessor::Key Key;
  typedef typename Accessor::Value Value;

  TopkScanner(const Accessor& accessor, Eval* eval, unsigned int k, float r = std::numeric_limits<float>::max()) : accessor_(accessor), eval_(eval), knn_k_(k), radius_(r)
  {
  }

  void reset(Query query)
  {
    query_ = query;
    accessor_.reset();
    topk_.reset(knn_k_, radius_);
    cnt_ = 0;
  }

  unsigned int cnt() const
  {
    return cnt_;
  }

  const Topk<Key>& topk() const
  {
    return topk_;
  }

  Topk<Key>& topk()
  {
    return topk_;
  }

  void operator () (Key key)
  {
    if(accessor_.mark(key))
    {
      ++cnt_;
      topk_ << typename Topk<Key>::Element(key, (*eval_)(accessor_(key)));
    }
  }

  Query query()
  {
    return query_;
  }

protected:
  Accessor accessor_;
  Eval* eval_;
  unsigned int knn_k_;
  float radius_;
  Topk<Key> topk_;
  Query query_;
  unsigned int cnt_;
};

template<typename Accessor, typename Eval>
class TopkScanner<Accessor, const float*, Eval>
{
public:
  typedef typename Accessor::Key Key;
  typedef const float* Value;

  TopkScanner(const Accessor& accessor, Eval* eval, unsigned int k, float r = std::numeric_limits<float>::max())
    : accessor_(accessor), knn_k_(k), radius_(r)
  {
    eval_ = eval;
    dim_ = eval_->dim();
  }

  void reset(const float* query)
  {
    query_ = query;
    accessor_.reset();
    topk_.reset(knn_k_, radius_);
    cnt_ = 0;
  }

  unsigned int cnt() const
  {
    return cnt_;
  }

  const Topk<Key>& topk() const
  {
    return topk_;
  }

  Topk<Key>& topk()
  {
    return topk_;
  }

  void operator () (Key key)
  {
    if(accessor_.mark(key))
    {
      ++cnt_;
      topk_ << typename Topk<Key>::Element(key, (*eval_)(accessor_(key)));
    }
  }

  const float* query()
  {
    return query_;
  }

private:
  Accessor accessor_;
  unsigned int dim_;
  unsigned int knn_k_;
  float radius_;
  Topk<Key> topk_;
  const float* query_;
  unsigned int cnt_;
  Eval* eval_;
};

}

#endif
