#ifndef CSPACE_LEARNING_MULTIPROBE_LSH_H
#define CSPACE_LEARNING_MULTIPROBE_LSH_H


#include "hash/common.h"
#include "hash/LSH.h"
#include "hash/composite.h"
#include "hash/metric.h"
#include "hash/LSH_index.h"
#include "hash/topk.h"

namespace cspace_learning
{


/**
 * \section ref Reference
 * Wei Dong, Zhe Wang, William Josephson, Moses Charikar, Kai Li. Modeling LSH
 * for Performance Tuning.. To appear in In Proceedings of ACM 17th Conference
 * on Information and Knowledge Management (CIKM). Napa Valley, CA, USA.
 * October 2008.
 *
 * Qin Lv, William Josephson, Zhe Wang, Moses Charikar, Kai Li. Multi-Probe LSH:
 * Efficient Indexing for High-Dimensional Similarity Search. Proceedings of the
 * 33rd International Conference on Very Large Data Bases (VLDB). Vienna,
 * Austria. September 2007.
 *
 */


class MultiProbeLSH : public RepeatHash<GaussianLSH>
{
  unsigned int H_;

public:
  typedef RepeatHash<GaussianLSH> Super;
  typedef Super::Domain Domain;

  struct Parameter : public Super::Parameter
  {
    unsigned int range;

    template<class Archive>
    void serialize(Archive& ar)
    {
      ar & range;
      ar & repeat;
      ar & dim;
      ar & W;
    }
  };

  MultiProbeLSH() {}

  template<typename RNG>
  void reset(const Parameter& param, RNG& rng)
  {
    H_ = param.range;
    Super::reset(param, rng);
  }

  template<typename RNG>
  MultiProbeLSH(const Parameter& param, RNG& rng)
  {
    reset(param, rng);
  }

  unsigned int getRange() const
  {
    return H_;
  }

  unsigned int operator () (Domain obj) const
  {
    return Super::operator ()(obj) % H_;
  }

  template<typename Archive>
  void serialize(Archive& ar)
  {
    Super::serialize(ar);
    ar & H_;
  }

  void generateProbeSequence(Domain obj, std::vector<unsigned int>& seq, unsigned int T) const;
};

template<typename HashTable, typename Key>
class MultiProbeLSHIndex : public LSHIndex<MultiProbeLSH, HashTable, Key>
{
public:
  typedef LSHIndex<MultiProbeLSH, HashTable, Key> Super;

  typedef typename Super::Parameter Parameter;

  typedef typename Super::Domain Domain;

protected:
  Parameter param_;

public:

  MultiProbeLSHIndex() {}

  template<typename RNG>
  void init(const Parameter& param, RNG& rng, unsigned int L)
  {
    Super::init(param, rng, L);
    param_ = param;
  }

  template<typename Scanner>
  void query(unsigned int T, Scanner& scanner)
  {
    std::vector<unsigned int> seq;
    for(unsigned int i = 0; i < Super::lshs_.size(); ++i)
    {
      Super::lshs_[i].generateProbeSequence(scanner.query(), seq, T);
      for(unsigned int j = 0; j < seq.size(); ++j)
        Super::tables_[i].query(j, scanner);
    }
  }

};

}

#endif
