#ifndef CSPACE_LEARNING_COMPOSITE_H
#define CSPACE_LEARNING_COMPOSITE_H

#include <assert.h>

namespace cspace_learning
{

/**
 *  The mod of an LSH function by some value N is usually still locality
 *  sensitive.  This can be used to limit the hash value of certain LSH,
 *  so that the hash value can be used to index a fixed-sized hash table.
 *
 *  Let LSH be the original class, and N be the divisor, then the parameter type
 *  is defined as
 *  \code
 *      struct Parameter {
 *          unsigned range;     // the divisor.
 *          ...                 // the original parameters are inherited.
 *      };
 * \endcode
 * The domain is the same as the original LSH and the range is always N.
 *
 */
template<typename LSH>
class Tail
{
protected:
  LSH lsh_;
  unsigned int range_;

public:

  typedef LSH Super;

  struct Parameter : public Super::Parameter
  {
    unsigned int range;
  };

  typedef typename Super::Domain Domain;

  Tail()
  {
  }

  template<typename RNG>
  void reset(const Parameter& param, RNG& rng)
  {
    range_ = param.range;
    lsh_.reset(param, rng);
  }

  template<typename RNG>
  Tail(const Parameter& param, RNG& rng)
  {
    reset(param, rng);
  }

  unsigned int getRange() const
  {
    return range_;
  }

  unsigned int operator () (Domain obj) const
  {
    return lsh_(obj) % range_;
  }

  template<typename Archive>
  void serialize(Archive& ar)
  {
    ar & lsh_;
    ar & range_;
  }
};

/**
 * This is the same as Tail, except the divisor is determined at compile time
 * and passed in as a template parameter.  Both the Domain and Parameter of the
 * original LSH are kept the same.
 */

template<typename LSH, unsigned int range>
class FixedTail
{
protected:
  LSH lsh_;

public:
  typedef LSH Super;

  struct Parameter : public Super::Parameter
  {
  };

  typedef typename Super::Domain Domain;

  FixedTail() {}

  template<typename RNG>
  void reset(const Parameter& param, RNG& rng)
  {
    lsh_.reset(param, rng);
  }

  template<typename RNG>
  FixedTail(const Parameter& param, RNG& rng)
  {
    reset(param, rng);
  }

  unsigned int getRange() const
  {
    return range;
  }

  unsigned int operator () (Domain obj) const
  {
    return lsh_(obj) % range;
  }

  const LSH& getLSH() const
  {
    return lsh_;
  }

  template<typename Archive>
  void serialize(Archive& ar)
  {
    ar & lsh_;
  }
};

/**
 * This is a special case of FixedTail, with the divisor being 2, with the effect
 * of taking only the least significant bit of the hash value.  This is mainly used
 * to generate sketches.
 */

template<typename LSH>
class LSB : public FixedTail<LSH, 2>
{
public:
  typedef FixedTail<LSH, 2> Super;

  struct Parameter : public Super::Parameter
  {
  };

  typedef typename Super::Domain Domain;

  LSB() {}

  template<typename RNG>
  LSB(const Parameter& param, RNG& rng) : FixedTail<LSH, 2>(param, rng) {}
};

/**
 * The concatenation of a number of LSHes of the same class is usually used as a
 * new LSH to augment the locality sensitivity.  The Repeat class is to
 * concatenate N independent LSH instances.
 *
 * The domain of LSH remains the same.  The new parameter is
 * defined as:
 *      \code
 *      struct Parameter {
 *          unsigned repeat;     // # of LSHes to concatenate.
 *          ...                  // all parameters of the base LSH are inherited.
 *      };
 *      \endcode
 *
 * Because the hash value is represented as unsigned int, which has only 32
 * bits, the range of the original LSH need to be small enough so that the
 * concatenated value does not overflow.  Specifically, we require that
 * \f[
 * LSH::getRange()^N <= 2^{32}.
 * \f]
 * We also require that the range of the base LSH only depends on the
 * parameter, so an array of such LSHes initialized with the same parameter but
 * independent random numbers have the same range.
 */

template<typename LSH>
class Repeat
{
public:
  typedef LSH Super;

  struct Parameter : public Super::Parameter
  {
    unsigned int repeat;
  };

  typedef typename Super::Domain Domain;

  Repeat()
  {
  }

  template<typename RNG>
  void reset(const Parameter& param, RNG& rng)
  {
    dup_ = param.repeat;
    assert(dup_ > 0);
    lsh_.resize(dup_);
    lsh_[0].reset(param, rng);
    range_ = unit_ = lsh_[0].getRange();
    assert(unit_ > 0);
    assert((unsigned int)(1 << (sizeof(unsigned int) * 8 / dup_)) >= unit_);
    for(unsigned int i = 1; i < dup_; ++i)
    {
      lsh_[i].reset(param, rng);
      assert(unit_ == lsh_[i].getRange());
      range_ *= unit_;
    }
  }

  template<typename RNG>
  Repeat(const Parameter& param, RNG& rng)
  {
    reset(param, rng);
  }

  unsigned int getRange() const
  {
    return range_;
  }

  unsigned int operator() (Domain obj) const
  {
    unsigned int ret = 0;
    for(unsigned int i = 0; i < dup_; ++i)
    {
      ret *= unit_;
      ret += lsh_[i](obj);
    }
    return ret;
  }

  template<typename Archive>
  void serialize(Archive& ar)
  {
      ar & dup_;
      ar & range_;
      ar & unit_;
      ar & lsh_;
  }
protected:
  std::vector<Super> lsh_;
  unsigned int dup_;
  unsigned int range_;
  unsigned int unit_;
};

/**
 * This composition is to workaround the case where the range of individual
 * LSHes are so large that the concatenation can not be held in a single
 * unsigned int.  The method is to further hash the concatenated value.
 * Specifically, if <h1, h2, ..., hN> are the original values, this composition
 * produces (a1*h1 + a2*h2 + aN*hN), with a1~aN being random unsigned integers.
 * The range of the produced LSH is 0 (the whole range of unsigned).
 *
 * The domain of the LSH remains the same.  The new parameter is
 * defined as:
 *      \code
 *      struct Parameter {
 *          unsigned repeat;     // # of LSHes to concatenate
 *          ...                  // all parameters of the base LSH are inherited.
 *      };
 *      \endcode
 */
template<typename LSH>
class RepeatHash
{
public:
  typedef LSH Super;

  struct Parameter : public Super::Parameter
  {
    unsigned int repeat;
  };

  typedef typename LSH::Domain Domain;

  RepeatHash() {}

  template<typename RNG>
  void reset(const Parameter& param, RNG& rng)
  {
    assert(param.repeat > 0);
    lsh_.resize(param.repeat);
    for(unsigned int i = 0; i < param.repeat; ++i)
    {
      lsh_[i].reset(param, rng);
    }

    a_.resize(param.repeat);

    for(unsigned int i = 0; i < param.repeat; ++i)
      a_[i] = rng();
  }

  template<typename RNG>
  RepeatHash(const Parameter& param, RNG& rng)
  {
    reset(param, rng);
  }

  unsigned int getRange() const
  {
    return 0;
  }

  unsigned int operator () (Domain obj) const
  {
    unsigned int ret = 0;
    for(unsigned int i = 0; i < lsh_.size(); ++i)
      ret += lsh_[i](obj) * a_[i];

    return ret;
  }

  template<typename Archive>
  void serialize(Archive& ar)
  {
    ar & lsh_;
    ar & a_;
    assert(a_.size() == lsh_.size());
  }


  protected:
      std::vector<Super> lsh_;
      std::vector<unsigned int> a_;
};

/*
 * The XOR of a number of 1-bit LSHes has higher locality sensitivity than
 * the original LSH.  This serves a similar purpose as RepeatHash.
 *
 * The new parameter is defined as:
 *      \code
 *      struct Parameter {
 *          unsigned xor_;     // # of LSHes to XOR
 *          ...                  // all parameters of the base LSH are inherited.
 *      };
 *      \endcode
 */
template<typename LSH>
class Xor
{
public:
  typedef LSH Super;

  struct Parameter: public Super::Parameter
  {
    unsigned xor_;
  };

  typedef typename Super::Domain Domain;

  Xor() {}

  template<typename RNG>
  void reset(const Parameter& param, RNG &rng)
  {
    lsh_.resize(param.xor_);
    for(unsigned int i = 0; i < lsh_.size(); ++i)
    {
        lsh_[i].reset(param, rng);
        assert(lsh_[i].getRange() == 2);
    }
  }

  template<typename RNG>
  Xor(const Parameter& param, RNG &rng)
  {
    reset(param, rng);
  }

  unsigned int getRange () const
  {
    return 2;
  }

  unsigned int operator () (Domain obj) const
  {
    unsigned int ret = 0;
    for(unsigned int i = 0; i < lsh_.size(); ++i)
    {
        ret = ret ^ lsh_[i](obj);
    }
    return ret;
  }

  template<typename Archive>
  void serialize(Archive & ar)
  {
    ar & lsh_;
  }

protected:
  std::vector<Super> lsh_;
};



}

#endif
