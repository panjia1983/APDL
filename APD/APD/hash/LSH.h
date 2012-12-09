#ifndef CSPACE_LEARNING_LSH_H
#define CSPACE_LEARNING_LSH_H

#include "hash/common.h"

namespace cspace_learning
{

class TrivialLSH
{
public:
  struct Parameter
  {
  };

  typedef const float* Domain;

  TrivialLSH()
  {
  }

  template<typename RNG>
  void reset(const Parameter& param, RNG& rng)
  {
  }

  template<typename RNG>
  TrivialLSH(const Parameter& param, RNG& rng)
  {
  }

  unsigned getRange() const
  {
    return 1;
  }

  unsigned operator () (Domain obj) const
  {
    return 0;
  }

  template<typename Archive>
  void serialize(Archive& ar)
  {
  }
};

/**
 * This LSH is defined on the D-dimensional vector space.  For a vector X, the
 * hash value is defined as
 * \f[
 *                h(X) = [b + a_1*X_1 + a_2*X_2 + ... + a_D*X_D] / W
 * \f]
 * where W is a positive value called the window size; b is sampled uniformly from
 * [0, W); a1 ~ aD are N random variables independently sampled from the so-called stable
 * distribution, which is specifiled by the template parameter DIST.
 *
 * The domain of the LSH is (const float *), and the parameter is defined as
 * \code
 *      struct Parameter {
 *          unsigned dim;
 *          float W;
 *      };
 * \endcode
 * The range of this LSH is 0.
 *
 * Following are two spacial cases:
 * \code
 *     typedef StableDistLsh<Cauchy> CauchyLsh;
 *     typedef StableDistLsh<Gaussian> GaussianLsh;
 * \endcode
 * Cauchy distribution is 1-stable and Gaussian distribution is 2-stable.  These
 * two LSHes can be used to approximate L1 and L2 distances respectively.
 *
 * For more information on stable distribution based LSH, see the following reference.
 *
 *     Mayur Datar , Nicole Immorlica , Piotr Indyk , Vahab S. Mirrokni,
 *     Locality-sensitive hashing scheme based on p-stable distributions,
 *     Proceedings of the twentieth annual symposium on Computational geometry, June
 *     08-11, 2004, Brooklyn, New York, USA.
 */
template<typename Distrib>
class StableDistLSH
{
  std::vector<float> a_;
  float b_;
  float W_;
  unsigned int dim_;

public:

  struct Parameter
  {
    unsigned int dim;
    float W;
  };

  typedef const float* Domain;

  StableDistLSH()
  {
  }

  template<typename RNG>
  void reset(const Parameter& param, RNG& rng)
  {
    a_.resize(param.dim);
    W_ = param.W;
    dim_ = param.dim;

    boost::variate_generator<RNG&, Distrib> gen(rng, Distrib());

    for(unsigned int i = 0; i < dim_; ++i) a_[i] = gen();

    b_ = boost::variate_generator<RNG&, Uniform>(rng, Uniform(0, W_))();

  }

  template<typename RNG>
  StableDistLSH(const Parameter& param, RNG& rng)
  {
    reset(param, rng);
  }

  unsigned int getRange() const
  {
    return 0;
  }

  unsigned int operator () (Domain obj) const
  {
    dotprod<float> dot(dim_);
    float ret = b_ + dot(&a_[0], obj);
    return unsigned(int(std::floor(ret / W_)));
  }

  unsigned int operator () (Domain obj, float* delta) const
  {
    dotprod<float> dot(dim_);
    float ret = b_ + dot(&a_[0], obj);
    ret /= W_;

    float flr = std::floor(ret);
    *delta = ret - flr;
    return unsigned(int(flr));
  }

  template<typename Archive>
  void serialize(Archive& ar)
  {
    ar & a_;
    ar & b_;
    ar & W_;
    ar & dim_;
  }
};

/** LSH for L1 metric */
typedef StableDistLSH<Cauchy> CauchyLSH;

/** LSH for L2 metric */
typedef StableDistLSH<Gaussian> GaussianLSH;


template<typename Distrib>
class SparseStableDistLSH
{
  std::vector<float> a_;
  float b_;
  float W_;
  unsigned int dim_;

public:

  struct Parameter
  {
    unsigned int dim;
    float W;
  };

  typedef SparseVector Domain;

  SparseStableDistLSH()
  {
  }

  template<typename RNG>
  void reset(const Parameter& param, RNG& rng)
  {
    a_.resize(param.dim);
    W_ = param.W;
    dim_ = param.dim;

    boost::variate_generator<RNG&, Distrib> gen(rng, Distrib());

    for(unsigned int i = 0; i < dim_; ++i) a_[i] = gen();

    b_ = boost::variate_generator<RNG&, Uniform>(rng, Uniform(0, W_))();

  }

  template<typename RNG>
  SparseStableDistLSH(const Parameter& param, RNG& rng)
  {
    reset(param, rng);
  }

  unsigned int getRange() const
  {
    return 0;
  }

  unsigned int operator () (const SparseVector& obj) const
  {
    dotprod<float> dot(dim_);
    float ret = b_ + dot(obj, &a_[0]);
    return unsigned(int(std::floor(ret / W_)));
  }

  unsigned int operator () (const SparseVector& obj, float* delta) const
  {
    dotprod<float> dot(dim_);
    float ret = b_ + dot(obj, &a_[0]);
    ret /= W_;

    float flr = std::floor(ret);
    *delta = ret - flr;
    return unsigned(int(flr));
  }

  template<typename Archive>
  void serialize(Archive& ar)
  {
    ar & a_;
    ar & b_;
    ar & W_;
    ar & dim_;
  }
};

/** LSH for L1 metric */
typedef SparseStableDistLSH<Cauchy> SparseCauchyLSH;

/** LSH for L2 metric */
typedef SparseStableDistLSH<Gaussian> SparseGaussianLSH;


/**
 * Random hyperplane based LSH can be used to approximate cosine similarity.
 * This LSH is defined on the D-dimensional vector space.  For a vector X, the
 * hash value is defined as
 * \f[
 *                h(X) = [a_1*X_1 + a_2*X_2 + ... + a_D*X_D] > 0 ? 1 : 0
 * \f]
 * where <a1,...,aD> is a random vector sampled from the unit hypersphere.
 *
 * The domain of the LSH is (const float *), and the parameter is defined as
 * \code
 *      struct Parameter {
 *          unsigned dim;
 *      };
 * \endcode
 * The range of this LSH is 0.
 *
 * For more information on stable distribution based LSH, see the following reference.
 *
 *      Charikar, M. S. 2002. Similarity estimation techniques from rounding
 *      algorithms. In Proceedings of the Thiry-Fourth Annual ACM Symposium on
 *      theory of Computing (Montreal, Quebec, Canada, May 19 - 21, 2002). STOC '02.
 *      ACM, New York, NY, 380-388. DOI= http://doi.acm.org/10.1145/509907.509965
 *
 */
class HyperPlaneLSH
{
  std::vector<float> a_;

public:

  struct Parameter
  {
    unsigned int dim;
  };

  typedef const float* Domain;

  HyperPlaneLSH()
  {
  }

  template<typename RNG>
  void reset(const Parameter& param, RNG& rng)
  {
    a_ = boost::variate_generator<RNG&, boost::uniform_on_sphere<float> >(rng, boost::uniform_on_sphere<float>(param.dim))();
  }

  template<typename RNG>
  HyperPlaneLSH(const Parameter& param, RNG& rng)
  {
    reset(param, rng);
  }

  unsigned int getRange() const
  {
    return 2;
  }

  unsigned int operator () (Domain obj) const
  {
    float ret = 0;
    for(unsigned int i = 0; i < a_.size(); ++i)
    {
      ret += a_[i] * obj[i];
    }
    return ret >= 0 ? 1 : 0;
  }

  unsigned int operator () (Domain obj, float* delta) const
  {
    float ret = 0;
    for(unsigned int i = 0 ; i < a_.size(); ++i)
    {
      ret += a_[i] * obj[i];
    }

    if(ret >= 0)
    {
      *delta = ret;
      return 1;
    }
    else
    {
      *delta = -ret;
      return 0;
    }
  }

  template<typename Archive>
  void serialize(Archive& ar)
  {
    ar & a_;
  }
};

/**
 *   This LSH can be used to approximate L1 distance for closed D-dimensional
 *   space [min, max]^D.  It hashes each input vector into a 0-1 value, so its
 *   range is 2.  The hash function is the following:  a random dimension is
 *   chosen, and a threshold T is sampled uniformly in [min, max].  For each input
 *   vector, the value at that dimension is check.  If the value is larger than T,
 *   1 is returned and otherwise 0 is returned.  Following is the parameter.
 *
 *   \code
 *   struct Parameter
 *   {
 *       unsigned dim;  // Dimension of domain.
 *       float min;     // Each dimension is limited in [min, max].
 *       float max;
 *   };
 *   \endcode
 *
 *   The method is discussed in the following papers:
 *
 *      Zhe Wang, Wei Dong, William Josephson, Qin Lv, Moses Charikar, Kai Li.
 *      Sizing Sketches: A Rank-Based Analysis for Similarity Search. In
 *      Proceedings of the 2007 ACM SIGMETRICS International Conference on
 *      Measurement and Modeling of Computer Systems . San Diego, CA, USA. June
 *      2007.
 *
 *      Qin Lv, Moses Charikar, Kai Li. Image Similarity Search with Compact
 *      Data Structures. In Proceedings of ACM 13th Conference on Information
 *      and Knowledge Management (CIKM), Washington D.C., USA. November 2004.
 *
 *   Note that the original method allows the range of each dimension to be
 *   different and also allow each dimensioin to carry a weight.  The
 *   implementation here is simplified.
 */
class ThresholdingLSH
{
  unsigned int dim_;
  float threshold_;

public:

  struct Parameter
  {
    unsigned int dim;
    float min;
    float max;
  };

  typedef const float* Domain;

  ThresholdingLSH()
  {
  }

  template<typename RNG>
  void reset(const Parameter& param, RNG& rng)
  {
    dim_ = boost::variate_generator<RNG&, UniformUnsigned>(rng, UniformUnsigned(0, param.dim - 1));
    threshold_ = boost::variate_generator<RNG &, Uniform>(rng, Uniform(param.min, param.max))();
  }

  template<typename RNG>
  ThresholdingLSH(const Parameter& param, RNG& rng)
  {
    reset(param, rng);
  }

  unsigned int getRange() const
  {
    return 2;
  }

  unsigned int operator () (Domain obj) const
  {
    return obj[dim_] >= threshold_ ? 1 : 0;
  }

  unsigned int operator () (Domain obj, float* delta) const
  {
    float ret = obj[dim_] - threshold_;
    if(ret > 0)
    {
      *delta = ret;
      return 1;
    }
    else
    {
      *delta = -ret;
      return 0;
    }
  }

  template<typename Archive>
  void serialize(Archive& ar)
  {
    ar & dim_;
    ar & threshold_;
  }
};

}

#endif
