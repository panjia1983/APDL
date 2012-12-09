#ifndef CSPACE_LEARNING_LSH_INDEX_H
#define CSPACE_LEARNING_LSH_INDEX_H


namespace cspace_learning
{

template<typename LSH, typename HashTable, typename Key>
class LSHIndex
{
public:
  typedef typename LSH::Parameter Parameter;
  typedef typename LSH::Domain Domain;

protected:

  std::vector<LSH> lshs_;
  std::vector<HashTable> tables_;


public:
  template<typename RNG>
  void init(const Parameter& param, RNG& rng, unsigned int L)
  {
    lshs_.resize(L);
    tables_.resize(L);

    for(unsigned int i = 0; i < L; ++i)
    {
      lshs_[i].reset(param, rng);
      tables_[i].init(lshs_[i].getRange());
    }
  }

  template<typename Scanner>
  void query(Scanner& scanner) const
  {
    for(unsigned int i = 0; i < lshs_.size(); ++i)
    {
      unsigned int index = lshs_[i](scanner.query());
      tables_[i].query(index, scanner);
    }
  }

  void insert(Key key, const Domain& value)
  {
    for(unsigned int i = 0; i < lshs_.size(); ++i)
    {
      unsigned int index = lshs_[i](value);
      tables_[i].insert(index, key);
    }
  }
};

}

#endif
