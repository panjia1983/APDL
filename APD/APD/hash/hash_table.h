#ifndef CSPACE_LEARNING_HASH_TABLE_H
#define CSPACE_LEARNING_HASH_TABLE_H

#include <vector>
#include <list>
#include <google_hash/sparse_hash_map>
#include <google_hash/dense_hash_map>
#include <hash/cuckoo.hpp>
#include <map>
#include <hash_map>
#include <stdexcept>

namespace cspace_learning
{

template<typename Key, typename Data, typename HashFcn, typename EqualKey = std::equal_to<Key> >
class LinearHashTable
{
  HashFcn h_;
public:

  void init(unsigned int range)
  {
    if(range == 0)
    {
      throw std::logic_error("LSH with unlimited range should not use LinearHashTable.");
    }

    table_.resize(range);
  }

  void insert(Key key, Data value)
  {
    unsigned int index = h_(key);
    table_[index].push_back(value);
  }

  template<typename Scanner>
  void query(Key key, Scanner& scanner) const
  {
    unsigned int index = h_(key);
    for(unsigned int i = 0; i < table_[index].size(); ++i)
      scanner(table_[index][i]);
  }

protected:
  typedef std::vector<Data> Bin;
  std::vector<Bin> table_;

};

template<typename Key, typename Data, typename HashFcn, typename EqualKey>
class GoogleSparseHashTable
{
protected:
  typedef std::list<Data> Bin;
  typedef google::sparse_hash_map<Key, Bin, HashFcn, EqualKey> HashTable;

  HashTable table_;

public:

  void init(unsigned int range)
  {
  }

  void insert(Key key, Data value)
  {
    table_[key].push_back(value);
  }

  template<typename Scanner>
  void query(Key key, Scanner& scanner) const
  {
    typename HashTable::const_iterator p = table_.find(key);
    if(p != table_.end())
    {
      for(typename Bin::const_iterator i = (*p).second.begin(); i != (*p).second.end(); ++i)
        scanner(*i);
    }
  }
};

template<typename Key, typename Data, typename HashFcn, typename EqualKey>
class GoogleDenseHashTable
{
protected:
  typedef std::list<Data> Bin;
  typedef google::dense_hash_map<Key, Bin, HashFcn, EqualKey> HashTable;

  HashTable table_;

public:

  void init(unsigned int range)
  {
  }

  GoogleDenseHashTable()
  {
    table_.set_empty_key(NULL);
  }

  void insert(Key key, Data value)
  {
    table_[key].push_back(value);
  }

  template<typename Scanner>
  void query(Key key, Scanner& scanner) const
  {
    typename HashTable::const_iterator p = table_.find(key);
    if(p != table_.end())
    {
      for(typename Bin::const_iterator i = (*p).second.begin(); i != (*p).second.end(); ++i)
        scanner(*i);
    }
  }
};

template<typename Key, typename Data, typename HashFcn, typename EqualKey>
class MultiMapHashTable
{
protected:
  typedef std::list<Data> Bin;
  typedef __gnu_cxx::hash_multimap<Key, Bin, HashFcn, EqualKey> HashTable;

  HashTable table_;

public:
  void init(unsigned int range)
  {
  }

  void insert(Key key, Data value)
  {
    table_.insert(HashTable::value_type(key, value));
  }

  template<typename Scanner>
  void query(Key key, Scanner& scanner) const
  {
    std::pair<typename HashTable::const_iterator, typename HashTable::const_iterator> p = table_.equal_range(key);

    for(typename HashTable::const_iterator i = p.first; i != p.second; ++i)
      scanner((*i).second);
  }
};

template<typename Key, typename Data, typename Compare>
class MultiMapTable
{
protected:
  typedef std::multimap<Key, Data, Compare> Table;
  Table table_;

public:
  void init(unsigned int range)
  {
  }

  void insert(Key key, Data value)
  {
    table_.insert(std::make_pair(key, value));
  }

  template<typename Scanner>
  void query(Key key, Scanner& scanner) const
  {
    std::pair<typename Table::const_iterator, typename Table::const_iterator> p = table_.equal_range(key);

    for(typename Table::const_iterator i = p.first; i != p.second; ++i)
      scanner((*i).second);
  }
};

template<typename Key, typename Data, template <typename> class HashFcns, typename EqualKey>
class CuckooHashTable
{
protected:
  typedef std::list<Data> Bin;
  typedef cuckoo<Key, Bin, HashFcns<Key>, EqualKey> HashTable;

  HashTable table_;

public:
  void init(unsigned int range)
  {
  }

  void insert(Key key, Data value)
  {
    table_[key].push_back(value);
  }

  template<typename Scanner>
  void query(Key key, Scanner& scanner) const
  {
    typename HashTable::const_iterator p = table_.find(key);
    if(p != table_.end())
    {
      for(typename Bin::const_iterator i = (*p).second.begin(); i != (*p).second.end(); ++i)
        scanner(*i);
    }
  }
};


}

#endif
