#ifndef CSPACE_LEARNING_HASH_H
#define CSPACE_LEARNING_HASH_H

namespace cspace_learning
{

struct robert_jenkins_hash
{
  unsigned int operator () (unsigned int a) const
  {
    a = (a+0x7ed55d16) + (a<<12);
    a = (a^0xc761c23c) ^ (a>>19);
    a = (a+0x165667b1) + (a<<5);
    a = (a+0xd3a2646c) ^ (a<<9);
    a = (a+0xfd7046c5) + (a<<3);
    a = (a^0xb55a4f09) ^ (a>>16);
    return a;
  }
};

struct murmur_hash
{
  murmur_hash(unsigned int seed)
  {
    seed_ = seed;
  }

  unsigned int operator () (unsigned int a) const
  {
    const unsigned int m = 0x5bd1e995;
    const int r = 24;

    unsigned int h = seed_ ^ 4;

    a *= m;
    a ^= a >> r;
    a *= m;

    h *= m;
    h ^= a;

    h ^= h >> 13;
    h *= m;
    h ^= h >> 15;
    return h;
  }

  unsigned int seed_;
};


}

#endif
