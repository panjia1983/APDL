#ifndef CSPACE_LEARNING_ARCHIVE_H
#define CSPACE_LEARNING_ARCHIVE_H

#include <vector>
#include <iostream>

namespace cspace_learning
{

static inline std::ostream& operator & (std::ostream& os, int i)
{
  return os.write((const char *)&i, sizeof(i));
}

static inline std::ostream& operator & (std::ostream& os, unsigned int i)
{
  return os.write((const char *)&i, sizeof(i));
}

static inline std::ostream& operator & (std::ostream& os, float i)
{
  return os.write((const char *)&i, sizeof(i));
}

static inline std::istream& operator & (std::istream& is, int& i)
{
  return is.read((char *)&i, sizeof(i));
}

static inline std::istream& operator & (std::istream& is, unsigned int& i)
{
  return is.read((char *)&i, sizeof(i));
}

static inline std::istream& operator & (std::istream& is, float& i)
{
  return is.read((char *)&i, sizeof(i));
}

static inline std::ostream& operator & (std::ostream& os, std::vector<float>& v)
{
  unsigned int L = v.size();
  os & L;
  os.write((const char *)&v[0], v.size() * sizeof(float));
  return os;
}

static inline std::istream& operator & (std::istream& is, std::vector<float>& v)
{
  unsigned int L;
  is & L;
  v.resize(L);
  is.read((char *)&v[0], v.size() * sizeof(float));
  return is;
}

static inline std::ostream& operator & (std::ostream& os, std::vector<unsigned int>& v)
{
  unsigned int L = v.size();
  os & L;
  os.write((const char *)&v[0], v.size() * sizeof(unsigned int));
  return os;
}

static inline std::istream& operator & (std::istream& is, std::vector<unsigned int>& v)
{
  unsigned int L;
  is & L;
  v.resize(L);
  is.read((char *)&v[0], v.size() * sizeof(unsigned int));
  return is;
}

template <typename D>
std::ostream& operator & (std::ostream& s, D& t)
{
  t.serialize(s);
  return s;
}

template <typename D>
std::istream& operator & (std::istream& s, D& t)
{
  t.serialize(s);
  return s;
}

template <typename D>
std::ostream& operator & (std::ostream& os, std::vector<D>& v)
{
  unsigned int L = v.size();
  os & L;
  for(unsigned int i = 0; i < L; ++i) {
      os & v[i];
  }
  return os;
}

template <typename D>
std::istream& operator & (std::istream& is, std::vector<D>& v)
{
  unsigned int L;
  is & L;
  v.resize(L);
  for(unsigned int i = 0; i < L; ++i) {
      is & v[i];
  }
  return is;
}

}

#endif
