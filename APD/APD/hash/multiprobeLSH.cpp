#include "hash/multiprobeLSH.h"
#include <assert.h>
#include <queue>

namespace cspace_learning
{

static inline unsigned long long leftshift(unsigned int N)
{
  return (unsigned long long) 1 << N;
}

struct Probe
{
  unsigned long long mask;
  unsigned long long shift;
  float score;
  unsigned int reserve;
  bool operator < (const Probe& p) const
  {
    return score < p.score;
  }

  Probe operator + (const Probe& p) const
  {
    Probe ret;
    ret.mask = mask | p.mask;
    ret.shift = shift | p.shift;
    ret.score = score + p.score;
    return ret;
  }

  bool conflict(const Probe& p) const
  {
    return (mask & p.mask) != 0;
  }

  static const unsigned int MAX_M = 64;
  static const unsigned int MAX_T = 200;
};

typedef std::vector<Probe> ProbeSequence;

class ProbeGT
{
public:
  bool operator () (const Probe& p1, const Probe& p2) const
  {
    return p2 < p1;
  }
};

void generateExpectScores(ProbeSequence &seq, unsigned int M)
{
  assert(M <= sizeof(seq[0].mask)* 8);
  seq.resize(2 * M);
  for(unsigned int l = 0; l < M; ++l)
  {
    unsigned int r = 2 * M - l - 1;
    seq[l].mask = seq[r].mask = seq[r].shift = leftshift(l);
    seq[l].shift = 0;
    seq[l].reserve = seq[r].reserve = 0;
    float delta = (l + 1.0) / (M + 1.0) * 0.5;
    seq[l].score = (l + 1.0) * (l + 2.0) / (M + 1.0) / (M + 2.0) * 0.25;
    seq[r].score = 1.0 - 2.0 * delta + seq[l].score;
  }
}

void generateProbeSequenceTemplate(ProbeSequence &seq, unsigned int M, unsigned int T)
{
  ProbeSequence scores;
  generateExpectScores(scores, M);
  assert(T > 0);

  std::priority_queue<Probe, std::vector<Probe>, ProbeGT> heap;
  Probe init;
  init.mask = init.shift = 0;
  init.score = 0;
  init.reserve = 0;
  heap.push(init);

  seq.clear();

  for(;;)
  {
    if(heap.empty()) break;
    seq.push_back(heap.top());
    if(seq.size() == T) break;

    Probe shift = heap.top();
    heap.pop();

    for(unsigned next = shift.reserve; next < 2 * M; ++next)
    {
      if(!shift.conflict(scores[next]))
      {
        Probe tmp = shift + scores[next];
        tmp.reserve = next + 1;
        heap.push(tmp);
      }
    }
  }
}

class ProbeSequenceTemplates: public std::vector<ProbeSequence>
{
public:
  ProbeSequenceTemplates(unsigned int max_M, unsigned int max_T)
      : std::vector<ProbeSequence>(max_M + 1)
  {
    for(unsigned int i = 1; i <= max_M; ++i)
    {
      generateProbeSequenceTemplate(at(i), i, max_T);
    }
  }
};

static ProbeSequenceTemplates probeSequenceTemplates(Probe::MAX_M, Probe::MAX_T);


void MultiProbeLSH::generateProbeSequence(Domain obj, std::vector<unsigned int>& seq, unsigned int T) const
{
  ProbeSequence scores;
  std::vector<unsigned int> base;
  scores.resize(2 * lsh_.size());
  base.resize(lsh_.size());

  for(unsigned int i = 0; i < lsh_.size(); ++i)
  {
    float delta;
    base[i] = Super::lsh_[i](obj, &delta);
    scores[2*i].mask = i;
    scores[2*i].reserve = 1;    // direction
    scores[2*i].score = delta;
    scores[2*i+1].mask = i;
    scores[2*i+1].reserve = unsigned(-1);
    scores[2*i+1].score = 1.0 - delta;
  }
  std::sort(scores.begin(), scores.end());

  ProbeSequence& tmpl = probeSequenceTemplates[lsh_.size()];

  seq.clear();
  for(ProbeSequence::const_iterator it = tmpl.begin(); it != tmpl.end(); ++it)
  {
    if(seq.size() == T) break;
    const Probe& probe = *it;
    unsigned int hash = 0;
    for(unsigned int i = 0; i < lsh_.size(); ++i)
    {
      unsigned int h = base[scores[i].mask];
      if(probe.mask & leftshift(i))
      {
        if(probe.shift & leftshift(i))
        {
          h += scores[i].reserve;
        }
        else
        {
          h += unsigned(-1) * scores[i].reserve;
        }
      }
      hash += h * a_[scores[i].mask];
    }
    seq.push_back(hash % H_);
  }
}

}
