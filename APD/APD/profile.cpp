#include "profile.h"


APDL::tools::Profiler& APDL::tools::Profiler::Instance(void)
{
  static Profiler p(true, false);
  return p;
}

#if ENABLE_PROFILING

#include <vector>
#include <algorithm>
#include <sstream>

void APDL::tools::Profiler::start(void)
{
  lock_.lock();
  if (!running_)
  {
    tinfo_.set();
    running_ = true;
  }
  lock_.unlock();
}

void APDL::tools::Profiler::stop(void)
{
  lock_.lock();
  if (running_)
  {
    tinfo_.update();
    running_ = false;
  }
  lock_.unlock();
}

void APDL::tools::Profiler::clear(void)
{
  lock_.lock();
  data_.clear();
  tinfo_ = TimeInfo();
  if (running_)
    tinfo_.set();
  lock_.unlock();
}

void APDL::tools::Profiler::event(const std::string &name, const unsigned int times)
{
  lock_.lock();
  data_[boost::this_thread::get_id()].events[name] += times;
  lock_.unlock();
}

void APDL::tools::Profiler::average(const std::string &name, const double value)
{
  lock_.lock();
  AvgInfo &a = data_[boost::this_thread::get_id()].avg[name];
  a.total += value;
  a.totalSqr += value*value;
  a.parts++;
  lock_.unlock();
}

void APDL::tools::Profiler::begin(const std::string &name)
{
  lock_.lock();
  data_[boost::this_thread::get_id()].time[name].set();
  lock_.unlock();
}

void APDL::tools::Profiler::end(const std::string &name)
{
  lock_.lock();
  data_[boost::this_thread::get_id()].time[name].update();
  lock_.unlock();
}

void APDL::tools::Profiler::status(std::ostream &out, bool merge)
{
  stop();
  lock_.lock();
  printOnDestroy_ = false;

  out << std::endl;
  out << " *** Profiling statistics. Total counted time : " << time::seconds(tinfo_.total) << " seconds" << std::endl;

  if (merge)
  {
    PerThread combined;
    for (std::map<boost::thread::id, PerThread>::const_iterator it = data_.begin() ; it != data_.end() ; ++it)
    {
      for (std::map<std::string, unsigned long int>::const_iterator iev = it->second.events.begin() ; iev != it->second.events.end(); ++iev)
        combined.events[iev->first] += iev->second;
      for (std::map<std::string, AvgInfo>::const_iterator iavg = it->second.avg.begin() ; iavg != it->second.avg.end(); ++iavg)
      {
        combined.avg[iavg->first].total += iavg->second.total;
        combined.avg[iavg->first].totalSqr += iavg->second.totalSqr;
        combined.avg[iavg->first].parts += iavg->second.parts;
      }
      for (std::map<std::string, TimeInfo>::const_iterator itm = it->second.time.begin() ; itm != it->second.time.end(); ++itm)
      {
        TimeInfo &tc = combined.time[itm->first];
        tc.total = tc.total + itm->second.total;
        tc.parts = tc.parts + itm->second.parts;
        if (tc.shortest > itm->second.shortest)
          tc.shortest = itm->second.shortest;
        if (tc.longest < itm->second.longest)
          tc.longest = itm->second.longest;
      }
    }
    printThreadInfo(out, combined);
  }
  else
    for (std::map<boost::thread::id, PerThread>::const_iterator it = data_.begin() ; it != data_.end() ; ++it)
    {
      out << "Thread " << it->first << ":" << std::endl;
      printThreadInfo(out, it->second);
    }
  lock_.unlock();
}


/// @cond IGNORE
namespace APDL
{

struct dataIntVal
{
  std::string name;
  unsigned long int value;
};

struct SortIntByValue
{
  bool operator()(const dataIntVal &a, const dataIntVal &b) const
  {
    return a.value > b.value;
  }
};

struct dataDoubleVal
{
  std::string name;
  double value;
};

struct SortDoubleByValue
{
  bool operator()(const dataDoubleVal &a, const dataDoubleVal &b) const
  {
    return a.value > b.value;
  }
};
}
/// @endcond

void APDL::tools::Profiler::printThreadInfo(std::ostream &out, const PerThread &data)
{
  double total = time::seconds(tinfo_.total);

  std::vector<dataIntVal> events;
  for (std::map<std::string, unsigned long int>::const_iterator iev = data.events.begin() ; iev != data.events.end() ; ++iev)
  {
    dataIntVal next = {iev->first, iev->second};
    events.push_back(next);
  }
  std::sort(events.begin(), events.end(), SortIntByValue());
  if (!events.empty())
    out << "Events:" << std::endl;
  for (unsigned int i = 0 ; i < events.size() ; ++i)
    out << events[i].name << ": " << events[i].value << std::endl;

  std::vector<dataDoubleVal> avg;
  for (std::map<std::string, AvgInfo>::const_iterator ia = data.avg.begin() ; ia != data.avg.end() ; ++ia)
  {
    dataDoubleVal next = {ia->first, ia->second.total / (double)ia->second.parts};
    avg.push_back(next);
  }
  std::sort(avg.begin(), avg.end(), SortDoubleByValue());
  if (!avg.empty())
    out << "Averages:" << std::endl;
  for (unsigned int i = 0 ; i < avg.size() ; ++i)
  {
    const AvgInfo &a = data.avg.find(avg[i].name)->second;
    out << avg[i].name << ": " << avg[i].value << " (stddev = " <<
      sqrt(fabs(a.totalSqr - (double)a.parts * avg[i].value * avg[i].value) / ((double)a.parts - 1.)) << ")" << std::endl;
  }

  std::vector<dataDoubleVal> time;

  for (std::map<std::string, TimeInfo>::const_iterator itm = data.time.begin() ; itm != data.time.end() ; ++itm)
  {
    dataDoubleVal next = {itm->first, time::seconds(itm->second.total)};
    time.push_back(next);
  }

  std::sort(time.begin(), time.end(), SortDoubleByValue());
  if (!time.empty())
    out << "Blocks of time:" << std::endl;

  double unaccounted = total;
  for (unsigned int i = 0 ; i < time.size() ; ++i)
  {
    const TimeInfo &d = data.time.find(time[i].name)->second;

    double tS = time::seconds(d.shortest);
    double tL = time::seconds(d.longest);
    out << time[i].name << ": " << time[i].value << "s (" << (100.0 * time[i].value/total) << "%), ["
        << tS << "s --> " << tL << " s], " << d.parts << " parts";
    if (d.parts > 0)
      out << ", " << (time::seconds(d.total) / (double)d.parts) << " s on average";
    out << std::endl;
    unaccounted -= time[i].value;
  }
  out << "Unaccounted time : " << unaccounted;
  if (total > 0.0)
    out << " (" << (100.0 * unaccounted / total) << " %)";
  out << std::endl;

  out << std::endl;
}

#endif