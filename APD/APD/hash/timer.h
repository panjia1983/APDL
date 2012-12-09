#ifndef CSPACE_LEARNING_TIMER_H
#define CSPACE_LEARNING_TIMER_H

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

namespace cspace_learning
{

class Timer
{
public:
        Timer();
        ~Timer();

        void   start();                             ///< start timer
        void   stop();                              ///< stop the timer
        double getElapsedTime();                    ///< get elapsed time in milli-second
        double getElapsedTimeInSec();               ///< get elapsed time in second (same as getElapsedTime)
        double getElapsedTimeInMilliSec();          ///< get elapsed time in milli-second
        double getElapsedTimeInMicroSec();          ///< get elapsed time in micro-second

private:
        double startTimeInMicroSec;                 ///< starting time in micro-second
        double endTimeInMicroSec;                   ///< ending time in micro-second
        int    stopped;                             ///< stop flag
#ifdef _WIN32
        LARGE_INTEGER frequency;                    ///< ticks per second
        LARGE_INTEGER startCount;
        LARGE_INTEGER endCount;
#else
        timeval startCount;
        timeval endCount;
#endif
};

}

#endif
