#include "osqp.h"

#ifndef TIMING_H
#define TIMING_H


#ifdef _WIN32

    #define WIN32_LEAN_AND_MEAN
    #include <Windows.h>

    typedef LARGE_INTEGER typeTime;

#elif defined(__linux__) || defined(__APPLE__)

    #include <sys/time.h> 

    typedef struct timeval typeTime;

#else

    #include <time.h>

    typedef clock_t typeTime;

#endif

/* Get current timestamp */
void timer_now(typeTime* time);

/* Get elapsed time between two timestamps in milliseconds  */
OSQPFloat timer_diff_ms(const typeTime* start, const typeTime* end);

#endif