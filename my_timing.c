
#include "my_timing.h"

#ifdef _WIN32

    void timer_now(typeTime* time)
    {
        QueryPerformanceCounter(time);
    }

    OSQPFloat timer_diff_ms(const typeTime* start, const typeTime* end)
    {
        LARGE_INTEGER frequency;
        QueryPerformanceFrequency(&frequency);
        return (end->QuadPart - start->QuadPart) * ((OSQPFloat)1000) / frequency.QuadPart;
    }

#elif defined(__linux__) || defined(__APPLE__)

    void timer_now(typeTime* time)
    {
        gettimeofday(time, 0);
    }

    OSQPFloat timer_diff_ms(const typeTime* start, const typeTime* end)
    {
        return (end->tv_sec - start->tv_sec) * ((OSQPFloat)1000)
            + (end->tv_usec - start->tv_usec) / ((OSQPFloat)1000);
    }

#else

    void timer_now(typeTime* time)
    {
        *time = clock();
    }

    OSQPFloat timer_diff_ms(const typeTime* start, const typeTime* end)
    {
        return (*end - *start) * ((OSQPFloat)1000) / CLOCKS_PER_SEC;
    }

#endif