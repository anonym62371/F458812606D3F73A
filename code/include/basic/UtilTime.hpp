#ifndef __TIMEUTIL_H_
#define __TIMEUTIL_H_


#include <stack>
#include <iostream>
#include <fstream>

#ifndef WIN32
#include <sys/resource.h>
#include <sys/times.h>
#endif

float calculateExecutionTime(struct rusage *myTimeStart, struct rusage *myTimeEnd)
{
    return (
        ((float) (myTimeEnd->ru_utime.tv_sec  - myTimeStart->ru_utime.tv_sec)) +
        ((float) (myTimeEnd->ru_utime.tv_usec - myTimeStart->ru_utime.tv_usec)) * 1e-6
    );
}

std::stack<rusage> timer_stack;

enum Time_Unit {
    MILLISECOND, SECOND, MINUTE, HOUR
};

void timer_start() {
    rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    timer_stack.push(ru);
}


float timer_end(Time_Unit unit = SECOND) {
    if (timer_stack.empty()) {
        return 0.0f;
    }
    rusage start = timer_stack.top();
    timer_stack.pop();
    rusage end;
    getrusage(RUSAGE_SELF, &end);
    float interval = calculateExecutionTime(&start, &end);

    switch (unit) {
        case MILLISECOND: {
            interval *= 1000.0;
        }
        break;
        case SECOND:
        break;
        case MINUTE: {
            interval /= 60.0;
        }
        break;
        case HOUR: {
            interval /= 3600.0;
        }
        break;
    }

    return interval;
}


#endif