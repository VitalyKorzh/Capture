#ifndef __TIME_PROFILER_H__
#define __TIME_PROFILER_H__

#include <chrono>
#include <string>
#include <map>
#include <iostream>
#include <ostream>
#include <iomanip>

class TimeProfiler {
    static std::map<std::string, double> timings;
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::string name;
    
public:
    TimeProfiler(const std::string& name) : name(name) {
        start = std::chrono::high_resolution_clock::now();
    }
    
    ~TimeProfiler() {
        auto end = std::chrono::high_resolution_clock::now();
        double duration = std::chrono::duration<double, std::milli>(end - start).count();
        timings[name] += duration;
    }
    
    static void print(std::ostream &os);
    
    static void reset() {
        timings.clear();
    }
};

#endif