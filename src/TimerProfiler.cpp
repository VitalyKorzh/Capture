#include "TimeProfiler.h"

std::map<std::string, double> TimeProfiler::timings;

void TimeProfiler::print(std::ostream &os)
{
        os << "\n# === Time Profiling Results (ms) ===\n";
        os << std::fixed << std::setprecision(3);  // 3 знака после запятой
        for (const auto& it : timings) {
            os << "# " << std::setw(25) << std::left << it.first
                      << ": " << std::setw(10) << std::right << it.second << " ms\n";
        }
        os << std::resetiosflags(std::ios::fixed);  // Сброс форматирования
}