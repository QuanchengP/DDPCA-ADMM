#ifndef _Utility_hpp
#define _Utility_hpp

#include <string>
#include <fstream>
#include <iostream>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <cstdint>

namespace Ddpca {

/****************************************************************************************************/
// Integer type aliases
using I32 = int32_t;
using I64 = int64_t;
using Real = double;

inline std::string Double2String(double inpuNumb) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6) << std::noshowpoint;
    oss << inpuNumb;
    return oss.str();
}

/****************************************************************************************************/
constexpr Real PI = 3.141592653589793238462;
//no false sharing
constexpr I64 nfsAlign = 64;//std::hardware_destructive_interference_size;

// Get corresponding Gaussian points and weights based on the number of integration points
template <I64 N> 
struct GaussData {
    static constexpr std::array<Real, N> GetPoints() { return {}; }
    static constexpr std::array<Real, N> GetWeights() { return {}; }
};

// Specialization for 2-point Gaussian quadrature data
template <> 
struct GaussData<2> {
    static constexpr std::array<Real, 2> GetPoints() { 
        //1.0 / std::sqrt(3.0)
        return {-0.5773502691896257645091, 0.5773502691896257645091}; 
    }
    static constexpr std::array<Real, 2> GetWeights() { 
        return {1.0, 1.0}; 
    }
};

// Specialization for 3-point Gaussian quadrature data
template <> 
struct GaussData<3> {
    static constexpr std::array<Real, 3> GetPoints() { 
        //std::sqrt(3.0 / 5.0)
        return {-0.7745966692414833770359, 0.0, 0.7745966692414833770359}; 
    }
    static constexpr std::array<Real, 3> GetWeights() { 
        return {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0}; 
    }
};

/****************************************************************************************************/
// Internal helper class for logging functionality
class TxtLog {
private:
    // File stream object
    std::ofstream logFile;
    
    // Private constructor to enforce singleton pattern
    TxtLog() {
        // Open log file in append mode
        logFile.open("testLogg.txt", std::ios::out | std::ios::app);
        if (!logFile.is_open()) {
            std::cerr << "Failed to open log file: resuLogg.txt" << "\n";
        }
    }
    
    // Get current time as string
    std::string GetCurrentTime() {
        auto now = std::chrono::system_clock::now();
        auto now_c = std::chrono::system_clock::to_time_t(now);
        std::stringstream ss;
        ss << std::put_time(std::localtime(&now_c), "%Y-%m-%d %H:%M:%S");
        return ss.str();
    }
    
public:
    // Destructor
    ~TxtLog() {
        if (logFile.is_open()) {
            logFile.close();
        }
    }
    
    // Log message method
    void Log(const std::string& message) {
        if (logFile.is_open()) {
            logFile << "[" << GetCurrentTime() << "] " << message << "\n";
            logFile.flush(); // Ensure log is written immediately
        }
    }
    
    // Get singleton instance (thread-safe in C++11 and later)
    static TxtLog& GetInstance() {
        static TxtLog instance;
        return instance;
    }
};

// Convenient logging function - this is what users will actually call
inline void Log(const std::string& message) {
    TxtLog::GetInstance().Log(message);
}

} // namespace Ddpca

#endif // _Utility_hpp