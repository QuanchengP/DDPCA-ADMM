#ifndef _Coordinate_hpp
#define _Coordinate_hpp

#include "../General/General.hpp"

#include <cassert>
#include <numeric>   // for std::transform_reduce
#include <execution> // for std::execution::unseq

namespace Ddpca {

template <typename Container>
inline Container Cross(const Container& x, const Container& y);

/****************************************************************************************************/
//3D coordinate class with comparison capabilities
class Coordinate {
public:
    static constexpr I64 dimension = 3;
    std::array<Real,dimension> data;

    static constexpr Real coorError = 1.0E-10;
    
public:
    // Default constructor - fully inline
    Coordinate() noexcept : data{{0.0, 0.0, 0.0}} {}
    
    // Constructor with coordinates - fully inline
    Coordinate(Real x, Real y, Real z) noexcept : data{{x, y, z}} {}
    
    // Copy constructor - fully inline
    Coordinate(const Coordinate& other) = default;
    
    // Move constructor
    Coordinate(Coordinate&& other) noexcept = default;
    
    // Copy assignment operator
    Coordinate& operator=(const Coordinate& other) = default;
    
    // Move assignment operator
    Coordinate& operator=(Coordinate&& other) noexcept = default;
    
    // Destructor
    ~Coordinate() noexcept = default;

public:
    
    // Overload operator[] for index access
    inline Real& operator[](I64 index) {
        assert(index >= 0 && index < dimension && "Vector index out of range");
        return data[index];
    }

    inline const Real& operator[](I64 index) const {
        assert(index >= 0 && index < dimension && "Vector index out of range");
        return data[index];
    }

    // Comparison operator
    bool operator<(const Coordinate &other) const {
        if (std::abs(data[0] - other.data[0]) > coorError) {
            return data[0] < other.data[0];
        }
        if (std::abs(data[1] - other.data[1]) > coorError) {
            return data[1] < other.data[1];
        }
        return data[2] < other.data[2] - coorError;
    }

    // Create cross product of 3D vectors (only applicable for 3-element vectors)
    inline Coordinate Cross(const Coordinate& other) const {
        Coordinate result;
        result.data = Ddpca::Cross(data, other.data);
        return result;
    }
};

template <typename Container>
inline Container Cross(const Container& x, const Container& y) {
    Container result;
    result[0] = x[1] * y[2] - x[2] * y[1];
    result[1] = x[2] * y[0] - x[0] * y[2];
    result[2] = x[0] * y[1] - x[1] * y[0];
    return result;
}

template <typename Container>
inline Real NRM2(const Container& x) noexcept {
    return std::sqrt(std::transform_reduce(std::execution::unseq, 
        x.begin(), x.end(),    // 输入范围
        0.0,                   // 初始值
        std::plus<>(),         // 累加操作
        [](const auto& val) { return val * val; } // 转换函数：计算平方
    ));
}

//x += y
template<typename Container>
inline void XPEY(Container& x, const Container& y) noexcept {
    std::transform(std::execution::unseq, 
        x.begin(), x.end(), 
        y.begin(), 
        x.begin(), 
        [](const auto& x_val, const auto& y_val) { return x_val + y_val; });
}

//y += a x
template<typename Container>
inline void AXPY(const Real& alpha, const Container& x, Container& y) noexcept {
    std::transform(std::execution::unseq, 
        x.begin(), x.end(), 
        y.begin(), 
        y.begin(), 
        [alpha](const auto& x_val, const auto& y_val) { return alpha * x_val + y_val; });
}

//z = a x + y
template<typename Container>
inline void AXPY(const Real& alpha, const Container& x, Container& y, Container& z) noexcept {
    std::transform(std::execution::unseq, 
        x.begin(), x.end(), 
        y.begin(), 
        z.begin(), 
        [alpha](const auto x_val, const auto y_val) { return alpha * x_val + y_val; });
}

//x *= alpha
template <typename Container>
inline void SCAL(const Real& alpha, Container& x) noexcept {
    std::transform(std::execution::unseq, 
        x.begin(), x.end(), x.begin(),
        [alpha](const auto& val) { return alpha * val; });
}

//x^T * y
template<typename Container>
inline Real DOT(const Container& x, const Container& y) noexcept {
    return transform_reduce(std::execution::unseq, 
        x.begin(), x.end(),
        y.begin(), 0.0,
        std::plus<>(), 
        [](const auto& x_val, const auto& y_val) { return x_val * y_val; });
}

template<>
inline Real DOT(const std::array<Real, 3>& x, const std::array<Real, 3>& y) noexcept {
    return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

//x /= ||x||
template <typename Container>
inline void NORMALIZE(Container& x) noexcept {
    Real norm = NRM2(x);
    std::transform(std::execution::unseq, 
        x.begin(), x.end(), x.begin(),
        [norm](const auto& val) { return val / norm; });
}

} // namespace Ddpca

#endif // _Coordinate_hpp