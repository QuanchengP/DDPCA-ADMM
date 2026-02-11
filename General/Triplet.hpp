#ifndef _Triplet_hpp
#define _Triplet_hpp

#include "General.hpp"

#include <iostream>

namespace Ddpca {

//POD 类型（Plain Old Data 类型）
class Triplet {
public:
    // 公共成员变量，避免不必要的访问器函数调用开销
    I64 row;    // 行索引
    I64 col;    // 列索引
    Real val;        // 值
};

// 非成员函数：用于排序（先按行，再按列）
inline bool operator<(const Triplet& lhs, const Triplet& rhs) noexcept {
    if (lhs.row != rhs.row) {
        return lhs.row < rhs.row;
    }
    return lhs.col < rhs.col;
}

// 非成员函数：相等比较
inline bool operator==(const Triplet& lhs, const Triplet& rhs) noexcept {
    return lhs.row == rhs.row && lhs.col == rhs.col && lhs.val == rhs.val;
}

// 非成员函数：不相等比较
inline bool operator!=(const Triplet& lhs, const Triplet& rhs) noexcept {
    return !(lhs == rhs);
}

// 非成员函数：交换两个三元组
inline void swap(Triplet& a, Triplet& b) noexcept {
    std::swap(a.row, b.row);
    std::swap(a.col, b.col);
    std::swap(a.val, b.val);
}

// 输出操作符重载
inline std::ostream& operator<<(std::ostream& os, const Triplet& t) {
    os << "(" << t.row << ", " << t.col << ", " << t.val << ")";
    return os;
}

} // namespace Ddpca

#endif // _Triplet_hpp