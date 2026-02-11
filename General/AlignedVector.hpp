#ifndef _AlignedVector_hpp
#define _AlignedVector_hpp
#include "General.hpp"
#include "AlignedAllocate.hpp"

#include <vector>

namespace Ddpca {

// 自定义对齐分配器
template <typename T, size_t Alignment = nfsAlign>
struct AlignedAllocator {
    using value_type = T;
    using size_type = I64;
    using difference_type = I64;
    using pointer = T*;
    using const_pointer = const T*;
    using reference = T&;
    using const_reference = const T&;

    // 模板别名，用于rebind
    template <typename U>
    struct rebind {
        using other = AlignedAllocator<U, Alignment>;
    };

    // 默认构造函数
    constexpr AlignedAllocator() noexcept = default;

    // 拷贝构造函数
    constexpr AlignedAllocator(const AlignedAllocator& other) noexcept = default;

    // 不同类型的拷贝构造函数
    template <typename U>
    constexpr AlignedAllocator(const AlignedAllocator<U, Alignment>& other) noexcept {}

    // 分配内存：使用AlignedAllocate.hpp中的AlignedAllocate函数
    T* allocate(size_t n) {
        return AlignedAllocate<T>(n);
    }

    // 释放内存：使用AlignedAllocate.hpp中的Deallocate函数
    void deallocate(T* ptr, size_t n) {
        Deallocate<T>(ptr, n);
    }

    // 相等性判断
    template <typename U>
    constexpr bool operator==(const AlignedAllocator<U, Alignment>&) const noexcept {
        return true;
    }

    // 不等性判断
    template <typename U>
    constexpr bool operator!=(const AlignedAllocator<U, Alignment>&) const noexcept {
        return false;
    }
};

// 使用对齐分配器的 vector（64 字节对齐）
using AlignedVectorRx = std::vector<Real, AlignedAllocator<Real, nfsAlign>>;

} // namespace Ddpca

#endif // _AlignedVector_hpp