#ifndef _AlignedAllocate_hpp
#define _AlignedAllocate_hpp
#include "Utility.hpp"

#include <sys/mman.h>
#include <unistd.h>
#include <cerrno>
#include <cstring>
#include <iostream>
#include <stdexcept>

namespace Ddpca {

// static constexpr size_t HUGEPAGE_SIZE = 2 * 1024 * 1024; // 2MB

template <typename T>
T* AlignedAllocate(I64 n) {
    // 计算需要的内存大小
    // const size_t size = n * sizeof(T);
    // // 向上对齐到 2MB
    // size_t aligned_size = ((size + HUGEPAGE_SIZE - 1) / HUGEPAGE_SIZE) * HUGEPAGE_SIZE;
    // void* ptr = mmap(
    //     nullptr,
    //     aligned_size,
    //     PROT_READ | PROT_WRITE,
    //     MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB,
    //     -1,
    //     0
    // );
    // if (ptr == MAP_FAILED) {
    //     std::cerr << "❌ mmap(MAP_HUGETLB) failed for "
    //                 << aligned_size / (1024*1024) << " MB: "
    //                 << strerror(errno) << std::endl;
    //     return nullptr;
    // }
    // // Linux guarantees MAP_HUGETLB returns HUGEPAGE_SIZE-aligned address
    // // 所以无需手动对齐 ptr
    // return static_cast<T*>(ptr);

    // return new T[n];
    
    const I64 size = n * sizeof(T);
    const I64 aligned_size = ((size + nfsAlign - 1) / nfsAlign) * nfsAlign;
    void* ptr = std::aligned_alloc(nfsAlign, aligned_size);
    if (!ptr) throw std::bad_alloc();
    return static_cast<T*>(ptr);
}

template <typename T>
void Deallocate(T* ptr, I64/* n*/) {
    // if (!ptr) return;
    // const size_t size = n * sizeof(T);
    // size_t aligned_size = ((size + HUGEPAGE_SIZE - 1) / HUGEPAGE_SIZE) * HUGEPAGE_SIZE;
    // if (munmap(ptr, aligned_size) != 0) {
    //     perror("munmap failed");
    // }

    // delete[] ptr;

    std::free(static_cast<void*>(ptr));
}

} // namespace Ddpca

#endif // _AlignedAllocate_hpp