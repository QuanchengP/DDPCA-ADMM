#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#include <string.h>

#define ONE_GB (1UL << 30)  // 1G字节 = 2^30
#define TEST_FILE "/hugepages_1g/hugepage_test_1g"
// 低内核版本替代MAP_HUGE_1GB的宏（1G大页的标识位）
#define MAP_HUGE_1GB (21 << 26)

int main() {
    int fd;
    void *addr;

    // 1. 打开/创建1G大页文件（必须在1G hugetlbfs目录下）
    fd = open(TEST_FILE, O_CREAT | O_RDWR, 0755);
    if (fd < 0) {
        perror("open TEST_FILE failed");
        return EXIT_FAILURE;
    }

    // 2. 扩展文件大小到1G（匹配1G大页尺寸）
    if (ftruncate(fd, ONE_GB) < 0) {
        perror("ftruncate failed");
        close(fd);
        return EXIT_FAILURE;
    }

    // 3. 映射1G大页（用(1<<21)替代MAP_HUGE_1GB）
    addr = mmap(NULL, ONE_GB,
                PROT_READ | PROT_WRITE,  // 读写权限
                MAP_SHARED | MAP_HUGETLB | MAP_HUGE_1GB,  // 核心标识
                fd, 0);
    if (addr == MAP_FAILED) {
        perror("mmap 1G hugepage failed");
        close(fd);
        return EXIT_FAILURE;
    }

    // 4. 验证：写入+读取数据，确认内存可用
    char *test_str = "Test 1GB HugePage Success! (Low Kernel Version)";
    memcpy(addr, test_str, strlen(test_str) + 1);  // 更安全的写入方式
    printf("✅ 1G大页映射成功！\n");
    printf("   内存地址：%p\n", addr);
    printf("   写入内容：%s\n", (char*)addr);

    // 5. 释放资源（避免内存泄漏）
    munmap(addr, ONE_GB);
    close(fd);
    unlink(TEST_FILE);  // 删除临时文件

    return EXIT_SUCCESS;
}