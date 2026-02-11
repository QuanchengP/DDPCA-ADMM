#include "../General/General.hpp"

// ========== 2. 计算任务数据（主线程管理的vector） ==========
class Compute{

public:

    static const Ddpca::I64 MAX_TASK_CAPACITY = 40; // 预分配vector容量（足够容纳所有任务）
    Ddpca::I64 custom_int[MAX_TASK_CAPACITY];          // 整型参数
    std::string custom_str[MAX_TASK_CAPACITY];  // 字符串参数
    double custom_double[MAX_TASK_CAPACITY];    // 浮点参数
    Ddpca::I64 result[MAX_TASK_CAPACITY];              // 任务结果

public:

    void Initialize(){
        //
        for(Ddpca::I64 tid = 0; tid < Compute::MAX_TASK_CAPACITY; ++ tid){
            custom_int[tid] = 1;
            custom_str[tid] = "测试";
            custom_double[tid] = 3.14159 + tid;
            result[tid] = 1;
        }
    }

    // Task1：直接访问全局vector下标，无临时TaskData，无锁
    void Task1(Ddpca::I64 taskIndex) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        // 直接操作全局数据（读）
        result[taskIndex] = result[taskIndex] * 2; // 直接写结果到全局vector
    }

    // Task2：逻辑同Task1，直接下标访问全局vector
    void Task2(Ddpca::I64 taskIndex) {
        std::this_thread::sleep_for(std::chrono::milliseconds(150));
        // 直接操作全局数据（读）
        result[taskIndex] = result[taskIndex] * 3;
    }

};

int main(int argc, char **argv) {
    //
    Ddpca::Initialize(argc, argv);

    Compute comp;
    comp.Initialize();

    std::vector<std::pair<Ddpca::I64, std::vector<Ddpca::I64>>> threadTask_0 = {
        {0, {0,1}}, {2, {2,3}}, {4, {4,5}}, {6, {6,7}}, {8, {8,9}},
        {10, {10,11}}, {12, {12,13}}, {14, {14,15}}, {16, {16,17}}, {18, {18,19}}
    };
    std::vector<std::pair<Ddpca::I64, std::vector<Ddpca::I64>>> threadTask_1 = {
        {0, {0,1}}, {1, {2,3}}, {2, {4,5}}, {3, {6,7}}, {4, {8,9}},
        {5, {10,11}}, {6, {12,13}}, {7, {14,15}}, {8, {16,17}}, {9, {18,19}},
        {10, {20,21}}, {11, {22,23}}, {12, {24,25}}, {13, {26,27}}, {14, {28,29}},
        {15, {30,31}}, {16, {32,33}}, {17, {34,35}}, {18, {36,37}}, {19, {38,39}}
    };
    // 步骤5：主线程循环执行任务（全程管理vector）
    for (Ddpca::I64 round = 1; round <= 5; ++ round) {
        std::cout << "\n==================== 第" << round << "轮 ====================\n";
        //
        Ddpca::threadManager.RunTask(0, threadTask_0, 
            std::bind(&Compute::Task1, &comp, std::placeholders::_1));
        std::cout << comp.result[0] << "\n";
        //
        Ddpca::threadManager.RunTask(0, threadTask_1, 
            std::bind(&Compute::Task2, &comp, std::placeholders::_1));
        std::cout << comp.result[0] << "\n";
    }
    
    std::cout << "\n所有轮次执行完毕！\n";
    for(Ddpca::I64 tid = 0; tid < Compute::MAX_TASK_CAPACITY; ++ tid){
        std::cout << comp.result[tid] << " ";
    }
    std::cout << "\n";

    Ddpca::Finalize();
    return 0;
}