#ifndef _ThreadManager_hpp
#define _ThreadManager_hpp

#include "Utility.hpp"

#include <vector>
#include <thread>
#include <mutex>
#include <functional>
#include <semaphore>
#include <sched.h>

namespace Ddpca {

// Thread manager designed only for Ddpca: threadsPerDomain, threadsPerInterface,
// Support nested parallelism: at most two levels.
class ThreadManager {

public:

    I64 numbThreads;
    constexpr static I64 maxNestLevel = 2;
    constexpr static I64 maxSemaphores = std::counting_semaphore<>::max();

    std::counting_semaphore<maxSemaphores> startedSemaphore{0};

    std::vector<std::unique_ptr<std::counting_semaphore<maxSemaphores>>> taskDoneSemaphore;

    std::vector<std::thread> threads;
    //std::unique_ptr: no need to delete
    std::vector<std::unique_ptr<std::counting_semaphore<maxSemaphores>>> enableSemaphore;

    bool stopFlag;
    std::vector<std::vector<I64>> taskSet; // no atomic: be careful!
    std::vector<I64> nestFirst;
    std::vector<std::function<void(I64)>> Task; // no atomic: be careful!

public:

    void Establish(I64 numbThreads_);

    // threadId must be less than the number of CPU cores, 
    // because each thread will be bound to a core
    void ThreadLoop(I64 threadId);

    void RunTask(
        const I64 nestLevel, 
        std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
        const std::function<void(I64)>& taskFunction);

    void Stop();

    std::mutex printMutex;
    void SafePrint(const std::string& message){
        std::lock_guard<std::mutex> lock(printMutex);
        std::cout << message << std::endl;
    }

public:

    I64 threadsPerDomain;
    I64 threadsPerInterface;
    I64 threadsPerLevel;

    //S2S: single task(domain/interface) to single thread
    std::vector<std::pair<I64, std::vector<I64>>> domainS2S;
    std::vector<std::pair<I64, std::vector<I64>>> interfaceS2S;
    std::vector<std::pair<I64, std::vector<I64>>> levelS2S; // only used for NDD problem
    //S2M: single task(domain/interface) to multiple threads
    std::vector<std::vector<std::pair<I64, std::vector<I64>>>> domainS2M;
    std::vector<std::vector<std::pair<I64, std::vector<I64>>>> interfaceS2M;
    std::vector<std::vector<std::pair<I64, std::vector<I64>>>> levelS2M; // only used for NDD problem
    std::vector<std::pair<I64, std::vector<I64>>> one2oneS2S;

    void ThreadDistribute(const I64 threadsPerTask, const I64 numbTasks, 
        std::vector<std::pair<I64, std::vector<I64>>>& threadTask, 
        std::vector<std::vector<std::pair<I64, std::vector<I64>>>>& domainS2M);

}; // class ThreadManager

} // namespace Ddpca

#endif // _ThreadManager_hpp