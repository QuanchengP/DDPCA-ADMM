#include "ThreadManager.hpp"

namespace Ddpca {

void ThreadManager::Establish(I64 numbThreads_){
    //
    numbThreads = numbThreads_;
    stopFlag = false;
    //
    one2oneS2S.clear();
    one2oneS2S.resize(numbThreads);
    for(I64 ti = 0; ti < numbThreads; ++ ti){
        one2oneS2S[ti].first = ti;
        one2oneS2S[ti].second.emplace_back(ti);
    }
    //
    Task.resize(maxNestLevel * numbThreads);
    taskDoneSemaphore.resize(maxNestLevel * numbThreads);
    taskSet.resize(maxNestLevel * numbThreads);
    nestFirst.resize(maxNestLevel * numbThreads);
    enableSemaphore.resize(maxNestLevel * numbThreads);
    for (I64 ti = 0; ti < maxNestLevel * numbThreads; ++ ti) {
        taskSet[ti].clear();
        nestFirst[ti] = -1;
        taskDoneSemaphore[ti] = std::make_unique<std::counting_semaphore<maxSemaphores>>(0);
        enableSemaphore[ti] = std::make_unique<std::counting_semaphore<maxSemaphores>>(0);
    }

    for (I64 ti = 0; ti < maxNestLevel * numbThreads; ++ ti) {
        threads.emplace_back([this, ti]() { (this->ThreadLoop)(ti); });
    }

    for (I64 ti = 0; ti < maxNestLevel * numbThreads; ++ ti) {
        startedSemaphore.acquire();
    }
    Log("    " + std::to_string(maxNestLevel * numbThreads) 
        + " threads are started, waiting for tasks...");
}

void ThreadManager::ThreadLoop(I64 threadId) {
    //
    startedSemaphore.release();
    //
    pthread_t tid = pthread_self();
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET((threadId < numbThreads) ? threadId : threadId - numbThreads, &cpuset);
    if (pthread_setaffinity_np(tid, sizeof(cpuset), &cpuset) != 0) {
        Log("    Thread binding to core " + std::to_string(threadId) + " failed!");
    }

    while (true) {
        (* enableSemaphore[threadId]).acquire();
        
        if (stopFlag) {
            Log("    Thread " + std::to_string(threadId) + " is exiting...");
            break;
        }

        for(auto& idId : taskSet[threadId]){
            Task[threadId](idId);
        }

        taskSet[threadId].clear();
        (* taskDoneSemaphore[nestFirst[threadId]]).release();
    }
}

void ThreadManager::RunTask(
    const I64 nestLevel, 
    std::vector<std::pair<I64, std::vector<I64>>> threadTask, 
    const std::function<void(I64)>& taskFunction) {
    //
    I64 numbTasks = threadTask.size();
    for(I64 ti = 0; ti < numbTasks; ++ ti){
        threadTask[ti].first += nestLevel * numbThreads;
    }
    //
    I64 firstThread = threadTask[0].first;
    for(I64 ti = 0; ti < numbTasks; ++ ti){
        I64 tid = threadTask[ti].first;
        taskSet[tid].assign(threadTask[ti].second.begin(), threadTask[ti].second.end());
        Task[tid] = taskFunction;
        nestFirst[tid] = firstThread;
        (* enableSemaphore[tid]).release();
    }

    for(I64 ti = 0; ti < numbTasks; ++ ti){
        (* taskDoneSemaphore[firstThread]).acquire();
    }

    if(numbTasks > 16){
        Log("        " + std::to_string(numbTasks) + " tasks are done.");
    }
}

void ThreadManager::Stop() {
    //
    stopFlag = true;
    I64 numbThreads_ = threads.size();
    for (I64 ti = 0; ti < numbThreads_; ti++) {
        (* enableSemaphore[ti]).release();
    }
    for (auto& t : threads) {
        if (t.joinable()) {
            t.join();
        }
    }
}

void ThreadManager::ThreadDistribute(const I64 threadsPerTask, const I64 numbTasks, 
    std::vector<std::pair<I64, std::vector<I64>>>& threadTask, 
    std::vector<std::vector<std::pair<I64, std::vector<I64>>>>& ttS2M){
    //
    //
    threadTask.clear();
    ttS2M.clear();
    if(threadsPerTask * numbTasks <= numbThreads){
        //
        threadTask.resize(numbTasks);
        ttS2M.resize(numbTasks);
        I64 tempThread = 0;
        for(I64 ti = 0; ti < numbTasks; ++ ti){
            threadTask[ti].first = tempThread;
            threadTask[ti].second.emplace_back(ti);
            //
            ttS2M[ti].resize(threadsPerTask);
            for(I64 tj = 0; tj < threadsPerTask; ++ tj){
                ttS2M[ti][tj].first = tempThread + tj;
                ttS2M[ti][tj].second.clear();
                ttS2M[ti][tj].second.emplace_back(tj);
            }
            tempThread += threadsPerTask;
        }
    }
    else{
        I64 numbBox = numbThreads / threadsPerTask;
        threadTask.resize(numbBox);
        ttS2M.resize(numbTasks);
        for(I64 ti = 0; ti < numbTasks; ++ ti){
            ttS2M[ti].resize(threadsPerTask);
        }
        I64 numbTaskPerBox = numbTasks / numbBox;
        I64 remainderBox = numbTasks % numbBox;
        I64 remainderIndex = remainderBox * (numbTaskPerBox + 1);
        for(I64 ti = 0; ti < numbBox; ++ ti){
            threadTask[ti].first = ti * threadsPerTask;
            I64 tempTask, maxTj;
            if(ti < remainderBox){
                maxTj = numbTaskPerBox + 1;
                tempTask = ti * maxTj;
            }
            else{
                tempTask = remainderIndex + (ti - remainderBox) * numbTaskPerBox;
                maxTj = numbTaskPerBox;
            }
            for(I64 tj = 0; tj < maxTj; ++ tj){
                threadTask[ti].second.emplace_back(tempTask);
                for(I64 tk = 0; tk < threadsPerTask; ++ tk){
                    ttS2M[tempTask][tk].first = threadTask[ti].first + tk;
                    ttS2M[tempTask][tk].second.clear();
                    ttS2M[tempTask][tk].second.emplace_back(tk);
                }
                ++ tempTask;
            }
        }
    }
}

} // namespace Ddpca