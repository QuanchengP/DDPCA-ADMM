#ifndef _General_hpp
#define _General_hpp

#include "ThreadManager.hpp"

#include <cmath>
#include <array>

namespace Ddpca{

/****************************************************************************************************/

inline ThreadManager threadManager;

inline void EvenlyDistribute(const I64 totalSize, const I64 numbPart, 
    std::vector<I64>& startIndex, std::vector<I64>& endIndex){
    //
    I64 partSize = totalSize / numbPart;
    I64 remainderPart = totalSize % numbPart;
    I64 remainderIndex = remainderPart * (partSize + 1);
    //
    for(I64 tp = 0; tp < numbPart; tp ++){
        if(tp < remainderPart){
            startIndex[tp] = tp * (partSize + 1);
            endIndex[tp] = startIndex[tp] + partSize + 1;
        }
        else{
            startIndex[tp] = remainderIndex + (tp - remainderPart) * partSize;
            endIndex[tp] = startIndex[tp] + partSize;
        }
    }
}

inline void Initialize(I32 argc, char **argv){
    //
    std::cout << "Console output is redirected to testLogg.txt\n";
    //
	std::cout << std::setiosflags(std::ios::scientific) << std::setprecision(20);
    // std::hardware_destructive_interference_size
    Log("    Cache line size: " + std::to_string(nfsAlign) + " bytes");
    //
    I64 numbThreads;
    if(argc == 4) {
        numbThreads = std::stoi(argv[1]);
        threadManager.threadsPerDomain = std::stoi(argv[2]);
        threadManager.threadsPerInterface = std::stoi(argv[3]);
    }
    else {
        numbThreads = std::thread::hardware_concurrency() / 2;
        threadManager.threadsPerDomain = 1;
        threadManager.threadsPerInterface = 1;
        Log("    Default number of threads: " + std::to_string(numbThreads));
    }
    threadManager.Establish(numbThreads);
}

inline void Finalize(){
    threadManager.Stop();
}

} // namespace Ddpca

#endif // _General_hpp