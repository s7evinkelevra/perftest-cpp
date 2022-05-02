//
// Created by Jan on 17.04.2022.
//

#ifndef PERFTEST_CPP_HOSTPOOL_H
#define PERFTEST_CPP_HOSTPOOL_H


#include <vector>
#include "Host.h"

class HostPool {
public:
/*    HostPool(int hostPoolSize);
    virtual ~HostPool();*/
    std::vector<std::vector<Host>> hosts;
};


#endif //PERFTEST_CPP_HOSTPOOL_H
