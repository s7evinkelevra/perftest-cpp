//
// Created by Jan on 20.04.2022.
//

#ifndef PERFTEST_CPP_MERITCACHE_H
#define PERFTEST_CPP_MERITCACHE_H

#include <deque>

class MeritCache {
public:
    std::deque<std::deque<int>> cache;
};


#endif //PERFTEST_CPP_MERITCACHE_H
