#pragma once
#ifndef UTILS_HPP
#define UTILS_HPP

#define MAX_PRINT_NEDGE             (100000)

#include <numeric>
#include <cmath>

using GraphElem = int64_t;
using GraphWeight = double;

struct Edge
{   
    GraphElem tail_;
    GraphWeight weight_;
    
    Edge(): tail_(-1), weight_(-1.0) {}
};

#endif // UTILS
