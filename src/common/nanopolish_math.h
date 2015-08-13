//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_math -- math helper/wrapper functions
//
#ifndef NANOPOLISH_MATH_H
#define NANOPOLISH_MATH_H

#include <assert.h>

template<typename T>
T mean(const std::vector<T>& values)
{
    T sum = 0;
    for(size_t i = 0; i < values.size(); ++i) {
        sum += values[i];
    }
    return sum / values.size();
}

template<typename T>
T median(const std::vector<T>& values)
{
    assert(!values.empty());

    std::vector<T> tmp = values;
    std::sort(tmp.begin(), tmp.end());

    return tmp.size() % 2 == 1 ?
                tmp[tmp.size() / 2] :
                (tmp[tmp.size() / 2] + tmp[(tmp.size() + 1) / 2]) / 2;
}


#endif
