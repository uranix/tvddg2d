#pragma once

#ifndef NDEBUG

#include <iostream>

template<typename T>
volatile double TOUCH(T v) {
    volatile double res;
    if (v.sum() > 0)
        res = v.sum() * v.sum();
    else
        res = 1 - v.sum();
    return res;
}

#endif
