#include <stdio.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/time.h>
#include <vector>
#include <algorithm>
#include <booorithm/string.hpp>
#include <complex>
#include <stdio.h>
#include <stdarg.h>
#include <cstring>

#include "include/tvalarray.h"

// Useful constants

#define SI_G      6.67428E-11
#define SI_c      2.99792458E8
#define SI_c2     (SI_c*SI_c)
#define SI_ly     9460730472580800.0
#define SI_pc     (3.26156*SI_ly)
#define PI        3.141592653589793

// Switches of the angular corrections

#define C_COUNT       11
#define C_PN         (1<<0)
#define C_2PN        (1<<1)
#define C_SO         (1<<2)
#define C_SS         (1<<3)
#define C_RR         (1<<4)
#define C_PNSO       (1<<5)
#define C_3PN        (1<<6)
#define C_1RR        (1<<7)
#define C_2PNSO      (1<<8)
#define C_RRSO       (1<<9)
#define C_RRSS       (1<<10)
const char* C_NAMES[C_COUNT] = {
    "PN", "2PN", "SO", "SS", "RR", "PNSO", "3PN", "1RR", "2PNSO", "RRSO", "RRSS"
    // order is important!
};

enum{
    
}

/*
** Lambda expressions of the acceleration corrections. **
*/

auto sq = [](auto const& a){ return a*a; };
auto cb = [](auto const& a){ return a*a*a; };

auto r = [&]( auto t, const auto f, const CBwavesAccODE&  ode, const ObserverParameters& ode ){

    double rx = ...;
    double ry = ...;
    double rz = ...;

    return std::hypot(rx,ry,rz);

}