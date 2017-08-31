#include <iostream>

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