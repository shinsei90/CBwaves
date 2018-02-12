#pragma once

#include <Config.hpp>
#include <InitParams.hpp>
#include <PDE.hpp>
#include <Vector.hpp>

// Type aliases
using state = PDE::StateVector<Vector<mass, 3>>;

enum Component : int
{
    Radius
};

struct dynamicalParams {

            dynamicalParams(const state& state_, initParams iparams_);

            Vector<mass, 3> rr;
            Vector<mass, 3> x1;
            Vector<mass, 3> x2;
            Vector<mass, 3> v;

            Vector<mass, 3> Spin1;
            Vector<mass, 3> Spin2;

            Vector<mass, 3> x;
            Vector<mass, 3> n;
            Vector<mass, 3> sigma;
            Vector<mass, 3> LN;
            Vector<mass, 3> Delta;
            Vector<mass, 3> Spin; 

            mass rdot;
            mass r;

        };