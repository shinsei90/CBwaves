#pragma once

#include <Config.hpp>
#include <Vector.hpp>


struct initParams{  //Constant parameters of the equations.

        initParams(mass m1_, mass m2_, mass r0_){
            //Mass ratio
        m1 = m1_;
        m2 = m2_;
        m = m1 + m2;

        //Initial separation
        r0 = r0_;
        r_init = {r0, 0.0, 0.0};


        dm = m1 - m2;
        mu = m1*m2/m;
        eta = mu/m;
        }
      
        //Mass ratio
        mass m1;
        mass m2;
        mass m;

        //Initial separation
        mass r0;
        Vector<mass, 3>  r_init;
        Vector<mass, 3> N;


        mass dm;
        mass mu;
        mass eta;

    };