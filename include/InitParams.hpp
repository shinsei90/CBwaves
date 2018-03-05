#pragma once

#include <Config.hpp>
#include <Vector.hpp>


struct initParams{  //Constant parameters of the equations.

        initParams(mass m1_, mass m2_, mass r0_, mass ecc_){
        
        //Mass ratio
        m1 = m1_*MSUN;
        m2 = m2_*MSUN;
        m = m1 + m2;

        //Initial separation
        r0 = r0_; //Initial Separation of the objects
        r_init = {r0, 0.0, 0.0}; //Initial Separation vector of the objects


        dm = m1 - m2;
        mu = m1*m2/m;
        eta = mu/m;

        //Excenricity
        ecc = ecc_;

        //Initial velocity
        v_init = {0. , m/r0_*(1. - ecc_)};

        }
      
        //Mass ratio
        mass m1;
        mass m2;
        mass m;

        //Initial separation
        mass r0;
        Vector<mass, 3> r_init;
        Vector<mass, 3> N;


        mass dm;
        mass mu;
        mass eta;

        mass ecc;
        Vector<mass, 3> v_init;

    };