#include <Corrections.hpp>
#include <DynamicalParams.hpp>

dynamicalParams::dynamicalParams(const state& state_, initParams ip_){

            rr = state_.get<Radius>() /*- ip_.r_init*/; //Relative Separation vector of the objects
            r1 = {ip_.m2/ip_.m*rr[0], ip_.m2/ip_.m*rr[1], ip_.m2/ip_.m*rr[2]};                  //position of m1
            r2 = {ip_.m1/ip_.m*rr[0], ip_.m1/ip_.m*rr[1], ip_.m1/ip_.m*rr[2]};                 //position of m2
            r = length(rr); //Relative Separation of the objects
            
            x = r1 - r2;
            n = {rr[0]/r, rr[1]/r, rr[2]/r};
            v = state_.get<Velocity>();
            rdot = dot(n, v);
            LN = ip_.mu*cross(x, v);

            Spin1 = spin1(*this, ip_);
            Spin2 = spin2(*this, ip_);
            Spin = Spin1 + Spin2;
            sigma = (ip_.m2/ip_.m1)*Spin1 + (ip_.m1/ip_.m2)*Spin2;
            Delta = ip_.m*(Spin2/ip_.m2 - Spin1/ip_.m1);
            }