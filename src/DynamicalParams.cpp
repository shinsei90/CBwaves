#include <Corrections.hpp>
#include <DynamicalParams.hpp>

dynamicalParams::dynamicalParams(const state& state_, initParams ip_){

            Vector<mass, 3> rr = state_.get<Radius>() - ip_.r_init;
            Vector<mass, 3> x1 = {ip_.r0/2, 0, 0};
            Vector<mass, 3> x2 = {-ip_.r0/2, 0, 0};
            Vector<mass, 3> v = {0, SI_c/3, 0};

            Vector<mass, 3> x = x1 - x2;
            Vector<mass, 3> n = x/r;
            Vector<mass, 3> sigma = (ip_.m2/ip_.m1)*Spin1 + (ip_.m1/ip_.m2)*Spin2;
            Vector<mass, 3> LN = ip_.mu*cross(x, v);
            Vector<mass, 3> Delta = ip_.m*(Spin2/ip_.m2 - Spin1/ip_.m1);

            mass rdot = dot(n, v);
            mass r = length(rr);


            Vector<mass, 3> Spin1 = spin1(*this, ip_);
            Vector<mass, 3> Spin2 = spin2(*this, ip_);
            Vector<mass, 3> Spin = Spin1 + Spin2; 
            }