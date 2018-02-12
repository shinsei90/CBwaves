#include <Config.hpp>
#include <InitParams.hpp>
#include <PDE.hpp>
#include <Vector.hpp>

// Type aliases
using state = PDE::StateVector<Vector<mass, 3>>;

struct dynamicalParams {

            dynamicalParams(const state& state_, initParams iparams_){

            Vector<mass, 3> rr = state_.get<Radius>() - iparams_.r_init;
            Vector<mass, 3> x1 = {iparams_.r0/2, 0, 0};
            Vector<mass, 3> x2 = {-iparams_.r0/2, 0, 0};
            Vector<mass, 3> v = {0, SI_c/3, 0};

            Vector<mass, 3> Spin1 = spin1(iparams_.m, state_.get<Radius>());
            Vector<mass, 3> Spin2 = spin2(iparams_.m, state_.get<Radius>());

            Vector<mass, 3> x = x1 - x2;
            Vector<mass, 3> n = x/r;
            Vector<mass, 3> sigma = (m2/m1)*Spin1 + (m1/m2)*Spin2;
            Vector<mass, 3> LN = mu*cross(x, v);
            Vector<mass, 3> Delta = m*(Spin2/m2 - Spin1/m1);
            Vector<mass, 3> Spin = Spin1 + Spin2; 

            mass rdot = dot(n, v);
            mass r = length(rr);
            }

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