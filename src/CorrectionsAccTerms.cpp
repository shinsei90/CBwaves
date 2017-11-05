#include <Corrections.hpp>

#include <iostream>
#include <cmath>

dm = m1 - m2;
mass m = m1 + m2;
mass mu = m1*m2/m;
mass eta = mu/m;
mass rdot = dot(n, v); 

Vector<mass, 3> x1 = {r/2, 0, 0};
Vector<mass, 3> x2 = {-r/2, 0, 0};
Vector<mass, 3> v = {0, SI_c/3, 0};

Vector<mass, 3> x = x1 - x2;
Vector<mass, 3> n = x/r;
Vector<mass, 3> sigma = (m2/m1)*Spin1 + (m1/m2)*Spin2;
Vector<mass, 3> LN = mu*cross(x, v);
Vector<mass, 3> Delta = m*(Spin2/m2 - Spin1/m1);


auto sq = [](auto const& a){ return a*a;

};
auto cb = [](auto const& a){ return a*a*a; };



Vector<mass, 3> c_Newtonian(mass m, mass r){

    return Vector<mass, 3> rN = -(SI_G*m)/sq(r) * n;

}
Vector<mass, 3> c_PostNewtonian(mass m, mass r){

    Vector<mass,3> rPN = -(SI_G*m)/(SI_c2*cb(r))*(n*((1 + 3*eta)*sq(length(v)) - 2*(2 + eta)*(SI_G*m)/r - 3/2*eta*rdot) - 2*(2 - eta)*rdot*v);
    return rPN;

}
Vector<mass, 3> c_2PostNewtonian(mass m, mass r){

    Vector<mass, 3> r2PN = -(SI_G*m)/(sq(SI_c2*r)) * (n*(3/4*(12 + 29*eta)*sq((SI_G*m)/r) + eta*(3 - 4*eta)*std::pow(length(v),4) 
                         + 15/8*eta*(1 - 3*eta)*satd::pow(rdot,4) - 3/2*eta*(3 - 4*eta)*sq(length(v)*rdot) 
                         - 1/2*eta*(13 - 4*eta)*(SI_G*m)/r*sq(length(v)) - (2 + 25*eta + 2*sq(eta))(SI_G*m)/r*sq(rdot)) 
                         - 1/2*rdot*v*(eta*(15 + 4*eta)*sq(length(v)) - (4 + 41*eta + 8*sq(eta))*(SI_G*m)/r -3*eta*(3 + 2*eta)*sq(rdot)));
    return r2PN;

}
Vector<mass, 3> c_3PostNewtonian(mass m, mass r){

    Vector<mass, 3> r3PN = (SI_G*m)/(cb(SI_c)*sq(r))*(n*((16 + (1399/12 - 41/16*sq(PI))*eta + 71/2*sq(eta))*cb((SI_G*m)/r) 
                         + eta*(20827/840 + 123/64*sq(PI) - sq(eta))*sq((SI_G*m)/r*length(v)) 
                         - (1 + (22717/168 + 615/64*sq(PI))*eta + 11/8*sq(eta) - 7*cb(eta))*sq((SI_G*m)/r*rdot) 
                         - 1/4*eta*(11 - 49*eta + 52*sq(eta))*std::pow(length(v),6) + 35/16*eta*(1 - 5*eta + 5*sq(eta))*std::pow(rdot,6)
                         - 1/4*eta*(75 + 32*eta - 40*sq(eta))*(SI_G*m)/r*std::pow(length(v),4)
                         - 1/2*eta*(158 - 69*eta - 60*sq(eta))*(SI_G*m)/r*std::pow(rdot,4) 
                         + eta*(121 - 16*eta - 20*sq(eta))*(SI_G*m)/r*sq(length(v)*rdot) 
                         + 3/8*eta*(20 - 79*eta + 60*sq(eta))*sq(sq(length(v))*rdot) 
                         - 15/8*eta*(4 - 18*eta + 17*sq(eta))*sq(length(v)*sq(rdot))) 
                         + rdot*v*((4 + (5849/840 + 123/32*sq(PI))*eta - 25*sq(eta) - 8*cb(eta))*sq((SI_G*m)/r) 
                         + 1/8*eta*(65 - 152*eta - 48*sq(eta))*std::pow(length(v),4) + 15/8*eta*(3 - 8*eta - 2*sq(eta))*std::pow(rdot,4)
                         + eta*(15 + 27*eta + 10*sq(eta))*(SI_G*m)/r*sq(length(v)) 
                         - 11/6*eta*(329 + 177*eta + 108*sq(eta))*(SI_G*m)/r*sq(rdot) 
                         - 3/4*eta*(16 - 37*eta - 16*sq(eta))*sq(length(v)*rdot)));
    return r3PN;

}
Vector<mass, 3> c_SpinOrbit(mass m, mass r, Vector<mass, 3> const& Spin, Vector<mass, 3> const& sigma){

    Vector<mass, 3> rSO = SI_G/(SI_c2*cb(r))*(6*n*dot(cross(n,v),(Spin + sigma)) - cross(v,(4*Spin + 3*sigma)) 
                        + 3*rdot*cross(n,(2*Spin + sigma)));
    return rSO;
}
Vector<mass, 3> c_SpinSpin(mass m, mass r, Vector<mass, 3> const& Spin1, Vector<mass, 3> const& Spin2){

    Vector<mass, 3> rSS = -3*SI_G/(SI_c2*mu*std::pow(r,4))*(n*dot(Spin1,Spin2) + Spin1*dot(n,Spin2) + Spin2*dot(n,Spin2)
                        - 5*n*dot(n,Spin1)*dot(n,Spin2));
    return rSS;

}
Vector<mass, 3> c_BT_RR(mass m, mass r){

    Vector<mass, 3> rBTRR = 8/5*eta*sq(SI_G*m)/(std::pow(SI_c,5)*cb(r))*(rdot*n*(18*sq(length(v)) + 2/3*(SI_G*m)/r - 25*sq(rdot)) 
                          - v*(6*sq(length(v)) - 2(SI_G*m)/r - 15*sq(rdot)));
    return rBTRR;

}
Vector<mass, 3> c_PostNewtonianSO(mass m, mass r, Vector<mass, 3> const& Spin, Vector<mass, 3> const& Delta){

    Vector<mass, 3> rPNSO = SI_G/(sq(SI_c2)*cb(r))*(n*(dot(cross(n, v), Spin)*(-30*eta*sq(rdot) + 24*eta*sq(length(v)) 
                          - (SI_G*m)/r*(44 + 25*eta)) + dm/m*dot(cross(n, v), Delta)*(-15*eta*sq(rdot) + 12*eta*sq(length(v)) 
                          - (SI_G*m)/r*(24 + 29/2*eta))) + rdot*v*(dot(cross(n, v),Spin)*(-9 + 9*eta) + dm/m*dot(cross(n,v), Delta)*(-3 - 6*eta)) 
                          + cross(n, v)*(3/2*rdot*dot(v, Spin)*(-1 + eta) - 8*(SI_G*m)/r*eta*dot(n, Spin) 
                          - dm/m*(4*(SI_G*m)/r*eta*dot(n, Delta) + 3/2*rdot*dot(v, Delta))) 
                          + rodt*cross(n, Spin)*(-45/2*eta*sq(rdot) + 21*eta*sq(length(v)) - (SI_G*m)/r*(28 + 21*eta)) 
                          + dm/m*rdot*cross(n, Delta)*(-15*eta*sq(rdot) + 12*eta*sq(length(v)) - (SI_G*m)/r*(12 + 23/2*eta)) 
                          + cross(v,Spin)*(33/2*eta*sq(rdot) + (SI_G*m)/r*(24 + 11*eta) - 14*eta*sq(length(v))) 
                          + dm/m*cross(v, Delta)*(9*eta*sq(rdot) - 7*eta*sq(length(v)) + (SI_G*m)/r*(12 + 11/2*eta)));
    return rPNSO;

}
Vector<mass, 3> c_RR1PostNewtonian(mass m, mass r){

    Vector<mass, 3> rRR1PN = 8/5*eta*sq(SI_G*m)/(std::pow(SI_c, 7)*std::pow(r, 5))*(rdot*n*((87/14 - 48*eta)*std::pow(length(v),4) 
                           - (5379/28 - 136/3*eta)*sq(length(v))*(SI_G*m)/r + 25/2*(1 + 5*eta)*sq(length(v)*rdot) 
                           + (1353/4 - 133*eta)*sq(rdot)*(SI_G*m)/r - 35/2*(1 - eta)*std::pow(rdot,4) + (166/7 + 55/3*eta)*sq((SI_G*m)/r)) 
                           - v*(-27/14*std::pow(length(v), 4) - (4861/84 + 58/3*eta)*sq(length(v))*(SI_G*m)/r 
                           + 3/2*(13 - 37*eta)*sq(length(v)*rdot) + (2591/12 + 97*eta)*sq(rdot)*(SI_G*m)/r - 25/2*(1 - 7*eta)*std::pow(rdot, 4) 
                           + 1/3*(776/7 + 55*eta)*sq((SI_G*m)/r)));
    return rRR1PN;

}
Vector<mass, 3> c_RRSO(mass m, mass r, Vector<mass, 3> const& Spin, Vector<mass, 3> const& sigma, Vector<mass, 3> const& LN){

    Vector<mass, 3> rRRSO = - (sq(SI_G)*eta*m)/(5*std::pow(SI_c, 7)*std::pow(r, 4))*((rdot*n)/(mu*r)*((120*sq(length(v)) + 280*sq(rdot) 
                          + 453*(SI_G*m)/r)*dot(LN, Spin) + (120*sq(length(v)) + 280*sq(rdot) + 458*(SI_G*m)/r)*dot(LN, sigma)) 
                          + v/(mu*r)*((87*sq(length(v)) - 675*sq(rdot) -  901/3*(SI_G*m)/r)*dot(LN, Spin) 
                          + 4*(18*sq(length(v)) - 150*sq(rdot) - 66*(SI_G*m)/r)*dot(LN, sigma)) 
                          - 2/3*rdot*cross(v, Spin)*(48*sq(length(v)) + 15*sq(rdot) + 364*(SI_G*m)/r) 
                          + 1/3*rdot*cross(v, sigma)*(291*sq(length(v)) - 705*sq(rdot) - 772*(SI_G*m)/r) 
                          + 1/2*cross(n, Spin)*(31*std::pow(length(v), 4) - 260*sq(length(v)*rdot) + 245*std::pow(rdot, 4) 
                          + 537*sq(rdot)*(SI_G*m)/r + 4/3*sq((SI_G*m)/r)) + 1/2*cross(n,sigma)*(115*std::pow(length(v4)) 
                          - 1130*sq(length(v)*rdot) + 1295*std::pow(rdot, 4) - 869/3*sq(length(v))*(SI_G*m)/r +849*sq(rdot)*(SI_G*m)/r 
                          + 44/3*sq((SI_G*m)/r)));
    return rRRSO;

}
Vector<mass, 3> c_RRSS(mass m, mass r, Vector<mass, 3> const& Spin1, Vector<mass, 3> const& Spin2){

    Vector<mass, 3> rRRSS = sq(SI_G)/(std::pow(SI_c, 7)*std::pow(r, 5))*(n*((287*sq(rdot) - 99*sq(length(v)) 
                          + 541/5*(SI_G*m)/r)*rdot*dot(Spin1, Spin2) 
                          - (2646*sq(rdot) - 714*sq(length(v)) + 1961/5*(SI_G*m)/r)*rdot*dot(n, Spin1)*dot(n, Spin2) 
                          + (1029*sq(rdot) - 123*sq(length(v)) + 629/10*(SI_G*m)/r)*(dot(n, Spin1)*dot(n, Spin2) + dot(v, Spin2)*dot(v, Spin1))
                          - 336*rdot*dot(v, Spin1)*dot(v, Spin2)) + v*((171/3*sq(length(v)) - 195*sq(rdot) - 67*(SI_G*m)/r)*dot(Spin1, Spin2)
                          - (174*sq(length(v)) - 1386*sq(rdot) - 1038/5*(SI_G*m)/r)*dot(n, Spin1)*dot(n, Spin2) 
                          - 438*rdot*dot(n, Spin1)*dot(n, Spin2) + dot(v, Spin2)*dot(v, Spin1) + 96*rdot*dot(v, Spin1)*dot(v, Spin2))
                          + (27/10*sq(length(v)) - 75/2*sq(rdot) - 509/30*(SI_G*m)/r)*(dot(v,Spin2)*Spin1 + dot(v, Spin1)*Spin2) 
                          + (174*sq(length(v)) + 1386*sq(rdot) + 1038/5*(SI_G*m)/r)*rdot*(dot(n,Spin2)*Spin1 + dot(n, Spin1)*Spin2));
    return rRRSS;

}

mass e_Newtonian(mass m, mass r){

    return mass eN = mu*(1/2*sq(length(v)) - (SI_G*m)/r);

}
mass e_PostNewtonian(mass m, mass r){

    mass ePN = mu/SI_c2*(3/8*(1 - 3*eta)*std::pow(length(v), 4) + 1/2*(3 + eta)*sq(length(v))*(SI_G*m)/r + 1/2*eta*(SI_G*m)/r*sq(rdot) + 1/2*sq((SI_G*m)/r));
    return ePN;

}
mass e_2PostNewtonian(mass m, mass r){

    mass e2PN = mu/std::pow(SI_c, 4)*(5/16*(1 - 17*eta + 13*sq(eta))*std::pow(length(v),6) -3/8*eta*(1 - 3*eta)*(SI_G*m)/r*std::pow(rdot,4)
              + 1/8*(21 - 23*eta - 27*sq(eta))*(SI_G*m)/r*std::pow(length(v),4) + 1/8*(14 - 55*eta + 4*sq(eta))*sq((SI_G*m)/r*length(v))
              + 1/4*(1 - 15*eta)*(SI_G*m)/r*sq(length(v)*rdot) - 1/4*(2 + 15*eta)*cb((SI_G*m)/r) +1/8*(4 + 69*eta + 4*sq(eta))*sq((SI_G*m)/r*rdot));
    return e2PN;

}
mass e_3PostNewtonian(mass m, mass r){

    mass e3PN = mu/std::pow(SI_c, 6)*((3/8 + 18469/840*eta)std::pow((SI_G*m)/r, 4) + (5/4 - (6747/280 - 41/64*sq(PI))*eta 
              - 21/4*sq(eta) + 1/2*cb(eta))*cb((SI_G*m)/r)*sq(length(v))
              + (3/2 + (2321/280 - 123/64*sq(PI))*eta + 51/4*sq(eta) + 7/2*cb(eta))*cb((SI_G*m)/r)*sq(rdot )
              + 1/128*(35 - 413*eta + 1666*sq(eta) - 2261*cb(eta))*std::pow(length(v), 8) 
              + 1/16*(135 - 194*eta + 406*sq(eta) - 108*cb(eta))*sq((SI_G*m)/r*sq(length(v))) 
              + 1/16*(12 + 248*eta - 815*sq(eta) - 324*cb(eta))*sq((SI_G*m)/r*length(v)*rdot) 
              - 1/48*eta*(731 - 492*eta - 288*sq(eta))*sq((SI_G*m)/r*sq(rdot)) 
              + 1/16*(55 - 215*eta + 116*sq(eta) + 325*cb(eta))*(SI_G*m)/r*std::pow(length(v), 6)
              + 1/16*(5 - 25*eta +25*sq(eta))*(SI_G*m)/r*std::pow(rdot, 6) - 1/16*(21 + 75*eta - 375*sq(eta))*(SI_G*m)/r*sq(sq(length(v))*rdot)
              - 1/16*eta*(9 - 84*eta + 165*sq(eta))*(SI_G*m)/r*sq(length(v)*sq(rdot)));
    return e3PN;

}
mass e_SpinOrbit(mass r, Vector<mass, 3> const& sigma, Vector<mass, 3> const& LN){

    return mass eSO = SI_G/(SI_c2*cb(r))*dot(LN, sigma);

}
mass e_SpinSpin(mass r, Vector<mass, 3> const& Spin1, Vector<mass, 3> const& Spin2){

    return mass eSS = SI_G/(SI_c2*cb(r))*(3*dot(n, Spin1)*dot(n, Spin2) - dot(Spin1, Spin2));

}
mass e_PostNewtonianSO(mass m, mass r, Vector<mass, 3> const& Spin, Vector<mass, 3> const& Delta){

    mass ePNSO = (SI_G*mu)/(2*cb(SI_c)*sq(r))*cross(n, v)*(Delta*dm/m*((1 - 5*eta)*sq(length(v)) + 3*eta*(SI_G*m)/r) 
               - 3*Spin*((1 + eta)*sq(length(v)) + eta*sq(rdot) - 4/3*eta*(SI_G*m)/r));
    return ePNSO;

}

Vector<mass, 3> l_PostNewtonian(){

    return Vectro<mass, 3> lPN= LN/SI_c2*(1/2*sq(length(v))*(1 - 3*eta) + (3 + eta)*(SI_G*m)/r); 

}
Vector<mass, 3> l_2PostNewtonian(){

    Vector<mass, 3> l2PN = LN/std::pow(SI_c,4)*(3/8*(1 - 7*eta + 13*sq(eta))*std::pow(length(v), 4) - 1/2*eta*(2 + 5*eta)*(SI_G*m)/r*sq(rdot) 
                         + 1/2*(7 - 10*eta  - 9*sq(eta))*(SI_G*m)/r*sq(length(v)) + 1/4*(14 - 41*eta + 4*sq(eta))*sq((SI_G*m)/r));
    return l2PN;

}
Vector<mass, 3> l_3PostNewtonian(){

    Vector<mass, 3> l3PN = LN/std::pow(SI_c, 6)*((5/2 - (5199/280 - 41/32*sq(PI))*eta - 7*sq(eta) + cb(eta))*cb((SI_G*m)/r)
                         + 1/16*(5 - 59*eta + 238*sq(eta) - 323*cb(eta))*std::pow(length(v),6)
                         + 1/12*(135 - 322*eta + 315*sq(eta) - 108*cb(eta))*sq((SI_G*m)/r*length(v)) 
                         + 1/24*(12 - 287*eta - 951*sq(eta) - 324*cb(eta))*sq((SI_G*m)/r*rdot)
                         + 1/8*(33 - 142*eta + 106*sq(eta) + 195*cb(eta))*(SI_G*m)/r*std::pow(length(v), 4)
                         - 1/4*eta*(12 - 7*eta - 75*sq(eta))*(SI_G*m)/r*sq(length(v)*rdot)
                         + 3/8*eta*(2 - 2*eta - 11*sq(eta))*(SI_G*m)/r*std::pow(rdot, 4));
    return l3PN;

}
Vector<mass, 3> l_SpinOrbit(){

    return lSO = mu/(SI_c2*m)*((SI_G*m)/r*cross(n, cross(n, (2*Spin + sigma))) - 1/2*cross(v, cross(v, sigma)));

}
