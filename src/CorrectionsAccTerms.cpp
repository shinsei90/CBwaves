#include <Corrections.hpp>
#include <Vector.hpp>

#include <iostream>
#include <cmath>

mass m1 = 10;
mass m2 = 1.2;
mass dm = m1 - m2;
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
Vector<mass, 3> Spin = Spin1 + Spin2; 



auto sq = [](auto const& a){ return a*a;};
auto cb = [](auto const& a){ return a*a*a; };



Vector<mass, 3> c_Newtonian(mass m, mass r){

    Vector<mass, 3> rN = -(SI_G*m)/sq(r) * n;
    return rN;

}
Vector<mass, 3> c_PostNewtonian(mass m, mass r){

    Vector<mass,3> rPN = -(SI_G*m)/(SI_c2*cb(r))*(n*((1 + 3*eta)*sq(length(v)) - 2*(2 + eta)*(SI_G*m)/r - 3/2*eta*rdot) - 2*(2 - eta)*rdot*v);
    return rPN;

}
Vector<mass, 3> c_2PostNewtonian(mass m, mass r){

    Vector<mass, 3> r2PN = -(SI_G*m)/(sq(SI_c2*r)) * (n*(3/4*(12 + 29*eta)*sq((SI_G*m)/r) + eta*(3 - 4*eta)*std::pow(length(v),4) 
                         + 15/8*eta*(1 - 3*eta)*std::pow(rdot,4) - 3/2*eta*(3 - 4*eta)*sq(length(v)*rdot) 
                         - 1/2*eta*(13 - 4*eta)*(SI_G*m)/r*sq(length(v)) - (2 + 25*eta + 2*sq(eta))*(SI_G*m)/r*sq(rdot)) 
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

Vector<mass, 3> c_4PostNewtonian(mass m, mass r){

    mass a0 = (315/128*eta - 2205/128*sq(eta) + 2205/64*cb(eta) - 2205/128*std::pow(eta, 4))*std::pow(rdot, 8)
                       + (-175/16*eta + 595/8*sq(eta) - 2415/15*cb(eta) + 735/8*std::pow(eta, 4))*sq(cb(rdot)*length(v))
                       + (135/8*eta - 1875/16*sq(eta) + 4035/16*cb(eta) - 1335/8*std::pow(eta, 4))*std::pow(rdot*length(v))
                       + (-21/2*eta + 1191/16*sq(eta) - 327/2*cb(eta) + 99*std::pow(eta, 4))*sq(rdot*cb(length(v)))
                       + (21/8*eta - 175/8*sq(eta) + 61*cb(eta) - 54*std::pow(eta, 4))*std::pow(length(v), 8);

    mass a1 = m/r*((2973/40*eta + 407*sq(eta) + 181/2*cb(eta) - 86*std::pow(eta, 4))*std::pow(rdot, 6)
                       + (1497/32*eta - 1627/2*sq(eta) - 81*cb(eta) + 228*std::pow(eta,4))*sq(length(v)*sq(rdot))
                       - (2583/16*eta + 1009/2*sq(eta) + 47*cb(eta) - 104*std::pow(eta, 4))*sq(rdot*sqlength(v))
                       + (1067/32*eta - 58*sq(eta) - 44*cb(eta) + 58*std::pow(eta, 4))*std::pow(length(v), 6));

    mass a2 = sq(m/r)*(2094751/960*eta*std::pow(rdot,4) + 45255/1024*sq(PI)*eta*std::pow(rdot, 4)
                       + 326101/91*sq(eta)*std::pow(rdot, 4) - 4305/128*sq(PI*eta)*std::pow(rdot, 4)
                       - 1959/32*cb(eta)*std::pow(rdot, 4) - 126*std::pow(eta*rdot, 4) - 1636681/1120*eta*sq(rdot*length(v))
                       - 12585/512*eta*sq(PI*rdot*length(v)) - 255461/112*sq(eta*rdot*length(v)) + 3075/128*sq(PI*eta*rdot*length(v))
                       - 309/4*cb(eta)*sq(rdot*length(v)) + 63*std::pow(eta, 4)*sq(rdot*length(v)) 
                       + (1096941/11200*eta + 1155/1024*sq(PI)*eta + 7263*sq(eta) - 123/64*sq(PI*eta) + 145/2*cb(eta) 
                       - 16*std::pow(eta, 4))*std::pow(length(v),4));

    mass a3 = cb(m/r)*((-2 + (1297943/8400 - 2969/16*sq(PI))*eta + (1255151/840 + 7095/16*sq(PI))*sq(eta) 
                       - 17*cb(eta) - 24*std::pow(eta, 4))* sq(rdot) ((1237279/25200  + 3835/96*sq(PI))*eta 
                       - (693947/2520 + 229/8*sq(PI))*sq(eta) + 19/2*cb(eta))*sqlength(v));

    mass a4 = std::pow(m/r, 4)*(25 + (6625537/12600 - 4543/96*sq(PI))*eta + (477763/720 + 3/4*sq(PI))*sq(eta));
    mass A_4PN = a0 + a1 + a2 + a3 + a4;

    mass b0 = (105/16*eta - 245/8*sq(eta) + 385/16*cb(eta) + 35/8*std::pow(eta, 4))*std::pow(rdot, 7)
                       + (-165/8*eta + 1665/16*sq(eta) - 1725/16*cb(eta) - 105/4*std::pow(eta, 4))*std::pow(rdot, 5)*sqlength(v)
                       + (45/2*eta - 1869/16*sq(eta) + 129*cb(eta) + 54*std::pow(eta, 4))*cb(rdot)*std::pow(length(v),4)
                       + (-157/16*eta + 54*sq(eta) - 69*cb(eta) - 24*std::pow(eta, 4))*rdot*std::pow(length(v), 6);

    mass b1 = m/r*(-(54319/160*eta + 901/8*sq(eta) - 60*cb(eta) - 30*std::pow(eta, 4))*std::pow(rdot, 5)
                       + (25943/48*eta + 1199/12*sq(eta) - 349/2*cb(eta) - 98*std::pow(eta, 4))*cb(rdot)*sqlength(v)
                       - (5725/32*eta + 389/8*sq(eta) - 118*cb(eta)) - 44*std::pow(eta, 4)*rdot*std::pow(length(v), 4));

    mass b2 = sq(m/r)*((-(9130111/3306 + 4695/256*sq(PI))*eta - (184613/112 - 1845/64*sq(PI))*sq(eta) 
                       + 209/2*cb(eta) + 74*std::pow(eta, 4))*cb(rdot) ((8692601/5600 + 1455/256*sq(PI))*eta + (58557/70 - 123/8*sq(PI))*sq(eta)
                       - 70*cb(eta) - 34*std::pow(eta, 4))*rdot*sqlength(v));
    mass b3 = cb(m/r)*(2 - (619267/525 - 791/16*sq(PI))*eta - (28406/45 + 2201/32*sq(PI))*sq(eta) + 66*cb(eta) + 16*std::pow(eta, 4));
    mass B_4PN = b0 + b1 + b2 + b3;

    Vector<mass, 3> r4PN = -(SI_G*m)/sq(r)*((1 + A_4PN)*n + B_4PN*v);
    return r4PN;

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
                          - v*(6*sq(length(v)) - 2*(SI_G*m)/r - 15*sq(rdot)));
    return rBTRR;

}
Vector<mass, 3> c_PostNewtonianSO(mass m, mass r, Vector<mass, 3> const& Spin, Vector<mass, 3> const& Delta){

    Vector<mass, 3> rPNSO = SI_G/(sq(SI_c2)*cb(r))*(n*(dot(cross(n, v), Spin)*(-30*eta*sq(rdot) + 24*eta*sq(length(v)) 
                          - (SI_G*m)/r*(44 + 25*eta)) + dm/m*dot(cross(n, v), Delta)*(-15*eta*sq(rdot) + 12*eta*sq(length(v)) 
                          - (SI_G*m)/r*(24 + 29/2*eta))) + rdot*v*(dot(cross(n, v),Spin)*(-9 + 9*eta) + dm/m*dot(cross(n,v), Delta)*(-3 - 6*eta)) 
                          + cross(n, v)*(3/2*rdot*dot(v, Spin)*(-1 + eta) - 8*(SI_G*m)/r*eta*dot(n, Spin) 
                          - dm/m*(4*(SI_G*m)/r*eta*dot(n, Delta) + 3/2*rdot*dot(v, Delta))) 
                          + rdot*cross(n, Spin)*(-45/2*eta*sq(rdot) + 21*eta*sq(length(v)) - (SI_G*m)/r*(28 + 21*eta)) 
                          + dm/m*rdot*cross(n, Delta)*(-15*eta*sq(rdot) + 12*eta*sq(length(v)) - (SI_G*m)/r*(12 + 23/2*eta)) 
                          + cross(v,Spin)*(33/2*eta*sq(rdot) + (SI_G*m)/r*(24 + 11*eta) - 14*eta*sq(length(v))) 
                          + dm/m*cross(v, Delta)*(9*eta*sq(rdot) - 7*eta*sq(length(v)) + (SI_G*m)/r*(12 + 11/2*eta)));
    return rPNSO;

}

Vector<mass, 3> c_2PostNewtonian(mass m, mass r, Vector<mass, 3> const& Spin, Vector<mass, 3> const& Delta){

    mass a1 = dm/m*dot(n, cross(Delta, v))*n*((-105/4*eta + 315/4*sq(eta))*std::pow(rdot, 4) 
                       + (30*eta - 75*sq(eta))*rdot*sqlength(v) + (-9*eta + 24*sq(eta))*std::pow(length(v), 4))
                       + dot(n, cross(Spin, v))*n*((-105/2*eta + 315/2*sq(eta))*std::pow(rdot, 4) 
                       + (60*eta - 150*sq(eta))*sq(rdot)*sqlength(v) + (-18*eta + 48*sq(eta))*std::pow(length(v), 4))
                       + dm/m*dot(n, cross(Delta, v))*v*((-15/2*eta - 105/2*sq(eta))*cb(rdot) 
                       + (3/8 + 15/4*eta + 141/8*sq(eta))*rdot*sqlength(v))
                       + dot(n, cross(Spin, v))*v*((-15/2*eta - 195/4*sq(eta))*cb(rdot)
                       + (3/8 + 27/8*eta + 249/8*sq(eta))*rdot*sqlength(v))
                       + cross(n, Spin)*((315/8*eta - 945/8*sq(eta))*std::pow(rdot, 5)
                       + (-105/2*eta + 585/4*sq(eta))*cb(length(n))*std::pow(length(v), 5) 
                       + (-3/8 + 57/4*eta - 237/8*sq(eta))*length(n)*std::pow(length(v), 5))
                       + dm/m*cross(n, Delta)*((105/4*eta - 525/8*sq(eta))*std::pow(rdot, 5)
                       + (-75/2*eta + 345/4*sq(eta))*cb(length(n))*std::pow(length(v), 5)
                       + (-3/8 + 57/4*eta - 237/8*sq(eta))length(n)*std::pow(length(v), 5))
                       + cross(Spin, v)*((225/8*eta - 585/8*sq(eta))*std::pow(rdot, 4)
                       + (-3/8 - 255/8*eta + 627/8*sq(eta))*rdot*sqlength(v)
                       + (21/2*eta - 28*sq(eta))*std::pow(length(v), 4))
                       + dm/m*cross(Delta, v)*((15*eta - 315/8*sq(eta))*std::pow(rdot, 4)
                       + (-3/8 - 69/4*eta + 351/8*sq(eta))*rdot*sqlength(v)
                       + (-11/2*eta - 14*sq(eta))*std::pow(length(v),4));

    mass a2 = dm/m*dot(n, cross(Delta, v))*n*((3147/8*eta + 255/4*sq(eta))*sq(rdot)
                       + (-131/8*eta - 19*sq(eta))*sqlength(v)) 
                       + dot(n, cross(Spin, v))*n*((1635/2*eta + 117*sq(eta))*sq(rdot) 
                       + (-217/4*eta - 28*sq(eta))*sqlength(v)) 
                       + dot(n, cross(Delta, v))*v*dm/m*(-381/2*eta - 25*sq(eta))*(rdot)
                       + dot(n, cross(Spin, v))*v*(-777/2*eta - 87/2*sq(eta))*(rdot)
                       + cross(n, Spin)*((-1215/2*eta - 105*sq(eta))*cb(rdot) 
                       + (1067/4*eta + 79/2*sq(eta))*rdot*sqlength(v))
                       + dm/m*cross(n, Delta)*(-2193/8*eta + 279/4*sq(eta)*cb(rdot) 
                       + (945/8*eta + 23*sq(eta))*rdot*sqlength(v))
                       + cross(Spin, v)*((-352*eta - 123/2*sq(eta))*sq(rdot)
                       + (197/4*eta + 14*sq(eta))*sqlength(v))
                       + dm/m*cross(Delta, v)*((-1325/8*eta - 147/4*sq(eta))*sq(rdot)
                       + (177/8*eta + 7*sq(eta))*sqlength(v));

    mass a3 = dot(n, cross(Delta, v))*n*dm/m*(-111/2* - 441/4*eta + 5*sq(eta))
                       + dot(n, cross(Spin, v))*n*(-195/2 - 749/4*eta + 8*sq(eta))
                       + cross(n, Spin)*(121/2+ 65*eta - 8*sq(eta))*(rdot)
                       + dm/m*cross(n, Delta)*(57/2 + 85/4*eta - 6*sq(eta))*(rdot)
                       + cross(Spin, v)*(105/2 + 137/2*eta) + cross(Delta, v)*dm/m*(57/2 + 65/2*eta);
    
    Vector<mass, 3> r2PNSO = SI_G/(std::pow(SI_c,7)*cb(r))*(a1 + (SI_G*m)/r*a2 +sq((SI_G*m)/r)*a3);
    return r2PNSO;

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
                          + 537*sq(rdot)*(SI_G*m)/r + 4/3*sq((SI_G*m)/r)) + 1/2*cross(n,sigma)*(115*std::pow(length(v),4) 
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

    mass eN = mu*(1/2*sq(length(v)) - (SI_G*m)/r);
    return eN;

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

    mass e3PN = mu/std::pow(SI_c, 6)*((3/8 + 18469/840*eta)*std::pow((SI_G*m)/r, 4) + (5/4 - (6747/280 - 41/64*sq(PI))*eta 
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

    mass eSO = SI_G/(SI_c2*cb(r))*dot(LN, sigma);
    return eSO;

}
mass e_SpinSpin(mass r, Vector<mass, 3> const& Spin1, Vector<mass, 3> const& Spin2){

    mass eSS = SI_G/(SI_c2*cb(r))*(3*dot(n, Spin1)*dot(n, Spin2) - dot(Spin1, Spin2));
    return eSS;

}
mass e_PostNewtonianSO(mass m, mass r, Vector<mass, 3> const& Spin, Vector<mass, 3> const& Delta){

    mass ePNSO = (SI_G*mu)/(2*cb(SI_c)*sq(r))*cross(n, v)*(Delta*dm/m*((1 - 5*eta)*sq(length(v)) + 3*eta*(SI_G*m)/r) 
               - 3*Spin*((1 + eta)*sq(length(v)) + eta*sq(rdot) - 4/3*eta*(SI_G*m)/r));
    return ePNSO;

}

Vector<mass, 3> l_PostNewtonian(mass m, mass r, Vector<mass, 3> const& LN){

    Vector<mass, 3> lPN= LN/SI_c2*(1/2*sq(length(v))*(1 - 3*eta) + (3 + eta)*(SI_G*m)/r);
    return  lPN;

}
Vector<mass, 3> l_2PostNewtonian(mass m, mass r, Vector<mass, 3> const& LN){

    Vector<mass, 3> l2PN = LN/std::pow(SI_c,4)*(3/8*(1 - 7*eta + 13*sq(eta))*std::pow(length(v), 4) - 1/2*eta*(2 + 5*eta)*(SI_G*m)/r*sq(rdot) 
                         + 1/2*(7 - 10*eta  - 9*sq(eta))*(SI_G*m)/r*sq(length(v)) + 1/4*(14 - 41*eta + 4*sq(eta))*sq((SI_G*m)/r));
    return l2PN;

}
Vector<mass, 3> l_3PostNewtonian(mass m, mass r, Vector<mass, 3> const& LN){

    Vector<mass, 3> l3PN = LN/std::pow(SI_c, 6)*((5/2 - (5199/280 - 41/32*sq(PI))*eta - 7*sq(eta) + cb(eta))*cb((SI_G*m)/r)
                         + 1/16*(5 - 59*eta + 238*sq(eta) - 323*cb(eta))*std::pow(length(v),6)
                         + 1/12*(135 - 322*eta + 315*sq(eta) - 108*cb(eta))*sq((SI_G*m)/r*length(v)) 
                         + 1/24*(12 - 287*eta - 951*sq(eta) - 324*cb(eta))*sq((SI_G*m)/r*rdot)
                         + 1/8*(33 - 142*eta + 106*sq(eta) + 195*cb(eta))*(SI_G*m)/r*std::pow(length(v), 4)
                         - 1/4*eta*(12 - 7*eta - 75*sq(eta))*(SI_G*m)/r*sq(length(v)*rdot)
                         + 3/8*eta*(2 - 2*eta - 11*sq(eta))*(SI_G*m)/r*std::pow(rdot, 4));
    return l3PN;

}
Vector<mass, 3> l_SpinOrbit(mass m, mass r, Vector<mass, 3> const& Spin, Vector<mass, 3> const& sigma){

    Vector<mass, 3> lSO = mu/(SI_c2*m)*((SI_G*m)/r*cross(n, cross(n, (2*Spin + sigma))) - 1/2*cross(v, cross(v, sigma)));
    return lSO;

}


Vector<mass, 3> spin1(mass m, mass r){

    mass alpha1PN;
    mass alpha2PN;
    mass alpha3PN;
    
    alpha1PN = (SI_G*m)/sq(r)*(3/4 + 1/2*eta - 3/4*dm/m);
    
    alpha2PN = (SI_G*m)/sq(r)*((-3/2*eta + 3/4*sq(eta) - 3/2*eta*dm/m)*sq(rdot) 
             + (1/16 + 1/18*eta - 3/8*sq(eta) + dm/m*(-1/16 + 1/2*eta))*sqlength(v)))
             + sq(SI_G*m)/cb(r)*(-1/4 - 3/8*eta + 1/2*sq(eta) + dm/m*(1/4 - 1/8*eta));
    
    alpha2PN = (SI_G*m)/sq(r)*((15/8*eta - 195/32*sq(eta) + 15/16*cb(eta) + dm/m*(15/8*eta - 75/32*sq(eta)))*std::pow(length(n)length(v), 4)
             + (-3*eta + 291/32*sq(eta) - 45/16*cb(eta) + dm/m*(-3*eta + 177/32*sq(eta)))*sq(rdot)*sqlength(v)
             + (1/32 + 19/16*eta - 31/8*sq(eta) + 17/16*cb(eta) + dm/m*(-1/32 + 3/4*eta - 11/8*sq(eta)))*std::pow(length(v),4))
             + sq(SI_G*m)/cb(r)*((1/4 - 525/32*eta - 159/16*sq(eta) + 13/4*cb(eta) + dm/m*(-1/4 - 75/32*eta - 87/16*sq(eta)))*sq(rdot
             + (3/16 + 27/4*eta + 75/32*sq(eta) - 9/8*cb(eta) + dm/m*(-3/16 + 9/8*eta + 35/32*sq(eta)))*sqlength(v))
             + cb(SI_G*m)/std::pow(r,4)*(7/16 - 9/4*eta - 9/8*sq(eta)*1/2*cb(eta) + dm/m*(-7/16 - 1/8*eta - 1/8*sq(eta)));

    Vector<mass, 3> Spin1 = cross(n,v)*(1/SI_c2*alpha1PN + 1/sq(SI_c2)*alpha2PN + 1/cb(SI_c2)*alpha3PN);

    return Spin1;

}

Vector<mass, 3> spin2(mass m, mass r){

    mass alpha1PN;
    mass alpha2PN;
    mass alpha3PN;
    
    alpha1PN = (SI_G*m)/sq(r)*(3/4 + 1/2*eta + 3/4*dm/m);
    
    alpha2PN = (SI_G*m)/sq(r)*((-3/2*eta + 3/4*sq(eta) + 3/2*eta*dm/m)*sq(rdot) 
             + (1/16 + 1/18*eta - 3/8*sq(eta) - dm/m*(-1/16 + 1/2*eta))*sqlength(v)))
             + sq(SI_G*m)/cb(r)*(-1/4 - 3/8*eta + 1/2*sq(eta) - dm/m*(1/4 - 1/8*eta));
    
    alpha2PN = (SI_G*m)/sq(r)*((15/8*eta - 195/32*sq(eta) + 15/16*cb(eta) - dm/m*(15/8*eta - 75/32*sq(eta)))*std::pow(rdot, 4)
             + (-3*eta + 291/32*sq(eta) - 45/16*cb(eta) - dm/m*(-3*eta + 177/32*sq(eta)))*sq(rdot)*sqlength(v)
             + (1/32 + 19/16*eta - 31/8*sq(eta) + 17/16*cb(eta) - dm/m*(-1/32 + 3/4*eta - 11/8*sq(eta)))*std::pow(length(v),4))
             + sq(SI_G*m)/cb(r)*((1/4 - 525/32*eta - 159/16*sq(eta) + 13/4*cb(eta) - dm/m*(-1/4 - 75/32*eta - 87/16*sq(eta)))*sq(rdot
             + (3/16 + 27/4*eta + 75/32*sq(eta) - 9/8*cb(eta) - dm/m*(-3/16 + 9/8*eta + 35/32*sq(eta)))*sqlength(v))
             + cb(SI_G*m)/std::pow(r,4)*(7/16 - 9/4*eta - 9/8*sq(eta)*1/2*cb(eta) - dm/m*(-7/16 - 1/8*eta - 1/8*sq(eta)));

    Vector<mass, 3> Spin2 = cross(n,v)*(1/SI_c2*alpha1PN + 1/sq(SI_c2)*alpha2PN + 1/cb(SI_c2)*alpha3PN);

    return Spin2;

}