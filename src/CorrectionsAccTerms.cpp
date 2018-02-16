#include <Config.hpp>
#include <Corrections.hpp>
#include <DynamicalParams.hpp>
#include <InitParams.hpp>
#include <Vector.hpp>

#include <iostream>
#include <cmath>

auto sq = [](auto const& a){ return a*a;};
auto cb = [](auto const& a){ return a*a*a; };



Vector<mass, 3> c_Newtonian(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& r = dp.r;
    Vector<mass, 3> const& n = dp.n;

    Vector<mass, 3> rN = -(SI_G*m)/sq(r) * n;
    return rN;

}

Vector<mass, 3> c_PostNewtonian(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& eta = ip.eta;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& n = dp.n;
    Vector<mass, 3> const& v = dp.v;

    Vector<mass,3> rPN = -(SI_G*m)/(SI_c2*cb(r))*(n*((1. + 3.*eta)*sq(length(v)) - 2.*(2. + eta)*(SI_G*m)/r - 3./2.*eta*rdot) - 2.*(2. - eta)*rdot*v);
    return rPN;

}

Vector<mass, 3> c_2PostNewtonian(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& eta = ip.eta;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& n = dp.n;
    Vector<mass, 3> const& v = dp.v;

    Vector<mass, 3> r2PN = -(SI_G*m)/(sq(SI_c2*r)) * (n*(3./4.*(12. + 29.*eta)*sq((SI_G*m)/r) + eta*(3. - 4.*eta)*std::pow(length(v),4) 
                         + 15./8.*eta*(1. - 3.*eta)*std::pow(rdot,4) - 3./2.*eta*(3. - 4.*eta)*sq(length(v)*rdot) 
                         - 1./2.*eta*(13. - 4.*eta)*(SI_G*m)/r*sq(length(v)) - (2. + 25.*eta + 2.*sq(eta))*(SI_G*m)/r*sq(rdot)) 
                         - 1./2.*rdot*v*(eta*(15. + 4.*eta)*sq(length(v)) - (4. + 41.*eta + 8.*sq(eta))*(SI_G*m)/r -3.*eta*(3. + 2.*eta)*sq(rdot)));
    return r2PN;

}

Vector<mass, 3> c_3PostNewtonian(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& eta = ip.eta;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& n = dp.n;
    Vector<mass, 3> const& v = dp.v;

    Vector<mass, 3> r3PN = (SI_G*m)/(cb(SI_c)*sq(r))*(n*((16. + (1399./12. - 41./16.*sq(PI))*eta + 71./2.*sq(eta))*cb((SI_G*m)/r) 
                         + eta*(20827./840. + 123./64.*sq(PI) - sq(eta))*sq((SI_G*m)/r*length(v)) 
                         - (1. + (22717./168. + 615./64.*sq(PI))*eta + 11./8.*sq(eta) - 7.*cb(eta))*sq((SI_G*m)/r*rdot) 
                         - 1./4.*eta*(11. - 49.*eta + 52.*sq(eta))*std::pow(length(v),6) + 35./16.*eta*(1. - 5.*eta + 5.*sq(eta))*std::pow(rdot,6)
                         - 1./4.*eta*(75. + 32.*eta - 40.*sq(eta))*(SI_G*m)/r*std::pow(length(v),4)
                         - 1./2.*eta*(158. - 69.*eta - 60.*sq(eta))*(SI_G*m)/r*std::pow(rdot,4) 
                         + eta*(121. - 16.*eta - 20.*sq(eta))*(SI_G*m)/r*sq(length(v)*rdot) 
                         + 3./8.*eta*(20. - 79.*eta + 60.*sq(eta))*sq(sq(length(v))*rdot) 
                         - 15./8.*eta*(4. - 18.*eta + 17.*sq(eta))*sq(length(v)*sq(rdot))) 
                         + rdot*v*((4. + (5849./840. + 123./32.*sq(PI))*eta - 25.*sq(eta) - 8.*cb(eta))*sq((SI_G*m)/r) 
                         + 1./8.*eta*(65. - 152.*eta - 48.*sq(eta))*std::pow(length(v),4) + 15./8.*eta*(3. - 8.*eta - 2.*sq(eta))*std::pow(rdot,4)
                         + eta*(15. + 27.*eta + 10.*sq(eta))*(SI_G*m)/r*sq(length(v)) 
                         - 11./6.*eta*(329. + 177.*eta + 108.*sq(eta))*(SI_G*m)/r*sq(rdot) 
                         - 3./4.*eta*(16. - 37.*eta - 16.*sq(eta))*sq(length(v)*rdot)));
    return r3PN;

}

Vector<mass, 3> c_4PostNewtonian(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& eta = ip.eta;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& n = dp.n;
    Vector<mass, 3> const& v = dp.v;

    mass a0 = (315./128.*eta - 2205./128.*sq(eta) + 2205./64.*cb(eta) - 2205./128.*std::pow(eta, 4))*std::pow(rdot, 8)
            + (-175./16.*eta + 595./8.*sq(eta) - 2415./15.*cb(eta) + 735./8.*std::pow(eta, 4))*sq(cb(rdot)*length(v))
            + (135./8.*eta - 1875./16.*sq(eta) + 4035./16.*cb(eta) - 1335./8.*std::pow(eta, 4))*std::pow(rdot*length(v), 4)
            + (-21./2.*eta + 1191./16.*sq(eta) - 327./2.*cb(eta) + 99.*std::pow(eta, 4))*sq(rdot*cb(length(v)))
            + (21./8.*eta - 175./8.*sq(eta) + 61.*cb(eta) - 54.*std::pow(eta, 4))*std::pow(length(v), 8);

    mass a1 = SI_G*m/r*((2973./40.*eta + 407.*sq(eta) + 181./2.*cb(eta) - 86.*std::pow(eta, 4))*std::pow(rdot, 6)
            + (1497./32.*eta - 1627./2.*sq(eta) - 81.*cb(eta) + 228.*std::pow(eta,4))*sq(length(v)*sq(rdot))
            - (2583./16.*eta + 1009./2.*sq(eta) + 47.*cb(eta) - 104.*std::pow(eta, 4))*sq(rdot*sqlength(v))
            + (1067./32.*eta - 58.*sq(eta) - 44.*cb(eta) + 58.*std::pow(eta, 4))*std::pow(length(v), 6));

    mass a2 = sq(SI_G*m/r)*(2094751./960.*eta*std::pow(rdot,4) + 45255./1024.*sq(PI)*eta*std::pow(rdot, 4)
            + 326101./91.*sq(eta)*std::pow(rdot, 4) - 4305./128.*sq(PI*eta)*std::pow(rdot, 4)
            - 1959./32.*cb(eta)*std::pow(rdot, 4) - 126.*std::pow(eta*rdot, 4) - 1636681./1120.*eta*sq(rdot*length(v))
            - 12585./512.*eta*sq(PI*rdot*length(v)) - 255461./112.*sq(eta*rdot*length(v)) + 3075./128.*sq(PI*eta*rdot*length(v))
            - 309./4.*cb(eta)*sq(rdot*length(v)) + 63.*std::pow(eta, 4)*sq(rdot*length(v)) 
            + (1096941./11200.*eta + 1155./1024.*sq(PI)*eta + 7263.*sq(eta) - 123./64.*sq(PI*eta) + 145./2.*cb(eta) 
            - 16.*std::pow(eta, 4))*std::pow(length(v),4));

    mass a3 = cb(SI_G*m/r)*((-2. + (1297943./8400. - 2969./16.*sq(PI))*eta + (1255151./840. + 7095./16.*sq(PI))*sq(eta) 
            - 17.*cb(eta) - 24.*std::pow(eta, 4))*sq(rdot) + ((1237279./25200.  + 3835./96.*sq(PI))*eta 
            - (693947./2520. + 229./8.*sq(PI))*sq(eta) + 19./2.*cb(eta))*sqlength(v));

    mass a4 = std::pow(SI_G*m/r, 4)*(25. + (6625537./12600. - 4543./96.*sq(PI))*eta + (477763./720. + 3./4.*sq(PI))*sq(eta));
    mass A_4PN = a0 + a1 + a2 + a3 + a4;

    mass b0 = (105./16.*eta - 245./8.*sq(eta) + 385./16.*cb(eta) + 35./8.*std::pow(eta, 4))*std::pow(rdot, 7)
            + (-165./8.*eta + 1665./16.*sq(eta) - 1725./16.*cb(eta) - 105./4.*std::pow(eta, 4))*std::pow(rdot, 5)*sqlength(v)
            + (45./2.*eta - 1869./16.*sq(eta) + 129.*cb(eta) + 54.*std::pow(eta, 4))*cb(rdot)*std::pow(length(v),4)
            + (-157./16.*eta + 54.*sq(eta) - 69.*cb(eta) - 24.*std::pow(eta, 4))*rdot*std::pow(length(v), 6);

    mass b1 = SI_G*m/r*(-(54319./160.*eta + 901./8.*sq(eta) - 60.*cb(eta) - 30.*std::pow(eta, 4))*std::pow(rdot, 5)
            + (25943./48.*eta + 1199./12.*sq(eta) - 349./2.*cb(eta) - 98.*std::pow(eta, 4))*cb(rdot)*sqlength(v)
            - (5725./32.*eta + 389./8.*sq(eta) - 118.*cb(eta)) - 44.*std::pow(eta, 4)*rdot*std::pow(length(v), 4));

    mass b2 = sq(SI_G*m/r)*((-(9130111./3306. + 4695./256.*sq(PI))*eta - (184613./112. - 1845./64.*sq(PI))*sq(eta) 
            + 209./2.*cb(eta) + 74.*std::pow(eta, 4))*cb(rdot) + ((8692601./5600. + 1455./256.*sq(PI))*eta + (58557./70. - 123./8.*sq(PI))*sq(eta)
            - 70.*cb(eta) - 34.*std::pow(eta, 4))*rdot*sqlength(v));
    
    mass b3 = cb(SI_G*m/r)*(2. - (619267./525. - 791./16.*sq(PI))*eta - (28406./45. + 2201./32.*sq(PI))*sq(eta) + 66.*cb(eta) + 16.*std::pow(eta, 4));
    mass B_4PN = b0 + b1 + b2 + b3;

    Vector<mass, 3> r4PN = -(SI_G*m)/(std::pow(SI_c, 8)*sq(r))*((1. + A_4PN)*n + B_4PN*v);
    return r4PN;

}

Vector<mass, 3> c_SpinOrbit(dynamicalParams const& dp, initParams const& ip){

    mass const& eta = ip.eta;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& n = dp.n;
    Vector<mass, 3> const& v = dp.v;
    Vector<mass, 3> const& Spin = dp.Spin;
    Vector<mass, 3> const& sigma = dp.sigma;

    Vector<mass, 3> rSO = SI_G/(SI_c2*cb(r))*(6.*n*dot(cross(n,v),(Spin + sigma)) - cross(v,(4.*Spin + 3.*sigma)) 
                        + 3.*rdot*cross(n,(2.*Spin + sigma)));
    return rSO;
}

Vector<mass, 3> c_SpinSpin(dynamicalParams const& dp, initParams const& ip){

    mass const& mu = ip.mu;
    mass const& eta = ip.eta;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& n = dp.n;
    Vector<mass, 3> const& v = dp.v;
    Vector<mass, 3> const& Spin1 = dp.Spin1;
    Vector<mass, 3> const& Spin2 = dp.Spin2;

    Vector<mass, 3> rSS = -3.*SI_G/(SI_c2*mu*std::pow(r,4))*(n*dot(Spin1,Spin2) + Spin1*dot(n,Spin2) + Spin2*dot(n,Spin2)
                        - 5.*n*dot(n,Spin1)*dot(n,Spin2));
    return rSS;

}

Vector<mass, 3> c_BT_RR(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& eta = ip.eta;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& n = dp.n;
    Vector<mass, 3> const& v = dp.v;

    Vector<mass, 3> rBTRR = 8./5.*eta*sq(SI_G*m)/(std::pow(SI_c,5)*cb(r))*(rdot*n*(18.*sq(length(v)) + 2./3.*(SI_G*m)/r - 25.*sq(rdot)) 
                          - v*(6.*sq(length(v)) - 2.*(SI_G*m)/r - 15.*sq(rdot)));
    return rBTRR;

}

Vector<mass, 3> c_PostNewtonianSO(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& eta = ip.eta;
    mass const& dm = ip.dm;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& n = dp.n;
    Vector<mass, 3> const& v = dp.v;
    Vector<mass, 3> const& Spin = dp.Spin;
    Vector<mass, 3> const& Delta = dp.Delta;

    Vector<mass, 3> rPNSO = SI_G/(sq(SI_c2)*cb(r))*(n*(dot(cross(n, v), Spin)*(-30.*eta*sq(rdot) + 24.*eta*sq(length(v)) 
                          - (SI_G*m)/r*(44. + 25.*eta)) + dm/m*dot(cross(n, v), Delta)*(-15.*eta*sq(rdot) + 12.*eta*sq(length(v)) 
                          - (SI_G*m)/r*(24. + 29./2.*eta))) + rdot*v*(dot(cross(n, v),Spin)*(-9. + 9.*eta) + dm/m*dot(cross(n,v), Delta)*(-3. - 6.*eta)) 
                          + cross(n, v)*(3./2.*rdot*dot(v, Spin)*(-1. + eta) - 8.*(SI_G*m)/r*eta*dot(n, Spin) 
                          - dm/m*(4.*(SI_G*m)/r*eta*dot(n, Delta) + 3./2.*rdot*dot(v, Delta))) 
                          + rdot*cross(n, Spin)*(-45./2.*eta*sq(rdot) + 21.*eta*sq(length(v)) - (SI_G*m)/r*(28. + 21.*eta)) 
                          + dm/m*rdot*cross(n, Delta)*(-15.*eta*sq(rdot) + 12.*eta*sq(length(v)) - (SI_G*m)/r*(12. + 23./2.*eta)) 
                          + cross(v,Spin)*(33./2.*eta*sq(rdot) + (SI_G*m)/r*(24. + 11.*eta) - 14.*eta*sq(length(v))) 
                          + dm/m*cross(v, Delta)*(9.*eta*sq(rdot) - 7.*eta*sq(length(v)) + (SI_G*m)/r*(12. + 11./2.*eta)));
    return rPNSO;

}

Vector<mass, 3> c_2PostNewtonianSO(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& eta = ip.eta;
    mass const& dm = ip.dm;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& n = dp.n;
    Vector<mass, 3> const& v = dp.v;
    Vector<mass, 3> const& Spin = dp.Spin;
    Vector<mass, 3> const& Delta = dp.Delta;

    Vector<mass, 3> a1 = dm/m*dot(n, cross(Delta, v))*n*((-105./4.*eta + 315./4.*sq(eta))*std::pow(rdot, 4) 
                       + (30.*eta - 75.*sq(eta))*rdot*sqlength(v) + (-9.*eta + 24.*sq(eta))*std::pow(length(v), 4))
                       + dot(n, cross(Spin, v))*n*((-105./2.*eta + 315./2.*sq(eta))*std::pow(rdot, 4) 
                       + (60.*eta - 150.*sq(eta))*sq(rdot)*sqlength(v) + (-18.*eta + 48.*sq(eta))*std::pow(length(v), 4))
                       + dm/m*dot(n, cross(Delta, v))*v*((-15./2.*eta - 105./2.*sq(eta))*cb(rdot) 
                       + (3./8. + 15./4.*eta + 141./8.*sq(eta))*rdot*sqlength(v))
                       + dot(n, cross(Spin, v))*v*((-15./2.*eta - 195./4.*sq(eta))*cb(rdot)
                       + (3./8. + 27./8.*eta + 249./8.*sq(eta))*rdot*sqlength(v))
                       + cross(n, Spin)*((315./8.*eta - 945./8.*sq(eta))*std::pow(rdot, 5)
                       + (-105./2.*eta + 585./4.*sq(eta))*cb(length(n))*std::pow(length(v), 5) 
                       + (-3./8. + 57./4.*eta - 237./8.*sq(eta))*length(n)*std::pow(length(v), 5))
                       + dm/m*cross(n, Delta)*((105./4.*eta - 525./8.*sq(eta))*std::pow(rdot, 5)
                       + (-75./2.*eta + 345./4.*sq(eta))*cb(length(n))*std::pow(length(v), 5)
                       + (-3/8. + 57./4.*eta - 237./8.*sq(eta))*length(n)*std::pow(length(v), 5))
                       + cross(Spin, v)*((225./8.*eta - 585./8.*sq(eta))*std::pow(rdot, 4)
                       + (-3./8. - 255./8.*eta + 627./8.*sq(eta))*rdot*sqlength(v)
                       + (21./2.*eta - 28.*sq(eta))*std::pow(length(v), 4))
                       + dm/m*cross(Delta, v)*((15.*eta - 315./8.*sq(eta))*std::pow(rdot, 4)
                       + (-3./8. - 69./4.*eta + 351./8.*sq(eta))*rdot*sqlength(v)
                       + (-11./2.*eta - 14.*sq(eta))*std::pow(length(v),4));

    Vector<mass, 3> a2 = dm/m*dot(n, cross(Delta, v))*n*((3147./8.*eta + 255./4.*sq(eta))*sq(rdot)
                       + (-131./8.*eta - 19.*sq(eta))*sqlength(v)) 
                       + dot(n, cross(Spin, v))*n*((1635./2.*eta + 117.*sq(eta))*sq(rdot) 
                       + (-217./4.*eta - 28.*sq(eta))*sqlength(v)) 
                       + dot(n, cross(Delta, v))*v*dm/m*(-381./2.*eta - 25.*sq(eta))*(rdot)
                       + dot(n, cross(Spin, v))*v*(-777./2.*eta - 87./2.*sq(eta))*(rdot)
                       + cross(n, Spin)*((-1215./2.*eta - 105.*sq(eta))*cb(rdot) 
                       + (1067./4.*eta + 79./2.*sq(eta))*rdot*sqlength(v))
                       + dm/m*cross(n, Delta)*(-2193./8.*eta + 279./4.*sq(eta)*cb(rdot) 
                       + (945./8.*eta + 23.*sq(eta))*rdot*sqlength(v))
                       + cross(Spin, v)*((-352.*eta - 123./2.*sq(eta))*sq(rdot)
                       + (197./4.*eta + 14.*sq(eta))*sqlength(v))
                       + dm/m*cross(Delta, v)*((-1325./8.*eta - 147./4.*sq(eta))*sq(rdot)
                       + (177./8.*eta + 7.*sq(eta))*sqlength(v));

    Vector<mass, 3> a3 = dot(n, cross(Delta, v))*n*dm/m*(-111./2.* - 441./4.*eta + 5.*sq(eta))
                       + dot(n, cross(Spin, v))*n*(-195./2. - 749./4.*eta + 8.*sq(eta))
                       + cross(n, Spin)*(121./2. + 65.*eta - 8.*sq(eta))*(rdot)
                       + dm/m*cross(n, Delta)*(57./2. + 85./4.*eta - 6.*sq(eta))*(rdot)
                       + cross(Spin, v)*(105./2. + 137./2.*eta) + cross(Delta, v)*dm/m*(57./2. + 65./2.*eta);
    
    Vector<mass, 3> r2PNSO = SI_G/(std::pow(SI_c,7)*cb(r))*(a1 + (SI_G*m)/r*a2 +sq((SI_G*m)/r)*a3);
    return r2PNSO;

}

Vector<mass, 3> c_RR1PostNewtonian(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& eta = ip.eta;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& n = dp.n;
    Vector<mass, 3> const& v = dp.v;

    Vector<mass, 3> rRR1PN = 8./5.*eta*sq(SI_G*m)/(std::pow(SI_c, 7)*std::pow(r, 5))*(rdot*n*((87./14. - 48.*eta)*std::pow(length(v),4) 
                           - (5379./28. - 136./3.*eta)*sq(length(v))*(SI_G*m)/r + 25./2.*(1. + 5.*eta)*sq(length(v)*rdot) 
                           + (1353./4. - 133.*eta)*sq(rdot)*(SI_G*m)/r - 35./2.*(1. - eta)*std::pow(rdot,4) + (166./7. + 55./3.*eta)*sq((SI_G*m)/r)) 
                           - v*(-27./14.*std::pow(length(v), 4) - (4861./84. + 58./3.*eta)*sq(length(v))*(SI_G*m)/r 
                           + 3/2.*(13. - 37.*eta)*sq(length(v)*rdot) + (2591./12. + 97.*eta)*sq(rdot)*(SI_G*m)/r - 25./2.*(1. - 7.*eta)*std::pow(rdot, 4) 
                           + 1./3.*(776./7. + 55.*eta)*sq((SI_G*m)/r)));
    return rRR1PN;

}

Vector<mass, 3> c_RRSO(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& mu = ip.mu;
    mass const& eta = ip.eta;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& n = dp.n;
    Vector<mass, 3> const& v = dp.v;
    Vector<mass, 3> const& Spin = dp.Spin;
    Vector<mass, 3> const& sigma = dp.sigma;
    Vector<mass, 3> const& LN = dp.LN;

    Vector<mass, 3> rRRSO = - (sq(SI_G)*eta*m)/(5.*std::pow(SI_c, 7)*std::pow(r, 4))*((rdot*n)/(mu*r)*((120.*sq(length(v)) + 280.*sq(rdot) 
                          + 453.*(SI_G*m)/r)*dot(LN, Spin) + (120.*sq(length(v)) + 280.*sq(rdot) + 458.*(SI_G*m)/r)*dot(LN, sigma)) 
                          + v/(mu*r)*((87.*sq(length(v)) - 675.*sq(rdot) -  901./3.*(SI_G*m)/r)*dot(LN, Spin) 
                          + 4.*(18.*sq(length(v)) - 150.*sq(rdot) - 66.*(SI_G*m)/r)*dot(LN, sigma)) 
                          - 2./3.*rdot*cross(v, Spin)*(48.*sq(length(v)) + 15.*sq(rdot) + 364.*(SI_G*m)/r) 
                          + 1./3.*rdot*cross(v, sigma)*(291.*sq(length(v)) - 705.*sq(rdot) - 772.*(SI_G*m)/r) 
                          + 1./2.*cross(n, Spin)*(31.*std::pow(length(v), 4) - 260.*sq(length(v)*rdot) + 245.*std::pow(rdot, 4) 
                          + 537.*sq(rdot)*(SI_G*m)/r + 4./3.*sq((SI_G*m)/r)) + 1./2.*cross(n,sigma)*(115.*std::pow(length(v),4) 
                          - 1130.*sq(length(v)*rdot) + 1295.*std::pow(rdot, 4.) - 869./3.*sq(length(v))*(SI_G*m)/r +849.*sq(rdot)*(SI_G*m)/r 
                          + 44./3.*sq((SI_G*m)/r)));
    return rRRSO;

}

Vector<mass, 3> c_RRSS(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& eta = ip.eta;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& n = dp.n;
    Vector<mass, 3> const& v = dp.v;
    Vector<mass, 3> const& Spin1 = dp.Spin1;
    Vector<mass, 3> const& Spin2 = dp.Spin2;

    Vector<mass, 3> rRRSS = sq(SI_G)/(std::pow(SI_c, 7)*std::pow(r, 5))*(n*((287.*sq(rdot) - 99.*sq(length(v)) 
                          + 541./5.*(SI_G*m)/r)*rdot*dot(Spin1, Spin2) 
                          - (2646.*sq(rdot) - 714.*sq(length(v)) + 1961./5.*(SI_G*m)/r)*rdot*dot(n, Spin1)*dot(n, Spin2) 
                          + (1029.*sq(rdot) - 123.*sq(length(v)) + 629./10.*(SI_G*m)/r)*(dot(n, Spin1)*dot(n, Spin2) + dot(v, Spin2)*dot(v, Spin1))
                          - 336.*rdot*dot(v, Spin1)*dot(v, Spin2)) + v*((171./3.*sq(length(v)) - 195.*sq(rdot) - 67.*(SI_G*m)/r)*dot(Spin1, Spin2)
                          - (174.*sq(length(v)) - 1386.*sq(rdot) - 1038./5.*(SI_G*m)/r)*dot(n, Spin1)*dot(n, Spin2) 
                          - 438.*rdot*dot(n, Spin1)*dot(n, Spin2) + dot(v, Spin2)*dot(v, Spin1) + 96.*rdot*dot(v, Spin1)*dot(v, Spin2))
                          + (27./10.*sq(length(v)) - 75./2.*sq(rdot) - 509./30.*(SI_G*m)/r)*(dot(v,Spin2)*Spin1 + dot(v, Spin1)*Spin2) 
                          + (174.*sq(length(v)) + 1386.*sq(rdot) + 1038./5.*(SI_G*m)/r)*rdot*(dot(n,Spin2)*Spin1 + dot(n, Spin1)*Spin2));
    return rRRSS;

}



mass e_Newtonian(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& mu = ip.mu;
    mass const& r = dp.r;
    Vector<mass, 3> const& v = dp.v;

    mass eN = mu*(1./2.*sq(length(v)) - (SI_G*m)/r);
    return eN;

}

mass e_PostNewtonian(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& mu = ip.mu;
    mass const& eta = ip.eta;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& v = dp.v;

    mass ePN = mu/SI_c2*(3./8.*(1. - 3.*eta)*std::pow(length(v), 4) + 1./2.*(3. + eta)*sq(length(v))*(SI_G*m)/r + 1./2.*eta*(SI_G*m)/r*sq(rdot) + 1./2.*sq((SI_G*m)/r));
    return ePN;

}

mass e_2PostNewtonian(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& mu = ip.mu;
    mass const& eta = ip.eta;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& v = dp.v;

    mass e2PN = mu/std::pow(SI_c, 4)*(5./16.*(1. - 17.*eta + 13.*sq(eta))*std::pow(length(v),6) -3./8.*eta*(1. - 3.*eta)*(SI_G*m)/r*std::pow(rdot,4)
              + 1./8.*(21. - 23.*eta - 27.*sq(eta))*(SI_G*m)/r*std::pow(length(v),4) + 1./8.*(14. - 55.*eta + 4.*sq(eta))*sq((SI_G*m)/r*length(v))
              + 1./4.*(1. - 15.*eta)*(SI_G*m)/r*sq(length(v)*rdot) - 1./4.*(2. + 15.*eta)*cb((SI_G*m)/r) +1./8.*(4. + 69.*eta + 4.*sq(eta))*sq((SI_G*m)/r*rdot));
    return e2PN;

}

mass e_3PostNewtonian(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& mu = ip.mu;
    mass const& eta = ip.eta;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& v = dp.v;

    mass e3PN = mu/std::pow(SI_c, 6)*((3/8. + 18469./840.*eta)*std::pow((SI_G*m)/r, 4.) + (5./4. - (6747./280. - 41./64.*sq(PI))*eta 
              - 21./4.*sq(eta) + 1./2.*cb(eta))*cb((SI_G*m)/r)*sq(length(v))
              + (3/2. + (2321/280. - 123./64.*sq(PI))*eta + 51./4.*sq(eta) + 7./2.*cb(eta))*cb((SI_G*m)/r)*sq(rdot )
              + 1./128.*(35. - 413.*eta + 1666.*sq(eta) - 2261.*cb(eta))*std::pow(length(v), 8) 
              + 1./16.*(135. - 194.*eta + 406.*sq(eta) - 108.*cb(eta))*sq((SI_G*m)/r*sq(length(v))) 
              + 1./16.*(12. + 248.*eta - 815.*sq(eta) - 324.*cb(eta))*sq((SI_G*m)/r*length(v)*rdot) 
              - 1./48.*eta*(731. - 492.*eta - 288.*sq(eta))*sq((SI_G*m)/r*sq(rdot)) 
              + 1./16.*(55. - 215.*eta + 116.*sq(eta) + 325.*cb(eta))*(SI_G*m)/r*std::pow(length(v), 6)
              + 1./16.*(5. - 25.*eta +25.*sq(eta))*(SI_G*m)/r*std::pow(rdot, 6) - 1./16.*(21. + 75.*eta - 375.*sq(eta))*(SI_G*m)/r*sq(sq(length(v))*rdot)
              - 1./16.*eta*(9. - 84.*eta + 165.*sq(eta))*(SI_G*m)/r*sq(length(v)*sq(rdot)));
    return e3PN;

}

mass e_4PostNewtonian(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& mu = ip.mu;
    mass const& eta = ip.eta;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& v = dp.v;

    mass E4PN0 = (63./256. - 1089./256.*eta + 7065./256.*sq(eta) - 10143./128.*cb(eta) + 21735./256.*std::pow(eta, 4))*std::pow(length(v), 10);

    mass E4PN1 = SI_G*m/r*((-35./128.*eta + 245./128.*sq(eta) - 245./64.*cb(eta) + 245./128.*std::pow(eta, 4))*std::pow(rdot, 8)
               + (25./32.*eta - 125./16.*sq(eta) + 185./8.*cb(eta) - 595./32.*std::pow(eta,4))*std::pow(rdot, 6)*sqlength(v)
               + (27./64.*eta + 243./32.*sq(eta) - 1683./32.*cb(eta) + 4851./64.*std::pow(eta, 4))*std::pow(rdot*length(v), 4)
               + (-147./32.*eta + 369./32.*sq(eta) + 423./8.*cb(eta) - 4655./32.*std::pow(eta, 4))*sq(rdot)*std::pow(length(v), 6)
               + (525./128. - 4011./128.*eta + 9507./128.*sq(eta) - 357./64.*cb(eta) + 15827./128.*std::pow(eta, 4))*std::pow(length(v), 8));

    mass E4PN2 = sq(SI_G*m/r)*((-4771./640.*eta - 461./8.*sq(eta) - 17./2.*cb(eta) + 15./2.*std::pow(eta, 4))*std::pow(rdot, 6)
               + (5347./384.*eta + 19465./96.*sq(eta) - 489./8.*cb(eta) - 135./2.*std::pow(eta, 4))*std::pow(rdot, 4)*sqlength(v)
               + (15./16. - 5893./128.*eta - 12995./64.*sq(eta) + 18511./64.*cb(eta) + 2845./16.*std::pow(eta, 4))*sq(rdot*sqlength(v))
               + (575./32. - 4489./128.*eta + 5129./64.*sq(eta) - 8289./64.*cb(eta) + 975./16.*std::pow(eta, 4))*std::pow(length(v), 6));

    mass E4PN3 = cb(SI_G*m/r)*(-(2599207./6720. + 6465./1024.*sq(PI)*eta - (103205./224. - 615./128.*sq(PI))*sq(eta) 
               + 69./32.*cb(eta) + 87./4.*std::pow(eta, 4))*std::pow(rdot, 4) 
               + (21./4. + (1086923./1680. + 333./512.*sq(PI))*eta + (206013./560. +123./64.*sq(PI))*sq(eta)
               - 2437./16.*cb(eta) - 141./2.*std::pow(eta, 4))*sq(rdot*length(v))
               + (273./16. - (22649399./10800. - 1071./1024.*sq(PI))*eta + (521063./10080. - 205./128.*sq(PI))*sq(eta)
               + 2373./32.*cb(eta) - 45./4.*std::pow(eta, 4))*std::pow(length(v), 4));

    mass E4PN4 = std::pow(SI_G*m/r, 4)*((9./4. - (1622437./12600. - 2645./96.*sq(PI))*eta - (289351./2520. + 1367./32.*sq(PI))*sq(eta)
               + 213./8.*cb(eta) + 15./2.*std::pow(eta, 4))*sq(rdot)
               + (15./16. + (1859363./16800. - 149./32.*sq(PI))*eta + (22963./5040. + 311./32.*sq(PI))*sq(eta) - 29./8.*cb(eta) + 1./2.*std::pow(eta, 4.))*sqlength(v));

    mass E4PN5 = std::pow(SI_G*m/r, 5)*(-3./8. - (1697177./25200. + 105./32.*sq(PI))*eta - (55111./720. - 11.*sq(PI))*sq(eta));

    mass e4PN = mu/std::pow(SI_c, 8)*(E4PN0 + E4PN1 + E4PN2 + E4PN3 + E4PN4 + E4PN5);
    return e4PN;

}

mass e_SpinOrbit(dynamicalParams const& dp, initParams const& ip){

    mass const& r = dp.r;
    Vector<mass, 3> const& sigma = dp.sigma;
    Vector<mass, 3> const& LN = dp.LN;

    mass eSO = SI_G/(SI_c2*cb(r))*dot(LN, sigma);
    return eSO;

}

mass e_SpinSpin(dynamicalParams const& dp, initParams const& ip){

    mass const& r = dp.r;
    Vector<mass, 3> const& n = dp.n;
    Vector<mass, 3> const& Spin1 = dp.Spin1;
    Vector<mass, 3> const& Spin2 = dp.Spin2;

    mass eSS = SI_G/(SI_c2*cb(r))*(3.*dot(n, Spin1)*dot(n, Spin2) - dot(Spin1, Spin2));
    return eSS;

}

mass e_PostNewtonianSO(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& dm = ip.dm;
    mass const& mu = ip.mu;
    mass const& eta = ip.eta;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& n = dp.n;
    Vector<mass, 3> const& v = dp.v;
    Vector<mass, 3> const& Spin = dp.Spin;
    Vector<mass, 3> const& Delta = dp.Delta;

    mass ePNSO = (SI_G*mu)/(2.*cb(SI_c)*sq(r))*dot(cross(n, v), (Delta*dm/m*((1. - 5.*eta)*sq(length(v)) + 3.*eta*(SI_G*m)/r) 
               - 3.*Spin*((1. + eta)*sq(length(v)) + eta*sq(rdot) - 4./3.*eta*(SI_G*m)/r)));
    return ePNSO;

}

mass edot_Newtonian(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& mu = ip.mu;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& v = dp.v;

    mass edotN = -8./15.*(cb(SI_G)*sq(m*mu))/(std::pow(SI_c, 5)*std::pow(r, 4))*(12.*sqlength(v) - 11.*sq(rdot));
    return edotN;

}

mass edot_PostNewtonian(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& mu = ip.mu;
    mass const& eta = ip.eta;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& v = dp.v;

    mass edotPN = - 2./105.*(cb(SI_G)*sq(m*mu))/(std::pow(SI_c, 7)*std::pow(r, 4))*((785. - 852.*eta)*std::pow(length(v), 4) 
                - 160.*(17. - eta)*(SI_G*m)/r*sqlength(v) + 8.*(367. - 15.*eta)*(SI_G*m)/r*sq(rdot) 
                - 2.*(1487. - 1392.*eta)*sq(rdot*length(v)) + 3.*(687. - 620.*eta)*std::pow(rdot, 4)
                + 16.*(1. - 4.*eta)*sq((SI_G*m)/r));
    return edotPN;

}

mass edot_2PostNewtonian(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& mu = ip.mu;
    mass const& eta = ip.eta;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& v = dp.v;

    mass edot2PN = -8./15.*(cb(SI_G)*sq(m*mu))/(std::pow(SI_c, 9)*std::pow(r, 4))*(1./42.*(1692. - 5497.*eta + 4430.*sq(eta))*std::pow(length(v), 6)
                 - 1./14.*(1719. - 10278.*eta + 6292.*sq(eta))*sq(sqlength(v)*rdot) - 1./21.*(4446. - 5237.*eta + 1393.*sq(eta))*(SI_G*m)/r*std::pow(rdot, 4)
                 + 1./14.*(2018. - 15207.*eta + 7572.*sq(eta))*sq(length(v)*sq(rdot)) + 1./7.*(4987. - 8513.*eta + 2165.*sq(eta))*(SI_G*m)/r*sq(rdot*length(v))
                 + 1./756.*(281473. + 81828.*eta + 4368.*sq(eta))*sq((SI_G*m)/r*length(v))
                 - 1./42.*(2501. - 20234.*eta + 8404.*sq(eta))*std::pow(rdot, 6) - 1./63.*(33510. + 60971.*eta + 14290.*sq(eta))*(SI_G*m)/r*std::pow(rdot, 4)
                 - 1./252.*(106319. + 9798.*eta + 5376.*sq(eta))*sq((SI_G*m)/r*rdot) + 2./63.*(-253. + 1026.*eta - 56.*sq(eta))*cb((SI_G*m)/r));
    return edot2PN;

}

mass edot_25PostNewtonian(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& mu = ip.mu;
    mass const& eta = ip.eta;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& v = dp.v;

    mass edot25PN = -32./5.*(cb(SI_G)*sq(m*mu))/(std::pow(SI_c, 10)*std::pow(r, 4))*rdot*eta*(-12349./210.*(SI_G*m)/r*std::pow(rdot, 4)
                  + 4524./35.*(SI_G*m)/r*sq(length(v)*rdot) + 2753./126.*sq((SI_G*m)/r*length(v)) - 985./14.*(SI_G*m)/r*std::pow(rdot, 4)
                  + 13981./630.*sq((SI_G*m)/r*rdot) - 1./315.*cb((SI_G*m)/r));
    return edot25PN;

}

mass edot_SpinOrbit(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& dm = ip.dm;
    mass const& mu = ip.mu;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& v = dp.v;
    Vector<mass, 3> const& Spin = dp.Spin;
    Vector<mass, 3> const& Delta = dp.Delta;
    Vector<mass, 3> const& LN = dp.LN;

    mass edotSO = -8./15.*(cb(SI_G)*sq(m*mu))/(std::pow(SI_c, 7)*std::pow(r, 6))*dot(LN, Spin*(78.*sq(rdot) - 80.*sqlength(v) - 8.*(SI_G*m)/r)
                + dm/m*Delta*(51.*sq(rdot) - 43.*sqlength(v) + 4.*(SI_G*m)/r));
    return edotSO;

}

mass edot_SpinSpin(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& mu = ip.mu;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& n = dp.n;
    Vector<mass, 3> const& v = dp.v;
    Vector<mass, 3> const& Spin1 = dp.Spin1;
    Vector<mass, 3> const& Spin2 = dp.Spin2;

    mass edotSS = -4./15.*(cb(SI_G)*sq(m*mu))/(std::pow(SI_c, 7)*std::pow(r, 6))*(-3.*dot(n, Spin1)*dot(n, Spin2)*(168.*sqlength(v) - 269.*sq(rdot))
                + 3*dot(Spin1, Spin2)*(47*sqlength(v) - 55*sq(rdot)) + 71*dot(v, Spin1)*dot(v, Spin2) 
                - 171.*rdot*(dot(v, Spin1)*dot(n, Spin2) + dot(n, Spin1)*dot(v, Spin2)));
    return edotSS;

}

mass edot_PostNewtonianSO(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& dm = ip.dm;
    mass const& mu = ip.mu;
    mass const& eta = ip.eta;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& n = dp.n;
    Vector<mass, 3> const& v = dp.v;
    Vector<mass, 3> const& Spin = dp.Spin;
    Vector<mass, 3> const& Delta = dp.Delta;

    mass edotPNSO = -8./105.*(cb(SI_G)*sq(m*mu))/(std::pow(SI_c, 10)*std::pow(r, 5))*(dot(cross(n, v), Spin)*(std::pow(rdot, 4)*(3144.*eta - 2244.)
                  + sq((SI_G*m)/r)*(944. + 390.*eta) + (SI_G*m)/r*sq(rdot)*(526.*eta - 3223.) + sq(rdot*length(v))*(3519. - 5004.*eta)
                  + (SI_G*m)/r*sqlength(v)*(3805. - 224.*eta) + std::pow(length(v), 4)*(1810.*eta - 1207.))
                  + dot(cross(n, v), Delta)*dm/m*(std::pow(rdot, 4)*(2676.*eta - 7941./4.) + sq((SI_G*m)/r)*(238.*eta - 137.) 
                  + (SI_G*m)/r*sq(rdot)*(1199.*eta - 7327.) + sq(rdot*length(v))*(2364. - 3621.*eta)
                  + (SI_G*m)/r*sqlength(v)*(5387./2. - 497.*eta) + std::pow(length(v), 4.)*(1040.*eta - 2603./4.)));
    return edotPNSO;

}




Vector<mass, 3> l_PostNewtonian(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& eta = ip.eta;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& v = dp.v;
    Vector<mass, 3> const& LN = dp.LN;

    Vector<mass, 3> lPN= LN/SI_c2*(1./2.*sq(length(v))*(1. - 3.*eta) + (3. + eta)*(SI_G*m)/r);
    return  lPN;

}

Vector<mass, 3> l_2PostNewtonian(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& eta = ip.eta;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& v = dp.v;
    Vector<mass, 3> const& LN = dp.LN;

    Vector<mass, 3> l2PN = LN/std::pow(SI_c,4)*(3./8.*(1. - 7.*eta + 13.*sq(eta))*std::pow(length(v), 4) - 1./2.*eta*(2. + 5.*eta)*(SI_G*m)/r*sq(rdot) 
                         + 1./2.*(7 - 10.*eta  - 9.*sq(eta))*(SI_G*m)/r*sq(length(v)) + 1./4.*(14. - 41.*eta + 4.*sq(eta))*sq((SI_G*m)/r));
    return l2PN;

}

Vector<mass, 3> l_3PostNewtonian(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& eta = ip.eta;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& v = dp.v;
    Vector<mass, 3> const& LN = dp.LN;

    Vector<mass, 3> l3PN = LN/std::pow(SI_c, 6)*((5./2. - (5199./280. - 41./32.*sq(PI))*eta - 7.*sq(eta) + cb(eta))*cb((SI_G*m)/r)
                         + 1./16.*(5 - 59.*eta + 238.*sq(eta) - 323.*cb(eta))*std::pow(length(v),6)
                         + 1./12.*(135. - 322.*eta + 315.*sq(eta) - 108.*cb(eta))*sq((SI_G*m)/r*length(v)) 
                         + 1./24.*(12. - 287.*eta - 951.*sq(eta) - 324.*cb(eta))*sq((SI_G*m)/r*rdot)
                         + 1./8.*(33. - 142.*eta + 106.*sq(eta) + 195.*cb(eta))*(SI_G*m)/r*std::pow(length(v), 4)
                         - 1./4.*eta*(12. - 7.*eta - 75.*sq(eta))*(SI_G*m)/r*sq(length(v)*rdot)
                         + 3./8.*eta*(2. - 2.*eta - 11.*sq(eta))*(SI_G*m)/r*std::pow(rdot, 4));
    return l3PN;

}

Vector<mass, 3> l_4PostNewtonian(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& eta = ip.eta;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& v = dp.v;
    Vector<mass, 3> const& LN = dp.LN;

    mass L4PN0 = (35./128. - 605./128.*eta + 3925./128.*sq(eta) - 5635./54.*cb(eta) + 12075./128.*std::pow(eta, 4))*std::pow(length(v), 8);

    mass L4PN1 = SI_G*m/r*((-5/8*eta + 15/8*sq(eta) + 45/16*cb(eta) - 85/16*std::pow(eta, 4))*std::pow(rdot, 6)
               + (3.*eta - 45./4.*sq(eta) - 135./16.*cb(eta) + 693./16.*std::pow(eta, 4))*sq(length(v)*sq(rdot))
               - (53./8.*eta - 423./16.*sq(eta) - 299./16.*cb(eta) + 1995./16.*std::pow(eta, 4))*sq(sqlength(v)*rdot)
               + (75./16. - 151./4.*eta + 1553./16.*sq(eta) - 425./16.*cb(eta) - 2261./16.*std::pow(eta, 4))*std::pow(length(v), 6));

    mass L4PN2 = sq(SI_G*m/r)*((14773./320.*eta + 3235./48.*sq(eta) - 155./4.*cb(eta) - 27.*std::pow(eta, 4))*std::pow(rdot, 4)
               + (3./4. - 5551./60.*eta - 256./3.*sq(eta) + 4459./16.*cb(eta) + 569./4.*std::pow(eta, 4))*sq(rdot*length(v))
               + (345./26. - 65491./960.*eta + 12427./96.*sq(eta) - 3845./32.*cb(eta) + 585./8.*std::pow(eta, 4))*std::pow(length(v), 4));

    mass L4PN3 = cb(SI_G*m/r)*((7./2. + (7775977./16800. + 447./256.*sq(PI))*eta + 121449./560.*sq(eta) - 1025./8.*cb(eta) - 47.*std::pow(eta, 4))*sq(rdot)
               + (91./4. - (13576009./50400. - 469./256.*sq(PI))*eta + (276433./5040. - 41./16.*sq(PI))*sq(eta) + 637./8.*cb(eta) - 15.*std::pow(eta, 4))*sqlength(v));

    mass L4PN4 = std::pow(SI_G*m/r, 4)*(15./8. + (3809041./25200. - 85./8.*sq(PI))*eta - (20131./420. - 663./32.*sq(PI))*sq(eta) - 15./4.*cb(eta) + std::pow(eta, 4));

    Vector<mass, 3> l4PN = LN/std::pow(SI_c, 8)*(L4PN0 + L4PN1 + L4PN2 + L4PN3 + L4PN4);
    return l4PN;

}

Vector<mass, 3> l_SpinOrbit(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& mu = ip.mu;
    mass const& r = dp.r;
    Vector<mass, 3> const& n = dp.n;
    Vector<mass, 3> const& v = dp.v;
    Vector<mass, 3> const& Spin = dp.Spin;
    Vector<mass, 3> const& sigma = dp.sigma;

    Vector<mass, 3> lSO = mu/(SI_c2*m)*((SI_G*m)/r*cross(n, cross(n, (2.*Spin + sigma))) - 1./2.*cross(v, cross(v, sigma)));
    return lSO;

}

//Spin equations: Bohé, Alejandro et al. NNSO - CQG30(13)075017.pdf eq3.4 

Vector<mass, 3> spin1(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& dm = ip.dm;
    mass const& eta = ip.eta;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& n = dp.n;
    Vector<mass, 3> const& v = dp.v;
    Vector<mass, 3> const& Spin2 = dp.Spin2;

    mass alpha1PN; //eq3.4b
    mass alpha2PN; //eq3.4c
    mass alpha3PN; //eq3.4d
    
    alpha1PN = (SI_G*m)/sq(r)*(3./4. + 1./2.*eta - 3./4.*dm/m);
    
    alpha2PN = (SI_G*m)/sq(r)*((-3/2*eta + 3/4*sq(eta) - 3/2*eta*dm/m)*sq(rdot) 
             + (1./16. + 1./18.*eta - 3./8.*sq(eta) + dm/m*(-1./16. + 1./2.*eta))*sqlength(v))
             + sq(SI_G*m)/cb(r)*(-1./4. - 3./8.*eta + 1./2.*sq(eta) + dm/m*(1./4. - 1./8.*eta));
    
    alpha2PN = (SI_G*m)/sq(r)*((15./8.*eta - 195./32.*sq(eta) + 15./16.*cb(eta) + dm/m*(15./8.*eta - 75./32.*sq(eta)))*std::pow(length(n)*length(v), 4)
             + (-3.*eta + 291./32.*sq(eta) - 45./16.*cb(eta) + dm/m*(-3.*eta + 177./32.*sq(eta)))*sq(rdot)*sqlength(v)
             + (1./32. + 19./16.*eta - 31./8.*sq(eta) + 17./16.*cb(eta) + dm/m*(-1./32. + 3./4.*eta - 11./8.*sq(eta)))*std::pow(length(v),4))
             + sq(SI_G*m)/cb(r)*((1./4. - 525./32.*eta - 159./16.*sq(eta) + 13./4.*cb(eta) + dm/m*(-1./4. - 75./32.*eta - 87./16.*sq(eta)))*sq(rdot)
             + (3./16. + 27./4.*eta + 75./32.*sq(eta) - 9./8.*cb(eta) + dm/m*(-3./16. + 9./8.*eta + 35./32.*sq(eta)))*sqlength(v))
             + cb(SI_G*m)/std::pow(r,4)*(7./16. - 9./4.*eta - 9./8.*sq(eta)*1./2.*cb(eta) + dm/m*(-7./16. - 1./8.*eta - 1./8.*sq(eta)));

    Vector<mass, 3> Spin1 = cross(n,v)*(1./SI_c2*alpha1PN + 1./sq(SI_c2)*alpha2PN + 1./cb(SI_c2)*alpha3PN);

    return Spin1;

}

Vector<mass, 3> spin2(dynamicalParams const& dp, initParams const& ip){

    mass const& m = ip.m;
    mass const& dm = ip.dm;
    mass const& eta = ip.eta;
    mass const& r = dp.r;
    mass const& rdot = dp.rdot;
    Vector<mass, 3> const& n = dp.n;
    Vector<mass, 3> const& v = dp.v;
    Vector<mass, 3> const& Spin1 = dp.Spin1;

    mass alpha1PN;
    mass alpha2PN;
    mass alpha3PN;
    
    alpha1PN = (SI_G*m)/sq(r)*(3./4. + 1./2.*eta + 3./4.*dm/m);
    
    alpha2PN = (SI_G*m)/sq(r)*((-3./2.*eta + 3/4.*sq(eta) + 3./2.*eta*dm/m)*sq(rdot) 
             + (1./16. + 1./18.*eta - 3./8.*sq(eta) - dm/m*(-1./16. + 1/2*eta))*sqlength(v))
             + sq(SI_G*m)/cb(r)*(-1./4. - 3./8.*eta + 1./2.*sq(eta) - dm/m*(1./4. - 1./8.*eta));
    
    alpha2PN = (SI_G*m)/sq(r)*((15./8.*eta - 195./32.*sq(eta) + 15./16.*cb(eta) - dm/m*(15./8.*eta - 75./32.*sq(eta)))*std::pow(rdot, 4)
             + (-3.*eta + 291./32.*sq(eta) - 45./16.*cb(eta) - dm/m*(-3.*eta + 177./32.*sq(eta)))*sq(rdot)*sqlength(v)
             + (1./32. + 19./16.*eta - 31./8.*sq(eta) + 17./16.*cb(eta) - dm/m*(-1./32. + 3./4.*eta - 11./8.*sq(eta)))*std::pow(length(v),4))
             + sq(SI_G*m)/cb(r)*((1./4. - 525./32.*eta - 159./16.*sq(eta) + 13./4.*cb(eta) - dm/m*(-1./4. - 75./32.*eta - 87./16.*sq(eta)))*sq(rdot
             + (3./16. + 27./4.*eta + 75./32.*sq(eta) - 9./8.*cb(eta) - dm/m*(-3./16. + 9./8.*eta + 35./32.*sq(eta)))*sqlength(v)))
             + cb(SI_G*m)/std::pow(r, 4)*(7./16. - 9./4.*eta - 9./8.*sq(eta)*1./2.*cb(eta) - dm/m*(-7./16. - 1./8.*eta - 1./8.*sq(eta)));

    Vector<mass, 3> Spin2 = cross(n,v)*(1./SI_c2*alpha1PN + 1./sq(SI_c2)*alpha2PN + 1./cb(SI_c2)*alpha3PN);

    return Spin2;

}