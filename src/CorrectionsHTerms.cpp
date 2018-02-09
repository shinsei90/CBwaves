#include <Corrections.hpp>
#include <Vector.hpp>

#include <iostream>
#include <cmath>

auto sq = [](auto const& a){ return a*a;};
auto cb = [](auto const& a){ return a*a*a; };

Vector<mass, 3> h_Q(mass m, mass r){
    
    Vector<mass, 3> hQ;

    for(int i = 0; i < 3; ++i){
        for(int j; j < 3; ++j){
            hQ = 2*(v[i]*v[j] - (SI_G*m)/r*n[i]*n[j]);
        }
    }

    return hQ;

}

Vector<mass, 3> h_P05Q(mass m, mass r){

    Vector<mass, 3> hP05Q;

    for(int i = 0; i < 3; ++i){
        for(int j; j < 3; ++j){
            hP05Q = dm/(SI_c*m)*(3*(SI_G*m)/r*(n[i]*v[j] + n[j]*v[i] - rdot*n[i]*n[j])*dot(N,n) + ((SI_G*m)/r*n[i]*n[j] - 2*v[i]*v[j])*dot(N,v));
        }
    }

    return hP05Q;
}

Vector<mass, 3> h_PQ(mass m, mass r){

    Vector<mass, 3> hPQ;

    for(int i = 0; i < 3, ++i){
        for(int j; j < 3; ++j){
            hPQ = 1/(3*SI_c2)*(1 - 3*eta)*(4*(SI_G*m)/r*(3*rdot*n[i]*n[j] - 8*(n[i]*v[j] + n[j]*v[i]))*dot(N,n)*sq(dot(N,v))
                + 2*(3*v[i]*v[j] - (SI_G*m)/r*n[i]*n[j]) 
                + (SI_G*m)/r*((3*sqlength(v) - 15*sq(rdot) + 7*(SI_G*m)/r)*n[i]*n[j] + 30*rdot*(n[i]*v[j] + n[j]*v[i]) - 14*v[i]*v[j])*sq(dot(N,n)))
                + 3/4*(SI_G*m)/r*rdot*(5 + 3*eta)*(n[i]*v[j] + n[j]*v[i]) + ((1 - 3*eta)*sqlength(v) - 2/3*(2 - 3*eta)*(SI_G*m)/r)*v[i]*v[j]
                + (SI_G*m)/r*((1 - 3*eta)*sq(rdot) - 1/3*(10 + 3*eta)*sqlength(v) + 29/3*(SI_G*m)/r)*n[i]*n[j];
        }
    }

    return hPQ;
}

Vector<mass, 3> h_P15Q(mass m, mass r){

    Vector<mass, 3> hP15Q;

    for(int i = 0; i < 3, ++i){
        for(int j; j < 3; ++j){
            hP15Q = dm/(m*cb(SI_c))*(1 - 2*eta)*(1/4*(SI_G*m)/r*((45*sq(rdot) - 9*sqlength(v) - 28*(SI_G*m)/r)*n[i]*n[j]
                  + 58*n[i]*n[j] - 108*rdot*(n[i]*v[j] + n[j]*v[i]))*sq(dot(N,n))*dot(N,v)
                  + 1/2*((SI_G*m)/r*n[i]*n[j] - 4*v[i]*v[j])*cb(dot(N,v)) + (SI_G*m)/r*(5/4*(3*sqlength(v) - 7*sq(rdot) + 6*(SI_G*m)/r)*rdot*n[i]*n[j]
                  - 1/6*(21*sqlength(v - 105*sq(rdot) + 44*(SI_G*m)/r))*(n[i]*v[j] + n[j]*v[i]) - 17/2*v[i]*v[j])*cb(dot(N,n))
                  + 3/2*(SI_G*m)/r*(10*(n[i]*v[j] + n[j]*v[i]) - 3*rdot*n[i]*n[j])*dot(N,n)*sq(dot(N,v)))
                  + dm/(m*cb(SI_c))*1/12*(SI_G*m)/r*dot(N,n)*(n[i]*n[j]*rdot*(sq(rdot)*(15 - 90*eta) - sqlength(v)*(63 - 54*eta) 
                  - (SI_G*m)/r*(242 - 24*eta)) - rdot*v[i]*v[j]*(86 + 24*eta)
                  + 2*(n[i]*v[j] + n[j]*v[i])*(sq(rdot)*(63 + 54*eta) - (SI_G*m)/r*(128 - 36*eta) - sqlength(v)*(33 - 18*eta)))
                  + dm/(m*cb(SI_c)*dot(N,v)*(1/2*v[i]*v[j]*((SI_G*m)/r*(3 - 8*eta) - 2*sqlength(v)*(1 - 5*eta)) - (n[i]*v[j] 
                  + n[j]*v[i])*(SI_G*m)/r*rdot*(7 + 4*eta) - n[i]*n[j]*(SI_G*m)/r*(3/4*(1 - 2*eta)*sq(rdot) + 1/3*(26 - 3*eta)*(SI_G*m)/r
                  - 1/4*(7 - 2*eta)*sqlength(v)));
        }
    }

    return hP15Q;

}

Vector<mass, 3> h_P2Q(mass m, mass r){

    Vector<mass, 3> hP2Q;

    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3, ++j){

            hP2Q = 1/sq(SI_c2)*(1/60*(1 - 5*eta + 5*sq(eta))*(24*std::pow(dot(N,v),4)*(5*v[i]*v[j] - m/r*n[i]*n[j])
                 + m/r*std::pow(dot(N,n),4)*(2*(175*m/r - 465*sq(rdot) + 93*sqlength(v))*v[i]*v[j]
                 + 30*rdot*(63*sq(rdot) - 50*m/r - 27*sqlength(v))*(n[i]*v[j] + n[j]*v[i])
                 + (1155*m/r*sq(rdot) - 172*sq(m/r) - 945*sq(sq(rdot,4)) - 159*m/r*sqlength(v) + 630*sq(rdot*length(v)) - 45*sq(sqlength(v)))*n[i]*n[j])
                 + 24*m/r*cb(dot(N,n))*dot(N,v)*(87*rdot*v[i]*v[j] + 5*rdot*(14*sq(rdot) - 15*m/r - 6*sqlength(v))*n[i]*n[j]
                 + 16*(5*m/r - 10*sq(edot) + sqlength(v))*(n[i]*v[j] + n[j]*v[i])) + 288*m/r*dot(N,n)*cb(dot(N,v))*(rdot*n[i]*n[j]
                 - 4*(n[i]*v[j] + n[j]*v[i])) + 24*sq(dot(N,n)*dot(N,v))*((35*m/r - 45*sq(rdot) + 9*sqlength(v))*n[i]*n[j] - 76*v[i]*v[j]
                 + 126*(n[i]*v[j] + n[j]*v[i]))) + 1/15*sq(dot(N,v))*((5*(25 - 78*eta + 12*sq(eta))*m/r - (18 - 65*eta + 45*sq(eta))*sqlength(v)
                 + 9*(1 - 5*eta + 5*sq(eta)))*m/r*n[i]*n[j] + 3*(5*(1 - 9*eta + 21*sq(eta))*sqlength(v) - 2*(4 - 25*eta + 45*sq(eta))*m/r)*v[i]*v[j]
                 + 18*(6 - 15*eta - 10*sq(eta))*m/r*rdot*(n[i]*v[j] + n[j]*v[i])) + 1/15*dot(N,n)*dot(N,v)*m/r*((3*(36 - 145*eta + 150*sq(eta))*sqlength(v)
                 - 5*(127 - 392*eta + 36*sq(eta))*m/r - 15*(2 - 15*eta + 30*sq(eta))*sq(rdot))*rdot*n[i]*n[j]
                 + 6*(98 - 295*eta - 30*sq(eta))*rdot*v[i]*v[j] 2*(5*(66 - 221*eta + 96*sq(eta))*m/r - 9*(18 - 45*eta - 40*sq(eta))*sq(rdot)
                 - (66 - 265*eta + 365*sq(eta))*sqlength(v))*(n[i]*v[j] + n[j]*v[i])) + 1/60*sq(dot(N,n))*m/r*((3*(33 - 130*eta + 150*sq(eta))*sq(sqlength(v))
                 + 105*(1 - 10*eta + 30*sq(eta))*std::pow(rdot,4) + 15*(181 - 572*eta + 84*sq(eta))*m/r*sq(rdot) - (131 - 770*eta + 930*sq(eta))*m/r*sqlength(v)
                 - 60*(9 - 40*eta + 60*sq(eta))*sq(length(v)*rdot) - 8*(131 - 390*eta + 30*sq(eta))*sq(m/r))*n[i]*n[j]
                 + 4*((12 + 5*eta - 315*sq(eta))*sqlength(v) - 9*(39 - 115*eta - 35*sq(eta))*sq(rdot) + 5*(29 - 104*eta + 84*sq(eta))*m/r)*v[i]*v[j]
                 + 4*(15*(18 - 40*eta - 75*sq(eta))*sq(rdot) - 5*(197 - 640*eta + 180*sq(eta))*m/r 
                 + 3*(21 - 130*eta + 375*sq(eta))*sqlength(v))*rdot*(n[i]*v[j] + n[j]*v[i]))
                 + 1/60*(((467 + 780*eta - 120*sq(eta))*m/r*sqlength(v) - 15*(61 - 96*eta + 48*sq(eta))*m/r*sq(rdot)
                 - (144 - 265*eta - 135*sq(eta))*std::pow(length(v),4) + 6*(24 - 95*eta + 75*sq(eta))*sq(length(v)*rdot) - 2*(642 + 545*eta)*sq(m/r)
                 - 45*(1 - 5*eta + 5*sq(eta))*std::pow(rdot,4))*m/r*n[i]*n[j] + (4*(69 + 10*eta - 135*sq(eta))*m/r*sqlength(v)
                 - 12*(3 + 60*eta + 25*sq(eta))*m/r*sq(rdot) + 45*(1 - 7*eta + 13*sq(eta))*std::pow(length(v),4) 
                 - 10*(56 + 165*eta - 12*sq(eta))*sq(m/r))*v[i]*v[j] + 4*(2*(36 + 5*eta - 75*sq(eta))*sqlength(v) - 6*(7 - 15*eta - 15*sq(eta))*sq(rdot)
                 + 5*(35 + 45*eta + 36*sq(eta))*m/r)*m/r*rdot*(n[i]*v[j] + n[j]*v[i])));

        }
    }

    return hP2Q;

}

Vector<mass, 3> h_P2QSS(mass r, Vector<mass,3> const& Spin1, Vector<mass,3> const& Spin2){

    Vector<mass, 3> hP2QSS;

    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){

            hP2QSS = -6*SI_G/(SI_c2*mu*cb(r))*(n[i]*n[j]*(dot(Spin1,Spin2) - 5*dot(n,Spin1)*dot(n,Spin2)) 
                   + 2*(n[i]*Spin1[j] + n[j]*Spin1[i])*dot(n,Spin2) + 2*(n[i]*Spin2[j] + n[j]*Spin2[i])*dot(n,Spin1));

        }
    }

    return hP2QSS;

}

Vector<mass, 3> h_P15QSO(mass r, Vector<mass,3> const& Delta, Vector<mass, 3> Spin){

    Vector<mass, 3> hP15QSO;

    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){
            hPQ15SO = 2*SI_G/(SI_c2*sq(r))*(n[i]*n[j]*dot(cross(n, v), (12*Spin + 6*dm/m*Delta))
                  - (n[i]*cross(v, (9*Spin + 5*dm/m*Delta))[j] + n[j]*cross(v, (9*Spin + 5*dm/m*Delta))[i])
                  + (3*rdot*dot(N,n) - 2*dot(N,v))*(cross((Spin + dm/m*Delta), N)[i]*n[j] + cross((Spin + dm/m*Delta), N)[j]*n[i])
                  - (v[i]*cross(n, (2*Spin + 2*dm/m*Delta))[j] + v[j]*cross(n, (2*Spin + 2*dm/m*Delta))[i])
                  + rdot*(n[i]*cross(n, (12*Spin + 6*dm/m*Delta))[j] + n[j]*cross(n, (12*Spin + 6*dm/m*Delta))[i])
                  - 2*dot(N,n)*(cross((Spin + dm/m*Delta), N)[i]*v[j] + cross((Spin + dm/m*Delta), N)[j]*v[i]));
        }
    }

    return hP15QSO;

}
Vector<mass, 3> h_PQSO(mass r, Vector<mass,3> const& Delta){

    Vector<mass, 3> hPQSO;

    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < 3; ++j){

            hPQS0 = 2*SI_G/sq(SI_c*r)*(cross(Delta, N)[i]*n[j] + cross(Delta, N)[j]*n[i]);
        }
    }

    return hPQSO;

}