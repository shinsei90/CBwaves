#include <Corrections.hpp>
#include <Vector.hpp>

#include <iostream>
#include <cmath>

auto sq = [](auto const& a){ return a*a;};
auto cb = [](auto const& a){ return a*a*a; };

Vector<mass, 3> h_Q(mass m, mass r){
    
    Vector<mass, 3> hQ;

    for(int i = 0; i < 3, ++i){
        for(int j; j < 3; ++j){
            hQ = 2*(v[i]*v[j] - (SI_G*m)/r*n[i]*n[j]);
        }
    }

    return hQ;

}

Vector<mass, 3> h_P05Q(mass m, mass r){

    Vector<mass, 3> hP05Q;

    for(int i = 0; i < 3, ++i){
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