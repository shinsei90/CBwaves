#pragma once

using mass = double;                    // real type for mass dimension variables
using solver_internal = double;         // scalars used inside the solver

// If we are using c = G = 1 
#define SI_G      1                         // The universal gravitational constant
#define SI_c      1                         // The speed of light
#define MSUN      1477.1                     // The mass of the sun in meters

// If we are not using c = G = 1 system
// #define SI_G      6.67428E-11               // The universal gravitational constant in [m^3 kg^-1 s^-2]
// #define SI_c      2.99792458E8              // The speed of light in [m/s]
// #define MSUN      1477. * SI_c2/SI_G        // The mass of the sun in [kg]

// Constants that are not effected by G or c
#define SI_c2     (SI_c*SI_c)               // The square of the speed of lignt
#define SI_ly     9460730472580800.0        // Light year in [m]
#define SI_pc     (3.26156*SI_ly)           // Prasec
#define PI        3.141592653589793         // π