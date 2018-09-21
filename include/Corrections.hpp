#pragma once

#include "Config.hpp"
#include <DynamicalParams.hpp>
#include <InitParams.hpp>
#include "Vector.hpp"



//std::array<double, 3> dr(const double r, const mass m);

//Corrections of the orbit

state c_Newtonian(dynamicalParams const& dp, initParams const& ip);
state c_PostNewtonian(dynamicalParams const& dp, initParams const& ip);
state c_2PostNewtonian(dynamicalParams const& dp, initParams const& ip);
state c_3PostNewtonian(dynamicalParams const& dp, initParams const& ip);
state c_4PostNewtonian(dynamicalParams const& dp, initParams const& ip);
state c_SpinOrbit(dynamicalParams const& dp);
state c_SpinSpin(dynamicalParams const& dp, initParams const& ip);
state c_BT_RR(dynamicalParams const& dp, initParams const& ip);
state c_PostNewtonianSO(dynamicalParams const& dp, initParams const& ip);
state c_2PostNewtonianSO(dynamicalParams const& dp, initParams const& ip); //Bohe et al. NNSO - CQG30(13)075017
state c_RR1PostNewtonian(dynamicalParams const& dp, initParams const& ip);
state c_RRSO(dynamicalParams const& dp, initParams const& ip);
state c_RRSS(dynamicalParams const& dp, initParams const& ip);

//Energy
mass e_Newtonian(dynamicalParams const& dp, initParams const& ip);
mass e_PostNewtonian(dynamicalParams const& dp, initParams const& ip);
mass e_2PostNewtonian(dynamicalParams const& dp, initParams const& ip);
mass e_3PostNewtonian(dynamicalParams const& dp, initParams const& ip);
mass e_4PostNewtonian(dynamicalParams const& dp, initParams const& ip);
mass e_SpinOrbit(dynamicalParams const& dp, initParams const& ip);
mass e_SpinSpin(dynamicalParams const& dp, initParams const& ip);
mass e_PostNewtonianSO(dynamicalParams const& dp, initParams const& ip);

mass edot_Newtonian(dynamicalParams const& dp, initParams const& ip);
mass edot_PostNewtonian(dynamicalParams const& dp, initParams const& ip);
mass edot_2PostNewtonian(dynamicalParams const& dp, initParams const& ip);
mass edot_25PostNewtonian(dynamicalParams const& dp, initParams const& ip);
mass edot_SpinOrbit(dynamicalParams const& dp, initParams const& ip);
mass edot_SpinSpin(dynamicalParams const& dp, initParams const& ip);
mass edot_PostNewtonianSO(dynamicalParams const& dp, initParams const& ip);

//Momentum
Vector<mass, 3> l_PostNewtonian(dynamicalParams const& dp, initParams const& ip);
Vector<mass, 3> l_2PostNewtonian(dynamicalParams const& dp, initParams const& ip);
Vector<mass, 3> l_3PostNewtonian(dynamicalParams const& dp, initParams const& ip);
Vector<mass, 3> l_4PostNewtonian(dynamicalParams const& dp, initParams const& ip);
Vector<mass, 3> l_SpinOrbit(dynamicalParams const& dp, initParams const& ip);

//Spin equations

Vector<mass, 3> spin1(dynamicalParams const& dp, initParams const& ip);
Vector<mass, 3> spin2(dynamicalParams const& dp, initParams const& ip);

//Corrections of the Waveform

mass h_Q(dynamicalParams const& dp, initParams const& ip);
mass h_P05Q(dynamicalParams const& dp, initParams const& ip);
mass h_PQ(dynamicalParams const& dp, initParams const& ip);
mass h_P15Q(dynamicalParams const& dp, initParams const& ip);
mass h_P2Q(dynamicalParams const& dp, initParams const& ip);
mass h_PQSO(dynamicalParams const& dp, initParams const& ip);
mass h_P15QSO(dynamicalParams const& dp, initParams const& ip);
mass h_P2QSS(dynamicalParams const& dp, initParams const& ip);