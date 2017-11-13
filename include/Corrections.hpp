#include "Config.hpp"
#include "Vector.hpp"

std::array<double, 3> dr(const double r, const mass m);

//Corrections of the orbit

Vector<mass, 3> c_Newtonian(mass m, mass r);
Vector<mass, 3> c_PostNewtonian(mass m, mass r);
Vector<mass, 3> c_2PostNewtonian(mass m, mass r);
Vector<mass, 3> c_3PostNewtonian(mass m, mass r);
Vector<mass, 3> c_SpinOrbit(mass m, mass r, Vector<mass, 3> const& Spin, Vector<mass, 3> const& sigma);
Vector<mass, 3> c_SpinSpin(mass m, mass r, Vector<mass, 3> const& Spin1, Vector<mass, 3> const& Spin2);
Vector<mass, 3> c_BT_RR(mass m, mass r);
Vector<mass, 3> c_PostNewtonianSO(mass m, mass r, Vector<mass, 3> const& Spin, Vector<mass, 3> const& Delta);
Vector<mass, 3> c_RR1PostNewtonian(mass m, mass r);
Vector<mass, 3> c_RRSO(mass m, mass r, Vector<mass, 3> const& Spin, Vector<mass, 3> const& sigma, Vector<mass, 3> const& LN);
Vector<mass, 3> c_RRSS(mass m, mass r, Vector<mass, 3> const& Spin1, Vector<mass, 3> const& Spin2);

mass e_Newtonian(mass m, mass r);
mass e_PostNewtonian(mass m, mass r);
mass e_2PostNewtonian(mass m, mass r);
mass e_3PostNewtonian(mass m, mass r);
mass e_SpinOrbit(mass r, Vector<mass, 3> const& sigma, Vector<mass, 3> const& LN);
mass e_SpinSpin(mass r, Vector<mass, 3> const& Spin1, Vector<mass, 3> const& Spin2);
mass e_PostNewtonianSO(mass m, mass r, Vector<mass, 3> const& Spin, Vector<mass, 3> const& Delta);

Vector<mass, 3> l_PostNewtonian(mass m, mass r, Vector<mass, 3> const& LN);
Vector<mass, 3> l_2PostNewtonian(mass m, mass r, Vector<mass, 3> const& LN);
Vector<mass, 3> l_3PostNewtonian(mass m, mass r, Vector<mass, 3> const& LN);
Vector<mass, 3> l_SpinOrbit(mass m, mass r, Vector<mass, 3> const& Spin, Vector<mass, 3> const& sigma);

//Spin equations

Vector<mass, 3> Spin1(mass m, mass r);
Vector<mass, 3> Spin2(mass m, mass r);

//Corrections of the Waveform

Vector<mass, 3> h_Q(mass m, mass r);
Vector<mass, 3> h_P05Q(mass m, mass r);
Vector<mass, 3> h_PQ(mass m, mass r);
Vector<mass, 3> h_P15Q(mass m, mass r);
Vector<mass, 3> h_P2Q(mass m, mass r);
Vector<mass, 3> h_PQSO(mass r, Vector<mass,3> const& Delta);
Vector<mass, 3> h_P15QSO(mass m, mass r, Vector<mass,3> const& Delta);
Vector<mass, 3> h_P2QSS(mass r, Vector<mass,3> const& Spin1, Vector<mass,3> const& Spin2);