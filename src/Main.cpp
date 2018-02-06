// Custom includes
#include <PDE.hpp>
#include <Corrections.hpp>
#include <Config.hpp>

// Standard C++ includes
#include <complex>
#include <iostream>

/// <note>If PDE::StateVector<...> has std::complex in it and
/// std::complex::value_type does not match PDE::RK4::Solver::floating_type,
/// that means that mixed type scalar * complex operations will manifest
/// inside the RK4 solver (resulting in criptic template error messages).
/// To account for such scenarios, the user has to provide the compiler an
/// operator it can use in such cases. (It is not evident which operand should
/// be promoted, however generally the scalar might be a good choice.</note>
///
template <typename S, typename C>
auto operator*(const S& s, const std::complex<C>& c)
{
    return static_cast<typename std::complex<C>::value_type>(s) * c;
}

enum Component : int
{
    Radius
};

int main(){
    // Type aliases
    using mass_type = double;               // mass_needs to be precise
    using phase_type = std::complex<double>; // phase calculation can be less precise
    using solver_internal = double;         // scalars used inside the solver

    using state = PDE::StateVector<mass_type, phase_type>;
    using solver = PDE::RK4::Solver<solver_internal, state>;
    
    // Model params (read from config file if you want)
    mass_type eq1_contirb1_param1 = 0.99,
              eq1_contirb2_param1 = 1.01,
              eq1_contirb2_param2 = 1.02;
    phase_type::value_type eq2_contrib1_param1 = 0.5f,
                           eq2_contrib2_param1 = 0.1f,
                           eq2_contrib2_param2 = 0.001f;
    solver_internal dt = 0.1;

    struct initParams{  //Constant parameters of the equations.

        initParams(mass m1_, mass m2_, mass r0_){
            //Mass ratio
        m1 = m1_;
        m2 = m2_;
        m = m1 + m2;

        //Initial separation
        r0 = r0_;
        r_init = {r0, 0.0, 0.0};


        dm = m1 - m2;
        mu = m1*m2/m;
        eta = mu/m;
        }
      
        //Mass ratio
        mass m1;
        mass m2;
        mass m;

        //Initial separation
        mass r0;
        Vector<mass, 3>  r_init;


        mass dm;
        mass mu;
        mass eta;

    };

    initParams iparams(10, 1.2, 100);

    // Model switches (Compile time constants)
    constexpr bool use_c_Newtonian = true,
                   use_c_PostNewtonian = true,
                   use_c_2PostNewtonian = true,
                   use_c_3PostNewtonian = true,
                   use_c_4PostNewtonian = true,
                   use_c_SpinOrbit = true,
                   use_c_SpinSpin = true,
                   use_c_BT_RR = true,
                   use_c_PostNewtonianSO = true,
                   use_c_2PostNewtonianSO = true,
                   use_c_RR1PostNewtonian = true,
                   use_c_RRSO = true,
                   use_c_RRSS = true;
    
    constexpr bool use_h_Q = true,
                   use_h_P05Q = true,
                   use_h_PQ = true,
                   use_h_P15Q = true,
                   use_h_P2Q = true,
                   use_h_P15QSO = true,
                   use_h_P2QSS = true;

    constexpr bool use_e_PostNewtonian = true,
                   use_e_2PostNewtonian = true,
                   use_e_3PostNewtonian = true,
                   use_e_4PostNewtonian = true,
                   use_e_SpinOrbit = true,
                   use_e_SpinSpin = true,
                   use_e_PostNewtonianSO = true;

    constexpr bool use_edot_Newtonian = true,
                   use_edot_PostNewtonian = true,
                   use_edot_2PostNewtonian = true,
                   use_edot_25PostNewtonian = true,
                   use_edot_SpinOrbit = true,
                   use_edot_SpinSpin = true,
                   use_edot_PostNewtonianSO = true;

    constexpr bool use_l_PostNewtonian = true,
                   use_l_2PostNewtonian = true,
                   use_l_3PostNewtonian = true,
                   use_l_4PostNewtonian = true,
                   use_l_SpinOrbit = true;


    //Example equations.

    // auto c_Newtonian = [=](const mass_type& mass)
    // {
    //     return mass * eq1_contirb1_param1;
    // };

    // auto c_PostNewtonian =ome_derivantity = std::pow(eq1_contirb2_param1,
    //                                                       eq1_contirb2_param2) ](const mass_type& mass)
    // {
    //     return mass * some_derived_quantity;
    // };

    // auto eq2_contrib1 = [=](const phase_type& phase)
    // {
    //     return phase_type{ phase.real(), phase.imag() * eq2_contrib1_param1 };
    // };

    // auto eq2_contrib2 = [=](const phase_type& phase)
    // {
    //     auto some_other_derived = std::cos(eq2_contrib2_param1) * phase.real();

    //     return some_other_derived * phase_type{ phase.real(), phase.imag() * eq2_contrib2_param2 };
    // };

    // Assembling equations from contributions
    
    // It is expected, that sensible C++11/14 optimizers eliminate branching on ( true ? val1 : val2 )
    // type expressions at compile time. Even more sensible optimizers should change (val + 0) and
    // (val * 1) type expressions to nop (no-operation).
    auto corrs = [=](const mass r){
        return (use_c_Newtonian ? c_Newtonian(m, r) : (mass)0 ) +
               (use_c_PostNewtonian ? c_PostNewtonian(m, r) : (mass)0 ) + 
               (use_c_2PostNewtonian ? c_2PostNewtonian(m, r) : (mass)0 ) +
               (use_c_3PostNewtonian ? c_3PostNewtonian(m, r) : (mass)0 ) +
               (use_c_4PostNewtonian ? c_4PostNewtonian(m, r) : (mass)0 ) +
               (use_c_SpinOrbit ? c_SpinOrbit(m, r, Spin, sigma) : (mass)0 ) +
               (use_c_SpinSpin ? c_SpinSpin(m, r, Spin1, Spin2) : (mass)0 ) +
               (use_c_BT_RR ? c_BT_RR(m, r) : (mass)0 ) +
               (use_c_PostNewtonianSO ? c_PostNewtonianSO(m, r, Spin, Delta) : (mass)0 ) +
               (use_c_2PostNewtonianSO ? c_2PostNewtonianSO(m, r, Spin, Delta) : (mass)0 ) +
               (use_c_RR1PostNewtonian ? c_RR1PostNewtonian(m, r) : (mass)0 ) +
               (use_c_RRSO ? c_RRSO(m, r, Spin, sigma, LN) : (mass)0) +
               (use_c_RRSS ? c_RRSS(m, r, Spin1, Spin2) : (mass)0);
    };

    auto hterms = [=](const mass m){
        return (use_h_Q ? h_Q(m) : (mass)0) +
               (use_h_P05Q ? h_P05Q(m) : (mass)0) +
               (use_h_PQ ? h_PQ(m) : (mass)0) +
               (use_h_P15Q ? h_P15Q(m) : (mass)0) +
               (use_h_P2Q ? h_P2Q(m) : (mass)0 ) +
               (use_h_P15QSO ? h_P15QSO(m) : (mass)0 ) +
               (use_h_P2QSS ? h_P2QSS(m) : (mass)0 );
    };

    auto eterms = [=](const mass m){
        return (use_e_Newtonian ? e_Newtonian(m) : (mass)0 ) +
               (use_e_PostNewtonian ? e_PostNewtonian(m) : (mass)0 ) +
               (use_e_2PostNewtonian ? e_2PostNewtonian(m) : (mass)0 ) +
               (use_e_3PostNewtonian ? e_3PostNewtonian(m) : (mass)0 ) +
               (use_e_4PostNewtonian ? e_4PostNewtonian(m) : (mass)0 ) +
               (use_e_SpinOrbit ? e_SpinOrbit(m) : (mass)0 ) +
               (use_e_SpinSpin ? e_SpinSpin(m) (mass)0 ) +
               (use_e_PostNewtonianSO ? e_PostNewtonianSO(m) : (mass)0);
    }

    auto edot = [=](const mass m){
        return (use_edot_Newtonian ? edot_Newtonian(m) : (mass)0) +
               (use_edot_PostNewtonian ? edot_PostNewtonian(m) : (mass)0 ) +
               (use_edot_2PostNewtonian ? edot_PostNewtonian(m) : (mass)0 ) +
               (use_edot_25PostNewtonian ? edot_25PostNewtonian(m) : (mass)0 ) +
               (use_edot_SpinOrbit ? edot_SpinOrbit(m) : (mass)0 ) +
               (use_edot_SpinSpin ? edot_SpinSpin(m) : (mass)0 ) +
               (use_edot_PostNewtonianSO ? edot_PostNewtonianSO(m) : (mass)0 );
    }

    auto lterms = [=](const mass m){
        return (use_l_PostNewtonian ? l_PostNewtonian(m) : (mass)0) +
               (use_l_2PostNewtonian ? l_2PostNewtonian(m) : (mass)0) +
               (use_l_3PostNewtonian ? l_3PostNewtonian(m) : (mass)0) +
               (use_l_4PostNewtonian ? l_4PostNewtonian(m) : (mass)0) +
               (use_l_SpinOrbit ? l_SpinOrbit(m) : (mass)0);
    }

    // Create default constructed solver. Allocates storage but is invalid state.
    solver rk4;

    // Set the initial left-hand side values.
    rk4.lhs() = state{ {r0, 0.0, 0.0} };

    // Set the equation for each 
    rk4.equation() = [=](state& result, const state& rhs)
    {
        //auto new_rho = rhs.get<Radius>() * 0.1;
        //auto separation = rhs.get<Radius>() - r_init;

        struct dynamicalParams {

            dynamicalParams(state& state_, struct initParams_){

            Vector<mass, 3> rr = rhs.get<Radius>() - r_init;
            Vector<mass, 3> x1 = {r/2, 0, 0};
            Vector<mass, 3> x2 = {-r/2, 0, 0};
            Vector<mass, 3> v = {0, SI_c/3, 0};

            Vector<mass, 3> Spin1 = spin1(m, r);
            Vector<mass, 3> Spin2 = spin2(m, r);

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

        }

        dynamicalParams dparams(rhs.get<Radius>(), iparams);

        result = PDE::make_equation( corrs(rhs.get<Mass>(), dparams, iparams),
                                     hterms(rhs.get<Mass>(), dparams, iparams),
                                     eterms(rhs.get<Mass>(), dparams, iparams),
                                     lterms(rhs.get<Mass>(), dparams, iparams),
                                     edot(rhs.get<Mass>()), dparams, iparams);
    };

    for (solver_internal t = 0; t < 10; t += dt)
    {
        rk4.iterate(dt);
        
        std::cout << t << "\t" << rk4.lhs().get<Mass>() << "\t" << rk4.lhs().get<Phase>() << std::endl;
    }

    return 0;
}
