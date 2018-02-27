// Custom includes
#include <PDE.hpp>
#include <DynamicalParams.hpp>
#include <Corrections.hpp>
#include <Config.hpp>

// Standard C++ includes
#include <complex>
#include <iostream>
#include <fstream>

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

int main(){
    // Type aliases
    using solver = PDE::RK4::Solver<solver_internal, state>;
    
    // Model params (read from config file if you want)
    // mass_type eq1_contirb1_param1 = 0.99,
    //           eq1_contirb2_param1 = 1.01,
    //           eq1_contirb2_param2 = 1.02;
    // phase_type::value_type eq2_contrib1_param1 = 0.5f,
    //                        eq2_contrib2_param1 = 0.1f,
    //                        eq2_contrib2_param2 = 0.001f;
    
    std::ofstream myfile;
    myfile.open("debug.dat");

    mass rmin = 6.;
    solver_internal dt = 0.1;
    initParams iparams(10., 1.2, 100., 0.99);

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
                   use_h_PQSO = true,
                   use_h_P15QSO = true,
                   use_h_P2QSS = true;

    constexpr bool use_e_Newtonian = true,
                   use_e_PostNewtonian = true,
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

    constexpr bool use_spin1 = true;
    constexpr bool use_spin2 = true;


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

    auto corrs = [&](dynamicalParams const& dp){ // capture clause could be reference
        return (use_c_Newtonian ? c_Newtonian(dp, iparams) : Vector<mass, 3>{{0.,0.,0.}}) +              
               (use_c_PostNewtonian ? c_PostNewtonian(dp, iparams) : Vector<mass, 3>{{0.,0.,0.}} ) + 
               (use_c_2PostNewtonian ? c_2PostNewtonian(dp, iparams) : Vector<mass, 3>{{0.,0.,0.}} ) +
               (use_c_3PostNewtonian ? c_3PostNewtonian(dp, iparams) : Vector<mass, 3>{{0.,0.,0.}} ) +
               (use_c_4PostNewtonian ? c_4PostNewtonian(dp, iparams) : Vector<mass, 3>{{0.,0.,0.}} ) +
               (use_c_SpinOrbit ? c_SpinOrbit(dp, iparams) : Vector<mass, 3>{{0.,0.,0.}} ) +
               (use_c_SpinSpin ? c_SpinSpin(dp, iparams) : Vector<mass, 3>{{0.,0.,0.}} ) +
               (use_c_BT_RR ? c_BT_RR(dp, iparams) : Vector<mass, 3>{{0.,0.,0.}} ) +
               (use_c_PostNewtonianSO ? c_PostNewtonianSO(dp, iparams) : Vector<mass, 3>{{0.,0.,0.}} ) +
               (use_c_2PostNewtonianSO ? c_2PostNewtonianSO(dp, iparams) : Vector<mass, 3>{{0.,0.,0.}} ) +
               (use_c_RR1PostNewtonian ? c_RR1PostNewtonian(dp, iparams) : Vector<mass, 3>{{0.,0.,0.}} ) +
               (use_c_RRSO ? c_RRSO(dp, iparams) : Vector<mass, 3>{{0.,0.,0.}} ) +
               (use_c_RRSS ? c_RRSS(dp, iparams) : Vector<mass, 3>{{0.,0.,0.}} );
    };

    auto hterms = [=](dynamicalParams const& dp){
        return (use_h_Q ? h_Q(dp, iparams) : (mass)0. ) +
               (use_h_P05Q ? h_P05Q(dp, iparams) : (mass)0.) +
               (use_h_PQ ? h_PQ(dp, iparams) : (mass)0.) +
               (use_h_P15Q ? h_P15Q(dp, iparams) : (mass)0.) +
               (use_h_P2Q ? h_P2Q(dp, iparams) : (mass)0. ) +
               (use_h_PQSO ? h_PQSO(dp, iparams) : (mass)0. ) +
               (use_h_P15QSO ? h_P15QSO(dp, iparams) : (mass)0. ) +
               (use_h_P2QSS ? h_P2QSS(dp, iparams) : (mass)0. );
    };

    auto eterms = [=](dynamicalParams const& dp){
        return (use_e_Newtonian ? e_Newtonian(dp, iparams) : (mass)0. ) +
               (use_e_PostNewtonian ? e_PostNewtonian(dp, iparams) : (mass)0. ) +
               (use_e_2PostNewtonian ? e_2PostNewtonian(dp, iparams) : (mass)0. ) +
               (use_e_3PostNewtonian ? e_3PostNewtonian(dp, iparams) : (mass)0. ) +
               (use_e_4PostNewtonian ? e_4PostNewtonian(dp, iparams) : (mass)0. ) +
               (use_e_SpinOrbit ? e_SpinOrbit(dp, iparams) : (mass)0. ) +
               (use_e_SpinSpin ? e_SpinSpin(dp, iparams) : (mass)0. ) +
               (use_e_PostNewtonianSO ? e_PostNewtonianSO(dp, iparams) : (mass)0.);
    };

    auto edot = [=](dynamicalParams const& dp){
        return (use_edot_Newtonian ? edot_Newtonian(dp, iparams) : (mass)0.) +
               (use_edot_PostNewtonian ? edot_PostNewtonian(dp, iparams) : (mass)0. ) +
               (use_edot_2PostNewtonian ? edot_PostNewtonian(dp, iparams) : (mass)0. ) +
               (use_edot_25PostNewtonian ? edot_25PostNewtonian(dp, iparams) : (mass)0. ) +
               (use_edot_SpinOrbit ? edot_SpinOrbit(dp, iparams) : (mass)0. ) +
               (use_edot_SpinSpin ? edot_SpinSpin(dp, iparams) : (mass)0. ) +
               (use_edot_PostNewtonianSO ? edot_PostNewtonianSO(dp, iparams) : (mass)0. );
    };

    auto lterms = [=](dynamicalParams const& dp){
        return (use_l_PostNewtonian ? l_PostNewtonian(dp, iparams) : Vector<mass, 3>{{0.,0.,0.}} ) +
               (use_l_2PostNewtonian ? l_2PostNewtonian(dp, iparams) : Vector<mass, 3>{{0.,0.,0.}} ) +
               (use_l_3PostNewtonian ? l_3PostNewtonian(dp, iparams) : Vector<mass, 3>{{0.,0.,0.}} ) +
               (use_l_4PostNewtonian ? l_4PostNewtonian(dp, iparams) : Vector<mass, 3>{{0.,0.,0.}} ) +
               (use_l_SpinOrbit ? l_SpinOrbit(dp, iparams) : Vector<mass, 3>{{0.,0.,0.}} );
    };

    auto S1 = [=](dynamicalParams const& dp){
        return (use_spin1 ? spin1(dp, iparams) : Vector<mass, 3>{{0.,0.,0.}} );
    };
    auto S2 = [=](dynamicalParams const& dp){
        return (use_spin2 ? spin2(dp, iparams) : Vector<mass, 3>{{0.,0.,0.}} );
    };

    // Create default constructed solver. Allocates storage but is invalid state.
    solver rk4;

    // Set the initial left-hand side values.
    rk4.lhs() = state{ iparams.r_init };

    // Set the equation for each 
    rk4.equation() = [=](state& result, const state& rhs){

        //dynamicalParams dparams(rhs, iparams);
        dynamicalParams dp(rhs, iparams);

        result = PDE::make_equation( corrs(dp) );
                                    //  hterms(dparams),
                                    //  eterms(dparams),
                                    //  lterms(dparams),
                                    //  edot(dparams) );

    };

    for (solver_internal t = 0; length(rk4.lhs().get<Radius>()) > rmin; t += dt){
        
        dynamicalParams dp(rk4.lhs().get<Radius>() ,iparams);
        rk4.iterate(dt);
        
        //std::cout << t << "\t" << rk4.lhs().get<Mass>() << "\t" << rk4.lhs().get<Phase>() << std::endl;
        myfile << t << "\t" << dp.r1[1] << "\t" <<  dp.r1[2] << "\t" << dp.r1[3] << "\t" << dp.r2[1] << "\t" <<  dp.r2[2] << "\t" << dp.r2[3] <<"\n";

    }

    return 0;
}
