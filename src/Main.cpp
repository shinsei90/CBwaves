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
    
    Vector<mass, 3> nullVector = {0., 0., 0.};
    state nullState = state{nullVector, nullVector};
    
    std::ofstream myfile;
    myfile.open("debug.dat");

    // mass risco = 44310;
    initParams iparams(2., 2., 20., 0.);
    mass rmin = (6.* SI_G * iparams.m)/SI_c2;
    double printstep = 1000.;
    // double T = 2.*PI*iparams.r0/(std::sqrt(SI_G*(iparams.m1 + iparams.m2)/iparams.r0));
    solver_internal dt = 1./std::pow(2., 11);
    // solver_internal dt = T/printstep;

    // Model switches (Compile time constants)
    // constexpr bool use_c_Newtonian = true,
    //                use_c_PostNewtonian = true,
    //                use_c_2PostNewtonian = false, 
    //                use_c_3PostNewtonian = false, 
    //                use_c_4PostNewtonian = false, 
    //                use_c_SpinOrbit = false,
    //                use_c_SpinSpin = false,
    //                use_c_BT_RR = false, 
    //                use_c_PostNewtonianSO = false, 
    //                use_c_2PostNewtonianSO = false, 
    //                use_c_RR1PostNewtonian = false,
    //                use_c_RRSO = false, 
    //                use_c_RRSS = false; 
    
    // constexpr bool use_h_Q = false,
    //                use_h_P05Q = false,
    //                use_h_PQ = false,
    //                use_h_P15Q = false,
    //                use_h_P2Q = false,
    //                use_h_PQSO = false,
    //                use_h_P15QSO = false,
    //                use_h_P2QSS = false;

    // constexpr bool use_e_Newtonian = false,
    //                use_e_PostNewtonian = false,
    //                use_e_2PostNewtonian = false,
    //                use_e_3PostNewtonian = false,
    //                use_e_4PostNewtonian = false,
    //                use_e_SpinOrbit = false,
    //                use_e_SpinSpin = false,
    //                use_e_PostNewtonianSO = false;

    // constexpr bool use_edot_Newtonian = false,
    //                use_edot_PostNewtonian = false,
    //                use_edot_2PostNewtonian = false,
    //                use_edot_25PostNewtonian = false,
    //                use_edot_SpinOrbit = false,
    //                use_edot_SpinSpin = false,
    //                use_edot_PostNewtonianSO = false;

    // constexpr bool use_l_PostNewtonian = false,
    //                use_l_2PostNewtonian = false,
    //                use_l_3PostNewtonian = false,
    //                use_l_4PostNewtonian = false,
    //                use_l_SpinOrbit = false;

    // constexpr bool use_spin1 = false;
    // constexpr bool use_spin2 = false;


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
    
    // It is expected, that sensible C++11/14 optimizers eliminate branching on ( false ? val1 : val2 )
    // type expressions at compile time. Even more sensible optimizers should change (val + 0) and
    // (val * 1) type expressions to nop (no-operation).

    auto corrs = [&](dynamicalParams const& dp) -> state { // capture clause could be reference
        return c_Newtonian(dp, iparams); //+ c_PostNewtonian(dp, iparams);
        // return (use_c_Newtonian ? c_Newtonian(dp, iparams) : nullState) +              
        //        (use_c_PostNewtonian ? c_PostNewtonian(dp, iparams) : nullState ) + 
        //        (use_c_2PostNewtonian ? c_2PostNewtonian(dp, iparams) : nullState ) +
        //        (use_c_3PostNewtonian ? c_3PostNewtonian(dp, iparams) : nullState ) +
        //        (use_c_4PostNewtonian ? c_4PostNewtonian(dp, iparams) : nullState ) +
        //        (use_c_SpinOrbit ? c_SpinOrbit(dp) : nullState ) +
        //        (use_c_SpinSpin ? c_SpinSpin(dp, iparams) : nullState ) +
        //        (use_c_BT_RR ? c_BT_RR(dp, iparams) : nullState ) +
        //        (use_c_PostNewtonianSO ? c_PostNewtonianSO(dp, iparams) : nullState ) +
        //        (use_c_2PostNewtonianSO ? c_2PostNewtonianSO(dp, iparams) : nullState ) +
        //        (use_c_RR1PostNewtonian ? c_RR1PostNewtonian(dp, iparams) : nullState ) +
        //        (use_c_RRSO ? c_RRSO(dp, iparams) : nullState ) +
        //        (use_c_RRSS ? c_RRSS(dp, iparams) : nullState );
    };

    // auto hterms = [=](dynamicalParams const& dp){
    //     return (use_h_Q ? h_Q(dp, iparams) : (mass)0. ) +
    //            (use_h_P05Q ? h_P05Q(dp, iparams) : (mass)0.) +
    //            (use_h_PQ ? h_PQ(dp, iparams) : (mass)0.) +
    //            (use_h_P15Q ? h_P15Q(dp, iparams) : (mass)0.) +
    //            (use_h_P2Q ? h_P2Q(dp, iparams) : (mass)0. ) +
    //            (use_h_PQSO ? h_PQSO(dp, iparams) : (mass)0. ) +
    //            (use_h_P15QSO ? h_P15QSO(dp, iparams) : (mass)0. ) +
    //            (use_h_P2QSS ? h_P2QSS(dp, iparams) : (mass)0. );
    // };
/*
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
        return (use_l_PostNewtonian ? l_PostNewtonian(dp, iparams) : nullVector ) +
               (use_l_2PostNewtonian ? l_2PostNewtonian(dp, iparams) : nullVector ) +
               (use_l_3PostNewtonian ? l_3PostNewtonian(dp, iparams) : nullVector ) +
               (use_l_4PostNewtonian ? l_4PostNewtonian(dp, iparams) : nullVector ) +
               (use_l_SpinOrbit ? l_SpinOrbit(dp, iparams) : nullVector );
    };

    auto S1 = [=](dynamicalParams const& dp){
        return (use_spin1 ? spin1(dp, iparams) : nullVector );
    };
    auto S2 = [=](dynamicalParams const& dp){
        return (use_spin2 ? spin2(dp, iparams) : nullVector );
    };
*/
    // Create default constructed solver. Allocates storage but is invalid state.
    solver rk4;

    // Set the initial left-hand side values.
    rk4.lhs() = state{ iparams.v_init , iparams.r_init };

    // Set the equation for each 
    rk4.equation() = [=](state& result, const state& rhs){

        dynamicalParams dp(rhs, iparams);
        auto temp = corrs(dp);

        result = PDE::make_equation( temp.get<Velocity>(), temp.get<Radius>() );
                                    //  hterms(dparams),
                                    //  eterms(dparams),
                                    //  lterms(dparams),
                                    //  edot(dparams) );

    };

    double step = 0;

    for (solver_internal t = 0; length(rk4.lhs().get<Radius>()) > rmin; t += dt){
        
        dynamicalParams dp(state{rk4.lhs().get<Velocity>(), rk4.lhs().get<Radius>()} ,iparams);
        rk4.iterate(dt);

        if(step++ == printstep){
            step = 0;
            //std::cout << t << "\t" << rk4.lhs().get<Mass>() << "\t" << rk4.lhs().get<Phase>() << std::endl;
            myfile << t << "\t" << dp.r << "\t" << dp.r1[0] << "\t" <<  dp.r1[1] << "\t" << dp.r1[2] << "\t" << dp.r2[0] << "\t" 
               <<  dp.r2[1] << "\t" << dp.r2[2] << "\t" << dp.v[0] << "\t" << dp.v[1] << "\t" << dp.v[2] << "\n" << std::endl;
        }
    }

    return 0;
}
