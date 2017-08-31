// Custom includes
#include <PDE.hpp>

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
    Mass = 0,
    Phase
};

int main()
{
    // Type aliases
    using mass_type = double;               // mass_needs to be precise
    using phase_type = std::complex<float>; // phase calculation can be less precise
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

    // Model switches (Compile time constants)
    constexpr bool use_eq1_contrib1 = true,
                   use_eq1_contrib2 = false,
                   use_eq2_contrib1 = false,
                   use_eq2_contrib2 = true;

    auto eq1_contrib1 = [=](const mass_type& mass)
    {
        return mass * eq1_contirb1_param1;
    };

    auto eq1_contrib2 = [some_derived_quantity = std::pow(eq1_contirb2_param1,
                                                          eq1_contirb2_param2) ](const mass_type& mass)
    {
        return mass * some_derived_quantity;
    };

    auto eq2_contrib1 = [=](const phase_type& phase)
    {
        return phase_type{ phase.real(), phase.imag() * eq2_contrib1_param1 };
    };

    auto eq2_contrib2 = [=](const phase_type& phase)
    {
        auto some_other_derived = std::cos(eq2_contrib2_param1) * phase.real();

        return some_other_derived * phase_type{ phase.real(), phase.imag() * eq2_contrib2_param2 };
    };

    // Assembling equations from contributions
    
    // It is expected, that sensible C++11/14 optimizers eliminate branching on ( true ? val1 : val2 )
    // type expressions at compile time. Even more sensible optimizers should change (val + 0) and
    // (val * 1) type expressions to nop (no-operation).
    auto eq1 = [=](const mass_type& mass)
    {
        return (use_eq1_contrib1 ? eq1_contrib1(mass) : (mass_type)0 ) +
               (use_eq1_contrib2 ? eq1_contrib2(mass) : (mass_type)0 );
    };

    auto eq2 = [=](const phase_type& phase)
    {
        return (use_eq2_contrib1 ? eq2_contrib1(phase) : (phase_type)1) *
               (use_eq2_contrib2 ? eq2_contrib2(phase) : (phase_type)1);
    };

    // Create default constructed solver. Allocates storage but is invalid state.
    solver rk4;

    // Set the initial left-hand side values.
    rk4.lhs() = state{ 1.0 ,                // real
                       { 0.0, 1.0 } };      // complex

    // Set the equation for each 
    rk4.equation() = [=](state& result, const state& rhs)
    {
        result = PDE::make_equation( eq1(rhs.get<Mass>()),
                                     eq2(rhs.get<Phase>()) );
    };

    for (solver_internal t = 0; t < 10; t += dt)
    {
        rk4.iterate(dt);
        
        std::cout << t << "\t" << rk4.lhs().get<Mass>() << "\t" << rk4.lhs().get<Phase>() << std::endl;
    }

    return 0;
}
