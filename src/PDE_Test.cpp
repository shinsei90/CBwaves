#include <PDE.hpp>
#include <fstream>

enum LV
{
    Rabbits = 0,
    Wolves = 1
};

int main()
{
    // Using directives
    using real = double;
    using solver_internal = real;
    using state = PDE::StateVector<real, real>;
    using solver = PDE::RK4::Solver<solver_internal, state>;

    // Simulation params
    real alpha = 0.2,
         beta = 0.6,
         gamma = 0.3,
         delta = beta;
    solver_internal dt = 0.05,
                    max_t = 200.0;

    solver rk4;

    rk4.lhs() = state{0.3, 0.5};

    rk4.equation() = [=](state& result, const state& rhs)
    {
        result = PDE::make_equation(alpha * rhs.get<Rabbits>() - beta * rhs.get<Rabbits>() * rhs.get<Wolves>(),
                                    delta * rhs.get<Rabbits>() * rhs.get<Wolves>() - gamma * rhs.get<Wolves>());
    };

    std::ofstream out("lv.dat");

    for(solver_internal t = 0.0 ; t < max_t ; t += dt)
    {
        out << t << '\t' << rk4.lhs().get<Rabbits>() << '\t' << rk4.lhs().get<Wolves>() << std::endl;

        rk4.iterate(dt);
    }

    return 0;
}