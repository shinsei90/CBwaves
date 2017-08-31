//////////////////////////////////////////////////////////////////////
//                                                                  //
// Unexported templated PDE solver to be used by client codebase    //
//                                                                  //
// Author: Nagy-Egri Máté Ferenc                                    //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#pragma once


// Standard C++ includes
#include <cstddef>
#include <cassert>
#include <tuple>
#include <array>
#include <functional>
#include <type_traits>      // std::result_of

namespace PDE
{
    /// <summary>Traits class consisting of type aliases describing heterogenous state vectors.</summary>
    ///
    template <typename... VT>
    struct StateVectorTraits
    {
        // STL type aliases

        using container_type = std::tuple<VT...>;
        using size_type = std::size_t;
    };


    /// <summary>Read-only Expression Template base class of heterogenous state vectors to be implemented statically.</summary>
    ///
    template <typename E, typename... VT>
    struct ConstExpression : public StateVectorTraits<VT...>
    {
        // STL type aliases

        using container_type = typename StateVectorTraits<VT...>::container_type;
        using size_type = typename StateVectorTraits<VT...>::size_type;

        // ConstStateVector interface

        /// <summary>Obtain the N-th element of the state vector.</summary>
        ///
        template <int N> const auto get() const { return static_cast<const E&>(*this).template get<N>(); }

        /// <summary>Obtan the number of elements in the state vector.</summary>
        ///
        constexpr size_type size() const { return sizeof...(VT); }

        // ConstExpression interface

        operator const E&() const { return static_cast<const E&>(*this); }
    };


    /// <summary>Read-write Expression Template base class of heterogenous state vectors to be implemented statically.</summary>
    ///
    template <typename E, typename... VT>
    struct Expression : public ConstExpression<E, VT...>
    {
        // STL type aliases

        using container_type = typename StateVectorTraits<VT...>::container_type;
        using size_type = typename StateVectorTraits<VT...>::size_type;

        // StateVector interface

        template <int N> auto& get() { return static_cast<E&>(*this).template get<N>(); }

        // Expression interface

        operator E&() { return static_cast<E&>(*this); }
    };


    /// <summary>Namespace that contains implementation details of StateVector operations.</summary>
    ///
    namespace impl
    {
        /// <summary>Recursive class template implementing the assignment operator.</summary>
        ///
        template <int N>
        struct Assign
        {
            template <typename ET, typename... T>
            static void assign(std::tuple<T...>& dst, const ET& src)
            {
                //std::cout << "pde::assign<" << N << "> dst = " << std::get<N>(dst).extent() << " src = " << src.template get<N>().extent() << std::endl;
                std::get<N>(dst) = src.template get<N>();
                Assign<N - 1>::assign(dst, src);
            }
        };


        /// <summary>Closing stage of the recursive class template implementing the assignment operator.</summary>
        ///
        template <>
        struct Assign<0>
        {
            template <typename ET, typename... T>
            static void assign(std::tuple<T...>& dst, const ET& src)
            {
                //std::cout << "pde::assign<0> dst = " << std::get<0>(dst).extent() << " src = " << src.template get<0>().extent() << std::endl;
                std::get<0>(dst) = src.template get<0>();
            }
        };
    }


    /// <summary>Class storing the states of a heterogenous vector of states.</summary>
    /// <remarks>The StateVector class template is not an Expression Template. There are seperate View classes for that.</remarks>
    ///
    template <typename... VT>
    class StateVector : public StateVectorTraits<VT...>
    {
    public:

        // STL type aliases

        using container_type = typename StateVectorTraits<VT...>::container_type;
        using size_type = typename StateVectorTraits<VT...>::size_type;

        // Constructors / Destructors / Assignment operators

        /// <summary>Default constructor.</summary>
        /// <remarks>Default constructed objects are in an invalid state.</remarks>
        ///
        StateVector() = default;

        /// <summary>Default copy constructor.</summary>
        ///
        StateVector(const StateVector& in) = default;

        /// <summary>Default move constructor.</summary>
        ///
        StateVector(StateVector&& in) = default;

        /// <summary>Default destructor.</summary>
        ///
        ~StateVector() = default;

        /// <summary>Default copy assignment operator.</summary>
        ///
        StateVector& operator=(const StateVector&) = default;

        /// <summary>Default move assignment operator.</summary>
        ///
        StateVector& operator=(StateVector&&) = default;

        ///<summary>Constructs a state vector from the expression <c>expr</c>.</summary>
        ///
        template <typename StateVecExpr, typename... ValueType>
        StateVector(const ConstExpression<StateVecExpr, ValueType...>& expr)
        {
            serial_evaluator(expr);
        }

        ///<summary>Constructs a state vector from the expression <c>expr</c>.</summary>
        ///
        template <typename StateVecExpr, typename... ValueType>
        StateVector& operator=(const ConstExpression<StateVecExpr, ValueType...>& expr)
        {
            serial_evaluator(expr);

            return *this;
        }

        ///<summary>Copy constructs a state vector from the states <c>values</c>.</summary>
        ///
        StateVector(const VT&... values) : m_values(values...) {}

        ///<summary>Move constructs a state vector from the states <c>values</c>.</summary>
        ///
        StateVector(VT&&... values) : m_values(values...) {}

        // ConstStateVector interface

        /// <summary>Obtain the N-th element of the state vector.</summary>
        ///
        template <int N> auto& get() const { return std::get<N>(m_values); }

        /// <summary>Obtan the number of elements in the state vector.</summary>
        ///
        constexpr size_type size() const { return std::tuple_size<VT...>::value; }

        // StateVector interface

        /// <summary>Obtain the N-th element of the state vector.</summary>
        ///
        template <int N> auto& get() { return std::get<N>(m_values); }

    private:

        /// <summary>Initializes internal states and evaluates elements of the expression <c>expr</c> in a serial manner.</summary>
        ///
        template <typename StateVecExpr, typename... ValueType>
        void serial_evaluator(const ConstExpression<StateVecExpr, ValueType...>& expr)
        {
            // Extract type from encapsulating expression
            const StateVecExpr& v = expr;

            impl::Assign<(sizeof...(ValueType)) - 1>::assign(m_values, v);
        }

        container_type m_values;
    };


    /// <summary>Expression Template providing read-only access with reference semantics to the elements of a state vector with storage.</summary>
    ///
    template <typename... VT>
    class ConstView : public ConstExpression<ConstView<VT...>, VT...>
    {
    public:

        // Expression type aliases

        using expression_type = ConstExpression<ConstView<VT...>, VT...>;

        // STL type aliases

        using container_type = typename expression_type::container_type;
        using size_type = typename expression_type::size_type;

    private:

        using vector_type = StateVector<VT...>;
        using reference_type = std::reference_wrapper<const vector_type>;

    public:

        // Constructors / Destructors / Assignment operators

        /// <summary>Default constructor.</summary>
        /// <remarks>Default constructor deleted, as underlying reference type cannot be default initialized.</remarks>
        ///
        ConstView() = delete;

        /// <summary>Default copy constructor.</summary>
        ///
        ConstView(const ConstView&) = default;

        /// <summary>Default move constructor.</summary>
        /// <remarks>Default move constructor deleted, as underlying reference type cannot be moved.</remarks>
        ///
        ConstView(ConstView&&) = delete;

        /// <summary>Default destructor.</summary>
        ///
        ~ConstView() = default;

        /// <summary>Default copy assign operator.</summary>
        ///
        ConstView& operator=(const ConstView&) = default;

        /// <summary>Default move assign operator.</summary>
        /// <remarks>Default move assign operator deleted, as underlying reference type cannot be moved.</remarks>
        ///
        ConstView& operator=(ConstView&&) = delete;

        /// <summary>Constructs a <c>ConstView</c> from a <c>StateVector</c>.</summary>
        ///
        ConstView(const vector_type& v) : _v(std::cref(v)) {}

        // ConstStateVector interface

        /// <summary>Obtain the N-th element of the state vector.</summary>
        ///
        template <int N> const auto& get() const { return _v.get().template get<N>(); }

        /// <summary>Obtan the number of elements in the state vector.</summary>
        ///
        constexpr size_type size() const { return std::tuple_size<VT...>::value; }

    private:

        const reference_type _v;
    };


    /// <summary>Expression Template providing read-write access with reference semantics to the elements of a state vector with storage.</summary>
    ///
    template <typename... VT>
    class View : public Expression<View<VT...>, VT...>
    {
    public:

        // Expression type aliases

        using expression_type = ConstView<VT...>;

        // STL type aliases

        using container_type = typename expression_type::container_type;
        using size_type = typename expression_type::size_type;

    private:

        using vector_type = StateVector<VT...>;
        using reference_type = std::reference_wrapper<vector_type>;

    public:

        // Constructors / Destructors / Assignment operators

        /// <summary>Default constructor.</summary>
        /// <remarks>Default constructor deleted, as underlying reference type cannot be default initialized.</remarks>
        ///
        View() = delete;

        /// <summary>Default copy constructor.</summary>
        ///
        View(const View&) = default;

        /// <summary>Default move constructor.</summary>
        /// <remarks>Default move constructor deleted, as underlying reference type cannot be moved.</remarks>
        ///
        View(View&&) = delete;

        /// <summary>Default destructor.</summary>
        ///
        ~View() = default;

        /// <summary>Default copy assign operator.</summary>
        ///
        View& operator=(const View&) = default;

        /// <summary>Default move assign operator.</summary>
        /// <remarks>Default move assign operator deleted, as underlying reference type cannot be moved.</remarks>
        ///
        View& operator=(View&&) = delete;

        /// <summary>Constructs a <c>ConstView</c> from a <c>StateVector</c>.</summary>
        ///
        View(vector_type& v) : _v(std::ref(v)) {}

        ///<summary>Constructs a state vector from the expression <c>expr</c>.</summary>
        ///
        template <typename StateVecExpr, typename... ValueType>
        View(const ConstExpression<StateVecExpr, ValueType...>& expr)
        {
            // Extract type from encapsulating expression
            const StateVecExpr& v = expr;

            _v.get() = v;
        }

        ///<summary>Constructs a state vector from the expression <c>expr</c>.</summary>
        ///
        template <typename StateVecExpr, typename... ValueType>
        View& operator=(const ConstExpression<StateVecExpr, ValueType...>& expr)
        {
            // Extract type from encapsulating expression
            const StateVecExpr& v = expr;

            _v.get() = v;

            return *this;
        }

        // ConstStateVector interface

        /// <summary>Obtain the N-th element of the state vector.</summary>
        ///
        template <int N> auto get() const { return _v.get().template get<N>(); }

        /// <summary>Obtan the number of elements in the state vector.</summary>
        ///
        constexpr size_type size() const { return std::tuple_size<VT...>::value; }

        // StateVector interface

        /// <summary>Obtain the N-th element of the state vector.</summary>
        ///
        template <int N> auto& get() { return _v.get().template get<N>(); }

    private:

        reference_type _v;
    };


    /// <summary>Expression Template that evaluates function objects to yield elements of a state vector.</summary>
    ///
    template <typename... P>
    class Proxy : public ConstExpression<Proxy<P...>, P...>
    {
    public:

        // Expression type aliases

        using expression_type = ConstExpression<Proxy<P...>, P...>;

        // STL type aliases

        using container_type = typename expression_type::container_type;
        using size_type = typename expression_type::size_type;

        // Constructors / Destructors / Assignment operators

        Proxy(P... values) : _p(values...) {}

        // ConstStateVector interface

        /// <summary>Obtain the N-th element of the state vector.</summary>
        ///
        template <int N> auto get() const { return std::get<N>(_p); }

        /// <summary>Obtan the number of elements in the state vector.</summary>
        ///
        constexpr size_type size() const { return sizeof...(P); }

    private:

        std::tuple<P...> _p;
    };


    /// <summary>Expression Template that applies a unary function object to every element of a state vector.</summary>
    ///
    template <typename E, typename F, typename... VT>
    class Map : public ConstExpression<Map<E, F, VT...>, VT...>
    {
    public:

        // Expression type aliases

        using expression_type = ConstExpression<Map<E, F, VT...>, VT...>;
        using outer_expression_type = E;

        // STL type aliases

        using container_type = typename expression_type::container_type;
        using size_type = typename expression_type::size_type;

        /// <summary>Constructs a <c>Map</c> from an expression <c>u</c> and a function object <c>f</c>.</summary>
        ///
        Map(const ConstExpression<E, VT...>& u, const F f)
            : _u(static_cast<const E&>(u))
            , _f(f)
        {}

        // ConstStateVector interface

        /// <summary>Obtain the N-th element of the state vector.</summary>
        ///
        template <int N> auto get() const { return _f(_u.template get<N>()); }

        /// <summary>Obtan the number of elements in the state vector.</summary>
        ///
        constexpr size_type size() const { return sizeof...(VT); }

    private:

        const E _u;
        const F _f;
    };


    /// <summary>Expression Template that applies a binary function object to every element of two state vectors to yield the result.</summary>
    ///
    template <typename E1, typename E2, typename F, typename... VT>
    class Zip : public ConstExpression<Zip<E1, E2, F, VT...>, VT...>
    {
    public:

        // Expression type aliases

        using expression_type = ConstExpression<Zip<E1, E2, F, VT...>, VT...>;
        using outer_expression_type1 = E1;
        using outer_expression_type2 = E2;

        // STL type aliases

        using container_type = typename expression_type::container_type;
        using size_type = typename expression_type::size_type;

        // Constructors / Destructors / Assignment operators

        /// <summary>Constructs a <c>Zip</c> from two possibly varying expressions <c>u</c> and <c>v</c> and a function object <c>f</c>.</summary>
        ///
        Zip(const ConstExpression<E1, VT...>& u, const ConstExpression<E2, VT...>& v, const F f)
            : _u(static_cast<const E1&>(u))
            , _v(static_cast<const E2&>(v))
            , _f(f)
        {}

        // ConstStateVector interface

        /// <summary>Obtain the N-th element of the state vector.</summary>
        ///
        template <int N> auto get() const { return _f(_u.template get<N>(), _v.template get<N>()); }

        /// <summary>Obtan the number of elements in the state vector.</summary>
        ///
        constexpr size_type size() const { return sizeof...(VT); }

    private:

        const E1 _u;
        const E2 _v;
        F _f;
    };


    namespace RK4
    {
        template <typename T, typename States>
        class Solver;

        template <typename T, typename... States>
        class Solver<T, StateVector<States...>>
        {
        public:

            // Solver type aliases

            using floating_type = T;
            using state_vector = StateVector<States...>;
            //using view_type = View<States...>;
            //using const_view_type = ConstView<States...>;
            //using equation_type = std::function<void(view_type, const const_view_type)>;
            using equation_type = std::function<void(state_vector&, const state_vector&)>;

            // Constructors / Destructors / Assignment operators

            /// <summary>Default constructor.</summary>
            /// <remarks>Default constructed objects are in an invalid state.</remarks>
            ///
            Solver() = default;

            /// <summary>Default copy constructor.</summary>
            ///
            Solver(const Solver& in) = default;

            /// <summary>Default move constructor.</summary>
            ///
            Solver(Solver&& in) = default;

            /// <summary>Default destructor.</summary>
            ///
            ~Solver() = default;

            /// <summary>Default copy assignment operator.</summary>
            ///
            Solver& operator=(const Solver&) = default;

            /// <summary>Default move assignment operator.</summary>
            ///
            Solver& operator=(Solver&&) = default;

            // Solver interface

            /// <summary>Obtain the left-hand side stored within the solver.</summary>
            ///
            state_vector& lhs() { return m_lhs; }

            /// <summary>Obtain the left-hand side stored within the solver.</summary>
            ///
            const state_vector& lhs() const { return m_lhs; }

            /// <summary>Obtain the equation of the right-hand side stored within the solver.</summary>
            ///
            equation_type& equation() { return m_func; }

            /// <summary>Obtain the equation of the right-hand side stored within the solver.</summary>
            ///
            const equation_type& equation() const { return m_func; }

            /// <summary>Iterate over the indepependent variable of size dt.</summary>
            ///
            floating_type iterate(const floating_type dt)
            {
                return nice_iteration(dt);
            }

        //private:
        public:

            /// <summary>Performs a RK4 step in human readable form.</summary>
            ///
            floating_type nice_iteration(const floating_type& dt)
            {
                m_func(m_k.at(k1), m_lhs);

                //std::cout << "eval k2" << std::endl;

                m_func(m_k.at(k2), tmp = m_lhs + (static_cast<floating_type>(0.5) * dt) * m_k.at(k1));

                //std::cout << "eval k3" << std::endl;

                m_func(m_k.at(k3), tmp = m_lhs + (static_cast<floating_type>(0.5) * dt) * m_k.at(k2));

                //std::cout << "eval k4" << std::endl;

                m_func(m_k.at(k4), tmp = m_lhs + (static_cast<floating_type>(1.0) * dt) * m_k.at(k3));

                //std::cout << "eval new lhs" << std::endl;

                m_lhs = m_lhs + (dt / static_cast<floating_type>(6.0))*(m_k.at(k1) + static_cast<floating_type>(2.0) * (m_k.at(k2) + m_k.at(k3)) + m_k.at(k4));

                return dt;
            }

            /// <summary>Performs a RK4 step in human non-readable form.</summary>
            ///
            floating_type ugly_iteration(const floating_type& dt)
            {
                /*
                m_k.at(k1) = m_func(m_lhs);

                tmp1 = static_cast<floating_type>(0.5) * dt * m_k.at(k1);
                tmp2 = m_lhs + tmp1;
                m_k.at(k2) = m_func(tmp2);

                tmp1 = static_cast<floating_type>(0.5) * dt * m_k.at(k2);
                tmp2 = m_lhs + tmp1;
                m_k.at(k3) = m_func(tmp2);

                tmp1 = static_cast<floating_type>(1.0) * dt * m_k.at(k3);
                tmp2 = m_lhs + tmp1;
                m_k.at(k4) = m_func(tmp2);

                tmp1 = static_cast<floating_type>(2.0) * m_k.at(k2);
                tmp2 = m_k.at(k1) + tmp1;
                tmp1 = static_cast<floating_type>(2.0) * m_k.at(k3);
                tmp3 = tmp2 + tmp1;
                tmp1 = tmp3 + m_k.at(k4);
                tmp2 = (dt / static_cast<floating_type>(6.0)) * tmp1;
                tmp3 = m_lhs + tmp2;

                m_lhs = tmp3;
                */
                return dt;
            }

            enum Step
            {
                k1 = 0,
                k2 = 1,
                k3 = 2,
                k4 = 3
            };

            state_vector m_lhs, tmp, tmp1, tmp2, tmp3;
            std::array<state_vector, 4> m_k;
            equation_type m_func;
        };

    } // namespace RK4

    template <typename E, typename F, typename... VT> auto map(const ConstExpression<E, VT...>& u, const F f) { return Map<E, F, VT...>(u, f); }

    template <typename E1, typename E2, typename F, typename... VT>
    auto zip(const ConstExpression<E1, VT...>& u,
             const ConstExpression<E2, VT...>& v,
             const F f)
    {
        return Zip<E1, E2, F, VT...>(u, v, f);
    }

    template <typename... P> auto make_equation(P... equations) { return Proxy<P...>(equations...); }
    
} // namespace PDE

///////////////////////////////////////////
// PDE::StateVector non-member operators //
///////////////////////////////////////////

// Unary
template <typename E, typename... T> auto operator+(const PDE::ConstExpression<E, T...>& v) { return PDE::map(v, [=](const auto& val) { return +val; }); }
template <typename E, typename... T> auto operator-(const PDE::ConstExpression<E, T...>& v) { return PDE::map(v, [=](const auto& val) { return -val; }); }
template <typename... T> auto operator+(const PDE::StateVector<T...>& v) { return PDE::map(PDE::ConstView<T...>(v), [=](const auto& val) { return +val; }); }
template <typename... T> auto operator-(const PDE::StateVector<T...>& v) { return PDE::map(PDE::ConstView<T...>(v), [=](const auto& val) { return -val; }); }

// Binary
template <typename E1, typename E2, typename... T> auto operator+(const PDE::ConstExpression<E1, T...>& u, const PDE::ConstExpression<E2, T...>& v) { return PDE::zip(u, v, [](const auto& lhs, const auto& rhs) { return lhs + rhs; }); }
template <typename E1, typename E2, typename... T> auto operator-(const PDE::ConstExpression<E1, T...>& u, const PDE::ConstExpression<E2, T...>& v) { return PDE::zip(u, v, [](const auto& lhs, const auto& rhs) { return lhs - rhs; }); }
template <typename E1, typename... T> auto operator+(const PDE::ConstExpression<E1, T...>& u, const PDE::StateVector<T...>& v) { return PDE::zip(u, PDE::ConstView<T...>(v), [](const auto& lhs, const auto& rhs) { return lhs + rhs; }); }
template <typename E1, typename... T> auto operator-(const PDE::ConstExpression<E1, T...>& u, const PDE::StateVector<T...>& v) { return PDE::zip(u, PDE::ConstView<T...>(v), [](const auto& lhs, const auto& rhs) { return lhs - rhs; }); }
template <typename E2, typename... T> auto operator+(const PDE::StateVector<T...>& u, const PDE::ConstExpression<E2, T...>& v) { return PDE::zip(PDE::ConstView<T...>(u), v, [](const auto& lhs, const auto& rhs) { return lhs + rhs; }); }
template <typename E2, typename... T> auto operator-(const PDE::StateVector<T...>& u, const PDE::ConstExpression<E2, T...>& v) { return PDE::zip(PDE::ConstView<T...>(u), v, [](const auto& lhs, const auto& rhs) { return lhs - rhs; }); }
template <typename... T> auto operator+(const PDE::StateVector<T...>& u, const PDE::StateVector<T...>& v) { return PDE::zip(PDE::ConstView<T...>(u), PDE::ConstView<T...>(v), [](const auto& lhs, const auto& rhs) { return lhs + rhs; }); }
template <typename... T> auto operator-(const PDE::StateVector<T...>& u, const PDE::StateVector<T...>& v) { return PDE::zip(PDE::ConstView<T...>(u), PDE::ConstView<T...>(v), [](const auto& lhs, const auto& rhs) { return lhs - rhs; }); }

template <typename Scalar, typename E, typename... T> auto operator*(const Scalar alpha, const PDE::ConstExpression<E, T...>& v) { return PDE::map(v, [=](const auto& val) { return alpha * val; }); }
template <typename Scalar, typename E, typename... T> auto operator*(const PDE::ConstExpression<E, T...>& v, const Scalar alpha) { return PDE::map(v, [=](const auto& val) { return alpha * val; }); }
template <typename Scalar, typename... T> auto operator*(const Scalar alpha, const PDE::StateVector<T...>& v) { return PDE::map(PDE::ConstView<T...>(v), [=](const auto& val) { return alpha * val; }); }
template <typename Scalar, typename... T> auto operator*(const PDE::StateVector<T...>& v, const Scalar alpha) { return PDE::map(PDE::ConstView<T...>(v), [=](const auto& val) { return val * alpha; }); }

template <typename Scalar, typename E, typename... T> auto operator/(const PDE::ConstExpression<E, T...>& v, const Scalar alpha) { return PDE::map(v, [=](const auto& val) { return val / alpha; }); }
template <typename Scalar, typename... T> auto operator/(const PDE::StateVector<T...>& v, const Scalar alpha) { return PDE::map(PDE::ConstView<T...>(v), [=](const auto& val) { return val / alpha; }); }
