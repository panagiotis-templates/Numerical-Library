#pragma once
#include <iostream>
#include <cmath>
#include<cassert>
#include<type_traits>
#include<utility>
#include<functional>
template<typename _Ty>

inline constexpr bool is_decimal_v = std::disjunction_v<std::is_same<_Ty, float>, std::is_same<_Ty, double>, std::is_same<_Ty, long double>>;

template<typename _Ty>
struct is_decimal :std::bool_constant<is_decimal_v<_Ty>>{};//tag dispatching must support it

template<typename _Ty>
requires(is_decimal_v<_Ty>)
[[nodiscard]] bool inline isEqual(const _Ty& a, const _Ty& b, const _Ty& epsilon = static_cast<_Ty>(10e-10))noexcept {

    return std::abs(a - b) < epsilon;
}

template<typename _Ty, typename u>
    requires(is_decimal_v<_Ty>)
[[nodiscard]] std::optional<_Ty> simpson(const _Ty& a, const  _Ty& b, const  _Ty& dx, u&& f) {
    static_assert(is_decimal_v<std::invoke_result_t<decltype(f), _Ty>>, "return type of f must be a floating point type");
    static_assert(std::is_invocable_r_v<_Ty, u, _Ty>, "4rd argument must be a callable that returns a floating point value and takes only one floating point value");
    if (b < a||isEqual<_Ty>(b,a) || dx < 0||isEqual<_Ty>(dx,0.0)) {
        std::cerr << "b>a &&dx>0" << '\n';
        return std::nullopt;
    }

    const size_t& iterations = static_cast<size_t>((b - a) / dx);
    _Ty result = 0;
    _Ty xi = a;

    for (size_t i = 0; i < iterations; i++)
    {
        result += (dx / 6) * (std::invoke(f, xi) + 4 * std::invoke(f, (xi + (xi + dx)) * 0.5) + std::invoke(f, xi + dx));
        std::cout << xi << "-" << xi + dx << " " << result << '\n';
        xi += dx;
    }

    return std::optional{ result };
}


/*#pragma once
#include<iostream>
#include<cmath>
#include<type_traits>
template<typename _Ty>
requires(std::disjunction_v<std::is_same<_Ty,float>,std::is_same<_Ty,double>,std::is_same<_Ty,long double>>)
[[nodiscard]]_Ty f(_Ty x)noexcept
{
    //return pow(2,x) + 5 * x +3;
    return std::exp(2 * x);
    //return sin(2 * M_PI *x); 
}

template<typename _Ty>
requires(std::disjunction_v<std::is_same<_Ty, float>, std::is_same<_Ty, double>, std::is_same<_Ty, long double>>)
[[nodiscard]]_Ty simpson(_Ty a, _Ty b, _Ty dx) {
    size_t  iterations =static_cast<size_t>( (b - a) / dx);
    _Ty result = 0;
    _Ty xi = a;

    for (size_t i = 0; i <iterations; i++)
    {
        result += (dx / 6) * (f(xi) + 4 * f((xi + (xi + dx)) * 0.5) + f(xi + dx));
        std::cout << xi << "-" << xi + dx << " " << result << '\n';
        xi += dx;
    }
   
    return result;
}
*/
