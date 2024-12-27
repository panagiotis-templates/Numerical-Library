#pragma once
#include <iostream>
#include <cmath>
#include<cassert>
#include<type_traits>
#include<utility>
#include<functional>
template<typename _Ty>
struct is_decimal :std::disjunction<std::is_same<_Ty, float>, std::is_same<_Ty, double>, std::is_same<_Ty, long double>> {};

template<typename _Ty>
inline constexpr bool is_decimal_v = is_decimal <_Ty>::value;

template<typename _Ty,typename u>
requires(is_decimal_v<_Ty>)
[[nodiscard]] inline _Ty trapezoid_integral(const _Ty& a, const  _Ty& b, const  _Ty& dx,u&&f) {
    static_assert(is_decimal_v<std::invoke_result_t<decltype(f), _Ty>>, "return type of f must be a floating point type");
    static_assert(std::is_invocable_r_v<_Ty, u, _Ty>, "3rd argument must be a callable that returns a floating point value and takes only one floating point value");
    assert(b > a);
    assert(dx > 0);
    const size_t& iterations = static_cast<size_t>((b - a) / dx);

    _Ty result = 0;
    _Ty xi = a;
    for (size_t i = 0; i < iterations; i++) {
        result += static_cast<_Ty>((std::invoke(f,xi+dx) + std::invoke(f,xi)) * (dx * 0.5));
        std::cout << xi << "-" << xi + dx << " " << result << '\n';
        xi += dx;
    }

    return result;
}

/*#pragma once
#include <iostream>
#include <cmath>
#include<cassert>
#include<type_traits>
template<typename _Ty>
requires(std::disjunction_v<std::is_same<_Ty,float>,std::is_same<_Ty,double>,std::is_same<_Ty,long double>>)
[[nodiscard]] inline _Ty f(const _Ty& x) noexcept
{
    //{for pi  #include<numbers> => std::numbers::pi}
    //return pow(2,x) + 5 * x +3;
    return std::exp(2 * x);
    //return sin(2 * M_PI *x);
}
template<typename _Ty>
requires(std::disjunction_v<std::is_same<_Ty,float>,std::is_same<_Ty,double>,std::is_same<_Ty,long double>>)
[[nodiscard]] inline _Ty trapezoid_integral(const _Ty& a,const  _Ty& b,const  _Ty& dx) {
    assert(b > a);
    assert(dx > 0);
    const size_t &iterations = static_cast<size_t>((b - a) / dx);

    _Ty result = 0;
    _Ty xi = a;
    for (size_t i = 0; i < iterations; i++) {
        result += static_cast<_Ty>((f(xi + dx) + f(xi)) * (dx * 0.5));
        std::cout << xi << "-" << xi + dx << " " << result << '\n';
        xi += dx;
    }

    return result;
}
*/







