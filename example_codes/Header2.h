#pragma once
#include<iostream>
#include<cmath>
#include<cassert>
#include<type_traits>
#include<utility>
#include<functional>
template<typename _Ty>

inline constexpr bool is_decimal_v = std::disjunction_v<std::is_same<_Ty, float>, std::is_same<_Ty, double>, std::is_same<_Ty, long double>>;

template<typename _Ty>
struct is_decimal :std::bool_constant<is_decimal_v<_Ty>>{};//tag dispatching must support it

template<typename _Ty, typename u, typename v>
requires(is_decimal_v<_Ty>)
inline void finite_difference(const _Ty& a, const _Ty& b, const  _Ty& h, u&& f, v&& ddf) {
    static_assert(std::is_same_v<std::invoke_result_t<decltype(f), _Ty>,_Ty>, "return type of f  must be the same with a,b,h");
    static_assert(std::is_invocable_r_v<_Ty, u, _Ty>, "4th argument must be a callable that returns a floating point value and takes only one floating point value");

    static_assert(std::is_same_v<std::invoke_result_t<decltype(ddf), _Ty>,_Ty>, "return type of ddf  must be the same with a,b,h");
    static_assert(std::is_invocable_r_v<_Ty, v, _Ty>, "5th argument must be a callable that returns a floating point value and takes only one floating point value");
    _Ty xi{}, res{}, hsq{ static_cast<_Ty>(std::pow(h, 2)) };
    if (b <= a || h <= 0) {
        std::cerr << "b>a &&h>0" << '\n';
        return;
    }
    
    const size_t& n{ static_cast<size_t>((b - a) / h + 1) };
    xi = a;
    
    for (size_t i = 0; i < n; i++) {
        res = (std::invoke(f, xi - h) - 2 * std::invoke(f, xi) + std::invoke(f, xi + h)) / hsq;
        std::cout << xi << " " << res << " " << std::invoke(ddf, xi) << " " << std::abs(std::invoke(ddf, xi) - res) << '\n';
        xi += h;
    }

}





/*#pragma once
#include<iostream>
#include<cmath>
#include<cassert>
#include<type_traits>
template<typename _Ty>
requires(std::disjunction_v<std::is_same<_Ty,float>,std::is_same<_Ty,double>,std::is_same<_Ty,long double>>)
[[nodiscard]]  inline _Ty f(const _Ty& x) noexcept
{
    return std::exp(2 * x);
}
template<typename _Ty>
requires(std::disjunction_v<std::is_same<_Ty,float>,std::is_same<_Ty,double>,std::is_same<_Ty,long double>>)
[[nodiscard]] inline  _Ty ddf(const _Ty& x) noexcept
{
    return 4 * std::exp(2 * x);
}
template<typename _Ty>
requires(std::disjunction_v<std::is_same<_Ty,float>,std::is_same<_Ty,double>,std::is_same<_Ty,long double>>)
inline void finite_difference(const _Ty& a, const _Ty& b, const  _Ty& h) {
   
    _Ty xi{}, res{}, hsq{ static_cast<_Ty>(std::pow(h, 2)) };
    const size_t& n{ static_cast<size_t>((b - a) / h + 1) };
    xi = a;
    assert(b>a);
    assert(h>0);
    for (size_t i = 0; i < n; i++) {
        res = (f(xi - h) - 2 * f(xi) + f(xi + h)) / hsq;
        std::cout << xi << " " << res << " " << ddf(xi) << " " << std::abs(ddf(xi) - res) << '\n';
        xi += h;
    }

}
*/
