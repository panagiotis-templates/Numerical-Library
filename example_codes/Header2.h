#pragma once
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
   
    _Ty xi{}, res{}, hsq{ std::pow(h, 2) };
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
