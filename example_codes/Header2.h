#pragma once

#include<concepts>
#include<vector>
#include<iostream>
#include<iomanip>
template<std::floating_point _Floating>
inline [[nodiscard]] _Floating  f(_Floating x) noexcept
{
    return exp(2 * x);
}
template<std::floating_point _Floating>
inline [[nodiscard]] _Floating ddf(_Floating x) noexcept
{
    return 4 * exp(2 * x);
}
template<std::floating_point _Floating>
inline void finite_difference(_Floating a, _Floating b, _Floating h) {
    _Floating xi{}, res{}, hsq{};
    size_t n{ static_cast<size_t>((b - a) / h + 1) };
    xi = a;
    hsq = std::pow(h, 2);
    for (size_t i = 0; i < n; i++) {
        res = (f(xi - h) - 2 * f(xi) + f(xi + h)) / hsq;
        std::cout << xi << " " << res << " " << ddf(xi) << " " << fabs(ddf(xi) - res) << '\n';
        xi += h;
    }

}