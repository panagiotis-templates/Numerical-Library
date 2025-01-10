#pragma once
#include<iostream>
#include<type_traits>
#include<optional>
#include<cmath>
template<typename _Ty>

inline constexpr bool is_decimal_v = std::disjunction_v<std::is_same<_Ty, float>, std::is_same<_Ty, double>, std::is_same<_Ty, long double>>;

template<typename _Ty>
struct is_decimal :std::bool_constant<is_decimal_v<_Ty>>{};//tag dispatching must support it

template<typename _Ty>
requires(is_decimal_v<_Ty>)
[[nodiscard]] bool inline isEqual(const _Ty& a, const _Ty& b, const _Ty& epsilon = static_cast<_Ty>(10e-10))noexcept {

    return std::abs(a - b) < epsilon;
}
template<typename _Ty>
requires(is_decimal_v<_Ty>)
[[nodiscard]] inline std::optional<_Ty> better_sum(const _Ty& n)noexcept
{
    if (n <= 0) {
        std::cerr << "b must be greater than zero\n";
        return std::nullopt;
    }

    _Ty sum = 1. / (std::pow(n, 2) + n);
    for (size_t i = 1; i < n; i++)
    {
        sum += 1 / ((n - i) * (n - i + 1));
    }
    sum += 1;
    return std::optional{ sum };
}
