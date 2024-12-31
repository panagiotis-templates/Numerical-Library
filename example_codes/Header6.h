#pragma once
#include<cmath>
#include<iostream>
#include<vector>
#include<type_traits>
#include<cassert>
template<typename _Ty>

inline constexpr bool is_decimal_v = std::disjunction_v<std::is_same<_Ty, float>, std::is_same<_Ty, double>, std::is_same<_Ty, long double>>;

template<typename _Ty>
struct is_decimal :std::bool_constant<is_decimal_v<_Ty>>{};//tag dispatching must support it

template<typename _Ty>
requires(is_decimal_v<_Ty>)
[[nodiscard]] inline  bool  isEqual(const _Ty& a, const _Ty& b,
    const _Ty& epsilon = static_cast<_Ty>(1.0E-8))noexcept {

    return std::abs(a - b) < epsilon;
}



template<typename _Ty>
requires(is_decimal_v<_Ty>)
[[nodiscard]]inline  std::optional<bool> strict_diagonal_dominace(const std::vector<std::vector<_Ty>>& A)noexcept //Pass it as refrence so you dont need copy,also add const if you dont goining to modify
{
    _Ty sum = 0;
    size_t n = A.size();
    if (A.size() <= 0) {
        std::cerr << "vector must be nxn with n>0" << '\n';
        return std::nullopt;
    }

    for (size_t i = 0; i < A.size(); i++)
    {
        sum = 0;//Reset the sum for new row
        if (A[i].size() != A.size()) {
            std::cerr << "vector must be nxn with n>0" << '\n';
            return std::nullopt;
        }

        for (size_t j = 0; j < A[i].size(); j++)
        {
            if (i != j) { sum += A[i][j]; } //Sum the values except the diagonal
            else { continue; }
        }
        if (std::abs(A[i][i]) < std::abs(sum) || isEqual(std::abs(A[i][i]), std::abs(sum))) //Check if one diagonal element is smaller than the sum
        {
            return std::optional{ false };
        }
    }
    return std::optional{ true };
}



