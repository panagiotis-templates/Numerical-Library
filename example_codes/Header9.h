#pragma once
#include<iostream>
#include<vector>
#include<type_traits>
#include<optional>
#include <iomanip>
#include <chrono>
#include <random>
template<typename _Ty>inline constexpr bool is_decimal_v = std::disjunction_v<std::is_same<_Ty, float>, std::is_same<_Ty, double>, std::is_same<_Ty, long double>>;

template<typename _Ty>
struct is_decimal :std::bool_constant<is_decimal_v<_Ty>>{};//tag dispatching must support it

template<typename _Ty>
requires(is_decimal_v<_Ty>)
[[nodiscard]] bool inline isEqual(const _Ty& a, const _Ty& b, const _Ty& epsilon = static_cast<_Ty>(10e-10))noexcept {

    return std::abs(a - b) < epsilon;
}

template<typename _Ty>
requires(is_decimal_v<_Ty>)
[[nodiscard]] inline std::optional<_Ty> fi(const _Ty& x, const std::vector<_Ty>& knot, const size_t& i)
{

    if (i == 0) //First function
    {
        if ((x > knot[0] || isEqual(x, knot[0])) && (x < knot[1] || isEqual(x, knot[1])))
        {
            //cout << "case1";
            if (isEqual<_Ty>(knot[1] - knot[0],0.0)){
                return std::nullopt;
            }
            return std::optional{ -(x - knot[1]) / (knot[1] - knot[0]) };
        }
        else
        {
            return std::optional<_Ty>{ 0 };
        }
    }
    //For every other function
    if (i != 0 || i != knot.size() - 1)
    {
        if ((x > knot[i - 1] || isEqual(x, knot[i - 1])) && (x < knot[i] || isEqual(x, knot[i])))
        {
            //cout << "case2";
            if (isEqual<_Ty>(knot[i] - knot[i-1],0.0)){
                return std::nullopt;
            }
            return std::optional{ (x - knot[i - 1]) / (knot[i] - knot[i - 1]) };
        }
        else if ((x > knot[i] || isEqual(x, knot[i])) && (x < knot[i + 1] || isEqual(x, knot[i + 1])))
        {
            //cout << "case3";
            if (isEqual<_Ty>(knot[i+1] - knot[i],0.0)){
                return std::nullopt;
            }
            return std::optional{ -(x - knot[i + 1]) / (knot[i + 1] - knot[i]) };
        }
        else
        {
            return std::optional< _Ty>{ 0 }; //No support at that x
        }
    }


    if (i == knot.size() - 1) //Last fucntion
    {
        if ((x > knot[knot.size() - 2] || isEqual(x, knot[knot.size() - 2])) && (x < knot[knot.size() - 1] || isEqual(x, knot[knot.size() - 1])))
        {
            //cout << "case4";
            if (isEqual<_Ty>(knot[knot.size() - 1] - knot[knot.size() - 2],0.0)){
                return std::nullopt;
            }
            return std::optional{ (x - knot[knot.size() - 2]) / (knot[knot.size() - 1] - knot[knot.size() - 2]) };
        }
        else
        {
            return std::optional< _Ty>{0};
        }
    }

    return std::nullopt;


}

//Faster implementation
template<typename _Ty>
requires(is_decimal_v< _Ty>)
[[nodiscard]] inline std::optional< _Ty> fi2(const _Ty& x, const std::vector<_Ty>& knot, const size_t& i)
{

    if (i > 0&& i<knot.size())
        if ((x > knot[i - 1] || isEqual(x, knot[i - 1])) && (x < knot[i] || isEqual(x, knot[i])))
        {
            if (isEqual<_Ty>(knot[i] - knot[i-1],0.0)){
                return std::nullopt;
            }
            return std::optional{ (x - knot[i - 1]) / (knot[i] - knot[i - 1]) };
        }
    if (i != knot.size() - 1)
        if ((x > knot[i] || isEqual(x, knot[i])) && (x < knot[i + 1] || isEqual(x, knot[i + 1])))
        {
            if (isEqual<_Ty>(knot[i+1] - knot[i],0.0)){
                return std::nullopt;
             }
            return std::optional{ -(x - knot[i + 1]) / (knot[i + 1] - knot[i]) };
        }

    if (i >= 0 && i < knot.size())
    {
        return std::optional<_Ty>{ 0 };
    }
    return std::nullopt;
}
