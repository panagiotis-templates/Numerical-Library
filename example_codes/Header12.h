#pragma once
#include<type_traits>
#include<exception>
#include<string>
#include<iostream>
#include<cmath>
#include<numbers>
#include<iomanip>
#include<utility>
#include<functional>

template<typename _Ty>

inline constexpr bool is_decimal_v = std::disjunction_v<std::is_same<_Ty, float>, std::is_same<_Ty, double>, std::is_same<_Ty, long double>>;

template<typename _Ty>
struct is_decimal :std::bool_constant<is_decimal_v<_Ty>>{};//tag dispatching must support it

class divisionWithZero :public std::exception {
private:
    std::string errorMessage; // To store the error message
public:
    // Constructor to initialize the error message
    explicit divisionWithZero(const std::string& message)
        : errorMessage(message) {
    }

    // Override the what() method
    const char* what() const noexcept override {
        return errorMessage.c_str();
    }
};
template<typename _Ty>
requires(is_decimal_v<_Ty>)
[[nodiscard]] bool inline isEqual(const _Ty& a, const _Ty& b, const _Ty& epsilon = static_cast<_Ty>(10e-10))noexcept {

    return std::abs(a - b) < epsilon;
}
template<typename _Ty, typename u>
requires(is_decimal_v<_Ty>)
[[nodiscard]] inline  _Ty derivative(u&& f, const  _Ty& x0, int order, const _Ty& delta = static_cast<_Ty>(1.0e-6)) //Numerical differation
{

    _Ty x1 = x0 - delta;
    _Ty x2 = x0 + delta;
    if (order == 1) {
        _Ty y1 = std::invoke(f, x1);
        _Ty y2 = std::invoke(f, x2);
        if (isEqual<_Ty>(x2 - x1, 0.0))throw divisionWithZero{ "division with zero\n" };
        return  (y2 - y1) / (x2 - x1);
    }
    else {
        _Ty  y1 = derivative(f, x1, order - 1);
        _Ty  y2 = derivative(f, x2, order - 1);
        if (isEqual<_Ty>(x2 - x1, 0.0))throw divisionWithZero{ "division with zero\n" };
        return  (y2 - y1) / (x2 - x1);
    }
}
template<typename _Ty, typename u>
requires(is_decimal_v<_Ty>)
inline void  finite_diff(const _Ty& a, const  _Ty& b, const _Ty& h, u&& f) {
    static_assert(std::is_same_v<std::invoke_result_t<decltype(f), _Ty>, _Ty>, "return type of f  must be the same with a,b,h");
    static_assert(std::is_invocable_r_v<_Ty, u, _Ty>, "4th argument must be a callable that returns a floating point value and takes only one floating point value");
    if (b - a < 0 ||isEqual<_Ty>(b,a))return;
    if (h < 0||isEqual<_Ty>(h,0.0))return;
    size_t n = static_cast<size_t>((b - a) / h + 1);
    _Ty xi = a;
    _Ty hsq = h * h;
    _Ty result{};
    for (size_t i = 0; i < n; i++)
    {
        result =static_cast<_Ty>( (std::invoke(f,xi - h) - 2 *std::invoke (f,xi) +std::invoke(f,xi + h)) / hsq); //Finite difference
        auto d1 = derivative(f, xi, 2);

        std::cout << xi << " " << result << " " << d1 << " " << std::abs(d1 - result) << '\n';

        xi += h; //Increment the xi for the next iteration 
    }
    return;

}
