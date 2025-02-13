#pragma once
#include<type_traits>
#include<exception>
#include<string>

#define _NODISCARD [[nodiscard]]

#define _CONSTEXPR20 constexpr
#define _PA_BEGIN namespace panagiotis{

#define _PA_END }
#define _PANAGIOTIS panagiotis

#if __cplusplus>201703L
template<typename _Ty>

inline constexpr bool is_double_or_long_double_v = std::disjunction_v<std::is_same<_Ty, double>, std::is_same<_Ty, long double>>;

template<typename _Ty>
struct is_double_or_long_double :std::bool_constant<is_double_or_long_double_v<_Ty>> {};//tag dispatching must support it


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



#endif
