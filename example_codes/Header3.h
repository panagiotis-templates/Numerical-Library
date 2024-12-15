/*#pragma once
this is the first version!!!
#include<iostream>
#include<cmath>
#include<optional>
#include<concepts>
#include<cassert>
template<std::floating_point _Floating>
[[nodiscard]] inline _Floating f(_Floating x) noexcept
{
    return static_cast<_Floating>(std::pow(x, 3) - 2 * x - 5);
    //return x * x + 2 * x;
}

template<std::floating_point _Floating>
[[nodiscard]] bool inline isEqual(const _Floating a, const _Floating b, const _Floating epsilon = static_cast<_Floating>(10e-10))noexcept {

    return std::abs(a - b) < epsilon;
}

template<std::floating_point _Floating>
[[nodiscard]] inline std::optional<_Floating> dicection(_Floating a, _Floating b, _Floating(*f)(_Floating), _Floating e = static_cast<_Floating>(10e-10))
{
    assert(b > a);
    assert(e > 0);
    _Floating d = b - a, c{}; //[a,b] interval ,d interval span ,e precision
    while (true)
    {
        d *= 0.5;
        if (std::signbit(f(a)) == std::signbit(f(b)) || d < e || isEqual(d, e))
        {
            std::cout << "Not found in interval: " << "[" << a << "," << b << "]";
            return std::nullopt;
        }
        c = static_cast<_Floating>(a + (b - a) * 0.5);
        std::cout << a << " " << b << " " << c << " " << d << " " << f(c) << '\n';
        if (std::abs(f(c) - 0) < 10e-6) //Smaller precision for the final result 
        {

            return c;
        }
        if (std::signbit(f(c)) == std::signbit(f(a))) //if(f(c) * f(a) > 0) is going to produce overflow if f(a) or f(c) is  small number
        {

            a = c;

        }
        else
        {

            b = c;
        }
    }
}
*/
//this is the second more flexible!!! gia tyxaies f
#pragma once
#include<iostream>
#include<cmath>
#include<optional>
#include<concepts>
#include<cassert>
#include<type_traits>
template<typename _Ty>
requires(std::disjunction_v<std::is_same<_Ty,float>,std::is_same<_Ty,double>,std::is_same<_Ty,long double>>)
[[nodiscard]] inline  _Ty  f(const _Ty& x) noexcept
{
    return static_cast<_Ty>(std::pow(x, 3) - 2 * x - 5);
    //return x * x + 2 * x;
}

template<typename _Ty>
requires(std::disjunction_v<std::is_same<_Ty,float>,std::is_same<_Ty,double>,std::is_same<_Ty,long double>>)
[[nodiscard]] bool inline isEqual(const _Ty& a, const _Ty& b, const _Ty& epsilon = static_cast<_Ty>(10e-10))noexcept {

    return std::abs(a - b) < epsilon;
}

template<typename _Ty,typename u>
requires(std::disjunction_v<std::is_same<_Ty,float>,std::is_same<_Ty,double>,std::is_same<_Ty,long double>>)
[[nodiscard]] inline std::optional<_Ty> dicection(_Ty a, _Ty b, u&& f, _Ty e = static_cast<_Ty>(10e-10))//u is the callable
{

    static_assert(std::is_same_v<decltype(f(a)), _Ty>, "return type of f must be a floating point type");
    static_assert(std::is_invocable_r_v<_Ty, u, _Ty>, "3rd argument must be a callable that returns a floating point value and takes only one floating point value");
   
    assert(b > a);
    assert(e > 0);
    _Ty d = b - a, c{}; //[a,b] interval ,d interval span ,e precision
    
    while (true)
    {
        d *= 0.5;
        if (std::signbit(f(a)) == std::signbit(f(b)) || d < e || isEqual(d, e))
        {
            std::cout << "Not found in interval: " << "[" << a << "," << b << "]";
            return std::nullopt;
        }
        c = static_cast<_Ty>(a + (b - a) * 0.5);
        std::cout << a << " " << b << " " << c << " " << d << " " << f(c) << '\n';
        if (std::abs(f(c) - 0) < 10e-6) //Smaller precision for the final result 
        {

            return std::optional{c};
        }
        if (std::signbit(f(c)) == std::signbit(f(a))) //if(f(c) * f(a) > 0) is going to produce overflow if f(a) or f(c) is  small number
        {

            a = c;

        }
        else
        {

            b = c;
        }
    }
}



