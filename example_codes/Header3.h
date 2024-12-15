#include<iostream>
#include<cmath>
#include<optional>
#include<concepts>

template<std::floating_point _Floating>
[[nodiscard]]inline _Floating f(_Floating x) noexcept
{
    return static_cast<_Floating>(std::pow(x, 3) - 2 * x - 5);
    //return x * x + 2 * x;
}

template<std::floating_point _Floating>
[[nodiscard]] bool inline isEqual(const _Floating a, const _Floating b,const _Floating epsilon= static_cast<_Floating>(10e-10))noexcept {

    return std::abs(a - b) < epsilon;
}

template<std::floating_point _Floating>
[[nodiscard]] inline std::optional<_Floating> dicection(_Floating a, _Floating b, _Floating(*f)(_Floating), _Floating e =static_cast<_Floating>( 10e-10))
{
    _Floating d = b - a, c{}; //[a,b] interval ,d interval span ,e precision
    while (true)
    {
        d *= 0.5;
        if (std::signbit(f(a)) == std::signbit(f(b)) || d < e || isEqual(d, e))
        {
            std::cout << "Not found in interval: " << "[" << a << "," << b << "]";
            return std::nullopt;
        }
        c =static_cast<_Floating>( a + (b - a) * 0.5);
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



