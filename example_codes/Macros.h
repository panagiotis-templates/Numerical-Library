#include<type_traits>
#include<exception>

//for our own customs things
#define _NODISCARD [[nodiscard]]

#define _CONSTEXPR20 constexpr
#define _PA_BEGIN namespace panagiotis{

#define _PA_END }
#define _PANAGIOTIS panagiotis

#if __cplusplus>201703L

template<typename _Ty>

inline constexpr bool is_decimal_v = std::disjunction_v<std::is_same<_Ty, float>, std::is_same<_Ty, double>, std::is_same<_Ty, long double>>;

template<typename _Ty>
struct is_decimal :std::bool_constant<is_decimal_v<_Ty>>{};//tag dispatching must support it




#endif
