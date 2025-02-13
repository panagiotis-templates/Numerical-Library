#include <iostream>
#include <vector>
#include<string>
#include <cmath>
#include <limits>
#include <iomanip>
#include<type_traits>
#include<exception>



template<typename _Ty>
inline constexpr bool is_decimal_v = std::disjunction_v<std::is_same<_Ty, float>, std::is_same<_Ty, double>, std::is_same<_Ty, long double>>;

template<typename _Ty>
struct is_decimal :std::bool_constant<is_decimal_v<_Ty>> {};//tag dispatching must support it

template<typename _Ty>
requires(is_decimal_v<_Ty>)
[[nodiscard]] bool inline isEqual(const _Ty& a, const _Ty& b, const _Ty& epsilon = static_cast<_Ty>(10E-9))noexcept {

    return std::abs(a - b) < epsilon;
}
// Helper function
template<typename _Ty>
requires(is_decimal_v<_Ty>)
[[nodiscard]] _Ty gdiv( _Ty a,  _Ty b)
{
    if (isEqual<_Ty>(a,0.0) && isEqual<_Ty>(b, 0.0))
    {
        return(0.0);
    }
    else
    {
        return(a / b);
    }
}


// Function to calculate the B-spline value
/*
    i is the index of the B-spline basis function in the B-spline basis set.
    This is a zero-based index, so the first B-spline basis function in the set has index 0.

    ord is the order of the B-spline basis function.
    This is also known as the degree of the B-spline basis function.
    The order is an integer greater than or equal to 1 that determines the number of
    knots used to define the B-spline basis function.

    x is the value at which the B-spline basis function should be evaluated.

    nk is the number of knots in the knot vector kns.

    kns is a vector of knots that define the B-spline basis function.
    The knots are sorted in ascending order,
    and the B-spline basis function is defined over the range between the first and last knot.
    The size of kns should be equal to nk.

*/
template<typename _Ty>
requires(is_decimal_v<_Ty>)
[[nodiscard]] _Ty bsp(int i, int ord, _Ty x, int nk,const std::vector<_Ty>& kns) {
    // Check for illegal value of i
    if (i<0 || i>nk - ord - 1) {
        std::cout << "illegal i value: i=" << i << "; nk-ord=" << nk << "-" << ord << "=" << nk - ord <<'\n';
        return std::numeric_limits< _Ty>::quiet_NaN();
    }

    // Return 0 if x is outside the interval defined by the knots
    if (x<kns[i] || x>kns[i + ord]) return 0.0;

    // Remove repeated knots
    int k = nk - 1;
    if (k < 0) {
        std::cout << "k greated or equal to zero\n";
    }
    while (isEqual<_Ty>(kns[k] , kns[k - 1])) k--;
    k--;

    // If ord is 1, return 1 if x is within the interval defined by the knots
    if (ord == 1) {
        if (i != k) {
            return static_cast<_Ty>(((kns[i] <= x && x < kns[i + 1]) ? 1.0 : 0.0));
        }
        else {
            if (i == k) {
                return static_cast<_Ty>(((kns[i] <= x && x <= kns[i + 1]) ? 1.0 : 0.0));
            }
            else return std::numeric_limits<_Ty>::quiet_NaN();
        }
    }
    // If ord is greater than 1, return the sum of two recursive calls
    else {
        return static_cast<_Ty>((gdiv((x - kns[i]) * bsp(i, ord - 1, x, nk, kns), kns[i + ord - 1] - kns[i]) +
            gdiv((kns[i + ord] - x) * bsp(i + 1, ord - 1, x, nk, kns), kns[i + ord] - kns[i + 1])
            ));
    }
}


template<typename _Ty>
requires(is_decimal_v<_Ty>)
inline void basic_calc_plot(const std::vector<_Ty> &kns,int order) {
    if (kns.empty() ) {
        // Handle the case where the vector is empty
        std::cout << "Error empty knot vector\n";
        return ;
    }
    if (order<1) {
        std::cout << "order must be greater or equal to 1\n";
        return ;
    }
    // Create an empty 2D vector to store the data
    std::vector<std::vector <_Ty>> data;

    // Create a vector of labels
    std::vector<std::string> labels = { "Xi" };

    // Create the labels
    for (int j = 0; j < kns.size() - order; j++)
    {

        labels.push_back("B_{" + std::to_string(j + 1) + "," + std::to_string(order - 1) + "}"); // order-1 and j+1 because the index is zero based
    }

    for (const std::string& s : labels)//string_view??
    {
        std::cout << s << "   ";
    }
    std::cout << '\n';

    for (_Ty i = 0.0; i <= kns.back(); i +=static_cast<_Ty>( 0.001))
    {
        // Create a new row vector
        std::vector< _Ty> row;

        row.push_back(i);
        // Starts from P=1 (Order)
        for (size_t j = 0; j < kns.size() - order; j++)
        {
            // Calculate the values for each column

            row.push_back(bsp(static_cast<int>(j), order, i, static_cast<int>(kns.size()), kns));
        }
        // Add the row to the data vector
        data.push_back(row);
    }

    // Print the contents of the data vector 
    for (const auto& row : data) {
        for (const auto& value : row)
        {
            std::cout << std::fixed << std::setprecision(4) << value << " ";
        }
        std::cout << '\n';
    }
    return;
}



