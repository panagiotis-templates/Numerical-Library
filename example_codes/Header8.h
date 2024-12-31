#pragma once
#include"Macros.h"
#include<iostream>
#include<cmath>
#include<numbers>
#include<vector>
#include<type_traits>
#include<optional>
#include<cassert>
#include<iomanip>
#include<utility>
#include<functional>
_PA_BEGIN
template<typename _Ty>
requires(is_decimal_v<_Ty>)
_NODISCARD bool inline isEqual(const _Ty& a, const _Ty& b, const _Ty& epsilon = static_cast<_Ty>(10e-10))noexcept {

    return std::abs(a - b) < epsilon;
}
template<typename _Ty>
requires(is_decimal_v<_Ty>)
_NODISCARD inline std::optional<std::vector<_Ty>> gauss_elim(std::vector<std::vector<_Ty>>& A, std::vector<_Ty>& b) {//tes reference edo giati kaneis copy 
    const size_t& n = A.size(); // Number of equations (rows)
    if (n <= 0)return std::nullopt;
    // Augment the coefficient matrix with the right-hand side vector
    for (size_t i = 0; i < n; i++) {
        A[i].emplace_back(b[i]); // Appending the corresponding element from b to each row of A
    }
    /* for(int i=0; i<n; i++)
    {
        for(int j=0; j<A[0].size(); j++)
        {
            cout << A[i][j] << " ";
        }
        cout << endl;
    } */
    std::cout << '\n';
    // Perform Gaussian elimination
    for (size_t i = 0; i < n; i++) { // Loop over each row (equation)
        // Find the row with the maximum absolute value in the ith column and swap rows
        size_t max_row = i;
        for (size_t j = i + 1; j < n; j++) {
            if (std::abs(A[j][i]) > std::abs(A[max_row][i])) {
                max_row = j; // Update the index of the row with the maximum absolute value
            }
        }
        std::swap(A[i], A[max_row]); // Swap the current row with the row with the maximum absolute value
        // Perform row operations to eliminate coefficients below the pivot element
        for (size_t j = i + 1; j < n; j++) { // Loop over rows below the pivot row
            double factor = A[j][i] / A[i][i]; // Compute the factor by which the pivot row will be multiplied
            for (size_t k = i; k < n + 1; k++) { // Loop over columns including the augmented column
                A[j][k] -= factor * A[i][k]; // Perform row operation to eliminate coefficients below the pivot element
            }
        }
    }
    if (isEqual(A[A.size() - 1][A.size() - 1], 0.0)) {
        return std::nullopt;
    }
    /*  for(int i=0; i<n; i++)
     {
         for(int j=0; j<A[0].size(); j++)
         {
             cout << A[i][j] << " ";
         }
         cout << endl;
     } */
     // Back-substitution to solve for x
    std::vector<_Ty> x(n, 0.0); // Initialize the solution vector x with zeros
    for (size_t i = n - 1; i > 0; i--) { // Start from the last equation and move upwards
        _Ty sum = 0.0;
        for (size_t j = i + 1; j < n; j++) { // Loop over elements to the right of the diagonal in the current row
            sum += A[i][j] * x[j]; // Compute the sum of products of coefficients and corresponding elements of x
        }
        if (isEqual(A[i][i], 0.0)) {
            return std::nullopt;
        }
        x[i] = (A[i][n] - sum) / A[i][i]; // Compute the value of x[i] using back-substitution
    }
    size_t i = 0;
    while (i == 0) {
        _Ty sum = 0.0;
        for (size_t j = i + 1; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        if (isEqual(A[i][i], 0.0)) {
            return std::nullopt;
        }
        x[i] = (A[i][n] - sum) / A[i][i];
        break;
    }
    return std::optional{ x }; // Return the solution vector x
}

template <typename _Ty>
requires(is_decimal_v<_Ty>)
inline void print2DVector(const std::vector<std::vector<_Ty>>& vec) {
    const size_t& numRows = vec.size();
    const size_t& numCols = vec[0].size();
   
    if (!(vec.size() > 0 && vec[0].size() > 0)) {
        std::cerr << "vec.size() > 0 && vec[0].size() > 0" << '\n';
        return;
    }
    // Print column headers
    std::cout << std::setw(12) << " ";
    for (size_t col = 0; col < numCols; ++col) {
        std::cout << std::setw(12) << col;
    }
    std::cout << '\n';
    // Print row headers and vector contents
    for (size_t row = 0; row < numRows; ++row) {
        std::cout << std::setw(12) << row;
        for (size_t col = 0; col < numCols; ++col) {
            std::cout << std::setw(12) << std::fixed << std::setprecision(6) << vec[row][col] << "("
                << row << "," << col << ")";
        }
        std::cout << '\n';
    }
}
template<typename _Ty>
requires(is_decimal_v<_Ty>)
inline void print_Vector(const std::vector<_Ty>& vec) {
    const size_t& numCols = vec.size();
   
    if (numCols <= 0) {
        std::cerr << "numCols>0" << '\n';
        return;
    }
    for (size_t col = 0; col < numCols; ++col) {
        std::cout << std::setw(12) << std::fixed << std::setprecision(6) << vec[col] << "(" << col
            << ")";
    }
    std::cout << '\n';
}
template<typename _Ty, typename u, typename v>
requires(is_decimal_v<_Ty>)
inline void finite_difference(const _Ty& a, const _Ty& b, const  _Ty& h, u&& f, v&& ddf) {
    static_assert(is_decimal_v<std::invoke_result_t<decltype(f), _Ty>>, "return type of f  must be a floating point type");
    static_assert(std::is_invocable_r_v<_Ty, u, _Ty>, "4th argument must be a callable that returns a floating point value and takes only one floating point value");

    static_assert(is_decimal_v<std::invoke_result_t<decltype(ddf), _Ty>>, "return type of ddf must be a floating point type");
    static_assert(std::is_invocable_r_v<_Ty, v, _Ty>, "5th argument must be a callable that returns a floating point value and takes only one floating point value");
    _Ty xi{}, res{}, hsq{ static_cast<_Ty>(std::pow(h, 2)) };
    if (b <= a || h <= 0) {
        std::cerr << "b>a &&h>0" << '\n';
        return;
    }
    
    const size_t& n{ static_cast<size_t>((b - a) / h + 1) };
    xi = a;
    
    for (size_t i = 0; i < n; i++) {
        res = (std::invoke(f, xi - h) - 2 * std::invoke(f, xi) + std::invoke(f, xi + h)) / hsq;
        std::cout << xi << " " << res << " " << std::invoke(ddf, xi) << " " << std::abs(std::invoke(ddf, xi) - res) << '\n';
        xi += h;
    }

}



template<typename _Ty, typename u>
requires(is_decimal_v<_Ty>)
_NODISCARD inline std::optional<_Ty> dicection(_Ty a, _Ty b, u&& f, _Ty e = static_cast<_Ty>(1.0E-10))//u is the callable
{
    static_assert(is_decimal_v<std::invoke_result_t<decltype(f), _Ty>>, "return type of f must be a floating point type");
    static_assert(std::is_invocable_r_v<_Ty, u, _Ty>, "3rd argument must be a callable that returns a floating point value and takes only one floating point value");

    if (b <= a || e <= 0) {
        return std::nullopt;
    }
    _Ty d = b - a, c{}; //[a,b] interval ,d interval span ,e precision

    while (true)
    {
        d *= 0.5;
        if (std::signbit(std::invoke(f, a)) == std::signbit(std::invoke(f, b)) || d < e || isEqual(d, e))
        {
            std::cout << "Not found in interval: " << "[" << a << "," << b << "]";
            return std::nullopt;
        }
        c = static_cast<_Ty>(a + (b - a) * 0.5);
        std::cout << a << " " << b << " " << c << " " << d << " " << f(c) << '\n';
        if (std::abs(std::invoke(f, c) - 0) < 10e-6) //Smaller precision for the final result 
        {

            return std::optional{ c };
        }
        if (std::signbit(std::invoke(f, c)) == std::signbit(std::invoke(f, a))) //if(f(c) * f(a) > 0) is going to produce overflow if f(a) or f(c) is  small number
        {

            a = c;

        }
        else
        {

            b = c;
        }
    }
}

template<typename _Ty, typename u>
requires(is_decimal_v<_Ty>)
_NODISCARD inline std::optional<_Ty> trapezoid_integral(const _Ty& a, const  _Ty& b, const  _Ty& dx, u&& f) {
    static_assert(is_decimal_v<std::invoke_result_t<decltype(f), _Ty>>, "return type of f must be a floating point type");
    static_assert(std::is_invocable_r_v<_Ty, u, _Ty>, "3rd argument must be a callable that returns a floating point value and takes only one floating point value");
    if (b <= a || dx <= 0) {
        std::cerr << "b>a &&dx>0" << '\n';
        return std::nullopt;
    }
    

const size_t& iterations = static_cast<size_t>((b - a) / dx);

_Ty result = 0;
_Ty xi = a;
for (size_t i = 0; i < iterations; i++) {
    result += static_cast<_Ty>((std::invoke(f, xi + dx) + std::invoke(f, xi)) * (dx * 0.5));
    std::cout << xi << "-" << xi + dx << " " << result << '\n';
    xi += dx;
}

return std::optional{ result };
}

template<typename _Ty, typename u>
    requires(is_decimal_v<_Ty>)
_NODISCARD std::optional<_Ty> simpson(const _Ty& a, const  _Ty& b, const  _Ty& dx, u&& f) {
    static_assert(is_decimal_v<std::invoke_result_t<decltype(f), _Ty>>, "return type of f must be a floating point type");
    static_assert(std::is_invocable_r_v<_Ty, u, _Ty>, "4rd argument must be a callable that returns a floating point value and takes only one floating point value");
    if (b <= a || dx <= 0) {
        std::cerr << "b>a &&dx>0" << '\n';
        return std::nullopt;
    }

    const size_t& iterations = static_cast<size_t>((b - a) / dx);
    _Ty result = 0;
    _Ty xi = a;

    for (size_t i = 0; i < iterations; i++)
    {
        result += (dx / 6) * (std::invoke(f, xi) + 4 * std::invoke(f, (xi + (xi + dx)) * 0.5) + std::invoke(f, xi + dx));
        std::cout << xi << "-" << xi + dx << " " << result << '\n';
        xi += dx;
    }

    return std::optional{ result };
}
template<typename _Ty>
    requires(is_decimal_v<_Ty>)
_NODISCARD inline  std::optional<bool> strict_diagonal_dominace(const std::vector<std::vector<_Ty>>& A)noexcept //Pass it as refrence so you dont need copy,also add const if you dont goining to modify
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
template<typename _Ty>
    requires(is_decimal_v<_Ty>)
_NODISCARD  inline std::optional<_Ty> better_sum(const _Ty& n)noexcept
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



//Algorithm from page 159 dougalis
//Slower implementation but more organized code
template<typename _Ty>
    requires(is_decimal_v<_Ty>)
_NODISCARD inline std::optional<_Ty> fi(const _Ty& x, const std::vector<_Ty>& knot, const size_t& i)noexcept
{

    if (i == 0) //First function
    {
        if ((x > knot[0] || isEqual(x, knot[0])) && (x < knot[1] || isEqual(x, knot[1])))
        {
            //cout << "case1";
            if (isEqual(knot[1] - knot[0],0.0)){
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
            if (isEqual(knot[i] - knot[i-1],0.0)){
                return std::nullopt;
            }
            return std::optional{ (x - knot[i - 1]) / (knot[i] - knot[i - 1]) };
        }
        else if ((x > knot[i] || isEqual(x, knot[i])) && (x < knot[i + 1] || isEqual(x, knot[i + 1])))
        {
            //cout << "case3";
            if (isEqual(knot[i+1] - knot[i],0.0)){
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
            if (isEqual(knot[knot.size() - 1] - knot[knot.size() - 2],0.0)){
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
_NODISCARD inline std::optional< _Ty> fi2(const _Ty& x, const std::vector<_Ty>& knot, const size_t& i)noexcept
{

    if (i > 0&& i<knot.size())
        if ((x > knot[i - 1] || isEqual(x, knot[i - 1])) && (x < knot[i] || isEqual(x, knot[i])))
        {
            if (isEqual(knot[i] - knot[i-1],0.0)){
                return std::nullopt;
            }
            return std::optional{ (x - knot[i - 1]) / (knot[i] - knot[i - 1]) };
        }
    if (i != knot.size() - 1)
        if ((x > knot[i] || isEqual(x, knot[i])) && (x < knot[i + 1] || isEqual(x, knot[i + 1])))
        {
            if (isEqual(knot[i+1] - knot[i],0.0)){
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

_PA_END



