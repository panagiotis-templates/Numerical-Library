#pragma once
#include<iostream>
#include<cmath>
#include<numbers>
#include<vector>
#include<type_traits>
#include<optional>
#include<iomanip>
#include<utility>
#include<functional>

template<typename _Ty>

inline constexpr bool is_decimal_v = std::disjunction_v<std::is_same<_Ty, float>, std::is_same<_Ty, double>, std::is_same<_Ty, long double>>;

template<typename _Ty>
struct is_decimal :std::bool_constant<is_decimal_v<_Ty>> {};//tag dispatching must support it



//Utility fucntions
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





template<typename _Ty>
    requires(is_decimal_v<_Ty>)
[[nodiscard]] bool inline isEqual(const _Ty& a, const _Ty& b, const _Ty& epsilon = static_cast<_Ty>(10e-10))noexcept {

    return std::abs(a - b) < epsilon;
}
template<typename _Ty>
    requires(is_decimal_v<_Ty>)
[[nodiscard]] inline std::optional<std::vector<_Ty>> gauss_elim(std::vector<std::vector<_Ty>>& A, std::vector<_Ty>& b) {//tes reference edo giati kaneis copy 
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
            _Ty factor = A[j][i] / A[i][i]; // Compute the factor by which the pivot row will be multiplied
            for (size_t k = i; k < n + 1; k++) { // Loop over columns including the augmented column
                A[j][k] -= factor * A[i][k]; // Perform row operation to eliminate coefficients below the pivot element
            }
        }
    }
    if (isEqual<_Ty>(A[A.size() - 1][A.size() - 1], 0.0)) {
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
        if (isEqual<_Ty>(A[i][i], 0.0)) {
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
        if (isEqual<_Ty>(A[i][i], 0.0)) {
            return std::nullopt;
        }
        x[i] = (A[i][n] - sum) / A[i][i];
        break;
    }
    return std::optional{ x }; // Return the solution vector x
}






template<typename _Ty,typename u,typename v>
requires(is_decimal_v<_Ty>)
std::optional<_Ty> finite_diff_central(const _Ty& a,const _Ty& b,const _Ty& h,u&&f,v&&y) {
    static_assert(std::is_same_v<std::invoke_result_t<decltype(f), _Ty,_Ty>,_Ty>, "return type of f  must be the same with a,b,h");
    static_assert(std::is_invocable_r_v<_Ty, u, _Ty,_Ty>, "4th argument must be a callable that returns a floating point value and takes only one floating point value");

    static_assert(is_decimal_v<std::invoke_result_t<decltype(y), _Ty>>, "return type of ddf must be a floating point type");
    static_assert(std::is_invocable_r_v<_Ty, v, _Ty>, "5th argument must be a callable that returns a floating point value and takes only one floating point value");

    _Ty xi, hsq, max = -1, sum = 0; //Utility variables
    _Ty  y_a = 0, y_b = static_cast<_Ty>(std::exp(-2.0)), q = 2;
    if (h <= 0)return std::nullopt;
    if (b - a <= 0)return std::nullopt;
    size_t n=static_cast<size_t>( (b - a) / h + 1);
    std::optional<std::vector<_Ty>> U;
    std::vector<_Ty>F(n); //Solution vector U,Right hand side F
    std::vector<std::vector<_Ty>> A(n, std::vector<_Ty>(n, 0));
    xi = a;
    hsq = h * h; //Single time calculation of h squared
    for (size_t i = 0; i < n; i++)
    {
        //Build the right hand side 
        F[i] = std::invoke(f,xi, q) * hsq;
        //Build the left hand side 
        for (size_t g = 0; g < n; g++)
        {
            if (i == g) { A[i][g] = 2 + hsq * q; }
            else if (i == g - 1 || i == g + 1) { A[i][g] = -1; }
        }
        xi += h; //Increment the xi for the next iteration 
    }
    F[0] = y_a;
    F[n - 1] = y_b;
    U = gauss_elim(A, F);
    if (!U.has_value())return std::nullopt;
    std::vector<_Ty>u = std::move(U.value());
    xi = a; //Reset
    for (size_t i = 0; i < n; i++)
    {
        sum +=static_cast<_Ty>( std::sqrt(std::pow(std::abs(u[i] - std::invoke(y,xi)), 2)));
        std::cout << xi << " " << std::invoke(y, xi) << " " << u[i] << '\n';
        if (std::abs(u[i] - std::invoke(y, xi)) > max) //Find max difference
        {
            max = std::abs(u[i] - std::invoke(y, xi));
        }
        xi += h; //Increment the xi for the next iteration 
    }
    return std::optional{ h * sum };
}
