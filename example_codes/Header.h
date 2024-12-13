#pragma once
#include<cmath>
#include<vector>
#include<iostream>
#include<cassert>
#include<concepts>
template<std::floating_point _Floating>
inline[[nodiscard]] _Floating  f(_Floating x) noexcept 
{
    return exp(2 * x);
}
template<std::floating_point _Floating>
inline [[nodiscard]] _Floating ddf(_Floating x) noexcept 
{
    return 4 * exp(2 * x);
}
template<std::floating_point _Floating>
inline void finite_difference(_Floating a, _Floating b, _Floating h) {
    _Floating xi{}, res{}, hsq{};
    size_t n{ static_cast<size_t>((b - a) / h + 1) };
    xi = a;
    hsq= std::pow(h, 2);
    for (size_t i = 0; i < n; i++) {
        res = (f(xi - h) - 2 * f(xi) + f(xi + h)) / hsq;
        std::cout << xi << " " << res << " " << ddf(xi) << " " << fabs(ddf(xi) - res) << '\n';
        xi +=  h;
    }
    
}
template<std::floating_point _Floating>
[[nodiscard]] std::vector<_Floating> gauss_elim(std::vector<std::vector<_Floating>>& A, std::vector<_Floating>& b) {//tes reference edo giati kaneis copy 
    const size_t &n = A.size(); // Number of equations (rows)
    assert(n > 0);
    // Augment the coefficient matrix with the right-hand side vector
    for (size_t  i = 0; i < n; i++) {
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
    for (size_t  i = 0; i < n; i++) { // Loop over each row (equation)
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
    if (A[A.size() - 1][A.size() - 1]==0) {
        std::cout << "error\n";
        std::exit(-1);
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
    std::vector<_Floating> x(n, 0.0); // Initialize the solution vector x with zeros
    for (size_t i = n - 1; i > 0; i--) { // Start from the last equation and move upwards
        _Floating sum = 0.0;
        for (size_t j = i + 1; j < n; j++) { // Loop over elements to the right of the diagonal in the current row
            sum += A[i][j] * x[j]; // Compute the sum of products of coefficients and corresponding elements of x
        }
        assert(A[i][i] != 0);
        x[i] = (A[i][n] - sum) / A[i][i]; // Compute the value of x[i] using back-substitution
    }
    size_t i = 0;
    while (i == 0) {
        _Floating sum = 0.0;
        for (size_t j = i + 1; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        x[i] = (A[i][n] - sum) / A[i][i];
        break;
    }
    return x; // Return the solution vector x
}
template<std::floating_point _Floating>
[[nodiscard]] _Floating  f(_Floating x)noexcept {
    return std::sin(x * std::M_PI);
    //return pow(x,2) + 3*x + 5 ;
}
template <std::floating_point _Floating>  void print2DVector(const std::vector<std::vector<_Floating>>& vec) {
    const size_t &numRows = vec.size();
    const size_t &numCols = vec[0].size();
    assert(vec.size() > 0&&vec[0].size());
    // Print column headers
    std::cout << std::setw(12) << " ";
    for (size_t col = 0; col <numCols; ++col) {
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
template<std::floating_point _Floating>
inline void print_Vector(const std::vector<_Floating>& vec) {
    const size_t&numCols = vec.size();
    for (size_t col = 0; col < numCols; ++col) {
        std::cout << std::setw(12) << std::fixed << std::setprecision(6) << vec[col] << "(" << col
            << ")";
    }
    std::cout << '\n';
}
//void polyonomial_Interpolation(std::ofstream &output,size_t size) {
//    double sum = 0;
//    double sum = 0;
//    std::vector<double> Pn(n);
//    for (int i = 0; i < n; i++)
//    {
//        sum = 0;
//        for (int j = 0; j < n; j++)
//        {
//            sum += c[j] * A[i][j];
//        }
//        Pn[i] = sum;
//    }
//    output << "x  " << " Pn(x) " << " f(x) " << '\n';
//    for (int i = 0; i < n; i++)
//    {
//        std::cout << "x:" << x[i] << " Pn(x):" << Pn[i] << " f(x)" << f(x[i]) << '\n';
//        output << x[i] << "  " << Pn[i] << "  " << f(x[i]) << '\n';
//    }
//
//
//    
//
//}


