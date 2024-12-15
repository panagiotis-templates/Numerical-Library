#pragma once
#include<cmath>
#include<vector>
#include<iostream>
#include<cassert>
#include<iomanip>
#include<numbers>
#include<type_traits>
template<typename _Ty>
requires(std::disjunction_v<std::is_same<_Ty,float>,std::is_same<_Ty,double>,std::is_same<_Ty,long double>>)
[[nodiscard]] inline std::vector<_Ty> gauss_elim(std::vector<std::vector<_Ty>>& A, std::vector<_Ty>& b) {//tes reference edo giati kaneis copy 
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
    std::vector<_Ty> x(n, 0.0); // Initialize the solution vector x with zeros
    for (size_t i = n - 1; i > 0; i--) { // Start from the last equation and move upwards
        _Ty sum = 0.0;
        for (size_t j = i + 1; j < n; j++) { // Loop over elements to the right of the diagonal in the current row
            sum += A[i][j] * x[j]; // Compute the sum of products of coefficients and corresponding elements of x
        }
        assert(A[i][i] != 0);
        x[i] = (A[i][n] - sum) / A[i][i]; // Compute the value of x[i] using back-substitution
    }
    size_t i = 0;
    while (i == 0) {
        _Ty sum = 0.0;
        for (size_t j = i + 1; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        x[i] = (A[i][n] - sum) / A[i][i];
        break;
    }
    return x; // Return the solution vector x
}
template<typename _Ty>
requires(std::disjunction_v<std::is_same<_Ty,float>,std::is_same<_Ty,double>,std::is_same<_Ty,long double>>)
[[nodiscard]] inline  _Ty  f(_Ty x)noexcept {
    return std::sin(x * std::numbers::pi);
    //return pow(x,2) + 3*x + 5 ;
}
template <typename _Ty> 
requires(std::disjunction_v<std::is_same<_Ty,float>,std::is_same<_Ty,double>,std::is_same<_Ty,long double>>)
inline    void print2DVector(const std::vector<std::vector<_Ty>>& vec) {
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
template<typename _Ty>
requires(std::disjunction_v<std::is_same<_Ty,float>,std::is_same<_Ty,double>,std::is_same<_Ty,long double>>)
inline void print_Vector(const std::vector<_Ty>& vec) {
    const size_t&numCols = vec.size();
    for (size_t col = 0; col < numCols; ++col) {
        std::cout << std::setw(12) << std::fixed << std::setprecision(6) << vec[col] << "(" << col
            << ")";
    }
    std::cout << '\n';
}




int main() {
    size_t  n ;// Dimension of the system
    std::cin >> n;
    assert(n > 0);
    //ofstream output;
    //output.open("output.txt", ios::out);
   /* if (!output.is_open())
    {
        cout << "Error opening file";
        exit(-1);
    }*/

    // Build the x,c,b vector
    std::vector< double> x(n), c(n), b(n);

    for (size_t i = 0; i < n; i++)
    {
        x[i] =static_cast<long double>( i);
        b[i] = f(static_cast<long double>(i));
    }


    //Build the left side (A)
    std::vector<std::vector<double>> A(n, std::vector< double>(n, 0));
    for (size_t i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A[i][j] = std::pow(x[i], j);
        }
    }

    c = gauss_elim(A, b);

    std::cout << "x:";
    print_Vector(x);
    std::cout << '\n';


    std::cout << "A matrix:" << '\n';
    print2DVector(A);


    std::cout << "c:" << '\n';
    print_Vector(c);
    std::cout << '\n';


    //Build the interpolation polynomial
    double sum = 0;
    std::vector<double> Pn(n);
    for (size_t i = 0; i < n; i++)
    {
        sum = 0;
        for (int j = 0; j < n; j++)
        {
            sum += c[j] * std::pow(x[i], j);
        }
        Pn[i] = sum;
    }
    std::cout << "x  " << " Pn(x) " << " f(x) " << '\n';
    for (size_t i = 0; i < n; i++)
    {
        std::cout << "x:" << x[i] << " Pn(x):" << Pn[i] << " f(x)" << f(x[i]) << '\n';

    }


    return 0;
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


