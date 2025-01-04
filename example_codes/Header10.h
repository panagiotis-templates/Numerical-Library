#include<iostream>
#include<type_traits>
#include<cmath>
#include<vector>
#include<optional>
#include<utility>
template<typename _Ty>

inline constexpr bool is_decimal_v = std::disjunction_v<std::is_same<_Ty, float>, std::is_same<_Ty, double>, std::is_same<_Ty, long double>>;

template<typename _Ty>
struct is_decimal :std::bool_constant<is_decimal_v<_Ty>>{};//tag dispatching must support it

template<typename _Ty>
requires(is_decimal_v<_Ty>)
[[nodiscard]] bool inline isEqual(const _Ty& a, const _Ty& b, const _Ty& epsilon = static_cast<_Ty>(10e-10))noexcept {

    return std::abs(a - b) < epsilon;
}

template<typename _Ty>
requires(is_decimal_v<_Ty>)
[[nodiscard]] inline  std::optional<std::vector<_Ty>> Cholesky_method(const std::vector<std::vector<_Ty>>&A, const std::vector<_Ty>& b) {
    if (A.size() <= 0||A.size()!=b.size())return std::nullopt;
    std::vector<std::vector<_Ty>> L(A.size(), std::vector<_Ty>(A.size())), Lt(A.size(), std::vector<_Ty>(A.size()));
    std::vector<_Ty> x(L.size()), y(L.size());
    for (size_t i = 0; i < A.size();i++) {
        if (A.size() != A[i].size()) {
            std::cout << "Not a square matrix";
            return std::nullopt;
        }
    }
   
    
    std::vector<std::vector<_Ty>> At(A.size(), std::vector<_Ty>(A.size()));
    //Transpose the A into At
    for (size_t i = 0; i < At.size(); i++) {
        for (size_t j = 0; j < At[0].size(); j++) {
            At[i][j] = A[j][i];
        }
    }
    if (A != At) {//Checks for symmetric A ,A=At
        std::cout << "Not symetric A";
        return std::nullopt;
    }
    //Decomposition of A into lower triangular based on the formula (3.24) page 95 Numerical Analysis Dougalis
    _Ty sum = 0, sum2 = 0, sum_x = 0, sum_k = 0;
    for (size_t i = 0; i < A.size(); i++) {
        for (size_t j = 0; j < i; j++) {
            for (size_t k = 0; k < j; k++) {
                sum += L[i][k] * L[j][k];
            }
           
            L[i][j] = (A[i][j] - sum) / L[j][j];
            sum = 0;
        }
        for (size_t g = 0; g < i; g++) {
            sum2 += L[i][g]*L[i][g];
        }
        if (A[i][i] - sum2 < 0)return std::nullopt;
        L[i][i] = std::sqrt(A[i][i] - sum2);
        if (L[i][i] < 0) {//Check for positive definite of A based on the decomposition of L Lt

            std::cout << "A is not positive definite";
            return std::nullopt;
        }
        sum2 = 0;
    }
    for (size_t i = 0; i < Lt.size(); i++) { //Build the L transpose Lt
        for (size_t j = 0; j < Lt[0].size(); j++) {
            Lt[i][j] = L[j][i];
        }
    }
    for (size_t i = 0; i < L.size(); i++) { //Foward substitution L y= b
        for (size_t j = 0; j < i; j++) {
            sum_x += L[i][j] * y[j];
        }
        if (isEqual(L[i][i], 0.0)) {
            return std::nullopt;
        }
        y[i] = (b[i] - sum_x) / L[i][i];
        sum_x = 0;
    }
    x[x.size() - 1] = y[y.size() - 1] / Lt[Lt.size() - 1][Lt.size() - 1];
    for (size_t k = L.size() - 2; k > 0; k--) {//Backwards substitution Lt x = y
        for (size_t j = k + 1; j <= Lt.size() - 1; j++) {
            sum_k += Lt[k][j] * x[j];
        }
        if (isEqual(Lt[k][k], 0.0))return std::nullopt;
        x[k] = (y[k] - sum_k) / Lt[k][k];
        sum_k = 0;
    }
    size_t k = 0;
    while (k == 0) {
        for (size_t j = k + 1; j <= Lt.size() - 1; j++) {
            sum_k += Lt[k][j] * x[j];
        }
        if (isEqual(Lt[k][k], 0.0))return std::nullopt;
        x[k] = (y[k] - sum_k) / Lt[k][k];
        sum_k = 0;
        break;
    }
    return std::optional{ x };
}
