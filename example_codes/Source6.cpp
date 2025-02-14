#include"Header6.h"


int main()
{
    //Tests
    std::vector<std::vector<double>> A_dominant = { {5,0,0},{0,5,0},{0,0,5} };
    std::vector<std::vector<double>> A_dominant_not = { {5,0,0},{0,5,0},{10,0,5} };
    std::cout << strict_diagonal_dominace(A_dominant).value() <<'\n';
    std::cout << strict_diagonal_dominace(A_dominant_not).value()<<'\n';


    return 0;
}
