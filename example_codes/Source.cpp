#include"Header.h"


int main() {
    size_t size;
    std::cin >> size; // Dimension of the system
    assert(size > 0);
    /*  std::ofstream output;
    output.open("output.txt", std::ios::out);
    if (!output.is_open())
    {
        std::cout << "Error opening file";
        std::exit(-1);
    }*/

    // Build the x,c,b vector
    std::vector<double> x(size), c(size), b(size);

    for (size_t i = 0; i < size; i++)
    {
        x[i] =static_cast<double>(i);
        b[i] = f(static_cast<double>(i));
    }
   

    //Build the left side (A)
    std::vector<std::vector<double>> A(size, std::vector<double>(size, 0));
    for (size_t i = 0; i < size; i++)
    {
        for (size_t j = 0; j < size; j++)
        {
            A[i][j] = std::pow(x[i], j);
        }
    }
     // std::sin(5.5);
     c=gauss_elim(A, b);
     
    //print_Vector(c);
    std::cout << "x:";
    print_Vector(x);
    std::cout << '\n';


    std::cout << "A matrix:" << '\n';
    print2DVector(A);
    //system("pause");

    std::cout << "c:" << '\n';
    print_Vector(c);
    std::cout << '\n';


    //Build the interpolation polynomial
    double sum = 0;
    std::vector<double> Pn(size);
    for (size_t i = 0; i < size; i++)
    {
        sum = 0;
        for (size_t j = 0; j < size; j++)
        {
            sum += c[j] * std::pow(x[i], j);
        }
        Pn[i] = sum;
    }
    std::cout<< "x  " << " Pn(x) " << " f(x) " << '\n';
    for (size_t i = 0; i < size; i++)
    {
        std::cout << "x:" << x[i] << " Pn(x):" << Pn[i] << " f(x)" << f(x[i]) << '\n';
        //std::cout << x[i] << "  " << Pn[i] << "  " << f(x[i]) << '\n';
    }


    return 0;
}