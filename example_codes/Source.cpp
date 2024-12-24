#include"Header.h"

int main() {
    size_t  n;// Dimension of the system
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
    std::vector<double> x(n), b(n);

    for (size_t i = 0; i < n; i++)
    {
        x[i] = static_cast<double>(i);
        b[i] = f(static_cast<double>(i));
    }


    //Build the left side (A)
    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0));
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            A[i][j] = std::pow(x[i], j);
        }
    }

    auto c = gauss_elim(A, b);

    std::cout << "x:";
    print_Vector(x);
    std::cout << '\n';


    std::cout << "A matrix:" << '\n';
    print2DVector(A);
    if (c.has_value()) {

        std::cout << "c:" << '\n';
        print_Vector(c.value());
        std::cout << '\n';


        //Build the interpolation polynomial
        double sum = 0;
        std::vector<double> Pn(n);
        for (size_t i = 0; i < n; i++)
        {
            sum = 0;
            for (size_t j = 0; j < n; j++)
            {
                sum += c->at(j) * std::pow(x[i], j);
            }
            Pn[i] = sum;
        }
        std::cout << "x  " << " Pn(x) " << " f(x) " << '\n';
        for (size_t i = 0; i < n; i++)
        {
            std::cout << "x:" << x[i] << " Pn(x):" << Pn[i] << " f(x)" << f(x[i]) << '\n';

        }

    }
    return 0;
}
    return 0;
}
