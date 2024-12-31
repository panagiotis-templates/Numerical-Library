#include "Header9.h"
int main()
{
   
    using namespace panagiotis;
    double a = 0, b = 1, dx = 0.1, n = (b - a) / dx + 1, result = 0, xi = a, h;
    std::vector<double> knot(static_cast<size_t>(n));


    /*  h=doubleDist(rnd); */
    h = 1;

    for (size_t i = 0; i < n; i++) //Build the knot vector
    {
        //cout <<xi << "-" << xi+dx << " "<<result << endl;
        knot[i] = xi;
        xi += dx;
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    xi = a;
    for (size_t i = 0; i < n; i++) //test 
    {
        //cout <<xi << "-" << xi+dx << " "<<result << endl;
        if(auto value= fi(xi, knot, i);value)
        std::cout << xi << " " << *value << '\n';
        xi += dx * h;
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    std::cout << "f1: " << duration.count() << "us\n";
    start = std::chrono::high_resolution_clock::now();
    xi = a;
    for (size_t i = 0; i < n; i++) //test 
    {
        //cout <<xi << "-" << xi+dx << " "<<result << endl;
        if(auto value= fi2(xi, knot, i);value)
        std::cout << xi << " " << *value << '\n';
        xi += dx * h;
    }
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "f2: " << duration.count() << "us\n";



    return 0;
}
