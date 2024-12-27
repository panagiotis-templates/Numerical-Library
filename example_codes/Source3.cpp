#include"Header3.h"

int main()
{
    double a = 2, b = 3, d = b - a, e = 10e-10; //[a,b] interval ,d interval span ,e precision 
 if (auto c = dicection(a, b, [](const double& x) {return std::pow(x, 3) - 2 * x - 5; }); c.has_value()) {
     std::cout << *c << '\n';
 }
 return 0;
    

}
