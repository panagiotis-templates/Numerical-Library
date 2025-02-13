#include"Header12.h"







int main()
{
   
    using namespace panagiotis;
    finite_diff(1.0, 3.0, 0.1, [](double x)->double {return std::sin(2 * std::numbers::pi * x); });


    return 0;
}
