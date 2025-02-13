#include"Header13.h"
int main()
{

    // Define the knots vector
    std::vector<double> kns = { 00,0,0,0.25,0.5,0.75,1,1,1 };
    //std::vector<double> kns = {0.0,  1.0, 2.0, 3.0};
   

    // Define the order of the B-spline
    int order = 3;
    basic_calc_plot(kns, order);
    

    return 0;
}
