#include"Header10.h"
int main()
{
   
    //using namespace panagiotis;
    //vector<vector<double>> A={{25,15,-5},{15,18,0},{-5,0,11}};
    std::vector<std::vector<double>> A = { {1,1,2},{1,2,2},{2,2,8} };
    std::vector<double> b = { 1,3,-2 };
    std::optional<std::vector<double>> x(b.size());
    //Transpose the A into At
    x = Cholesky_method(A, b);
    
    
    if (x.has_value()) {

        std::cout << "x:" << std::endl;
        
        for (double i : x.value())
        {
            std::cout << i << " " << std::endl;
        }
        print2DVector(A);
    }

  
    return 0;
}
