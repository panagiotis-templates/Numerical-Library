#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;



double f(double x)
{
    //return exp(2*x);
    return sin(2*M_PI *x);
}
/* 
double ddf(double x) // DEBUG
{
    return 4 * exp(2 * x);
}
 */
double derivative(double (*f)(double), double x0, int order) //Numerical differation
{
         const double delta = 1.0e-6;
         double x1 = x0 - delta;
         double x2 = x0 + delta;
         if (order == 1) {
                  double y1 = f(x1);
                  double y2 = f(x2);

                  return (y2 - y1) / (x2 - x1);
          } else {
                  double y1 = derivative(f, x1, order - 1);
                  double y2 = derivative(f, x2, order - 1);

                  return (y2 - y1) / (x2 - x1);
          }
 }


int main()
{
    ofstream output;
    char filename[] = "0output.txt";
    vector<double> h={0.1,0.01,0.001}; 
    //vector <double> h={0.1};
    double a=1,b=3,n,xi,result,hsq;
    for(int j=0; j<h.size(); j++)
    {
        filename[0]=48+j; //Change the file name based on number(dont forget it is char)
        output.open(filename); 
        n=(b-a)/h[j] + 1; //Calculate the number of elements inside the partition based on the current h
        xi=a; //Reset the partition count to the start of the interval
        //cout << "For step size of: " << h[j] << endl;
        output << "xi  aprox real error " << endl;
        hsq=h[j] * h[j]; //Single time calculation of h squared
        for(int i=0; i<n; i++)
        {
           result=(f(xi-h[j])- 2* f(xi)+f(xi+h[j]))/hsq; //Finite difference
           output << xi << " " << result << " " << derivative(f,xi,2) << " "<<fabs(derivative(f,xi,2) -result) << endl;
           xi+=h[j]; //Increment the xi for the next iteration 
        }
        output.close();
    }
    return 0;
}