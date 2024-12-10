#include <cmath>
#include <iostream>
#include <vector>


using namespace std;



double f(double x)
{
    return exp(2*x);
}

double ddf(double x)
{
    return 4 * exp(2 * x);
}

int main()
{
     vector<double> h={0.1,0.01,0.001}; 
    //vector <double> h={0.1};
    double a=1,b=3,n,xi,result,hsq;
    for(int j=0; j<h.size(); j++)
    {
        n=(b-a)/h[j] + 1; //Calculate the number of elements inside the partition based on the current h
        xi=a; //Reset the partition count to the start of the interval
        cout << "For step size of: " << h[j] << endl;
        hsq=h[j] * h[j]; //Single time calculation of h squared
        for(int i=0; i<n; i++)
        {
           result=(f(xi-h[j])- 2* f(xi)+f(xi+h[j]))/hsq; //Finite difference
           cout << xi << " " << result << " " << ddf(xi) << " "<<fabs(ddf(xi)-result) << endl;
           xi+=h[j]; //Increment the ti for the next iteration 
        }
    }
    return 0;
}