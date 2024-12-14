#include <iostream>
#include <cmath>

using namespace std;

double f(double x)
{
    //return pow(2,x) + 5 * x +3;
    return exp(2*x);
    //return sin(2 * M_PI *x);
}


int main()
{
    double a=0,b=1,dx=0.1,n=(b-a)/dx ,result=0,xi=a;
    for(int i=0; i<n; i++)
    {
        result+=(f(xi+dx) + f(xi)) * (dx * 0.5);  // More calculation but better for debug
        //result+=(f(xi+dx) + f(xi)) //Alternative way less calculations
        cout <<xi << "-" << xi+dx << " "<<result << endl;
        xi+=dx;
    }
    //result *=(dx * 0.5); //Alternative way less calculations
    cout << "result: "<<result; 
    return 0;
}