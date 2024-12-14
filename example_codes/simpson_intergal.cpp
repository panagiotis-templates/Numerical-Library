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
        result+= (dx/6) * (f(xi)+ 4 * f((xi+(xi+dx)) * 0.5) +  f(xi+dx)); 
        cout <<xi << "-" << xi+dx << " "<<result << endl;
        xi+=dx;
    }
    cout << "result: "<<result; 
    return 0;
}