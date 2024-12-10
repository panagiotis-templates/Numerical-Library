#include <iostream>
#include <cmath>
#include <limits>



using namespace std;


 double f(double x)
{
    return pow(x,3) -  2*x - 5;
    //return x * x + 2 * x;
}


bool inline isEqual(const double a, const double b,
                    const double epsilon = 10e-10) {

  return fabs(a - b) < epsilon;
}


double dicection(double a,double b,double (*f)(double),double e=10e-10)
{
    double d=b-a,c; //[a,b] interval ,d interval span ,e precision
    while(true)
    {
        d *= 0.5;
        if(signbit(f(a)) == signbit(f(b)) || d < e || isEqual(d, e))
        {
            cout << "Not found in interval: "<< "[" << a << "," << b <<"]";
            return numeric_limits<double>::quiet_NaN();
        }
        c=a+(b-a) * 0.5;
        cout << a << " " << b << " " << c << " " << d << " "<< f(c)<<endl;
        if(fabs(f(c)- 0) < 10e-6) //Smaller precision for the final result 
        {
            cout << "Found the solution: " << c;
            return c;
        }
        if(signbit(f(c)) == signbit(f(a))) //if(f(c) * f(a) > 0) is going to produce overflow if f(a) or f(c) is  small number
        {
            a = c;
        }
        else
        {
            b = c;
        }
    }
}


//Algorithm based on the page 35 from dougalis
int main()
{
    double a=2,b=3,d=b-a,e=10e-10,c; //[a,b] interval ,d interval span ,e precision 
    dicection(a, b,&f) ;
    return 0; 
}