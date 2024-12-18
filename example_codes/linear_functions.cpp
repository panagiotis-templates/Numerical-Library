#include <iostream>
#include <cmath>
#include <vector>
#include <limits>
#include <chrono>
#include <random>

using namespace std::chrono;

using namespace std;


bool inline isEqual(const double a, const double b,
                    const double epsilon = 10e-8) {

  return fabs(a - b) < epsilon;
}

//Algorithm from page 159 dougalis
//Slower implementation but more organized code
double fi(double x,const vector <double> &knot,const int &i) 
{
    //For every other function
    if(i != 0 || i != knot.size()-1)
    {
        if((x  > knot[i-1] || isEqual(x, knot[i-1])) && (x  < knot[i] || isEqual(x, knot[i])))
        {
                //cout << "case2";
                return (x-knot[i-1])/(knot[i]-knot[i-1]);
        }
        else if((x  > knot[i] || isEqual(x, knot[i])) && (x  < knot[i+1] || isEqual(x, knot[i+1])))
        {
                //cout << "case3";
                return -(x-knot[i+1])/(knot[i+1]-knot[i]);
        } 
        else 
        {
            return 0; //No support at that x
        }
    }
    if(i == 0) //First function
    {
        if((x  > knot[0] || isEqual(x, knot[0])) && (x  < knot[1] || isEqual(x, knot[1])))
        {
            //cout << "case1";
            return -(x-knot[1])/(knot[1]-knot[0]);
        }
        else 
        {
            return 0;
        }
    }
    if(i == knot.size()-1) //Last fucntion
    {
        if((x  > knot[knot.size()-2] || isEqual(x, knot[knot.size()-2])) && (x  < knot[knot.size()-1] || isEqual(x, knot[knot.size()-1])))
        {
            //cout << "case4";
            return (x-knot[knot.size()-2])/(knot[knot.size()-1]-knot[knot.size()-2]);
        }
        else 
        {
            return 0;
        }
    }
    
    return numeric_limits<double>::quiet_NaN();
    
    
}

//Faster implementation
double fi2(double x,const vector <double> &knot,const int &i) 
{
    if((x  > knot[i-1] || isEqual(x, knot[i-1])) && (x  < knot[i] || isEqual(x, knot[i])))
    {
            return (x-knot[i-1])/(knot[i]-knot[i-1]);
    }
    else if((x  > knot[i] || isEqual(x, knot[i])) && (x  < knot[i+1] || isEqual(x, knot[i+1])))
    {
            return -(x-knot[i+1])/(knot[i+1]-knot[i]);
    } 
    else 
    {
        return 0; 
    }
}



int main()
{
    
    double a=0,b=1,dx=0.1,n=(b-a)/dx + 1 ,result=0,xi=a,h;
    vector<double> knot(n);
    

   /*  h=doubleDist(rnd); */
    h=1;

    for(int i=0; i<n; i++) //Build the knot vector
    {
        //cout <<xi << "-" << xi+dx << " "<<result << endl;
        knot[i]=xi;
        xi+=dx;
    }
    auto start = high_resolution_clock::now();
    xi=a;
    for(int i=0; i<n; i++) //test 
    {
        //cout <<xi << "-" << xi+dx << " "<<result << endl;
        cout << xi <<" "<< fi(xi,knot,i) << endl;
        xi+=dx * h;
    }
    auto stop = high_resolution_clock::now();
     auto duration = duration_cast<microseconds>(stop - start);

    std::cout <<"f1: "<<duration.count() << std::endl;
     start = high_resolution_clock::now();
    xi=a;
    for(int i=0; i<n; i++) //test 
    {
        //cout <<xi << "-" << xi+dx << " "<<result << endl;
        cout << xi <<" "<< fi2(xi,knot,i) << endl;
        xi+=dx * h;
    }
    stop = high_resolution_clock::now();
      duration = duration_cast<microseconds>(stop - start);

    std::cout <<"f2: "<<duration.count() << std::endl;


    
    return 0;
}