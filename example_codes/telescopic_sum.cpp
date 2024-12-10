#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;


double better_sum(double n)
{
    double sum=1./(pow(n,2)+n);
    for(int i=1; i<n; i++)
    {
        sum+=1/((n-i)*(n-i+1));
    }
    sum+=1;
    return sum;
}


int main()
{
    double sum=1,n=9999;
    for(int i=1; i<=n; i++)
    {
        sum+= 1/ (pow(i,2) + i);
        //cout << i << " " << sum << endl; 
    }
    double sum_b=better_sum(n);
    cout << setprecision(30)<< fixed <<sum  << " better sum :" << sum_b << endl;
    /* cout << log(pow(640320,3)+ 744) /sqrt(163) << endl; */
    return 0;
}