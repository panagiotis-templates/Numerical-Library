#include"Header7.h"
int main()
{
    
    double sum = 1, n = 9999;
    for (int i = 1; i <= n; i++)
    {
        sum += 1 / (pow(i, 2) + i);
        //cout << i << " " << sum << endl; 
    }
    auto sum_b = better_sum(n);
    if(sum_b){
    std::cout << std::setprecision(30) << std::fixed << sum << " better sum :" << sum_b.value() << '\n';
    }
    /* cout << log(pow(640320,3)+ 744) /sqrt(163) << endl; */
    return 0;
}
