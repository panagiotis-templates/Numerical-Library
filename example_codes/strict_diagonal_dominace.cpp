#include <cmath>
#include <iostream>
#include <vector>

using namespace std;


bool inline isEqual(const double a, const double b,
                    const double epsilon = 10e-8) {

  return fabs(a - b) < epsilon;
}




bool strict_diagonal_dominace(const vector<vector<double>>  &A) //Pass it as refrence so you dont need copy,also add const if you dont goining to modify
{
    double sum=0;
    for(int i=0; i<A.size(); i++)
    {
        sum=0;//Reset the sum for new row
        for(int j=0; j<A[0].size(); j++)
        {
            if(i!=j){sum+= A[i][j];} //Sum the values except the diagonal
            else{continue;}
        }
        if(fabs(A[i][i]) < fabs(sum) || isEqual(fabs(A[i][i]) , fabs(sum))) //Check if one diagonal element is smaller than the sum
        {
            return false;
        }
    }
    return true;
}


int main()
{
    //Tests
    vector<vector<double>> A_dominant={{5,0,0},{0,5,0},{0,0,5}};
    vector<vector<double>> A_dominant_not ={{5,0,0},{0,5,0},{10,0,5}};
    cout << strict_diagonal_dominace(A_dominant) << endl;
    cout << strict_diagonal_dominace(A_dominant_not);


    return 0;
}