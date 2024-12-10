#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>



using namespace std;


bool inline isEqual(const double a, const double b,
                    const double epsilon = 10e-8) {

  return fabs(a - b) < epsilon;
}

void print2DVector(const vector<vector<double>>& vec) {
    

    const int numRows = vec.size();
    const int numCols = vec[0].size();

    // Print column headers
    cout << setw(12) << " ";
    for (int col = 0; col < numCols; ++col) {
        cout << setw(12) << col;
    }
    cout << '\n';

    // Print row headers and vector contents
    for (int row = 0; row < numRows; ++row) {
        cout << setw(12) << row;
        for (int col = 0; col < numCols; ++col) {
            cout << setw(12) << fixed << setprecision(6) << vec[row][col] << "(" << row << "," << col << ")";
        }
        cout << '\n';
    }
}




vector<double> Cholesky_method(const vector<vector<double>> A,const vector<double> b){
    vector<vector<double>> L(A.size(),vector<double> (A.size())),Lt(A.size(),vector<double> (A.size()));
    vector<double> x(L.size(),numeric_limits<double>::quiet_NaN()),y(L.size());
    if(A.size() != A[0].size()){ //Check the dimension of the matrix A and b
        cout << "Not a square matrix";
        return x;}
    vector<vector<double>> At(A.size(),vector<double>(A.size()));
    //Transpose the A into At
    for(int i=0; i<At.size(); i++){
        for(int j=0; j<At[0].size(); j++){
            At[i][j] = A[j][i];}}
    if(A != At) {//Checks for symmetric A ,A=At
        cout << "Not symetric A";
        return x;}
     //Decomposition of A into lower triangular based on the formula (3.24) page 95 Numerical Analysis Dougalis
    double sum=0,sum2=0,sum_x=0,sum_k=0;
    for(int i=0; i<A.size(); i++) {
        for(int j=0; j<i; j++){
            for(int k=0; k<j; k++){
                sum+=L[i][k] * L[j][k];}
            L[i][j]=(A[i][j]-sum)/L[j][j];
            sum=0;}
        for(int g=0; g<i; g++){
            sum2+=pow(L[i][g],2);}
        L[i][i]=sqrt(A[i][i] - sum2);
        if(L[i][i] < 0) {//Check for positive definite of A based on the decomposition of L Lt
        
            cout << "A is not positive definite" ;
            return x;}
        sum2=0;}
    for(int i=0; i<Lt.size(); i++){ //Build the L transpose Lt
        for(int j=0; j<Lt[0].size(); j++){
            Lt[i][j] = L[j][i]; }}
     for(int i=0; i<L.size(); i++) { //Foward substitution L y= b
        for(int j=0; j<i; j++){
            sum_x+= L[i][j] * y[j];}
        y[i]= (b[i]-sum_x) / L[i][i];
        sum_x=0;}
    x[x.size()- 1] = y[y.size()- 1]/Lt[Lt.size()- 1][Lt.size()- 1];
    for(int k=L.size()-2; k>=0; k--) {//Backwards substitution Lt x = y
        for(int j=k+1; j<=Lt.size()-1; j++){
            sum_k+= Lt[k][j] * x[j];}
        x[k]= (y[k]-sum_k) / Lt[k][k];
        sum_k=0;}
    return x;}


int main()
{
    //vector<vector<double>> A={{25,15,-5},{15,18,0},{-5,0,11}};
    vector<vector<double>> A={{1,1,2},{1,2,2},{2,2,8}};
    vector<double> b={1,3,-2};
     vector<double> x(b.size());
    //Transpose the A into At
    x=Cholesky_method(A,b);

    
    cout <<"x:" << endl;
    for(double i:x)
    {
        cout << i << " " << endl ;
    } 
   
    print2DVector(A);
    return 0;
}