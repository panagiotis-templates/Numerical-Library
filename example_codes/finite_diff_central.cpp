#include <iostream>
#include <cmath>
#include <iomanip>
#include <limits>
#include <vector>
#include <fstream>

using namespace std;

//Utility fucntions
template <typename T>  void print2DVector(const vector<vector<T>> &vec,ofstream &out) {

  const int numRows = vec.size();
  const int numCols = vec[0].size();

  // Print column headers
  out << setw(12) << " ";
  for (int col = 0; col < numCols; ++col) {
    out << setw(12) << col;
  }
  out << '\n';

  // Print row headers and vector contents
  for (int row = 0; row < numRows; ++row) {
    out << setw(12) << row;
    for (int col = 0; col < numCols; ++col) {
      out << setw(12) << fixed << setprecision(6) << vec[row][col] << "("
           << row << "," << col << ")";
    }
    out << '\n';
  }
}

void print_Vector(const vector<double> &vec,ofstream &out) {
  const int numCols = vec.size();
  for (int col = 0; col < numCols; ++col) {
    out << setw(12) << fixed << setprecision(6) << vec[col] << "(" << col
         << ")" << endl;
  }
  out << endl;
}





vector<double> gauss_elim(vector<vector<double>> A, vector<double> b) {
    int n = A.size(); // Number of equations (rows)

    // Augment the coefficient matrix with the right-hand side vector
    for (int i = 0; i < n; i++) {
        A[i].push_back(b[i]); // Appending the corresponding element from b to each row of A
    }

    /* for(int i=0; i<n; i++)
    {
        for(int j=0; j<A[0].size(); j++)
        {
            cout << A[i][j] << " ";
        }
        cout << endl;
    } */

    cout << endl;

    // Perform Gaussian elimination
    for (int i = 0; i < n; i++) { // Loop over each row (equation)
        // Find the row with the maximum absolute value in the ith column and swap rows
        int max_row = i;
        for (int j = i + 1; j < n; j++) {
            if (abs(A[j][i]) > abs(A[max_row][i])) {
                max_row = j; // Update the index of the row with the maximum absolute value
            }
        }
        
        swap(A[i], A[max_row]); // Swap the current row with the row with the maximum absolute value

        // Perform row operations to eliminate coefficients below the pivot element
        for (int j = i + 1; j < n; j++) { // Loop over rows below the pivot row
            double factor = A[j][i] / A[i][i]; // Compute the factor by which the pivot row will be multiplied
            for (int k = i; k < n + 1; k++) { // Loop over columns including the augmented column
                A[j][k] -= factor * A[i][k]; // Perform row operation to eliminate coefficients below the pivot element
            }
        }
    }

   /*  for(int i=0; i<n; i++)
    {
        for(int j=0; j<A[0].size(); j++)
        {
            cout << A[i][j] << " ";
        }
        cout << endl;
    } */

    // Back-substitution to solve for x
    vector<double> x(n, 0.0); // Initialize the solution vector x with zeros
    for (int i = n - 1; i >= 0; i--) { // Start from the last equation and move upwards
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) { // Loop over elements to the right of the diagonal in the current row
            sum += A[i][j] * x[j]; // Compute the sum of products of coefficients and corresponding elements of x
        }
        x[i] = (A[i][n] - sum) / A[i][i]; // Compute the value of x[i] using back-substitution
    }

    return x; // Return the solution vector x
}






double f(double t,double q)
{
    return (-0.25 * exp(-2 * t)) * (-(pow(M_PI,2)-16) *sin(M_PI_2 * t)- 8 * M_PI  * cos(M_PI_2 * t))+ q * exp(-2 *t) * sin(M_PI_2 * t);
    //return pow(M_PI,2) * sin(M_PI * t) + sin(M_PI * t); //Promblem 2
}
double y(double t)
{
    return exp(-2 * t) * sin(M_PI_2 * t);
    //return sin(M_PI * t); //Problem 2
}



int main()
{
    ofstream output;
    char filename[] = "0output.txt";
    vector<double> h={0.1,0.01,0.05,0.02,0.01}; 
    //vector <double> h={0.1};
    double xi,result,hsq,max= -1,sum=0; //Utility variables
    double a=0,b=1,y_a=0,y_b=exp(-2),q=2;
    //double a=0,b=1,y_a=0,y_b=0,q=1; //Problem 2
    int n;
    
    for(int j=0; j<h.size(); j++)
    {
        max= -1,sum=0;
        filename[0]=48+j; //Change the file name based on number(dont forget it is char)
        output.open(filename); 
        n=(b-a)/h[j] + 1; //Calculate the number of elements inside the partition based on the current h
        vector<double> U(n,numeric_limits<double>::quiet_NaN()),F(n); //Solution vector U,Right hand side F
        vector<vector<double>> A(n,vector<double> (n,0)); 
        xi=a; //Reset the partition count to the start of the interval
        //cout << "For step size of: " << h[j] << endl;
        
        hsq=h[j] * h[j]; //Single time calculation of h squared
        for(int i=0; i<n; i++)
        {
           //Build the right hand side 
           F[i]=f(xi,q) * hsq;
           //Build the left hand side 
           for(int g=0; g<n; g++)
           {
            if(i == g){A[i][g]= 2 + hsq * q;}
            else if (i == g-1 || i == g+1) {A[i][g]=-1;}
           }
           xi+=h[j]; //Increment the xi for the next iteration 
        }
        //Boundary Conditions

        F[0]=y_a  ;
        F[n-1]=y_b  ; 
    
        //Solve the system 
        U=gauss_elim(A, F);
        xi=a; //Reset
        
        for(int i=0; i<n; i++) 
        {
            sum += sqrt(pow(fabs(U[i]-y(xi)),2));
            output << xi << " " << y(xi) << " "<<U[i] << endl;
            if(fabs(U[i]-y(xi)) > max) //Find max difference
            {
                max= fabs(U[i]-y(xi));
            }
            xi+=h[j]; //Increment the xi for the next iteration 
        }
        //output << "Max difference:" << max<<endl;
        output << "Error Estimate " << h[j] * sum << endl; 

        output.close();
    }
    return 0;
}