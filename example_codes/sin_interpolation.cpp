#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;


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



double f(double x) { return sin(x * M_PI_2); 
//return pow(x,2) + 3*x + 5 ;
}

template <typename T>  void print2DVector(const vector<vector<T>> &vec) {

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
      cout << setw(12) << fixed << setprecision(6) << vec[row][col] << "("
           << row << "," << col << ")";
    }
    cout << '\n';
  }
}

void print_Vector(const vector<double> &vec) {
  const int numCols = vec.size();
  for (int col = 0; col < numCols; ++col) {
    cout << setw(12) << fixed << setprecision(6) << vec[col] << "(" << col
         << ")";
  }
  cout << '\n';
}

int main() {
  int n = 3; // Dimension of the system
  ofstream output;
  output.open("output.txt",ios::out);
  if(!output.is_open())
  {
    cout << "Error opening file";
    exit(-1);
  }
  
  // Build the x,c,b vector
  vector<double> x(n),c(n),b(n);

   for (int i = 0; i < n; i++) 
   {
    x[i]=i;
    b[i]=f(i);
   }  


 //Build the left side (A)
 vector<vector<double>> A(n,vector<double> (n,0)); 
  for(int i=0; i<n; i++)
  {
    for(int j=0; j<n; j++)
    {
        A[i][j]=pow(x[i],j);
    }
  }

  c=gauss_elim(A, b);

  cout << "x:";
  print_Vector(x);
  cout << endl;


  cout << "A matrix:" << endl;
  print2DVector(A);


   cout << "c:" << endl;
  print_Vector(c);
  cout << endl;


  //Build the interpolation polynomial
  double sum=0;
  vector<double> Pn(n);
  for(int i=0; i<n; i++)
  {
    sum=0;
    for(int j=0; j<n; j++)
    {
        sum+= c[j] * A[i][j]; 
    }
    Pn[i]=sum;
  }
   output<< "x  " << " Pn(x) " << " f(x) " <<endl;
  for(int i=0; i<n; i++)
  {
    cout << "x:" <<x[i] << " Pn(x):" << Pn[i] << " f(x)" << f(x[i])<<endl;
    output <<x[i] << "  " << Pn[i] << "  " << f(x[i])<<endl;
  }


  return 0;
}