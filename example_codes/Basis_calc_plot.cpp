#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <fstream>
#include <iomanip>
#include <string>

using namespace std;


// Helper function 
double gdiv(double a,double b)
{
    if(a==0.0&&b==0.0)
    {
        return(0.0);
    }    
    else
    { 
        return(a/b);
    }
}


// Function to calculate the B-spline value
/*
    i is the index of the B-spline basis function in the B-spline basis set. 
    This is a zero-based index, so the first B-spline basis function in the set has index 0.

    ord is the order of the B-spline basis function. 
    This is also known as the degree of the B-spline basis function. 
    The order is an integer greater than or equal to 1 that determines the number of 
    knots used to define the B-spline basis function.

    x is the value at which the B-spline basis function should be evaluated.

    nk is the number of knots in the knot vector kns.

    kns is a vector of knots that define the B-spline basis function. 
    The knots are sorted in ascending order, 
    and the B-spline basis function is defined over the range between the first and last knot. 
    The size of kns should be equal to nk.

*/
double bsp(int i, int ord, double x, int nk, vector<double> kns){
  // Check for illegal value of i
  if(i<0 || i>nk-ord-1){
    cout<<"illegal i value: i="<<i<<"; nk-ord="<<nk<<"-"<<ord<<"="<<nk-ord<<endl;
    return numeric_limits<double>::quiet_NaN();
  }

  // Return 0 if x is outside the interval defined by the knots
  if(x<kns[i] || x>kns[i+ord]) return 0.0;

  // Remove repeated knots
  int k=nk-1; 
  while(kns[k]==kns[k-1]) k--; 
  k--; 

  // If ord is 1, return 1 if x is within the interval defined by the knots
  if(ord == 1){
    if(i != k){
      return((kns[i]<=x && x<kns[i+1]) ? 1.0 : 0.0);
    }else{
      if(i == k){
        return((kns[i]<=x && x<=kns[i+1]) ? 1.0 : 0.0);
      }else return numeric_limits<double>::quiet_NaN();;
    }
  }
  // If ord is greater than 1, return the sum of two recursive calls
  else{
    return(gdiv((x-kns[i])*bsp(i, ord-1, x, nk, kns), kns[i+ord-1]-kns[i]) +
           gdiv((kns[i+ord]-x)*bsp(i+1, ord-1, x, nk, kns), kns[i+ord]-kns[i+1])
    );
  }
}



int main()
{
    // Define the knots vector
    vector<double> kns = {00,0,0,0.25,0.5,0.75,1,1,1};
    //vector<double> kns = {0.0,  1.0, 2.0, 3.0};

    // Check if the vector is empty(Used for the .back() method)
    if (kns.empty()) {
      // Handle the case where the vector is empty
      cout << "Error empty knot vector";
      return 1 ;
    } 

    // Define the order of the B-spline
    int order = 3;
   
    
    // Create an output file stream object
    ofstream outfile;

    // Open a file for writing
    outfile.open("data.txt", ios::out);

    // Check if the file opened successfully
    if (!outfile.is_open()) {
      // Handle the error
      cout << "Error opening the file";
      return 1;
    }

    // Create an empty 2D vector to store the data
    vector<vector<double>> data;

    // Create a vector of labels
    vector<string> labels={"Xi"};

    // Create the labels
    for(int j=0; j<kns.size()-order; j++)
    {
     
      labels.push_back("B_{" + to_string(j+1) + "," + to_string(order-1)+"}"); // order-1 and j+1 because the index is zero based
    }

    // Iterate over the vector of labels and prints them to the file
    for (const string& s : labels) 
    {
       outfile << s << "   "; 
    }
    outfile << endl;

    for(double i=0.0; i<= kns.back(); i+=0.001)
    {
      // Create a new row vector
      vector<double> row;
      
      row.push_back(i);
      // Starts from P=1 (Order)
      for(int j=0; j<kns.size()-order; j++)
      {
        // Calculate the values for each column
        
        row.push_back(bsp(j,order,i,kns.size(),kns));
      }  
      // Add the row to the data vector
      data.push_back(row);
    }

    // Print the contents of the data vector into the file
    for (const auto& row : data) {
      for (const auto& value : row) 
      {
        outfile <<fixed << setprecision(4) << value << " ";
      }
      outfile << endl;
    }

    // Close the file
    outfile.close();

    
    return 0;
}
