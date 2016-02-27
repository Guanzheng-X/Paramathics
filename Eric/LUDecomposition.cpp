#include<iostream>
#include<Eigen/Core>
#include<Eigen/Sparse>
#include<Eigen/SparseLU>
#include<Eigen/OrderingMethods>
#include<fstream>
#include<vector>

using namespace Eigen;
using namespace std;

int main(int argc, char *argv[])
{
   /******************************Construct a Sparse Matrix***********************/
   int rows, cols, nonezero;
   
   typedef Eigen::Triplet<double> T;

   vector<T> tripletList;

   ifstream data;

   data.open("b1_ss.mtx");   //"testMatrix.txt"; the condition num of b1_ss is 102.6863

   int row = 0, col = 0;
   double value = 0;
   
   data >> rows;
   data >> cols;
   data >> nonezero;

   cout << "The matrix has " << rows << "rows, " << cols << "cols, " << nonezero <<
	" nonezero values." << endl;
   cout << endl;

   SparseMatrix<double, ColMajor> mat(rows,cols);
   tripletList.reserve(nonezero);         //reserve the space to store nonezero value

   while (data >> row && data >> col && data >> value)
   {
      tripletList.push_back(T(row - 1, col - 1, value));
   }   

   mat.setFromTriplets(tripletList.begin(),tripletList.end());  //create the sparse matrix
   mat.makeCompressed();
   //cout << mat << endl;

   /*****************************LU Decomposition***********************/
   VectorXd x(rows), b(rows);

   //b(0) = 3;
   //b(1) = 4;

   b(0) = 2;
   b(1) = 4;
   b(2) = 0;
   b(3) = 2;
   b(4) = 2;
   b(5) = 4;
   b(6) = 0;

   SparseLU<SparseMatrix<double, ColMajor>, COLAMDOrdering<int> > solver;
   solver.analyzePattern(mat);
   solver.factorize(mat);
   x = solver.solve(b);
   cout << x << endl;     //the solution to the linear algebra function
}


/*
test matrix is:
A = [1 -1; 3 -8]
b = [3; 4]

the solution to the pair linear algebra functions is:
x = [4; 1]
*/
