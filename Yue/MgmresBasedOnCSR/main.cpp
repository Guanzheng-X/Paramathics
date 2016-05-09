#include <iostream>
#include <stdio.h>
#include <vector>
#include "SparseMatrix.hpp"
#include <fstream>
#include "mgmres.hpp"
using namespace std;
using namespace myg;



int rows; int cols; int nonezero; //get the size of input matrix
const double kEpsilon = 1e-8;
bool isConvergence = false;



int main(int argc, char* argv[]) 
{
    
    fstream read;
    
    read.open("/Users/wuyue/Desktop/dataset/bcsstm05.mtx", ios::in);
    
    int row = 0; int col = 0; float value = 0; //get the value
    read >> rows; read >> cols; read >> nonezero;
    
    std::cout << "matrix has " << rows << " rows " << cols << " cols " << nonezero << " nonezero values" << endl;
    
    //SparseMatrix mat(rows, cols);
    myg::SparseMatrix<int, double> mat(rows); // rows x rows identity matrix
    
    while (!read.eof()) {
        read >> row;
        read >> col;
        read >> value;
        
        mat(row - 1, col - 1, value);
        
        if (row != col)
        {
            mat(col - 1, row - 1, value);
        }
    }		// end of file
    
    read.close();
    
    /*
     rows = 2;
     myg::SparseMatrix<unsigned int, double> A(rows); // rows x rows identity matrix
     
     A(0, 0, 4);
     A(0, 1, 1);
     A(1, 0, 1);
     A(1, 1, 3);
     */
//    cout<<"test"<<endl;
//    int size =(int) mat().size();
//    std::cout<<size<<endl;
////    for(int i = 0; i < temp.size(); i++){
////        cout<<temp[i];
////    }
//    cout<<"size "<<mat.Dim()<<endl;
//    cout<<mat;
    
    vector<double> b;
    for (int i = 0; i < rows; i++) {
        b.push_back(1);
    }
    
    
    vector<double> x(rows);
    
//    x[0] = 2;
//    x[1] = 1;
    
    pmgmres_ilu_cr(mat, &x[0], &b[0], 10, 5, 1.0E-08, 1.0E-08);
    for(int i = 0; i < x.size(); i++){
        cout<<x[i]<<" ";
    }
    cout<<endl;
    // return
    return 0;
}