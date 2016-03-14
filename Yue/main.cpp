#include "math.hpp"
#include <fstream>
#include "sparse_matrix.hpp"
#include <iostream>
#include <cstdio>

#include <vector>
using namespace std;
void conjugate_gradient_test() {
    int rows; int cols ; int nonezero; //get the size of input matrix
    fstream read;
    
    read.open("/Users/wuyue/Desktop/dataset/bcsstm22.mtx", ios::in);
    
    int row=0 ; int col=0; float value=0; //get the value
    read>>rows; read >> cols; read>>nonezero;
    
    std::cout<<"matrix has "<<rows<<" rows "<< cols<<" cols "<<nonezero<<" nonezero values"<<endl;
    SparseMatrix mat(rows, cols);
    while (!read.eof()){
        read>>row;
        read>>col;
        read>>value;
        mat.set_element(row-1, col-1, value);
        //cout<<mat.get_element(row-1, col-1);
       
        
    }
    
    std::vector<double> vec;
    for(int i=0;i<rows;i++){
        vec.push_back(1);
    }
    
    std::vector<double> x(rows);
//    for(int i=0;i<rows;i++){
//        x.push_back(1);
//    }
    Math::conjugate_gradient(mat, &vec[0], 2000, &x[0]);
    for (std::vector<double>::iterator it = x.begin() ; it != x.end(); ++it)
        std::cout << ' ' << *it;
    std::cout<<endl;

}

int main() {
    conjugate_gradient_test();
    
    return 0;
}

