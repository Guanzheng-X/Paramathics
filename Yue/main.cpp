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
     //Test Case 1
//    SparseMatrix a(4, 2);
//    a.set_element(0, 0, 1.0);
//    a.set_element(0, 1, 2.0);
//    
//    a.set_element(1, 0, 2.0);
//    a.set_element(1, 1, -3.0);
//    
//    a.set_element(2, 0, 4.0);
//    a.set_element(2, 1, -1.0);
//    
//    a.set_element(3, 0, -5.0);
//    a.set_element(3, 1, 2.0);
//    
//    // Suppose solution is (2, 3).
//    std::vector<double> b(4);
//    b[0] = 3.0;  // 8.0
//    b[1] = 2.0;  // -5.0
//    b[2] = -4.0;
//    b[3] = 5.0;
//    
//    std::vector<double> x(2);
//    x[0] = x[1] = 0.0;
//    
//    Math::conjugate_gradient(a, &b[0], 2, &x[0]);
//        for (std::vector<double>::iterator it = x.begin() ; it != x.end(); ++it)
//            std::cout << ' ' << *it;
//        std::cout<<endl;
//    std::vector<double> c(4);
//    a.multiply_column(&x[0], &c[0]);
//    
//    printf("x: %lf %lf\n", x[0], x[1]);
//    
//    for (int i = 0; i < 4; i++) {
//        printf(" %lf(%lf)", c[i], b[i]);
//    }
//    printf("\n");
//    
//    // Test Case 2
//    a = SparseMatrix(1, 2);
//    a.set_element(0, 0, 1.0);
//    a.set_element(0, 1, 2.0);
//    
//    b.resize(1);
//    b[0] = 3.0;
//    
//    // x[0] = x[1] = 0.0;
//    x[0] = 7.0;
//    x[1] = -3.0;
//    
//    Math::conjugate_gradient(a, &b[0], 2, &x[0]);
//    
//    c.resize(1);
//    a.multiply_column(&x[0], &c[0]);
//    
//    printf("x: %lf %lf\n", x[0], x[1]);
//    
//    for (int i = 0; i < 1; i++) {
//        printf(" %lf(%lf)", c[i], b[i]);
//    }
//    printf("\n");
//    
//    printf("} conjugate_gradient_test\n");
}

int main() {
    conjugate_gradient_test();
    
    return 0;
}

