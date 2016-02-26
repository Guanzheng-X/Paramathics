//
//  main.cpp
//  CG2
//
//  Created by 吴越 on 16/2/19.
//  Copyright © 2016年 Yue Wu. All rights reserved.
//
#include <iostream>

#include <Core>
#include <Sparse>
#include <fstream>
#include <vector>
#include <ctime>
using namespace Eigen;
using namespace std;


int main(int argc, const char * argv[]) {
    
    int rows; int cols ; int nonezero; //get the size of input matrix
    //int n = 3;
   
    typedef Eigen::Triplet<double> T;//use to store row, col, v_ij
    
    vector<T> tripletList;	//use to store all the row
    
    /***************beging to construct the matrix****************/
    fstream read;
    
    read.open("/Users/wuyue/Desktop/test/Andrews1.mtx", ios::in);
    
    int row=0 ; int col=0; float value=0; //get the value
    read>>rows; read >> cols; read>>nonezero;
    
    cout<<"matrix has "<<rows<<" rows "<< cols<<" cols "<<nonezero<<" nonezero values"<<endl;
    
    SparseMatrix<double,RowMajor> mat(rows,cols);//declare a row major matrix
    tripletList.reserve(nonezero); //reserve the space to store nonezero value
    
    while (!read.eof()){
        read>>row;
        read>>col;
        read>>value;
        tripletList.push_back(T(row-1,col-1,value));
        
        
    }
    
    mat.setFromTriplets(tripletList.begin(), tripletList.end()); //crate the matrix
    mat.makeCompressed();
   // cout<<mat<<endl;
    
    VectorXd x(rows);
    VectorXd b(rows);
    srand((unsigned)time(0));
    for (int i=0; i<rows; i++) {
        b(i)=1;
       
    }
    ConjugateGradient<SparseMatrix<double,RowMajor>, Upper >cg;
    //LeastSquaresConjugateGradient<SparseMatrix<double> >cg;
  
     //ConjugateGradient<SparseMatrix<double,RowMajor>,Eigen::Lower|Eigen::Upper,Eigen::IdentityPreconditioner>cg;

    //BiCGSTAB<SparseMatrix<double,RowMajor> >cg;
    cg.compute(mat);
    //cg.analyzePattern(mat);
    //cg.factorize(mat);
    if(cg.info()!=Success){
        std::cout<<"decomposition failed"<<endl;
        return -1;
    }
   
    cg.setMaxIterations(1000);
    cg.setTolerance(0.0000001);
    //std::cout<<"#maxiteration: "<<cg.maxIterations()<<endl;
   
   // while(i<cg.maxIterations()) {
        

    asm("#it begins here!");
    long startTime = time(NULL);
    x=cg.solve(b);
    long endTime = time(NULL);
    std::cout << "#iterations: " << cg.iterations() << std::endl;
    std::cout << "estimated error: " << cg.error()      << std::endl;
    cout<<"solving time: "<<endTime-startTime<<" seconds"<<endl;
    asm("#it ends here!");
    if(cg.info()!=Success){
        std::cout<<"solving failed"<<endl;
        std::cout<<cg.info()<<endl;
        return -2;
    }
  
    //}
  //  VectorXd res = cg.solve(b);
   
        std::cout<<x<<endl;
}
