#include "math.hpp"
#include <fstream>
#include "sparse_matrix.hpp"
#include <iostream>
#include <cstdio>
#include <string>
#include <ctime>

#include <vector>
using namespace std;
void conjugate_gradient_test(int val1, int val2, int val3) {
	int rows; int cols; int nonezero; //get the size of input matrix
	fstream read;

	read.open("F:\\bcsstk27\\bcsstk27.mtx", ios::in);

	int row = 0; int col = 0; float value = 0; //get the value
	read >> rows; read >> cols; read >> nonezero;

	std::cout<<"matrix has "<<rows<<" rows "<< cols<<" cols "<<nonezero<<" nonezero values"<<endl;
	SparseMatrix mat(rows, cols);
	while (!read.eof()) {
		read >> row;
		read >> col;
		read >> value;
		mat.set_element(row - 1, col - 1, value);
		if (row != col)
			mat.set_element(col - 1, row - 1, value);
		//cout<<mat.get_element(row-1, col-1);       
	}

	read.close();

	std::vector<double> vec;
	for (int i = 0; i < rows; i++) {
		vec.push_back(1);
	}

	std::vector<double> x(rows);
	//    for(int i=0;i<rows;i++){
	//        x.push_back(1);
	//    }

	/*
	// x original answer vector
	std::vector<double> x0(rows);
	fstream xans;
	xans.open("F:\\XanswerBcsstk27.mtx");
	double xansValue = 0;
	while (!xans.eof())
	{
		xans >> xansValue;
		x0.push_back(xansValue);
	}
	xans.close();
	*/
	
	Math::conjugate_gradient(mat, &vec[0], 2000, &x[0], val1, val2, val3);

	/*
	// print the values of x
	std::cout << "Print the value of x" << std::endl;
	for (int i = 0; i < x.size(); i++)
	{
		std::cout << x[i] << endl;
	}
	*/

	
}
    //for (std::vector<double>::iterator it = x.begin() ; it != x.end(); ++it){
        // std::cout << ' ' << *it<<std::endl;
        //std::cout << ' ' << *it;
    //}

     //std::cout<<endl;
    


int main(int argc, char* argv[]) {
	
	string s1 = argv[1];
	string s2 = argv[2];
	string s3 = argv[3];

	int val1 = std::stoi(s1);
	int val2 = std::stoi(s2);
	int val3 = std::stoi(s3);

	clock_t start;
	double duration;

	start = std::clock();

	/* Your algorithm here */
	conjugate_gradient_test(val1, val2, val3);

	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;

	std::cout << duration << '\n';

	std::getchar();
    
    return 0;
}

