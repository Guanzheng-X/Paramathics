#include <iostream>
#include <stdio.h>
#include <vector>
#include "SparseMatrix.h"
#include <fstream>
#include <numeric>
#include <algorithm>
#include <chrono>

using namespace std;
using namespace myg;


struct blockCheck {
	vector<vector<double>> vec;
	vector<pair<int, int>> idx;

};


// checksum matirx A

int rows; int cols; int nonezero; //get the size of input matrix
const double kEpsilon = 0.00000001;
bool isConvergence = false;

bool run = false;
int cnt = 0;

//void conjugate_gradient(SparseMatrix<unsigned int, double> &a, vector<double>&b, int num_iterations, vector<double>&x, vector<double>& xans);
void conjugate_gradient(SparseMatrix<unsigned int, double> &a, vector<double>&b, int num_iterations, vector<double>&x, int idxFault1, int idxFault2);
double inner_product(vector<double>&a, vector<double>&b);


double inner_product(vector<double>&a, vector<double>&b) {
	double result = 0.0;
	int n = a.size();
	for (int c = 0; c < n; c++) {
		result += a[c] * b[c];
	}
	return result;
}

//void conjugate_gradient( SparseMatrix<unsigned int, double> &a, vector<double>&b, int num_iterations, vector<double>&x, vector<double>&xans)
void conjugate_gradient(SparseMatrix<unsigned int, double> &a, vector<double>&b, int num_iterations, vector<double>&x, int idxFault1, int idxFault2)
{
	
	//this is the original cg code

	// compute r=b-Ax0
	vector<double> ax(rows);
	vector<double> r(rows);

	//tmp
	//unsigned int corruptIt = 0;
	// use this declaration for fault injection in A matrix value
	// A multiply x0
	ax = a * x;


	// residual r = b - Ax0
	for (int i = 0; i < rows; i++)
	{
		r[i] = b[i] - ax[i];
	}

	// intitially p0 = r
	vector <double> p = r;

	int iteration;
	double converge = 1;

	// start the timer
	auto time_1 = std::chrono::high_resolution_clock::now();

	for (iteration = 0; iteration < num_iterations + 500; iteration++)
	{

		converge = inner_product(r, r);

		//std::cout << sqrt(converge) << std::endl;

		if (sqrt(converge)<kEpsilon)
		{
			//cout << "residual " << sqrt(converge);
			std::cout << "iteration number: " << iteration << std::endl;
			isConvergence = true;
			break;
		}


		// temp = A * p
		vector<double> q(rows);
		q = a * p;

		double alpha = inner_product(r, r) / inner_product(p, q);

		vector<double> xi(rows);

		for (int i = 0; i < rows; i++)
		{
			xi[i] = (x[i] + (alpha * p[i]));
		}

		x = xi;

		vector<double> ri(rows);

		for (int i = 0; i < rows; i++)
		{
			ri[i] = (r[i] - (alpha * q[i]));
		}


		double beta = inner_product(ri, ri) / inner_product(r, r);

		r = ri;

		vector<double> p1(rows);

		for (int i = 0; i < rows; i++)
		{
			p1[i] = (r[i] + (beta * p[i]));
		}

		p = p1;

		//std::cout << sqrt(converge) << endl;
		converge = inner_product(r, r);

	}

	// end timer
	auto time_2 = std::chrono::high_resolution_clock::now();
	cout << "Time for execution in ms " << chrono::duration_cast<std::chrono::milliseconds>(time_2 - time_1).count() << endl;

	if (!isConvergence) {
		std::cout << "Does not converge after " << iteration << std::endl;
		//std::cout << "residual is " << sqrt(converge) << std::endl;
	}

}





int main(int argc, char* argv[])
{

	fstream read;

	read.open(R"(gr_30_30.mtx)", ios::in);

	int row = 0; int col = 0; float value = 0; //get the value
	read >> rows; read >> cols; read >> nonezero;

	std::cout << "matrix has " << rows << " rows " << cols << " cols " << nonezero << " nonezero values" << endl;

	myg::SparseMatrix<unsigned int, double> A(rows); // rows x rows identity matrix

	while (!read.eof()) {
		read >> row;
		read >> col;
		read >> value;

		A(row - 1, col - 1, value);
		if (row != col)
		{
			A(col - 1, row - 1, value);
		}
	}		// end of file

	read.close();

	vector<double> b(rows);
	for (int i = 0; i < rows; i++) {
		b[i] = 1;
	}

	vector<double> x(rows);

	conjugate_gradient(A, b, 2000, x);

	//system("pause");
	return 0;

}
