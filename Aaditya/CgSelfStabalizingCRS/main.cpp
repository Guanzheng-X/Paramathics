#include <iostream>
#include <stdio.h>
#include <vector>
#include "SparseMatrix.h"
#include <fstream>

using namespace std;
using namespace myg;

void conjugate_gradient(SparseMatrix<unsigned int, double> &a, vector<double>&b, int num_iterations, vector<double>&x);
double inner_product(vector<double>a, vector<double>b);

int rows; int cols; int nonezero; //get the size of input matrix
const double kEpsilon = 1e-8;
bool isConvergence = false;



double inner_product(vector<double>a, vector<double>b) {
	double result = 0.0;
	int n = a.size();
	for (int c = 0; c < n; c++) {
		result += a[c] * b[c];
	}
	return result;
}

void conjugate_gradient( SparseMatrix<unsigned int, double> &a, vector<double>&b, int num_iterations, vector<double>&x)
{
	// compute r=b-Ax0
	vector<double> ax(rows);
	vector<double> r(rows);
	
	// A multiply x0
	ax = a * x;

	
	// residual r = b - Ax0
	for (int i = 0; i < rows; i++)
	{
		r[i] = b[i] - ax[i];
	}

	// intitially p0 = r

	vector <double> p = r;

	// compute scalar alpha0
	// alpha0 = (r0T * r) / (p0T * A * p0)


	// using self stabalizing algorithm

	int correction_frequency = 10;
	int iteration;
	double converge = 1;
	int jk = 0;

	for (iteration = 0; iteration < num_iterations + 2000; iteration++) 
	{

		converge = inner_product(r, r);
		//std::cout<<sqrt(converge)<<std::endl;
		jk++;
		/*
		if (jk < 39)
		{
		correction_frequency = 5 * iter_val;
		}
		else if (jk < 117)
		{
		correction_frequency = 15 * iter_val;
		}
		else if (jk < 234)
		{
		correction_frequency = 25* iter_val;
		}
		else if (jk < 312)
		{
		correction_frequency = 35* iter_val;
		}
		else if (jk < 389)
		{
		correction_frequency = 45* iter_val;
		}
		*/
		//std::cout << sqrt(converge) << std::endl;

		if (sqrt(converge)<kEpsilon) 
		{
			cout << "residual " << sqrt(converge);
			std::cout << "iteration number: " << iteration << std::endl;
			isConvergence = true;
			break;
		}

		
		if (iteration%correction_frequency == 0 % correction_frequency) {
			
			//r=A*x, temp=A*p
			r = a*x;

			vector<double> q = a * p;

			//ri=b-ri
			vector<double> tmp(rows);

			for (int i = 0; i < rows; i++)
			{
				tmp[i] = b[i] - r[i];
			}
			
			r = tmp;

			//alpha = ri^tpi/pi^tqi
			double alpha = inner_product(r, p) / inner_product(p, q);


			// std::cout<<alpha<<std::endl;

			//xi+1=xi+alpha*pi
			vector<double> xi(rows);

			for (int i = 0; i < rows; i++)
			{
				xi[i] = (x[i] + (p[i] * alpha));
			}

			x = xi;

			//ri+1=ri-alpha*q
			vector <double> ri(rows);

			for (int i = 0; i < rows; i++)
			{
				ri[i] = (r[i] - (q[i] * alpha));
			}

			r = ri;

			//beta=-r^t*q/p^t*q
			double temp1 = inner_product(r, q);
			double temp2 = inner_product(p, q);
			double beta = -1 * (temp1 / temp2);

			
			//pi+1=ri+1+beta*pi
			vector<double> pi(rows);

			for (int i = 0; i < rows; i++)
			{
				pi[i] = (r[i] +  (beta * p[i]));
			}

			p = pi;

		}
		else
		
		{
			
			// temp = A * p
			vector<double> q(rows);
			q = a * p;

			// alpha = r^2 / (p' * temp)
			double alpha = inner_product(r, r) / inner_product(p, q);

			// Fault injection code
			if (iteration == 7)
			{
				// alpha = rs_old / (p' * temp)
				std::cout << "Look at the alpha before flip " << alpha << std::endl;
				double old_alpha = alpha;
				long long te = *(long long*)&alpha;
				for (int i = 63; i >= 0; i--)
				{
					//std::cout << ((te >> i)&1);
					//if (i == 63) std::cout << " ";
					//else if (i == 52) std::cout << " ";
				}

				// create a temp number and then shift it to correct position
				long long num1 = 1;
				long long vv = num1 << 61;
				//std::cout << std::endl;
				// Now xor the alpha value
				long long ans = te ^ vv;
				for (int i = 63; i >= 0; i--)
				{
					//std::cout << ((ans >> i) & 1);
					//if (i == 63) std::cout << " ";
					//else if (i == 52) std::cout << " ";
				}
				//std::cout << std::endl;
				alpha = *(double*)&ans;
				long long tea = *(long long*)&alpha;
				for (int i = 63; i >= 0; i--)
				{
					//std::cout << ((tea >> i) & 1);
					//if (i == 63) std::cout << " ";
					//else if (i == 52) std::cout << " ";
				}
				//std::cout << "Number after bit " << flipBit + 1 << " is flipped" << " in iteration " << (itt) << "= " << alpha << std::endl;
				std::cout << "Change in alpha: " << alpha - old_alpha << std::endl;
			}



			//std::cout<<alpha<<std::endl;
			vector<double> xi(rows);

			// xi = x + alpha * p
			for (int i = 0; i < rows; i++)
			{
				xi[i] = (x[i] + (alpha * p[i]) );
			}

			
			x = xi;



			// ri+1 = ri - alpha * temp
			

			vector<double> ri(rows);

			for (int i = 0; i < rows; i++)
			{
				ri[i] = (r[i] - (alpha * q[i]) );
			}

			

			//beta=ri+1 ^2 / ri^2
			double beta = inner_product(ri, ri) / inner_product(r, r);
			
			r = ri;

			// pi = ri + (beta) * p
			vector<double> p1(rows);

			for (int i = 0; i < rows; i++)
			{
				p1[i] = (r[i] + (beta * p[i]));
			}

			p = p1;

			std::cout << converge << " " << sqrt(converge) << endl;
			converge = inner_product(r, r);

			
		}
	}

	if (!isConvergence) {
		std::cout << "not converge at " << iteration << std::endl;
		std::cout << "residual is " << sqrt(converge) << std::endl;
	}

}







int main(int argc, char* argv[]) 
{
	
	fstream read;

	read.open("F:\\gr_30_30\\gr_30_30.mtx", ios::in);

	int row = 0; int col = 0; float value = 0; //get the value
	read >> rows; read >> cols; read >> nonezero;

	std::cout << "matrix has " << rows << " rows " << cols << " cols " << nonezero << " nonezero values" << endl;

	//SparseMatrix mat(rows, cols);
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
	
	/*
	rows = 2;
	myg::SparseMatrix<unsigned int, double> A(rows); // rows x rows identity matrix
	
	A(0, 0, 4);
	A(0, 1, 1);
	A(1, 0, 1);
	A(1, 1, 3);
	*/
	vector<double> b;
	for (int i = 0; i < rows; i++) {
		b.push_back(i + 1);
	}


	vector<double> x(rows);
	x[0] = 2;
	x[1] = 1;


	conjugate_gradient(A, b, 2000, x);







	// return
	return 0;
}