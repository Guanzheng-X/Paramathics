#include <iostream>
#include <stdio.h>
#include <vector>
#include "SparseMatrix.h"
#include <fstream>

using namespace std;
using namespace myg;

//void conjugate_gradient(SparseMatrix<unsigned int, double> &a, vector<double>&b, int num_iterations, vector<double>&x, vector<double>& xans);
void conjugate_gradient(SparseMatrix<unsigned int, double> &a, vector<double>&b, int num_iterations, vector<double>&x);
double inner_product(vector<double>&a, vector<double>&b);

int rows; int cols; int nonezero; //get the size of input matrix
const double kEpsilon = 0.00000001;
bool isConvergence = false;



double inner_product(vector<double>&a, vector<double>&b) {
	double result = 0.0;
	int n = a.size();
	for (int c = 0; c < n; c++) {
		result += a[c] * b[c];
	}
	return result;
}

//void conjugate_gradient( SparseMatrix<unsigned int, double> &a, vector<double>&b, int num_iterations, vector<double>&x, vector<double>&xans)
void conjugate_gradient(SparseMatrix<unsigned int, double> &a, vector<double>&b, int num_iterations, vector<double>&x)
{
	// compute r=b-Ax0
	vector<double> ax(rows);
	vector<double> r(rows);

	//tmp
	unsigned int corruptIt = 0;
	
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
		// jk is used to get control over the iteration fault injection
		jk++;
		/*
		// code to inject fault in the same run at various interval 
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

		/*
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
		*/ // remove this if you want self stabalizing algorithm	

			// This code  snippet allows us inject fault in row, column, value 
			
			if (iteration == 100000)
			{
				//vector<unsigned int> asdf = a.IA();
				//vector<unsigned int> basdf = a.JA();
				/*
				vector<double> aaa = a();

				//unsigned int corruptIt = aaa[0];
				 corruptIt = aaa[0];

				 //double old_alpha = corruptIt;
				 long long te = *(long long*)&corruptIt;
				 //for (int i = 63; i >= 0; i--)
				// {
					 //std::cout << ((te >> i)&1);
					 //if (i == 63) std::cout << " ";
					 //else if (i == 52) std::cout << " ";
				// }

				 // create a temp number and then shift it to correct position
				 long long num1 = 1;
				 long long vv = num1 << 35;
				 //std::cout << std::endl;
				 // Now xor the alpha value
				 long long ans = te ^ vv;
				// for (int i = 63; i >= 0; i--)
				 //{
					 //std::cout << ((ans >> i) & 1);
					 //if (i == 63) std::cout << " ";
					 //else if (i == 52) std::cout << " ";
				 //}
				 //std::cout << std::endl;
				 corruptIt = *(double*)&ans;

				//unsigned int flip = 16;

				//unsigned int reVal = corruptIt ^ flip;

				aaa[0] = corruptIt;
				const vector<double> aaaaa = aaa;
				a.operator()( aaaaa);
				*/

				// row vector fault injection
				vector<unsigned int> asdf = a.IA();
				corruptIt = asdf[1];
				asdf[1] = 25;
				a.IA(asdf);
			}
			

			// temp = A * p
			vector<double> q(rows);
			q = a * p;

			// alpha = r^2 / (p' * temp)
			double alpha = inner_product(r, r) / inner_product(p, q);

			// This is the orginal code snippet which is used to get control over fault injection
			// in any iteration , here the fault is injected in alpha (step length of cg method)
			/*
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
				*/
			//} (def remoce it)

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

			//std::cout << sqrt(converge) << endl;
			converge = inner_product(r, r);

			// norm code
			// calculate relative A norm
			/*
			vector<double> ans1(rows);
			for (unsigned int i = 0; i < xans.size(); i++)
			{
				ans1[i] = xans[i] - x[i];
			}
			vector<double> norm1 = a * ans1;
			double norm = inner_product(norm1, ans1);
			double normAns = sqrt(norm);
			if (iteration % 50 == 0)
			{
				cout << normAns << endl;
			}
			*/
			// This is code snippet helps us to mimic the transient error since we flip the corrupted
			// values back to orginal value
			if (iteration == 100000)
			{
				/*
				vector<double> aaaa = a();				
				aaaa[0] = corruptIt;
				const vector<double>aa = aaaa;
				a.operator()(aa);
				*/
				// restoring corrupted row index value
				vector<unsigned int> asdf = a.IA();
				asdf[1] = corruptIt;
				a.IA(asdf);
			}
			
//remove this if you want to inlcude self stabalizing algorithm		} 			
			
	}

	if (!isConvergence) {
		std::cout << "Does not converge after " << iteration << std::endl;
		std::cout << "residual is " << sqrt(converge) << std::endl;
	}

}







int main(int argc, char* argv[]) 
{
	
	fstream read;

	read.open(R"(F:\FEM_3D_thermal1\FEM_3D_thermal1.mtx)", ios::in);

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
	vector<double> xans(rows);

	fstream read1;

	read1.open(R"(F:\XanswerBcsstk27.mtx)");

	int p = 0;
	double ansValue = 0;
	while (!read1.eof())
	{
		read1 >> ansValue;
		xans[p] = ansValue;
		p++;
	}
	*/
	/*
	rows = 2;
	myg::SparseMatrix<unsigned int, double> A(rows); // rows x rows identity matrix
	
	A(0, 0, 4);
	A(0, 1, 1);
	A(1, 0, 1);
	A(1, 1, 3);
	*/
	vector<double> b(rows);
	for (int i = 0; i < rows; i++) {
		b[i]=1;
	}

	vector<double> x(rows);

	
	//vector<unsigned int> asdf = A.IA();
	//vector<unsigned int> basdf = A.JA();

	/*
	for (int i = 0; i < asdf.size(); i++)
	{
		cout << asdf[i] << endl;
	}
	cout << "New Line" << endl;
	/*
	for (int i = 0; i < basdf.size(); i++)
	{
		cout << basdf[i] << "";
	}
	*/
	conjugate_gradient(A, b, 2000, x);
	//conjugate_gradient(A, b, 2000, x, xans);
	
	for (int i = 0; i < x.size(); i++)
	{
		cout << x[i] << endl;
	}
	

	
	return 0;
	
}