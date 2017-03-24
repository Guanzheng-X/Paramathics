#include <iostream>
#include <stdio.h>
#include <vector>
#include "SparseMatrix.h"
#include <fstream>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <string>

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
vector<vector<double>> matrixA(900, vector<double>(900, 0));


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
	// below code snippet will allow the preprocessing code to run and form a block checksum matrix
	// the idea here is to iterate over matrix A and based on number of nnz elements form a checksum matrix

	auto rowPointer = a.IA();
	auto colPointer = a.JA();
	auto valPointer = a();

	/*
    ofstream valFile("Position.txt");
	for (int i = 0; i < rowPointer.size(); i++) {
	valFile << rowPointer[i] << endl;
	}
	*/

	/*
	// create a vector which stores the index of value 8
	vector<int> idxVal8;

	// print the valPointer so that we know where to inject faults
	// now write the csv file into txt file
	ofstream valFile("valPointerOnly8Idx.txt");
	for (int i = 0; i < valPointer.size(); i++) {
	if (valPointer[i] == 8)
	idxVal8.push_back(i);
	//valFile << i + 1 << " " << valPointer[i] << endl;
	}

	for (int i = 0; i < 1001; i++) {
	auto idx1 = rand() % idxVal8.size();
	auto idx2 = rand() % idxVal8.size();
	valFile << idxVal8[idx1] << " " << idxVal8[idx2] << endl;
	}
	*/



	// This is recursive block checksum approach idea to implement checksum matrix

	// make one check sum matrix
	vector<double> check1(cols, 0);
	for (int i = 0; i < matrixA.size(); i++) {
		for (int j = 0; j < cols; j++) {
			//if (j == 0 && matrixA[i][j] != 0) {
			//cout << i << " " << j << endl;
			//cout << matrixA[i][j] << endl;
			//}
			check1[j] += matrixA[i][j];
		}
	}

	vector<vector<double>> check2(2, vector<double>(cols, 0));
	int indexOfChecksum = 0;
	bool flag = true;
	// iterate over the checksum matrix and insert the appropriate checksum vector in check2 2d vector
	for (int i = 0; i < matrixA.size(); i++) {
		for (int j = 0; j < cols; j++) {

			if (i > 455) {
				if (flag) {
					indexOfChecksum++;
					flag = false;
				}
			}
			check2[indexOfChecksum][j] += matrixA[i][j];
		}
	}

	// re-initialize the indexOfChecksum variable to zero
	indexOfChecksum = 0;
	//flag = true;
	vector<vector<double>> check4(4, vector<double>(cols, 0));
	vector<int> check4Range{ 230,455,670,900 };
	int check4It = 0;
	for (int i = 0; i < matrixA.size(); i++) {
		for (int j = 0; j < cols; j++) {
			if (i > check4Range[check4It] ) {
				if (check4It != 3) {
					check4It++;
					indexOfChecksum++;
				}
			}
			check4[indexOfChecksum][j] += matrixA[i][j];
		}
	}

	// re-initialize the indesOfChecksum variable to zero
	indexOfChecksum = 0;
	check4Range.clear();
	vector<vector<double>> check8(8, vector<double>(cols, 0));
	check4Range = { 120,230,340,455,560,670,780,900 };
	check4It = 0;
	for (int i = 0; i < matrixA.size(); i++) {
		for (int j = 0; j < cols; j++) {
			if (i > check4Range[check4It]) {
				if (check4It != 7) {
					check4It++;
					indexOfChecksum++;
				}
			}
			check8[indexOfChecksum][j] += matrixA[i][j];
		}
	}




	vector<blockCheck> checkVector;
	blockCheck bc1;
	bc1.vec = check2;
	bc1.idx.push_back(make_pair(0, 454));
	bc1.idx.push_back(make_pair(455, 899));
	checkVector.push_back(bc1);

	blockCheck bc2;
	bc2.vec = check4;
	bc2.idx.push_back(make_pair(0, 229));
	bc2.idx.push_back(make_pair(230, 454));
	bc2.idx.push_back(make_pair(455, 669));
	bc2.idx.push_back(make_pair(670, 899));
	checkVector.push_back(bc2);

	blockCheck bc3;
	bc3.vec = check8;
	bc3.idx.push_back(make_pair(0, 119));
	bc3.idx.push_back(make_pair(120, 229));
	bc3.idx.push_back(make_pair(230, 339));
	bc3.idx.push_back(make_pair(340, 454));
	bc3.idx.push_back(make_pair(455, 559));
	bc3.idx.push_back(make_pair(560, 669));
	bc3.idx.push_back(make_pair(670, 779));
	bc3.idx.push_back(make_pair(780, 899));
	checkVector.push_back(bc3);



	// ------------The preprocessing code snippet ends here

	//this is the original cg code

	// compute r=b-Ax0
	vector<double> ax(rows);
	vector<double> r(rows);

	// use this declaration for fault injection in A matrix value
	double corruptIt = 0;
	double corruptIt2 = 0;
	// A multiply x0
	ax = a * x;


	// residual r = b - Ax0
	for (int i = 0; i < rows; i++)
	{
		r[i] = b[i] - ax[i];
	}

	// intitially p0 = r
	vector <double> p = r;

	// using self stabalizing algorithm

	int iteration;
	double converge = 1;

	// start the timer
	auto time_1 = std::chrono::high_resolution_clock::now();

	for (iteration = 0; iteration < num_iterations + 500; iteration++)
	{
		converge = inner_product(r, r);
		//std::cout<<sqrt(converge)<<std::endl;

		if (sqrt(converge)<kEpsilon)
		{
			//cout << "residual " << sqrt(converge);
			std::cout << "iteration number: " << iteration << std::endl;
			isConvergence = true;
			break;
		}

        //cout << iteration << endl;
		// This code  snippet allows us inject fault in row, column, value 
		int ite = 4;

		if (iteration == ite)
		{
			vector<double> aaa = a();
			cout << "Iteration during which fault was injected " << iteration << endl;

			corruptIt = aaa[idxFault1];
			cout << "Value of in matrix A before corruption " << corruptIt << endl;
			double corr = corruptIt;

			long long te = *(long long*)&corr;
			long long num1 = 1;
			long long vv = num1 << 61;

			long long ans = te ^ vv;
			corr = *(double*)&ans;
			aaa[idxFault1] = corr;
			cout << "Value after corruption " << corr << endl;


			corruptIt2 = aaa[idxFault2];
			cout << "Value in matrix A before corruption " << corruptIt2 << endl;
			double corr2 = corruptIt2;

			// corrupt the second value
			long long te2 = *(long long*)&corr2;
			long long num2 = 1;
			long long vv2 = num2 << 61;
			long long ans2 = te2 ^ vv2;

			corr2 = *(double*)&ans2;
			aaa[idxFault2] = corr2;
			cout << "Value after corruption of second value " << corr2 << endl;

			const vector<double> aaaaa = aaa;
			a.operator()(aaaaa);
		}


		// temp = A * p
		vector<double> q(rows);
		q = a * p;

        
		// below is recursive block checksum mathod
		// use abft check here to detect if there was any fault while SpMv was done
		auto cSum = inner_product(p, check1);
		auto cqSum = std::accumulate(q.begin(), q.end(), 0);

		//auto ft = sqrt(inner_product(x, x));

        


		bool error = false;
		bool moreErrors = false;
		bool noErrorFlag = false;


        
		// error is detected if
		if (abs(cSum - cqSum) > 100) {
			error = true;
			cout << "error has been detected" << endl;
            //auto time_1 = std::chrono::high_resolution_clock::now();
			vector<pair<int, int>> extraError;
			if (error) {
				int k = 0;
				int j = 0;
				vector<bool> ans;
				while (k < 3) {
					auto f = checkVector[k].vec;

					for (int i = 0; i < 2; i++) {

						auto r1 = checkVector[k].idx[j + i].first;
						auto r2 = checkVector[k].idx[j + i].second;
						
						auto a2 = inner_product(checkVector[k].vec[j + i], p);
						auto a3 = accumulate(q.begin() + r1, q.begin() + r2 + 1, 0);
						if (abs(a2 - a3) > 100) {
							ans.push_back(true);
							// if this is last iteration then fix the value
							if (k == 2) {
								// recompute the corrupt part
								auto r3 = r1;
								for (; r3 < r2 + 1; r3++) {
									q[r3] = inner_product(matrixA[r3], p);
								}
								error = false;
								cout << "Error has been fixed. It was in  " << j + i << " block " << endl;
							
							}
						}
						else {
							ans.push_back(false);
						}
					}
					k++;
					int id = 0;
					for (auto const &i : ans) {
						if (i) {
							j += id;
							int kk = j;
							j = j * 2;
							// check if more error exists
							for (int o = id + 1; o < ans.size(); o++) {
								if (ans[o]) {
									// store j and id
									moreErrors = true;
									kk += o;
									auto tmpJ = kk * 2;
									extraError.push_back(make_pair(k, tmpJ));
								}
							}
							break;
						}
						id++;
					}
					ans.clear();
				}
			}

			cout << endl;

			auto cSum = inner_product(p, check1);
		auto cqSum = std::accumulate(q.begin(), q.end(), 0);

			// check if more errors exits
			if (moreErrors) {
				cout << "error has been detected" << endl;

				int k = extraError[0].first;
				int j = extraError[0].second;
				vector<bool> ans;

				while (k < 3) {
					auto f = checkVector[k].vec;
					for (int i = 0; i < 2; i++) {

						auto r1 = checkVector[k].idx[j + i].first;
						auto r2 = checkVector[k].idx[j + i].second;
						
						auto a2 = inner_product(checkVector[k].vec[j + i], p);
						auto a3 = accumulate(q.begin() + r1, q.begin() + r2 + 1, 0);
						if (abs(a2 - a3) > 100) {
							ans.push_back(true);
							// if this is last iteration then fix the value
							if (k == 2) {
								// recompute the corrupt part
								auto r3 = r1;
								
								for (; r3 < r2 + 1; r3++) {
									q[r3] = inner_product(matrixA[r3], p);
								}
								cout << "Error has been fixed. It was in  " << j + i << " block " << endl;
								moreErrors = false;
								break;
							}
						}
						else {
							ans.push_back(false);
						}
					}

					k++;
					int id = 0;
					for (auto const &i : ans) {
						if (i) {
							j += id;
							j = j * 2;
							break;
						}
						id++;
					}
					ans.clear();
				}


			}
    	//auto time_2 = std::chrono::high_resolution_clock::now() - time_1;
        //long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(time_2).count();
        //cout << "Time for error detection in ms " << microseconds << endl;

		}

       
		// recurisve block checksum approach ends here


		double alpha = inner_product(r, r) / inner_product(p, q);

		vector<double> xi(rows);

		// xi = x + alpha * p
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

		// This is code snippet helps us to mimic the transient error since we flip the corrupted
		// values back to orginal value
		if (iteration == ite)
		{
			vector<double> aaaa = a();

			aaaa[idxFault1] = corruptIt;
			cout << "Value has been restored back " << corruptIt << endl;

			// insert second fault
			aaaa[idxFault2] = corruptIt2;
			cout << "Second value has been restored back " << corruptIt2 << endl;

			const vector<double>aa = aaaa;
			a.operator()(aa);
		}

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
	cout << endl;

	fstream read;

	read.open(R"(gr_30_30.mtx)", ios::in);

	int row = 0; int col = 0; float value = 0; //get the value
	read >> rows; read >> cols; read >> nonezero;

	//std::cout << "matrix has " << rows << " rows " << cols << " cols " << nonezero << " nonezero values" << endl;

	//SparseMatrix mat(rows, cols);
	myg::SparseMatrix<unsigned int, double> A(rows); // rows x rows identity matrix



	while (!read.eof()) {
		read >> row;
		read >> col;
		read >> value;

		A(row - 1, col - 1, value);
		matrixA[row - 1][col - 1] = value;
		if (row != col)
		{
			A(col - 1, row - 1, value);
			matrixA[col - 1][row - 1] = value;
		}
	}		// end of file

	read.close();

	vector<double> b(rows);
	for (int i = 0; i < rows; i++) {
		b[i] = 1;
	}

	vector<double> x(rows);

	cout << "arg1: " << argv[1] << " " << "arg2: " << argv[2] << endl;
	conjugate_gradient(A, b, 2000, x, stoi(argv[1]), stoi(argv[2]));

	/*
	for (int i = 0; i < x.size(); i++)
	{
	cout << x[i] << endl;
	}
	*/


	//system("pause");
	return 0;
}
