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

	//int ere = 0;

	/*
	// let's start with the block size of 10
	vector<vector<double>> blockCheckSum;
	// get the count of the number of elements you want to add in each row of checksum matrix
	auto blockCount = (a().size())/10;
	int cntValue = 10;
	int rowIdx = 1;

	while (cntValue > 0) {
	auto tmpCount = blockCount;
	// add the row into blockCheckSum matrix and decrement the count of tmpCount
	while (tmpCount > 0) {
	if (tmpCount > rowPointer[rowIdx - 1] - rowPointer[rowIdx]) {
	int l = rowPointer[rowIdx - 1] - rowPointer[rowIdx];
	// get the rowSize of the matrix A
	vector<double> ins(rowPointer.size()/10, 0);
	for (; l > 0; l--) {
	ins[colPointer[l]] = valPointer[l];
	}
	}
	}
	cntValue--;
	}


	// Below code is for block checksum approach

	// create a 2d array to store the checksum of blocks and then convert it into sparse matrix
	vector<vector<double>> denseCheckMatrix(90, vector<double>(900,0));

	int denseRowIt = 0;
	for (int i = 0; i < matrixA.size(); i++) {
	if (i == 0) {
	}
	else {
	if (i % 10 == 0) {
	denseRowIt++;
	}
	}
	// now run through the matrixA to construct dense block checksum matrix
	for (int j = 0; j < cols; j++) {
	denseCheckMatrix[denseRowIt][j] += matrixA[i][j];
	}
	}

	// now write the csv file into txt file
	ofstream denseMat("densechecksummatrix.txt");
	denseMat << 90 << " " << 900 << " " << 2992 << endl;;

	for (int i = 0; i < denseCheckMatrix.size(); i++) {
	for (int j = 0; j < cols; j++) {
	if (denseCheckMatrix[i][j] != 0)
	denseMat << i + 1 << " " << j + 1 << " " << denseCheckMatrix[i][j] << endl;
	}
	}
	*/

	// now read the txt file for dense checksum matrix into a sparse matrix

	fstream readDenseMat;
	readDenseMat.open("densechecksummatrix.txt");
	int rr = 0; int cc = 0; double vv = 0;
	readDenseMat >> rr >> cc >> vv;
	myg::SparseMatrix<unsigned int, double> sparseCheckSum(cc); // rows x rows identity matrix

	int rSparse;
	int cSparse;
	double vSparse;

	while (!readDenseMat.eof()) {
		readDenseMat >> rSparse;
		readDenseMat >> cSparse;
		readDenseMat >> vSparse;
		sparseCheckSum(rSparse - 1, cSparse - 1, vSparse);
	}		// end of file

	//cout << "sparse checksum matrix size is : " << sparseCheckSum.Dim() << endl;

	readDenseMat.close();

	// block checksum code ends here

	// ------------The preprocessing code snippet ends here

	//this is the original cg code

	// compute r=b-Ax0
	vector<double> ax(rows);
	vector<double> r(rows);

	//tmp
	//unsigned int corruptIt = 0;
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

		int ite = 4;

		if (iteration == ite)
		{
			vector<double> aaa = a();
			//cout << "Iteration during which fault was injected " << iteration << endl;

			corruptIt = aaa[idxFault1];
			//cout << "Value of in matrix A before corruption " << corruptIt << endl;
			double corr = corruptIt;

			long long te = *(long long*)&corr;
			long long num1 = 1;
			long long vv = num1 << 61;

			long long ans = te ^ vv;
			corr = *(double*)&ans;
			aaa[idxFault1] = corr;
			//cout << "Value after corruption " << corr << endl;

			corruptIt2 = aaa[idxFault2];
			//cout << "Value in matrix A before corruption " << corruptIt2 << endl;
			double corr2 = corruptIt2;

			// corrupt the second value
			long long te2 = *(long long*)&corr2;
			long long num2 = 1;
			long long vv2 = num2 << 61;
			long long ans2 = te2 ^ vv2;

			corr2 = *(double*)&ans2;
			aaa[idxFault2] = corr2;
			//cout << "Value after corruption of second value " << corr2 << endl;

			const vector<double> aaaaa = aaa;
			a.operator()(aaaaa);
		}


		// temp = A * p
		vector<double> q(rows);
		q = a * p;

		// This is code snippet helps us to mimic the transient error since we flip the corrupted
		// values back to orginal value
		if (iteration == ite)
		{
			vector<double> aaaa = a();

			aaaa[idxFault1] = corruptIt;
			//cout << "Value has been restored back " << corruptIt << endl;

			// insert second fault
			aaaa[idxFault2] = corruptIt2;
			//cout << "Second value has been restored back " << corruptIt2 << endl;

			const vector<double>aa = aaaa;
			a.operator()(aa);
		}


		// below is the code for block checksum method to detect if there is any fault in Spmv
		auto sparseErrorScalarVector = sparseCheckSum * p;
		//cout << "size of sparseErrorScalarVector is : " << sparseErrorScalarVector.size() << endl;
		double sumAn = 0;
		int sumAnIt = 0;
		vector<double> sparseErrorVec(90, 0);


		// shrink q to appropriate block size 
		for (int i = 0; i < q.size(); i++) {

			if ((i+1) % 10 == 0) {
				sparseErrorVec[sumAnIt] = sumAn;
				sumAn = q[i];
				sumAnIt++;
			}
			else {
				sumAn += q[i];
			}
		}

		// update the last summed value
		sparseErrorVec[sparseErrorVec.size() - 1] += sumAn;


		vector<int> errorIndexVec;
		// now compare the two result
		for (int i = 0; i < 90; i++) {
			auto tempErrV = sparseErrorVec[i] - sparseErrorScalarVector[i];
			if (abs(tempErrV) > 1000) {
				// record the error index
				errorIndexVec.push_back(i);
				//cout << "caught the error at index " << i << endl;
			}
		}

		// now check for error and re-execute the required block
		for (auto&i : errorIndexVec) {
			int startRowPointer = i * 10;
			auto corAnsBlock = a.multiplySpecific(p, startRowPointer);
			// update the ans
			std::copy_n(corAnsBlock.begin(), corAnsBlock.size(), q.begin() + startRowPointer);
		}


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

	//std::cout << "matrix has " << rows << " rows " << cols << " cols " << nonezero << " nonezero values" << endl;

	//SparseMatrix mat(rows, cols);
	myg::SparseMatrix<unsigned int, double> A(rows); // rows x rows identity matrix

	//cout << "argv1: " << argv[1] << " argv2: " << argv[2] << endl;

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
