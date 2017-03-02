#include <iostream>
#include <stdio.h>
#include <vector>
#include "SparseMatrix.h"
#include <fstream>
#include <numeric>
#include <algorithm>

using namespace std;
using namespace myg;


struct blockCheck {
	vector<vector<double>> vec;
	vector<pair<int,int>> idx;

};


// checksum matirx A

int rows; int cols; int nonezero; //get the size of input matrix
const double kEpsilon = 0.00000001;
bool isConvergence = false;

bool run = false;
int cnt = 0;
vector<vector<double>> matrixA(900, vector<double>(900, 0));


//void conjugate_gradient(SparseMatrix<unsigned int, double> &a, vector<double>&b, int num_iterations, vector<double>&x, vector<double>& xans);
void conjugate_gradient(SparseMatrix<unsigned int, double> &a, vector<double>&b, int num_iterations, vector<double>&x);
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
void conjugate_gradient(SparseMatrix<unsigned int, double> &a, vector<double>&b, int num_iterations, vector<double>&x)
{
	// below code snippet will allow the preprocessing code to run and form a block checksum matrix
	// the idea here is to iterate over matrix A and based on number of nnz elements form a checksum matrix
	
	auto rowPointer = a.IA();
	auto colPointer = a.JA();
	auto valPointer = a();
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

	// check for the exit condition
	// there might be a chance thet last few rows of the matrix a might not be in checksum matrix
	*/

	/*
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

	readDenseMat.close();

	// block checksum code ends here
	*/

	
	// This is recursive block checksum approach idea to implement checksum matrix
	
	// make one check sum matrix
	vector<double> check1(cols,0);
	for (int i = 0; i < matrixA.size();i++) {
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
			
			if (i > 255) {
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
	vector<int> check4Range{133,255,378};
	int check4It = 0;
	for (int i = 0; i < matrixA.size(); i++) {
		for (int j = 0; j < cols; j++) {
			if (i == check4Range[check4It]+1) {
				if (check4It != 2) {
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
	check4Range = { 71,133,194,255,317,378,439 };
	check4It = 0;
	for (int i = 0; i < matrixA.size(); i++) {
		for (int j = 0; j < cols; j++) {
			if (i == check4Range[check4It] + 1) {
				if (check4It != 6) {
					check4It++;
					indexOfChecksum++;
				}
			}
			check8[indexOfChecksum][j] += matrixA[i][j];
		}
	}

	// re-initialize the indesOfChecksum variable to zero
	indexOfChecksum = 0;
	check4Range.clear();
	vector<vector<double>> check16(16, vector<double>(cols, 0));
	check4Range = { 40,71,102,133,163,194,224,255,286,317,348,378,408,439,470 };
	check4It = 0;
	for (int i = 0; i < matrixA.size(); i++) {
		for (int j = 0; j < cols; j++) {
			if (i == check4Range[check4It]+1) {
				if (check4It != 14) {
					check4It++;
					indexOfChecksum++;
				}
			}
			check16[indexOfChecksum][j] += matrixA[i][j];
		}
	}

	// re-initialize the indesOfChecksum variable to zerp
	indexOfChecksum = 0;
	check4Range.clear();
	vector<vector<double>> check32(32, vector<double>(cols, 0));
	check4Range = { 22,40,55,71,86,102,117,133,147,163,169,194,209,224,237,255,271,286,302,317,333,348,364,378,394,408,424,439,456,470,486 };
	check4It = 0;
	for (int i = 0; i < matrixA.size(); i++) {
		for (int j = 0; j < cols; j++) {
			if (i == check4Range[check4It]+1) {
				if (check4It != 30) {
					check4It++;
					indexOfChecksum++;
				}
			}
			check32[indexOfChecksum][j] += matrixA[i][j];
		}
	}
	
	
	vector<blockCheck> checkVector;
	blockCheck bc1;
	bc1.vec = check2;
	bc1.idx.push_back(make_pair(0,255));
	bc1.idx.push_back(make_pair(256, 900));
	checkVector.push_back(bc1);

	blockCheck bc2;
	bc2.vec = check4;
	bc2.idx.push_back(make_pair(0,135));
	bc2.idx.push_back(make_pair(136, 255));
	bc2.idx.push_back(make_pair(256, 378));
	bc2.idx.push_back(make_pair(379, 900));
	checkVector.push_back(bc2);

	blockCheck bc3;
	bc3.vec = check8;
	bc3.idx.push_back(make_pair(0, 71));
	bc3.idx.push_back(make_pair(72, 135));
	bc3.idx.push_back(make_pair(136, 194));
	bc3.idx.push_back(make_pair(195, 255));
	bc3.idx.push_back(make_pair(256, 317));
	bc3.idx.push_back(make_pair(318, 378));
	bc3.idx.push_back(make_pair(379, 439));
	bc3.idx.push_back(make_pair(440, 900));
	checkVector.push_back(bc3);

	blockCheck bc4;
	bc4.vec = check16;
	bc4.idx.push_back(make_pair(0, 40));
	bc4.idx.push_back(make_pair(41, 71));
	bc4.idx.push_back(make_pair(72, 102));
	bc4.idx.push_back(make_pair(103, 135));
	bc4.idx.push_back(make_pair(136, 163));
	bc4.idx.push_back(make_pair(164, 194));
	bc4.idx.push_back(make_pair(195, 224));
	bc4.idx.push_back(make_pair(225, 255));
	bc4.idx.push_back(make_pair(256, 286));
	bc4.idx.push_back(make_pair(287, 317));
	bc4.idx.push_back(make_pair(318, 348));
	bc4.idx.push_back(make_pair(349, 378));
	bc4.idx.push_back(make_pair(379, 408));
	bc4.idx.push_back(make_pair(409, 439));
	bc4.idx.push_back(make_pair(440, 470));
	bc4.idx.push_back(make_pair(471, 900));
	checkVector.push_back(bc4);

	blockCheck bc5;
	bc5.vec = check32;
	bc5.idx.push_back(make_pair(0, 22));
	bc5.idx.push_back(make_pair(23, 40));
	bc5.idx.push_back(make_pair(41, 55));
	bc5.idx.push_back(make_pair(56, 71));
	bc5.idx.push_back(make_pair(72, 86));
	bc5.idx.push_back(make_pair(87, 102));
	bc5.idx.push_back(make_pair(103, 117));
	bc5.idx.push_back(make_pair(118, 133));
	bc5.idx.push_back(make_pair(134, 147));
	bc5.idx.push_back(make_pair(148, 165));
	bc5.idx.push_back(make_pair(164, 169));
	bc5.idx.push_back(make_pair(170, 194));
	bc5.idx.push_back(make_pair(195, 209));
	bc5.idx.push_back(make_pair(210, 224));
	bc5.idx.push_back(make_pair(225, 237));
	bc5.idx.push_back(make_pair(238, 255));
	bc5.idx.push_back(make_pair(256, 271));
	bc5.idx.push_back(make_pair(272, 286));
	bc5.idx.push_back(make_pair(287, 302));
	bc5.idx.push_back(make_pair(303, 317));
	bc5.idx.push_back(make_pair(318, 333));
	bc5.idx.push_back(make_pair(334, 348));
	bc5.idx.push_back(make_pair(349, 364));
	bc5.idx.push_back(make_pair(365, 378));
	bc5.idx.push_back(make_pair(379, 394));
	bc5.idx.push_back(make_pair(395, 408));
	bc5.idx.push_back(make_pair(409, 424));
	bc5.idx.push_back(make_pair(425, 439));
	bc5.idx.push_back(make_pair(440, 456));
	bc5.idx.push_back(make_pair(457, 470));
	bc5.idx.push_back(make_pair(471, 486));
	bc5.idx.push_back(make_pair(487, 900));
	checkVector.push_back(bc5);

	// now let's create a random block size based on rows



			
	// ------------The preprocessing code snippet ends here
	
	//this is the original cg code

	// compute r=b-Ax0
	vector<double> ax(rows);
	vector<double> r(rows);

	//tmp
	//unsigned int corruptIt = 0;
	// use this declaration for fault injection in A matrix value
	double corruptIt = 0;
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
	int jk = -1;

	for (iteration = 0; iteration < num_iterations + 500; iteration++) 
	{
		
		converge = inner_product(r, r);
		//std::cout<<sqrt(converge)<<std::endl;
		// jk is used to get control over the iteration fault injection
		jk++;
		/*
		// code to inject fault in the same run at various interval 
		if (jk < 423)
		{
		correction_frequency = 10;
		}
		else if (jk < 846)
		{
		correction_frequency = 15 ;
		}
		else if (jk > 847)
		{
			correction_frequency = 20;
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
		// remove below if condition to include correction freq
		//if(iteration == 5) {
		//if (run) {
			// comment out the self stabilizing loop
			
			run = false;
			//r=A*x, temp=A*p
			r = a*x;
			cout << "Running self stabilizing loop in iteration: " << iteration << endl;
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
		*/
		 // remove this if you want self stabalizing algorithm	

			// This code  snippet allows us inject fault in row, column, value 
			int ite = 4;
			if (iteration == ite)
			{
				//vector<unsigned int> asdf = a.IA();
				//vector<unsigned int> basdf = a.JA();
				
				vector<double> aaa = a();
				cout << "Iteration during which fault was injected " << iteration << endl;
				
				//unsigned int corruptIt = aaa[0];
				 corruptIt = aaa[0];
				 cout << "Value of in matrix A before corruption " << corruptIt << endl;
				 double corr = corruptIt;

				 //double old_alpha = corruptIt;
				 long long te = *(long long*)&corr;
				 //for (int i = 63; i >= 0; i--)
				// {
					 //std::cout << ((te >> i)&1);
					 //if (i == 63) std::cout << " ";
					 //else if (i == 52) std::cout << " ";
				// }

				 // create a temp number and then shift it to correct position
				 long long num1 = 1;
				 long long vv = num1 << 61;
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
				 corr = *(double*)&ans;

				//unsigned int flip = 16;

				//unsigned int reVal = corruptIt ^ flip;

				aaa[0] = corr;
				cout << "Value after corruption " << corr << endl;
				const vector<double> aaaaa = aaa;
				a.operator()(aaaaa);
				

				/*
				// row vector fault injection
				vector<unsigned int> asdf = a.IA();
				corruptIt = asdf[1];
				asdf[1] = 25;
				a.IA(asdf);
				*/
				/*
				// column vector fault injection
				vector<unsigned int> asdf = a.JA();
				corruptIt = asdf[1];
				asdf[1] = 25;
				a.JA(asdf);
				*/
			}
			

			// temp = A * p
			vector<double> q(rows);
			q = a * p;

			/*
			// below is the code for block checksum method to detect if there is any fault in Spmv
			auto sparseErrorScalar = sparseCheckSum * p;
			double sumAn = 0;
			int sumAnIt = 0;
			vector<double> sparseErrorVec(90,0);
			for (int i = 1; i <= q.size(); i++) {
				if (i % 10 == 0) {
					sparseErrorVec[sumAnIt] = sumAn;
					sumAn = q[i-1];
					sumAnIt++;
				}
				else {
					sumAn += q[i - 1];
				}
			}

			vector<int> errorIndexVec;
			// now compare the two result
			for (int i = 0; i < 90; i++) {
				auto tempErrV = sparseErrorVec[i] - sparseErrorScalar[i];
				if (abs(tempErrV) > 1000) {
					// record the error index
					errorIndexVec.push_back(i);
					cout << "caught the error at index " << i << endl;
				}

			}

			// now check for error and re-execute the required block
			for (auto&i : errorIndexVec) {
				int startRowPointer = i * 10;
				auto corAnsBlock = sparseCheckSum.multiplySpecific(p, startRowPointer);
				// update the ans
				std::copy_n(corAnsBlock.begin(), corAnsBlock.size(), q.begin()+startRowPointer);
			}
			*/

			//int aaaaaaa = 0;

			
			// below is recursive block checksum mathod
			// use abft check here to detect if there was any fault while SpMv was done
			auto cSum = inner_product(p, check1);
			auto cqSum = std::accumulate(q.begin(), q.end(), 0);

			//auto ft = sqrt(inner_product(x, x));

			bool error = false;
			// error is detected if
			if (abs(cSum - cqSum) > 100) {
				error = true;
				cout << "error has been detected" << endl;
				// call the error correction subroutine
				//auto a1 = inner_product(check2[0], p);
				//auto a2 = inner_product(check2[1], p);
				//int expI = 1;
				while (error) {
					//int k = pow(2, expI);
					//vector<double> ans(k,0);
					int k = 0;
					int j = 0;
					while (k < 5){
						auto f = checkVector[k].vec;
						vector<bool> ans;
						for (int i = 0; i < 2; i++) {
							auto r1 = checkVector[k].idx[ (2*j) + i].first;
							auto r2 = checkVector[k].idx[ (2 *j) + i].second;
							// check for the block where error exits
							//auto endp = vector<double>(p.begin()+ r1, p.begin()+r2);
							auto a2 = inner_product(checkVector[k].vec[(2 * j) + i], p);
							auto a3 = accumulate(q.begin() + r1, q.begin() + r2, 0);
							if (abs(a2 - a3) > 100) {
								ans.push_back(true);
								// if this is last iteration then fix the value
								if (k == 4) {
									// recompute the lost part
									auto r3 = r1;
									for (; r3 < r2; r3++ ) {
										q[r3] = inner_product(matrixA[r3], p);
									}
									error = false;
								}
							}
							else {
								ans.push_back(false);
							}	
						}
						k++;
						int id = 0;
						for (auto &i : ans) {
							if (i) {
								j = id;
								break;
							}
							id++;
						}
						ans.clear();
					}
				}
			}
			// recurisve block checksum approach ends here
			 

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
			if (iteration == ite)
			{
				
				vector<double> aaaa = a();				
				aaaa[0] = corruptIt;
				cout << "Value has been flipped back " << corruptIt << endl;
				const vector<double>aa = aaaa;
				a.operator()(aa);
				
				/*
				// restoring corrupted row index value
				vector<unsigned int> asdf = a.IA();
				asdf[1] = corruptIt;
				a.IA(asdf);
				*/
				/*
				// resoring corrupted column index value
				vector<unsigned int> asdf = a.JA();
				asdf[1] = corruptIt;
				a.JA(asdf);
				*/

			}
			/*/
			//// in-variant condition
			auto inValue = inner_product(p, q);
			if (inValue > 0.1) {
				run = true;
				cnt++;
				std::cout << "Error was catched " << cnt << endl;
			}
			*/
			
//remove this if you want to inlcude self stabalizing algorithm		
// Add { if you want self stabilizng loop 			
			
	}

	if (!isConvergence) {
		std::cout << "Does not converge after " << iteration << std::endl;
		std::cout << "residual is " << sqrt(converge) << std::endl;
		//std::cout << "Total number of times error was catched " << cnt << endl;
	}

	

}







int main(int argc, char* argv[]) 
{
	
	fstream read;

	read.open(R"(C:\Users\aadit\Documents\gr_30_30.mtx)", ios::in);

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
		matrixA[row - 1][col - 1] = value;
		if (row != col)
		{
			A(col - 1, row - 1, value);
			matrixA[col - 1][row - 1] = value;
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
	/*
	for (int i = 0; i < x.size(); i++)
	{
		cout << x[i] << endl;
	}
	*/
	system("pause");
	return 0;
	
}
