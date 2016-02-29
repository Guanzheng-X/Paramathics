#include <Core>
#include <Dense>
#include <iostream>
#include <Eigen>
#include "CSVRow.h"


using namespace Eigen;
using namespace std;


int main()
{
	
	double b = 11.1;
	double a = 0.001;
	
	vector<double> target;  // target vector
	vector<int> input;     // input vector

	// Read the CSV file for input Data
	// It is kept blank for now
	ifstream file("");
	CSVRow row; // csv parser object


	double predict = 0;			// input value for which target value has to be predicted

	cout << "Enter the input value for which you want to run Bayesian Curve fitting algorithm: ";
	cin >> predict;

	// parsing csv and collecting target data in vector target
	while (file >> row)
	{
		if (row[0] == " ")		// give the name symbol which is present in csv
		{

			target.push_back(stod(row[ ]));			// push the particular index where the value of symbol is present
		}
	}


	int size = target.size();

	// declaring input vector (x)
	for (int k = 0; k < size; k++)
	{
		input.push_back(k + 1);
	}


	// Transpose the input set
	// Now define only Phi

	Matrix <double, 5, 1> Phi;

	Phi(0, 0) = 1;
	Phi(1, 0) = predict;
	Phi(2, 0) = pow(predict, 2);
	Phi(3, 0) = pow(predict, 3);
	Phi(4, 0) = pow(predict, 4);

	// Now transpose of Phi
	Matrix <double, 1, 5> PhiT = Phi.transpose();

	// define result Phi
	Matrix < double, 5, 5 > PhiR;
	PhiR.fill(0);

	// Start the algorithm
	for (int i = 0; i < size; i++)
	{
		Matrix < double, 5, 1 > PhiN;
		PhiN(0, 0) = 1;
		PhiN(1, 0) = input[i];
		PhiN(2, 0) = pow(PhiN(1, 0), 2);
		PhiN(3, 0) = pow(PhiN(1, 0), 3);
		PhiN(4, 0) = pow(PhiN(1, 0), 4);

		

		Matrix < double, 5, 5 > PhTemp = PhiN * PhiT;

		PhiR = PhiR + PhTemp;
	}

	// Matrix Identity creating identity matrix 
	Matrix<double, Dynamic, Dynamic> I = Matrix <double, 5, 5>::Identity();

	// Formula 1.72 S Inverse

	Matrix <double, Dynamic, Dynamic> sInverse = (a*I) + (b*PhiR);

	

	// Inverse of Inverse

	Matrix <double, Dynamic, Dynamic> s = sInverse.inverse();


	// Calculate variance

	Matrix <double, Dynamic, Dynamic> Var1 = PhiT * s * Phi;

	double var = Var1(0, 0) + (1 / b);

	// calculate mean

	Matrix < double, 5, 1 > MeanR;
	MeanR.fill(0);

	for (int i = 0; i < size; i++)
	{
		Matrix < double, 5, 1 > Mean1;
		Matrix < double, Dynamic, Dynamic > Mean3;

		Mean1(0, 0) = 1;
		Mean1(1, 0) = input[i];
		Mean1(2, 0) = pow(Mean1(1, 0), 2);
		Mean1(3, 0) = pow(Mean1(1, 0), 3);
		Mean1(4, 0) = pow(Mean1(1, 0), 4);
		
		// Multiply by scalar (t)
		Mean3 = Mean1 * target[i];

		MeanR += Mean3;
	}

	Matrix <double, Dynamic, Dynamic> Mean;

	Mean = (b)* (PhiT)* (s)* (MeanR);

	double mean = Mean(0, 0);

	cout << "Mean is " << mean << endl << "Var is " << var << endl;

	return(0);
}