/*
	Created by Eric Xu

	This file defines a mpc struct that contains all the data in a case file.
*/

#ifndef MPC_H
#define MPC_H

#include <iostream>
#include <Eigen/Core>

using namespace Eigen;

using std::cout;
using std::endl;

typedef Matrix<std::complex<double>, Dynamic, Dynamic> MatrixComXd;
typedef Matrix<std::complex<double>, Dynamic, 1> VectorComXd;

class mpc
{

public:
//----------- member data -------------
	int baseMVA;
	int nb;
	int nl;
	MatrixComXd Bus;
	MatrixComXd Gen;
	MatrixComXd Branch;

//----------- member function ----------
	mpc()
	{
		baseMVA = 0;
		nb = 0;
		nl = 0;
	}

	mpc &operator=(mpc rightSide)
	{
		baseMVA = rightSide.baseMVA;
		nb = rightSide.nb;
		nl = rightSide.nl;
		Bus = rightSide.Bus;
		Gen = rightSide.Gen;
		Branch = rightSide.Branch;

		return *this;
	}

	~mpc()
	{
	}

	void print()
	{
		cout << "The baseMVA is:\n" << baseMVA << endl;
		cout << "\nThe Bus matrix is:\n" << Bus << endl;
		cout << "\nThe Gen matrix is:\n" << Gen << endl;
		cout << "\nThe Branch matrix is:\n" << Branch << endl;
	}
};


#endif  