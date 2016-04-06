/*
	Created by Guanzheng Xu

	This file defines a mpc struct that contains all the data in a case file.
*/

#ifndef MPC_H
#define MPC_H

#include <iostream>
#include <Eigen/Core>

using namespace Eigen;

using std::cout;
using std::endl;


class mpc
{

public:
//----------- member data -------------
	int baseMVA;
	int nb;
	int nl;
	MatrixXd Bus;
	MatrixXd Gen;
	MatrixXd Branch;

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