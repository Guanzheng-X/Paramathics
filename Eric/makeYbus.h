#ifndef MAKEYBUS_H
#define MAKEYBUS_H

#include "index.h"
#include "mpc.h"

typedef Matrix<std::complex<double>, Dynamic, Dynamic> MatrixComXd;
typedef Matrix<std::complex<double>, Dynamic, 1> VectorComXd;

VectorComXd getYff(mpc MPC)
{
	VectorComXd stat(MPC.nb);
	stat = MPC.Branch.col(10);
	//cout << stat.real() << endl;

	return stat;
}


#endif