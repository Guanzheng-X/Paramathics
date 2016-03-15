/*
	Created by Eric Xu
*/

#include "loadCase.h"
#include "makeYbus.h"

typedef Matrix<std::complex<double>, Dynamic, Dynamic> MatrixComXd;
typedef Matrix<std::complex<double>, Dynamic, 1> VectorComXd;

int main()
{
	mpc mycase = loadcase();
	getYff(mycase);

	return 0;
}