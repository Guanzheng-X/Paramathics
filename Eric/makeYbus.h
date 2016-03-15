/*
	Created by Eric Xu
*/

#ifndef MAKEYBUS_H
#define MAKEYBUS_H

#include "index.h"
#include "mpc.h"


int getYff(mpc MPC)
{
	VectorXd stat_part(MPC.nl);
	VectorXcd stat(MPC.nl);

	stat_part = MPC.Branch.col(BR_STATUS);
	stat.real() = stat_part;
	stat.imag() = VectorXd::Zero(MPC.nl);
	cout << stat << endl;

	VectorXcd Ys(MPC.nl);
	Ys.real() = MPC.Branch.col(BR_R);
	Ys.imag() = MPC.Branch.col(BR_X);
	cout << Ys << endl;

	Ys = stat.cwiseQuotient(Ys);
	cout << Ys << endl;



	return 0;
}


#endif