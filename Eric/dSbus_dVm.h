/*
	Created by Guanzheng Xu

	This file performs partial derivative.
*/

#ifndef DSBUS_DVM_H
#define DSBUS_DVM_H

#include <math.h>    /* sqrt */
#include "index.h"
#include "mpc.h"

MatrixXcd dSbus_dVm(MatrixXcd &Ybus, VectorXcd &V)     // V is the complex voltage of each bus
{
	int n = V.rows();

	VectorXcd Ibus(n);
	Ibus = Ybus * V;

	// in matlab, Ybus should be checked for its storage class (sparse or dense, but here this procedure is ignored and will be added soon)
	MatrixXcd diagV(n, n), diagIbus(n, n), diagVnorm(n, n);
	diagV = MatrixXcd::Zero(n, n);    
	diagIbus = MatrixXcd::Zero(n, n);
	diagVnorm = MatrixXcd::Zero(n, n);

	for(int i = 0; i < n; i++)
	{
		diagV(i, i) = V(i);
		diagIbus(i, i) = Ibus(i);
		double absolute = sqrt( V(i).real() * V(i).real() + V(i).imag() * V(i).imag() );
		diagVnorm(i, i) = V(i) / absolute;   // this line needs to be checked further
	}

	MatrixXcd dSbus_dVm(V.rows(), V.rows()), dSbus_dVm_1(V.rows(), V.rows()), dSbus_dVm_2(V.rows(), V.rows());
	dSbus_dVm_1 = (Ybus * diagVnorm).conjugate();
	dSbus_dVm_1 = diagV * dSbus_dVm_1;
	dSbus_dVm_2 = diagIbus.conjugate();
	dSbus_dVm_2 = dSbus_dVm_2 * diagVnorm;
	dSbus_dVm = dSbus_dVm_1 + dSbus_dVm_2;

	return dSbus_dVm;
}

#endif