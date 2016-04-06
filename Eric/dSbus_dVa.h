/*
	Created by Guanzheng Xu

	This file performs partial derivative.
*/

#ifndef DSBUS_DVA_H
#define DSBUS_DVA_H

#include <math.h>    /* sqrt */
#include "index.h"
#include "mpc.h"

MatrixXcd dSbus_dVa(MatrixXcd &Ybus, VectorXcd &V)     // V is the complex voltage of each bus
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

	MatrixXcd dSbus_dVa(V.rows(), V.rows());
	dSbus_dVa = (diagIbus - Ybus * diagV).conjugate();
	dSbus_dVa = diagV * dSbus_dVa;

	MatrixXd dSbus_dVa_temp(V.rows(), V.rows());
	dSbus_dVa_temp = - dSbus_dVa.imag();
	dSbus_dVa.imag() = dSbus_dVa.real();
	dSbus_dVa.real() = dSbus_dVa_temp;

	return dSbus_dVa;
}

#endif