/*
	Created by Eric Xu

	This file creats admittance matrix from case files.
*/

#ifndef MAKEYBUS_H
#define MAKEYBUS_H

#include "index.h"
#include "mpc.h"


MatrixXcd makeYbus(mpc MPC)
{
	VectorXd stat_part(MPC.nl);
	VectorXcd stat(MPC.nl);

	stat_part = MPC.Branch.col(BR_STATUS);
	stat.real() = stat_part;
	stat.imag() = VectorXd::Zero(MPC.nl);

	VectorXcd Ys(MPC.nl);
	Ys.real() = MPC.Branch.col(BR_R);
	Ys.imag() = MPC.Branch.col(BR_X);

	Ys = stat.cwiseQuotient(Ys);

	VectorXd Bc(MPC.nl);
	Bc = MPC.Branch.col(BR_B);
	Bc = stat_part.cwiseProduct(Bc);

	// tap ratio here is assumed to be 1

	VectorXcd Ytt(MPC.nl), Yff(MPC.nl), Yft(MPC.nl), Ytf(MPC.nl);
	Ytt.imag() = Bc / 2;
	Ytt.real() = Ys.real();
	Ytt.imag() = Ys.imag() + Ytt.imag();
	Yff = Ytt;
	Yft = -Ys;
	Ytf = -Ys;

	// Ysh = Psh + Qsh is assumed to be 0. Shunt admittance.

	VectorXd f(MPC.nl);    // list of "from" buses
	f = MPC.Branch.col(F_BUS);

	VectorXd t(MPC.nl);    // list of "to" buses
	t = MPC.Branch.col(T_BUS);

	MatrixXd Cf(MPC.nl, MPC.nb);   // column-major
	Cf = MatrixXd::Zero(MPC.nl, MPC.nb);

	for(int i = 0; i < MPC.nl; i++)
	{
		int j = int(MPC.nl * (*(f.data() + i) - 1));
		*(Cf.data() + j + i ) = 1;
	}

	MatrixXd Ct(MPC.nl, MPC.nb);   // column-major
	Ct = MatrixXd::Zero(MPC.nl, MPC.nb);

	for(int i = 0; i < MPC.nl; i++)
	{
		int j = int(MPC.nl * (*(t.data() + i) - 1));
		*(Ct.data() + j + i ) = 1;
	}

	MatrixXcd Yf(MPC.nl, MPC.nb), Yt(MPC.nl, MPC.nb);
	Yf = MatrixXcd::Zero(MPC.nl, MPC.nb);
	Yt = MatrixXcd::Zero(MPC.nl, MPC.nb);

	MatrixXd Yf_real(MPC.nl, MPC.nb), Yf_imag(MPC.nl, MPC.nb), 
				Yt_real(MPC.nl, MPC.nb), Yt_imag(MPC.nl, MPC.nb);
	Yf_real = MatrixXd::Zero(MPC.nl, MPC.nb);   
	Yf_imag = MatrixXd::Zero(MPC.nl, MPC.nb);
	Yt_real = MatrixXd::Zero(MPC.nl, MPC.nb);
	Yt_imag = MatrixXd::Zero(MPC.nl, MPC.nb);

	VectorXd Yff_real(MPC.nl), Yff_imag(MPC.nl), Yft_real(MPC.nl), Yft_imag(MPC.nl),
				Ytt_real(MPC.nl), Ytt_imag(MPC.nl), Ytf_real(MPC.nl), Ytf_imag(MPC.nl);
	Yff_real = Yff.real();		Yff_imag = Yff.imag();
	Yft_real = Yft.real();		Yft_imag = Yft.imag();
	Ytt_real = Ytt.real();		Ytt_imag = Ytt.imag();
	Ytf_real = Ytf.real();		Ytf_imag = Ytf.imag();

	for(int i = 0; i < MPC.nl; i++)
	{
		int j = int(MPC.nl * (*(f.data() + i) - 1));
		*(Yf_real.data() + j + i ) += *(Yff_real.data() + i);
		*(Yf_imag.data() + j + i ) += *(Yff_imag.data() + i);

		int k = int(MPC.nl * (*(t.data() + i) - 1));
		*(Yf_real.data() + k + i ) += *(Yft_real.data() + i);
		*(Yf_imag.data() + k + i ) += *(Yft_imag.data() + i);

		int u = int(MPC.nl * (*(f.data() + i) - 1));
		*(Yt_real.data() + u + i ) += *(Ytf_real.data() + i);
		*(Yt_imag.data() + u + i ) += *(Ytf_imag.data() + i);

		int v = int(MPC.nl * (*(t.data() + i) - 1));
		*(Yt_real.data() + v + i ) += *(Ytt_real.data() + i);
		*(Yt_imag.data() + v + i ) += *(Ytt_imag.data() + i);
	}

	Yf.real() = Yf_real;		Yf.imag() = Yf_imag;
	Yt.real() = Yt_real;		Yt.imag() = Yt_imag;

	MatrixXcd Ybus(MPC.nb, MPC.nb);
	Ybus = Cf.transpose() * Yf + Ct.transpose() * Yt;

	return Ybus;
}


#endif