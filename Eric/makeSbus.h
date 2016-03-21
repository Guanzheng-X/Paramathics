/*
	Created by Eric Xu

	This file creats S bus from case files.
*/

#ifndef MAKESBUS_H
#define MAKESBUS_H

#include <vector>
#include "index.h"
#include "mpc.h"

using std::vector;

VectorXcd makeSbus(mpc MPC)
{
	//---------------------- construct gbus matrix -----------------------

	// gubus matrix is a matrix answering the following questions:
	// 1. which generators are on?    2. what buses are they at?

	vector<int> temp;
	
	for(int i = 0; i < MPC.Gen.rows(); i++)
	{
		if(MPC.Gen(i, GEN_STATUS) > 0)
			temp.push_back(MPC.Gen(i, GEN_BUS));
	}

	VectorXd gbus(temp.size());

	for(int i = 0; i < temp.size(); i++)
	{
		gbus(i) = temp.at(i);
	}

	//----------------------------- end ------------------------------------

	//------------ construct connection matrix, Cg -------------------------

	// element i, j is 1 is 1 if gen on (j) at bus i is ON.

	MatrixXd Cg(MPC.nb, gbus.rows());  // nb * (the number of ON generators)
	Cg = MatrixXd::Zero(MPC.nb, gbus.rows());

	for(int i = 0; i < gbus.rows(); i++)
	{
		Cg(gbus(i) - 1, i) = 1;
	}

	//---------------------------- end -------------------------------------

	//----------------------- construct S bus ------------------------------

	vector<int> temp_1, temp_2;

	for(int i = 0; i < MPC.Gen.rows(); i++)
	{
		if(MPC.Gen(i, GEN_STATUS) > 0)
		{
			temp_1.push_back(MPC.Gen(i, PG));
			temp_2.push_back(MPC.Gen(i, QG));
		}
	}

	VectorXd gen_on_PG(temp_1.size()), gen_on_QG(temp_2.size());

	for(int i = 0; i < temp.size(); i++)
	{
		gen_on_PG(i) = temp_1.at(i);
		gen_on_QG(i) = temp_2.at(i);
	}

	MatrixXcd Cg_Complex(MPC.nb, gbus.rows());
	Cg_Complex = MatrixXcd::Zero(MPC.nb, gbus.rows());
	Cg_Complex.real() = Cg;

	VectorXcd gen_on_power(gbus.rows());
	gen_on_power.real() = gen_on_PG;
	gen_on_power.imag() = gen_on_QG;

	VectorXcd Sbus(MPC.nb);
	Sbus = Cg_Complex * gen_on_power;

	VectorXcd load_power(MPC.nb);
	load_power.real() = MPC.Bus.col(PD);
	load_power.imag() = MPC.Bus.col(QD);

	Sbus -= load_power;
	Sbus = Sbus / MPC.baseMVA;

	//---------------------------- end -------------------------------------

	return Sbus;
}

#endif