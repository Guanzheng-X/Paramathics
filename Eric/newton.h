/*
	Created by Guanzheng Xu

	Newton Raphson Solver.
*/

#ifndef NEWTON_H
#define NEWTON_H

#include <Eigen/Dense>
#include "makeYbus.h"
#include "makeSbus.h"
#include "dSbus_dVm.h"
#include "dSbus_dVa.h"

const double PI = 3.1415926;
const double tol = 0.00000001;   // used to terminate NR iteration
const int max_it = 10;       // the max number of iterations

VectorXcd newton(mpc &MPC)
{
	bool converged = 0; // indicate whether iteration is converged or not
	int counter = 0;     // the current number of iterations

	// initialize Vm and Va
	VectorXd Va(MPC.nb), Vm(MPC.nb);
	Va = MPC.Bus.col(VA);    // measured in degrees
	Va = Va * PI / 180.0;    // measured in radians
	Vm = MPC.Bus.col(VM);

	// calculate complex vector V
	VectorXcd V(MPC.nb);

	V.real() = Vm.array().cwiseProduct(Va.array().cos());
	V.imag() = Vm.array().cwiseProduct(Va.array().sin());

	// calculate Ybus & Sbus
	MatrixXcd Ybus(makeYbus(MPC));
	VectorXcd Sbus(makeSbus(MPC));

	// store the indecies of PV buses into vector pv
	vector<int> temp;
	
	for(int i = 0; i < MPC.Bus.rows(); i++)
	{
		if(MPC.Bus(i, BUS_TYPE) == 2)  // find PV buses
			temp.push_back(MPC.Bus(i, BUS_I));
	}

	VectorXd pv(temp.size());

	for(int i = 0; i < temp.size(); i++)
	{
		pv(i) = temp.at(i);
	}

	// store the indecies of PQ buses into vector pq
	temp.resize(0);

	for(int i = 0; i < MPC.Bus.rows(); i++)
	{
		if(MPC.Bus(i, BUS_TYPE) == 1)  // find PQ buses
			temp.push_back(MPC.Bus(i, BUS_I));
	}

	VectorXd pq(temp.size());

	for(int i = 0; i < temp.size(); i++)
	{
		pq(i) = temp.at(i);
	}

	// store the index of reference bus into vector ref
	temp.resize(0);

	for(int i = 0; i < MPC.Bus.rows(); i++)
	{
		if(MPC.Bus(i, BUS_TYPE) == 3)  // find PQ buses
			temp.push_back(MPC.Bus(i, BUS_I));
	}

	VectorXd ref(temp.size());

	for(int i = 0; i < temp.size(); i++)
	{
		ref(i) = temp.at(i);
	}

	// establish the concatenation of pv and pq for further use
	VectorXd pv_pq(pv.rows() + pq.rows());
	pv_pq << pv, pq;

	// set up indexing for updating V
	int npv = pv.rows();
	int npq = pq.rows();

	int j1 = 1;			int j2 = npv;
	int j3 = j2 + 1;	int j4 = j2 + npq;
	int j5 = j4 + 1;	int j6 = j4 + npq;

	// evaluate F(x0)
	VectorXcd mis(MPC.nb);
	mis = (V.array().cwiseProduct((Ybus * V).array().conjugate())).matrix() - Sbus; //mis = V .* conj(Ybus * V) - Sbus

	temp.resize(0);   // release memory
	vector<double> F_temp(0);

	for(int i = 0; i < pv.rows(); i++)
	{
		F_temp.push_back(mis.real()[pv(i) - 1]);
	}

	for(int i = 0; i < pq.rows(); i++)
	{
		F_temp.push_back(mis.real()[pq(i) - 1]);
	}

	for(int i = 0; i < pq.rows(); i++)
	{
		F_temp.push_back(mis.imag()[pq(i) - 1]);
	}

	VectorXd F(2 * npq + npv);
	for(int i = 0; i < F_temp.size(); i++)
	{
		F(i) = F_temp.at(i);
	}

	F_temp.resize(0);   // release memory

	// check tolerance
	double normF = 0;
	normF = F.array().cwiseAbs().maxCoeff();

	if(normF < tol)
	{
		converged = 1;    // converged
		cout << "Converged in " << counter << " iterations." << endl << endl;
		cout << "Vm is the following: " << endl;
		cout << Vm << endl << endl;
		cout << "Va is the following: " << endl;
		cout << Va * 180.0 / PI << endl << endl;
	}

	// do Newton iterations
	while(!converged && (counter < max_it))
	{
		// update iteration counter
		counter++;

		// evaluate Jacobian
		MatrixXcd Sbus_Vm(dSbus_dVm(Ybus, V));  // matpower formula
		MatrixXcd Sbus_Va(dSbus_dVa(Ybus, V));  // matpower formula

		MatrixXd j11(npv + npq, npv + npq), j12(npv + npq, npq);
		MatrixXd j21(npq, npv + npq), j22(npq, npq);

		// construct matrix j11
		for(int i = 0; i < pv_pq.rows(); i++)
			for(int j = 0; j < pv_pq.rows(); j++)
			{
				j11(i, j) = Sbus_Va.real()(pv_pq(i) - 1, pv_pq(j) - 1);  
			}

		// construct matrix j12
		for(int i = 0; i < pv_pq.rows(); i++)
			for(int j = 0; j < pq.rows(); j++)
			{
				j12(i, j) = Sbus_Vm.real()(pv_pq(i) - 1, pq(j) - 1);  
			}

		// construct matrix j21
		for(int i = 0; i < pq.rows(); i++)
			for(int j = 0; j < pv_pq.rows(); j++)
			{
				j21(i, j) = Sbus_Va.imag()(pq(i) - 1, pv_pq(j) - 1);  
			}

		// construct matrix j22
		for(int i = 0; i < pq.rows(); i++)
			for(int j = 0; j < pq.rows(); j++)
			{
				j22(i, j) = Sbus_Vm.imag()(pq(i) - 1, pq(j) - 1);  
			}
		
		// construct Jacobian
		MatrixXd J(j11.rows() + j21.rows(), j11.cols() + j12.cols());
		J << j11, j12, j21, j22;

		// computer update step by using LU Decomposition
		MatrixXd dx = (-J).lu().solve(F);

		// update voltage
		for(int i = 0; i < npv; i++)
		{
			Va(pv(i) - 1) = Va(pv(i) - 1) + dx(i);
		}

		for(int i = 0; i < npq; i++)
		{
			Va(pq(i) - 1) = Va(pq(i) - 1) + dx(npv + i);
			Vm(pq(i) - 1) = Vm(pq(i) - 1) + dx(npv + npq + i);
		}

		V.real() = Vm.array().cwiseProduct(Va.array().cos());
		V.imag() = Vm.array().cwiseProduct(Va.array().sin());

		// evaluate F(x)
		mis = (V.array().cwiseProduct((Ybus * V).array().conjugate())).matrix() - Sbus; //mis = V .* conj(Ybus * V) - Sbus

		for(int i = 0; i < pv.rows(); i++)
		{
			F_temp.push_back(mis.real()[pv(i) - 1]);
		}

		for(int i = 0; i < pq.rows(); i++)
		{
			F_temp.push_back(mis.real()[pq(i) - 1]);
		}

		for(int i = 0; i < pq.rows(); i++)
		{
			F_temp.push_back(mis.imag()[pq(i) - 1]);
		}

		for(int i = 0; i < F_temp.size(); i++)
		{
			F(i) = F_temp.at(i);
		}

		F_temp.resize(0);   // release memory

		// check for tolerance
		normF = F.array().cwiseAbs().maxCoeff();

		if(normF < tol)
		{
			converged = 1;    // converged
			cout << "Converged in " << counter << " iterations." << endl << endl;
			cout << "Vm is the following: " << endl;
			cout << Vm << endl << endl;
			cout << "Va is the following: " << endl;
			cout << Va * 180.0 / PI  << endl << endl;
		}

	}

	if(!converged)
		cout << "Not Converged." << endl;

	return V;
}

#endif