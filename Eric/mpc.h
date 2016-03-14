/*
	Created by Eric Xu

	This file defines a mpc struct that contains all the data in a case file.
*/

#ifndef MPC_H
#define MPC_H

#include <Eigen/Core>
#include <iostream>

using namespace Eigen;

class mpc
{
public:
	MatrixXf Bus;
	MatrixXf Gen;
	MatrixXf Branch;

	void print()
	{
		std::cout << "The Bus matrix is:" << std::endl << Bus << std::endl;
		std::cout << "The Gen matrix is:" << std::endl << Gen << std::endl;
		std::cout << "The Branch matrix is:" << std::endl << Branch << std::endl;
	}
};


#endif  