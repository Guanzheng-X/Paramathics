/*
	Created by Guanzheng Xu
*/

#include "loadCase.h"
#include "newton.h"

int main()
{
	mpc mycase = loadcase();
	//mycase.print();

	newton(mycase);

	return 0;
}

// case9 -- 4 iterations (matpower 4 iterations)
// case39 -- 4 iterations (matpower 1 iterations)
// case57 -- 4 iterations (matpower 3 iterations)
// case118 -- 3 iterations (matpower 3 iterations)