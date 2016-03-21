/*
	Created by Eric Xu
*/

#include "loadCase.h"
#include "makeYbus.h"
#include "makeSbus.h"


int main()
{
	mpc mycase = loadcase();
	MatrixXcd Ybus(makeYbus(mycase));
	VectorXcd Sbus(makeSbus(mycase));
	//cout << Sbus << endl;

	return 0;
}