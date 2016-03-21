/*
	Created by Eric Xu
*/

#include "loadCase.h"
#include "makeYbus.h"


int main()
{
	mpc mycase = loadcase();
	MatrixXcd Ybus(getYbus(mycase));
	cout << Ybus << endl;

	return 0;
}