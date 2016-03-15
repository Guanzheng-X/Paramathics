/*
	Created by Eric Xu
*/

#include "loadCase.h"
#include "makeYbus.h"


int main()
{
	mpc mycase = loadcase();
	getYff(mycase);

	return 0;
}