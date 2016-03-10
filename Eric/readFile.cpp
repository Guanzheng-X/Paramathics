/*
	Created by Eric Xu
*/

#include <iostream>
#include <string.h>    //string object
#include <cstdlib>     //atoi function
#include <fstream>     //ifstream definition
#include <sstream>     //stringstream function
#include <vector>
#include <Eigen/Core>  //all data read from case file are stored in Eigen defined matrices
#include "mpc.h"

using namespace std;

mpc loadcase()
{
	mpc MP_Matrix;
	vector<vector<float> > busMatrix;
	vector<vector<float> > genMatrix;
	vector<vector<float> > branchMatrix;

	string case_line;
	ifstream casefile;
	casefile.open("case9.m");

	bool readBus = false;
	bool readGen = false;
	bool readBranch = false;

	if(casefile.is_open())
	{
		while(getline(casefile, case_line))
		{
			//cout << case_line << endl;
			if(case_line.find("mpc.baseMVA") != string::npos)
			{
				size_t index = case_line.find("=");
				string temp = case_line.substr(index + 1, string::npos);
				int baseMVA = atoi(temp.c_str());
				//cout << "baseMVA is " << baseMVA << endl;
			}

			if(case_line.find("mpc.bus ") != string::npos)
				readBus = true;

			if(case_line.find("mpc.gen ") != string::npos)
				readGen = true;

			if(case_line.find("mpc.branch ") != string::npos)
				readBranch = true;

			if((case_line.find("]") != string::npos) || (case_line.length() == 0))
			{
				readBus = false;
				readBranch = false;
				readGen = false;
			}

			if(readBus == true && case_line.find("mpc.bus ") == string::npos)     //read each line of the mpc.bus matrix in the case file
			{
				stringstream case_line_stream(case_line);
				vector<float> lineData;
				float value;

				while(case_line_stream >> value)
					lineData.push_back(value);

				busMatrix.push_back(lineData);
			}

			if(readGen == true && case_line.find("mpc.gen ") == string::npos)
			{
				stringstream case_line_stream(case_line);
				vector<float> lineData;
				float value;

				while(case_line_stream >> value)
					lineData.push_back(value);

				genMatrix.push_back(lineData);
			}

			if(readBranch == true && case_line.find("mpc.branch ") == string::npos)
			{
				stringstream case_line_stream(case_line);
				vector<float> lineData;
				float value;

				while(case_line_stream >> value)
					lineData.push_back(value);

				branchMatrix.push_back(lineData);
			}
		}

		MP_Matrix.Bus.resize(busMatrix.size(), (busMatrix.at(0)).size());
		MP_Matrix.Gen.resize(genMatrix.size(), (genMatrix.at(0)).size());
		MP_Matrix.Branch.resize(branchMatrix.size(), (branchMatrix.at(0)).size());

		for(int i = 0; i < busMatrix.size(); i++)
			for(int j = 0; j < (busMatrix.at(i)).size(); j++)
			{
				MP_Matrix.Bus(i, j) = (busMatrix.at(i)).at(j);
			}

		for(int i = 0; i < genMatrix.size(); i++)
			for(int j = 0; j < (genMatrix.at(i)).size(); j++)
			{
				MP_Matrix.Gen(i, j) = (genMatrix.at(i)).at(j);
			}

		for(int i = 0; i < branchMatrix.size(); i++)
			for(int j = 0; j < (branchMatrix.at(i)).size(); j++)
			{
				MP_Matrix.Branch(i, j) = (branchMatrix.at(i)).at(j);
			}

		MP_Matrix.print();
	}

	return MP_Matrix;
}

int main()
{
	loadcase();
}