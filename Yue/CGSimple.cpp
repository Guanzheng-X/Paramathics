#include <iostream>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <fstream>
#include <vector>
using namespace Eigen;
using namespace std;


int main()
{ 	
	Eigen::setNbThreads(4);
	Eigen::initParallel();
	cout<<"thread: "<< Eigen::nbThreads()<<endl;
	int rows; int cols ; int nonezero; //get the size of input matrix

	std::vector<int > mulResult;
	//int mulResult[4]; // use to store the mul of Matrix and vector

	typedef Eigen::Triplet<double> T;//use to store row, col, v_ij

	vector<T> tripletList;	//use to store all the row

	/***************beging to construct the matrix****************/
	fstream read;

	read.open("test2.txt", ios::in);

	int row=0 ; int col=0; float value=0; //get the value
	read>>rows; read >> cols; read>>nonezero;
	
	cout<<"matrix has "<<rows<<" rows "<< cols<<" cols "<<nonezero<<" nonezero values"<<endl;
	
	SparseMatrix<double,RowMajor> mat(rows,cols);//declare a row major matrix
	tripletList.reserve(nonezero); //reserve the space to store nonezero value

	while (!read.eof()){
			read>>row;
			read>>col;
			read>>value;
			tripletList.push_back(T(row-1,col-1,value));

		
	}

mat.setFromTriplets(tripletList.begin(), tripletList.end()); //crate the matrix
mat.makeCompressed();
cout<<mat<<endl;
/*****************matrix construct end*****************/

/*****************matrx mul vector begins**************/
int * outerPtr=mat.outerIndexPtr();//point to the outerstart
int * innerPtr=mat.innerIndexPtr(); // pointer to the innerstart

SparseVector<double,RowMajor> vec(cols);//declare a sparse vector, row major
for(int i=0;i<cols;i++){
	vec.insert(i)=i;
}
 

if(mat.isCompressed()){  //computing the matrix multiplication
#pragma omp parallel for
	for(int i=0;i<rows;i=i+1){
		int temp=0;
		for(int j=*(outerPtr+i);j<*(outerPtr+i+1);j++){
		
			temp+=*(mat.valuePtr()+j) * *(innerPtr+j); 
		// cout<<"j "<<j<<" ";
		// cout<<"mulResult "<<mulResult[i]<<" ";
		// cout<<"value "<<*(mat.valuePtr()+j)<<" ";
		// cout<<"innerPtr "<<*(innerPtr+j)<<endl;
	}
	//cout<<mulResult[i]<<" ";
	mulResult.push_back(temp);
	
}
	cout<<endl;

 }
cout<<"result from CSR computation "<<endl;
for(vector<int>::iterator iter=mulResult.begin();iter!=mulResult.end();++iter){
	cout<<*iter<<endl; //使用 * 访问迭代器所指向的元素
}
  	

 cout<<"result from Eigen "<<endl<<mat*vec.transpose();

/*****************matrx mul vector end****************/



}