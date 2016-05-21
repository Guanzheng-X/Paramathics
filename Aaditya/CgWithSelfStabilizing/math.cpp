#include "math.hpp"

#include "sparse_matrix.hpp"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <algorithm>

namespace CGMethod{
    
    const double kEpsilon = 0.00000001;
    bool isConvergence=false;
    // Multiply A^tA with b.
    // void multiply_ata_b(
    //                     const SparseMatrix &a_trans, const SparseMatrix &a,
    //                     const double *b, double *c) {
    //     double *temp = new double[a.get_num_rows()];
        
    //     a.multiply_column(b, temp);
    //     a_trans.multiply_column(temp, c);
        
    //     delete [] temp;
    // }
    
    double inner_product(const double *a, const double *b, int n) {
        double result = 0.0;
        for (int c = 0; c < n; c++) {
            result += a[c] * b[c];
        }
        return result;
    }
    
    void add(const double *a, const double *b, double beta, int n, double *c) {
        for (int i = 0; i < n; i++) {
            c[i] = a[i] + b[i] * beta;
            //std::cout<<"add result "<<c[i]<<std::endl;
        }
    }
    
}
using namespace CGMethod;
void Math::conjugate_gradient(
                              const SparseMatrix &a, const double *b, int num_iterations, double *x, int flipBit, int itt, int iter_val) {
    SparseMatrix a_trans = a.transpose();
    
    double *a_trans_b = new double[a_trans.get_num_rows()];
    
   // a_trans.multiply_column(b, a_trans_b);
	//double old_x = 0;

    // Solve a_trans * a * x = a_trans_b.
    double *r = new double[a_trans.get_num_rows()];
    double *p = new double[a_trans.get_num_rows()];    
    double *q = new double[a_trans.get_num_rows()];    
    double *ri1 = new double[a_trans.get_num_rows()];
    
   /*
   double *Anorm0 = new double[a_trans.get_num_rows()];
   double *Anorm1 = new double[a_trans.get_num_rows()];
   double *Anorm2 = new double[a_trans.get_num_rows()];
   */

    //r=b-Ax0
    double *temp3= new double [a_trans.get_num_rows()];
    a.multiply_column(x, temp3);//A*r0
    add(b, temp3, -1.0, a_trans.get_num_rows(), r);
    
    // p = r
    std::copy(r, r + a_trans.get_num_rows(), p);
    
    // rs_old = r' * r
    //double rs_old = inner_product(r, r, a_trans.get_num_rows());
    
//    if (sqrt(rs_old) < kEpsilon) {
//        printf("The initial solution is good.\n");
//        return;
//    }
    
    int correction_frequency=iter_val;
    int iteration;
    double converge = 0.0;
	int jk = 0;
    for (iteration = 0; iteration < num_iterations+2000;iteration++) {
        converge = inner_product(r, r, a_trans.get_num_rows());
        //std::cout<<sqrt(converge)<<std::endl;
		
		/*
		if (iteration == 1)
		{
			double old_x = x[591];
			long long te = *(long long*)&x[591];

			// create a temp number and then shift it to correct position
			long long num1 = 1;
			long long vv = num1 << flipBit;
			//std::cout << std::endl;
			// Now xor the alpha value
			long long ans = te ^ vv;
			x[591] = *(double*)&ans;
		}
		*/
		jk++;
		
		if (jk < 141)
		{
			correction_frequency = 5 * iter_val;
		}
		else if (jk < 423)
		{
			correction_frequency = 10 * iter_val;
		}
		else if (jk < 846)
		{
			correction_frequency = 15* iter_val;
		}
		else if (jk < 1128)
		{
			correction_frequency = 20* iter_val;
		}
		else if (jk < 1408)
		{
			correction_frequency = 25* iter_val;
		}
		
		//std::cout << sqrt(converge) << std::endl;
        if (sqrt(converge)<kEpsilon) {
            printf("residual = %.20lf\n", sqrt(converge));
            std::cout<<"iteration number: "<<iteration<<std::endl;
            isConvergence=true;
            break;
        }
		
		
       // std::cout<<rs_old<<std::endl;
        if(iteration%correction_frequency==0%correction_frequency){
            //r=A*x, temp=A*p
            a.multiply_column(x, r);
            a.multiply_column(p, q);
            
            //ri=b-ri
            add(b, r, -1.0, a_trans.get_num_rows(), r);
            //alpha = ri^tpi/pi^tqi
            double alpha = inner_product(r, p, a_trans.get_num_rows())/inner_product(p, q, a_trans.get_num_rows());
            
           // std::cout<<alpha<<std::endl;
            
            //xi+1=xi+alpha*pi
            add(x, p, alpha, a_trans.get_num_rows(), x);
            //ri+1=ri-alpha*q
            add(r, q, -alpha, a_trans.get_num_rows(), ri1);
            
            //beta=-r^t*q/p^t*q
            double temp1 = inner_product(ri1, q, a_trans.get_num_rows());
            double temp2 = inner_product(p, q, a_trans.get_num_rows());
            double beta=-1*temp1/temp2;
            
           // std::cout<<"temp1 "<<temp1<<" temp2 "<<temp2<<std::endl;
            //pi+1=ri+1+beta*pi
            add(ri1, p, beta, a_trans.get_num_rows(), p);
            std::copy(ri1, ri1 + a_trans.get_num_rows(), r);
           // rs_old = beta;

        }else
		
		{
		
        // temp = A'A * p
        //multiply_ata_b(a_trans, a, p, temp);
        a.multiply_column(p, q);
            
        // alpha = r^2 / (p' * temp)
        double alpha = inner_product(r, r, a_trans.get_num_rows()) / inner_product(p, q, a_trans.get_num_rows());

		
		// Fault injection code
		if (iteration == itt)
		{
			// alpha = rs_old / (p' * temp)
			std::cout << "Look at the alpha before flip " << alpha << std::endl;
			double old_alpha = alpha;
			long long te =*( long long*) &alpha;
			for (int i = 63; i >= 0; i--)
			{
				//std::cout << ((te >> i)&1);
				//if (i == 63) std::cout << " ";
				//else if (i == 52) std::cout << " ";
			}
			
			// create a temp number and then shift it to correct position
			long long num1 = 1;
			long long vv = num1 << flipBit;
			//std::cout << std::endl;
			// Now xor the alpha value
			long long ans = te ^ vv;
			for (int i = 63; i >= 0; i--)
			{
				//std::cout << ((ans >> i) & 1);
				//if (i == 63) std::cout << " ";
				//else if (i == 52) std::cout << " ";
			}
			//std::cout << std::endl;
			alpha = *(double*)&ans;
			long long tea = *(long long*)&alpha;
			for (int i = 63; i >= 0; i--)
			{
				//std::cout << ((tea >> i) & 1);
				//if (i == 63) std::cout << " ";
				//else if (i == 52) std::cout << " ";
			}
			std::cout << "Number after bit "<< flipBit+1 <<" is flipped" << " in iteration " <<(itt) <<"= " << alpha << std::endl;
			std::cout << "Change in alpha: " << alpha - old_alpha << std::endl;
		}
		
		
            
            
        //std::cout<<alpha<<std::endl;
            
            
        // x = x + alpha * p
        add(x, p, alpha, a_trans.get_num_rows(), x);
        
        // ri+1 = ri - alpha * temp
        add(r, q, -alpha, a_trans.get_num_rows(), ri1);
        
        //beta=ri+1 ^2 / ri^2
        double beta = inner_product(ri1, ri1, a_trans.get_num_rows())/inner_product(r, r, a_trans.get_num_rows());
        //std::cout<<sqrt(rs_new)<<std::endl;
        //std::cout<<iteration<<std::endl;
        // Traditionally, if norm(rs_new) is small enough, the iteration can stop.
//        if (sqrt(rs_new) < kEpsilon) {
//            /// DEBUG ///
//            printf("rs_new = %.20lf\n", rs_new);
//            std::cout<<"iteration number: "<<iteration<<std::endl;
//            isConvergence=true;
//            break;
//        }
        
        // p = r + (rs_new / rs_old) * p
        add(ri1, p, beta, a_trans.get_num_rows(), p);
        
        // rs_old = rs_new
        //rs_old = rs_new;
        //ri=ri+1
        std::copy(ri1, ri1 + a_trans.get_num_rows(), r);
        
        /// DEBUG ///
        //printf("  rs_old = %.20lf\n", rs_old);

		/*
		// calculate the A norm 
		add(x0,x, -1, a_trans.get_num_rows(), Anorm0);
		a.multiply_column(Anorm0, Anorm1);
		double vvv = inner_product(Anorm1, Anorm0, a_trans.get_num_rows());
		double vvvv = sqrt(vvv);
		std::cout << vvvv << std::endl;		
		if (iteration == 1)
		{
			x[591] = old_x;
		}
		*/
        }
		
    }
    
    /// DEBUG ///
    //printf("rs_old = %.20lf\n", rs_old);
    if (!isConvergence) {
        std::cout<<"not converge at "<<iteration<<std::endl;
        std::cout<<"residual is "<<sqrt(converge)<<std::endl;
    }
//    if (sqrt(rs_old)<=kEpsilon){
//        std::cout<<"converged at "<<iteration<<std::endl;
//        std::cout<<"residual is "<<rs_old<<std::endl;
//    }
    delete [] r;
    delete [] p;
    delete [] q;
    delete [] a_trans_b;
}