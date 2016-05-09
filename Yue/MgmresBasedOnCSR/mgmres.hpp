//
//  mgmres.hpp
//  MgmresDemo
//
//  Created by 吴越 on 16/5/7.
//  Copyright © 2016年 Yue Wu. All rights reserved.
//

#ifndef mgmres_hpp
#define mgmres_hpp
#include "SparseMatrix.hpp"
#include <stdio.h>
using namespace myg;
//  mgmres.hpp
//
//  Thanks to Philip Sakievich for correcting mismatches between this HPP
//  file and the corresponding CPP file, 23 March 2016!
//
void atx_cr ( int n, int nz_num, int ia[], int ja[], double a[], double x[],
             double w[] );
void ax_cr ( int n, int nz_num, int ia[], int ja[], double a[], double x[],
            double w[] );
void diagonal_pointer_cr ( int n, int nz_num, int ia[], int ja[], int ua[] );
void ilu_cr ( int n, int nz_num, int ia[], int ja[], double a[], int ua[],
             double l[] );
void lus_cr ( int n, int nz_num, int ia[], int ja[], double l[], int ua[],
             double r[], double z[] );
void mult_givens ( double c, double s, int k, double g[] );
void pmgmres_ilu_cr (SparseMatrix<int, double> &mat,
                     double *x, double *rhs, int itr_max, int mr, double tol_abs,
                     double tol_rel );
double r8vec_dot ( int n, double a1[], double a2[] );
double *r8vec_uniform_01 ( int n, int &seed );
void rearrange_cr ( int n, int nz_num, int ia[], int ja[], double a[] );
void timestamp ( );

#endif /* mgmres_hpp */
