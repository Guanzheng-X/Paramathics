/* File:     mpi_conjugate_gradient.c
 * Purpose:  Implement the iterative method known as 
 *           "method of conjugate gradients" with C and MPI. 
 *           The method of cojugate gradients is used to solve linear systems 
 *           when the coefficient matrix is symmetric and positive definite.
 *
 * Compile:  mpicc -g -Wall -o mpi_conjugate_gradient mpi_conjugate_gradient.c -lm
 *
 * Run:      mpiexec -n <number of processes> ./mpi_conjugate_gradient \
 *           <the row number of matrix> <tolerance> <max number of iterations> \
 *
 * Input:    n: the order of the system
 *           tolerance
 *           max_iter: the maximum number of iterations
 *           A: the symmetric positive definite n x n matrix
 *           b: the right-hand side b. 
 * Output:   The number of iterations
 *           The time used by the solver (not including I/O)
 *           The solution to the linear system (maybe suppressed)
 *           The norm of the residual calculated by the conjugate gradient method
 *           The norm of the residual calculated directly from the definition
 8           of residual
 *
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

/* Local functions */
int Get_n(int argc, char* argv[]);
double Get_tolerance(int argc, char* argv[]);
int Get_max_iter(int argc, char* argv[]);
void Usage(char* prog_name);
void Read_matrix_vector(double matrix[], double vector[], int n, int p);
void Print_matrix(double matrix[], int n);
void Print_vector(double vector[], int n);
/* Linear Algebra */
double Norm(double vector[], int n);
double Dot_product(double a[], double b[], int n);
void daxpy(double a, double x[], double y[], int n);
void daypx(double a, double x[], double y[], int n);
void Matrix_vector(double A[], double x[], double b[], int n);

/* Linear Algebra Functions involving communication */
double Parallel_norm(double vector[], int n, MPI_Comm comm);
double Parallel_dot_product(double a[], double b[], int n, MPI_Comm comm);
void Parallel_daypx(double a, double loc_x[], double loc_y[], double y[],
                    int n, MPI_Comm comm);
void Parallel_matrix_vector(double A[], double x[], double loc_b[],
                            int m, int n, MPI_Comm comm);
void Parallel_conjugate_gradient(double A[], double x[], int n,
                                 int p, double tolerance, int max_iter,
                                 int my_rank, MPI_Comm comm, double* residual,
                                 double* elapsed, int* total_iter);

/*---------------------------------------------------------------------------*/
int main(int argc, char* argv[]) {

    int my_rank, m, n, p, max_iter, suppress, i, total_iter, iter;
    double tolerance, residual, elapsed;
    double* temp_A;
    double* local_A;
    double* b;
    double* temp_b;
    double* x;
    double* xx;
    double* xxx;
    double* y;
    double* yy;
    double* yyy;
    double* residuals;
    MPI_Comm comm;
    MPI_Init(&argc, &argv);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &my_rank);
    
    if (my_rank == 0) {
        n = Get_n(argc, argv);
        tolerance = Get_tolerance(argc, argv);
        max_iter = Get_max_iter(argc, argv);
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, comm);
    MPI_Bcast(&tolerance, 1, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&max_iter, 1, MPI_INT, 0, comm);
    MPI_Bcast(&suppress, 1, MPI_INT, 0, comm);
    if (n < 2 || tolerance == -1.0 || max_iter == -1 || suppress == -1) {
        if (my_rank == 0)
            Usage(argv[0]);
        MPI_Finalize();
        exit(0);  
    }   
    
    b = malloc(n*sizeof(double));
    temp_b = malloc(n*sizeof(double));
    x = malloc(n*sizeof(double));
    xx = malloc(n*sizeof(double));
    xxx = malloc(n*sizeof(double));
    y = malloc(n*sizeof(double));
    yy = malloc(n*sizeof(double));
    yyy = malloc(n*sizeof(double));
    local_A = malloc((n*n/p)*sizeof(double));
    
    if (my_rank == 0) {
        residuals = malloc(n*sizeof(double));
        temp_A = malloc((n*n)*sizeof(double));
        Read_matrix_vector(temp_A, b, n, p);
        printf("\n");
    }
     
    MPI_Scatter(temp_A, n*n/p, MPI_DOUBLE, local_A, n*n/p, MPI_DOUBLE, 0, comm);
    MPI_Bcast(b, n, MPI_DOUBLE, 0, comm);
    MPI_Bcast(x, n, MPI_DOUBLE, 0, comm);
    MPI_Bcast(xx, n, MPI_DOUBLE, 0, comm);
    MPI_Bcast(y, n, MPI_DOUBLE, 0, comm);
    MPI_Bcast(yy, n, MPI_DOUBLE, 0, comm);
    for (i = 0; i < n; i++)
        temp_b[i] = b[i];


    Parallel_conjugate_gradient(local_A, x, n, p, tolerance, max_iter,
                                my_rank, comm, &residual, &elapsed,
                                &total_iter);
    printf("Parallel Conjugate Gradient Method was %d.\n", total_iter);
/*
    Parallel_cg_error(local_A, y, c, n, p, tolerance, max_iter,
                                my_rank, comm, &residual, &elapsed,
                                &total_iter);     
*/ 
    iter = 50;
    while ( iter <= 1000)
	{
	   Parallel_conjugate_gradient(local_A, xx, n, p, tolerance, iter,
                                my_rank, comm, &residual, &elapsed,
                                &total_iter);

           Parallel_cg_error(local_A, yy, n, p, tolerance, iter,
                                my_rank, comm);
 
    for (i = 0; i < n; i++)
    {
	xxx[i] = x[i] - xx[i];
	yyy[i] = x[i] -yy[i];
    }
    Matrix_vector(temp_A, xxx, xx, n); 
    Matrix_vector(temp_A, yyy, yy, n);
    double a1 = sqrt(Dot_product(xxx, xx, n));
    double a2 = sqrt(Dot_product(yyy, yy, n));
    double a3 = (a2 -a1) / a1;
    //printf( "At iteration %d, diff is %lf\n", iter, a1 );
    //printf( "At iteration %d, diff is %lf\n", iter, a2 );
    printf( "%lf\n", a3 );
    //printf("The a norm error definition is %lf.\n", a3);
    iter+=50;

	}
    
      	


    if (my_rank == 0) {
        printf("The total number of iterations for the ");
        printf("Parallel Conjugate Gradient Method was %d.\n", total_iter);
        printf("The Parallel Conjugate Gradient Method completed in");
        printf(" %.14e seconds.\n", elapsed);
        if (!suppress) {
           // printf("The solution to the linear system is:\n");
           // Print_vector(x, n);    
        }    
        printf("The norm of the residual calculated by the ");
        printf("Conjugate Gradient Method was %lf.\n", residual);
        Matrix_vector(temp_A, x, residuals, n);
        for (i = 0; i < n; i++) 
            residuals[i] = b[i] - residuals[i];
        printf("The norm of the residual calculated directly by the ");
        printf("definition of residual is %lf.\n", sqrt(Norm(residuals, n)));
           
    }

    free(b);
    free(x);
    free(local_A);
    if (my_rank == 0) {
        free(residuals);
        free(temp_A);
    }
    MPI_Finalize();
    return 0;
}  /* main */
    
/*-----------------------------------------------------------------------------
 * Function:    Get_n
 * Purpose:     Get the input value n
 * Input args:  argc:  number of command line args
 *              argv:  array of command line args
 */
int Get_n(int argc, char* argv[]) {
    if (argc == 4 || argc == 5 || argc == 6) 
        return strtol(argv[1], NULL, 10);
    return 1;    
}  /* Get_n */    
    
/*-----------------------------------------------------------------------------
 * Function:    Get_tolerance
 * Purpose:     Get the input value tolerance
 * Input args:  argc:  number of command line args
 *              argv:  array of command line args
 * Output:      tolerance or -1.0 if the value given isn't greater than zero.
 */
double Get_tolerance(int argc, char* argv[]) {
    double tolerance; 
    if ((tolerance = atof(argv[2])) > 0.0)
        return tolerance;
    return -1.0;    
}  /* Get_tolerance */    
    
/*-----------------------------------------------------------------------------
 * Function:    Get_max_iter
 * Purpose:     Get the input value max_iter
 * Input args:  argc:  number of command line args
 *              argv:  array of command line args
 * Output:      max_iter or -1 if the value given isn't greater than zero.
 */
int Get_max_iter(int argc, char* argv[]) {
    int max_iter;
    if ((max_iter = strtol(argv[3], NULL, 10)) > 0)
        return max_iter;
    return -1;    
}  /* Get_max_iter */        

      
/*-----------------------------------------------------------------------------
 * Function:  Usage
 * Purpose:   Print a brief message explaining how the program is run.
 *            Then quit.
 * In arg:    prog:  name of the executable
  */
void Usage(char prog[]) {

   fprintf(stderr, "usage: %s <order of the system> <tolerance>\n", prog);
   fprintf(stderr, "       <max number of iterations> <suppress output>\n");
   fprintf(stderr, "\"suppress output\" is an optional argument entered as\n");
   fprintf(stderr, "the character 'n' to suppress output of the ");
   fprintf(stderr, "solution vector.\n");
   exit(0);
}  /* Usage */

/*---------------------------------------------------------------------------
 * Function:  Read_matrix_vector
 * Purpose:   Read in the matrix A & vector b
 * In arg:    matrix, vector, n
 * Out arg:   matrix, vector
 */
void Read_matrix_vector(double matrix[], double vector[], int n, int p) {
    FILE* file;
    int i, row, col, nnz;
    double val;
    file = fopen("1138_bus.mtx", "r");
    fscanf(file,"%d %d %d\n",&n,&n,&nnz);
    printf ("ROWS = %d, COLUMNS = %d, NO OF NON-ZEROS = %d\n",n,n,nnz);
    for (i=0; i<nnz; i++) 
    {
	fscanf(file,"%d %d %le\n",&row,&col,&val);
        /*In matrix market format, rows and columns starts from 1 */
	row = row-1; col = col-1 ;
        matrix[row*n+col] = val;
	matrix[col*n+row] = val;
     }
    fclose(file);
    for (i = 0; i < n; i++)
        vector[i] = 1;
}  /* Read_matrix_vector */

/*---------------------------------------------------------------------------
 * Function:  Print_matrix
 * Purpose:   Print the contents of the matrix
 * In args:   matrix, n
 */
void Print_matrix(double matrix[], int n) {
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            printf("%lf ", matrix[i*n+j]);
        printf("\n");
    }
}  /* Print_matrix */

/*---------------------------------------------------------------------------
 * Function:  Print_vector
 * Purpose:   Print the contents of the vector
 * In args:   vector, n
 */
void Print_vector(double vector[], int n) {
    int i;

    for (i = 0; i < n; i++) 
        printf("%lf ", vector[i]);
    printf("\n");
}  /* Print_vector */

/*---------------------------------------------------------------------------
 * Function:  Norm
 * Purpose:   Calculate norm of a vector
 * In args:   vector, n
 * Out args:  norm
 * NOTE: This does NOT sqrt the final result
 */
double Norm(double vector[], int n) {
    double norm;
    int i;
    norm = 0.0;
    for (i = 0; i < n; i++)
        norm += vector[i] * vector[i];
    return norm;
}   /* Norm */

/*---------------------------------------------------------------------------
 * Function:  Parallel_norm
 * Purpose:   Calculate norm of a vector
 * In args:   vector, n, comm
 * Out args:  norm
 * NOTE: This does NOT sqrt the final result
 */
double Parallel_norm(double vector[], int n, MPI_Comm comm) {
    double norm, local_norm;
    norm = 0.0;
    local_norm = Norm(vector, n);
    MPI_Allreduce(&local_norm, &norm, 1, MPI_DOUBLE, MPI_SUM, comm);
    return norm;
}   /* Parallel_norm */

/*---------------------------------------------------------------------------
 * Function:  Dot_product
 * Purpose:   Calculate dot product of two vectors
 * In args:   a, b, n
 * Out args:  dot_product
 */
double Dot_product(double a[], double b[], int n) {
    double dot_product;
    int i;
    
    dot_product = 0.0;
    for (i = 0; i < n; i++)
        dot_product += a[i] * b[i];
    return dot_product;
}   /* Dot_product */

/*---------------------------------------------------------------------------
 * Function:  Parallel_dot_product
 * Purpose:   Calculate dot product of two vectors and collects the sum on
 *            each processor.
 * In args:   local_a, local_b, n, comm
 * Out args:  dot_product
 */
double Parallel_dot_product(double local_a[], double local_b[], int n, 
                            MPI_Comm comm) {
    double dot_product, local_dot;
    
    dot_product = 0.0;
    local_dot = Dot_product(local_a, local_b, n);
    MPI_Allreduce(&local_dot, &dot_product, 1, MPI_DOUBLE, MPI_SUM, comm);
    return dot_product;
}   /* Parallel_dot_product */

/*---------------------------------------------------------------------------
 * Function:  daxpy
 * Purpose:   scalar a * x + y = y
 * In args:   a, x, y, n
 * Out args:  y
 */
void daxpy(double a, double x[], double y[], int n) {
    int i;
    
    for (i = 0; i < n; i++) 
        y[i] += a * x[i];
}   /* daxpy */

/*-----------------------------------------------------------------------------
 * Function:  daypx
 * Purpose:   scalar a * y + x = y
 * In args:   a, x, y, n
 * Out args:  y
 */
void daypx(double a, double x[], double y[], int n) {
    int i;
    
    for (i = 0; i < n; i++) 
        y[i] = a * y[i] + x[i];
}   /* daypx */

/*-----------------------------------------------------------------------------
 * Function:  Parallel_daypx
 * Purpose:   scalar a * y + x = y and gathers the complete solution vector
 *            onto each processor
 * In args:   a, loc_x, loc_y, y, n, comm
 * Out args:  loc_y, y
 */
void Parallel_daypx(double a, double loc_x[], double loc_y[], double y[],
                    int n, MPI_Comm comm) {

    daypx(a, loc_x, loc_y, n);
    MPI_Allgather(loc_y, n, MPI_DOUBLE, y, n, MPI_DOUBLE, comm);
}   /* Parallel_daypx */

/*-----------------------------------------------------------------------------
 * Function:  Matrix_vector
 * Purpose:   Calculate Matrix Vector multiplication
 * In args:   A, x, b, n
 * Out args:  b
 */
void Matrix_vector(double A[], double x[], double b[], int n) {
    int i, j;
    
    for (i = 0; i < n; i++) {
        b[i] = 0.0;
        for (j = 0; j < n; j++)
            b[i] += A[i*n+j] * x[j];
    }
}   /* Matrix_vector */


/*-----------------------------------------------------------------------------
 * Function:  Parallel_matrix_vector
 * Purpose:   Calculate Matrix Vector multiplication
 * In args:   A, x, loc_b, m, n, comm
 * Out args:  loc_b
 */
void Parallel_matrix_vector(double A[], double x[], double loc_b[],
                            int m, int n, MPI_Comm comm) {
    int i, j;
    
    for (i = 0; i < m; i++) {
        loc_b[i] = 0.0;
        for (j = 0; j < n; j++)
            loc_b[i] += A[i*n+j] * x[j];
    }
}   /* Parallel_matrix_vector */

/*-----------------------------------------------------------------------------
 * Function:  Parallel_conjugate_gradient
 * Purpose:   Calculate conjugate gradient
 * In args:   A, x, b, n, p, tolerance, max_iter, my_rank, comm, residual,
 *            elapsed, total_iter
 * Out args: x(on process zero), residual, elapsed, total_iter 
 */
void Parallel_conjugate_gradient(double A[], double x[], int n,
                                 int p, double tolerance, int max_iter,
                                 int my_rank, MPI_Comm comm, double* residual,
                                 double* elapsed, int* total_iter) {
                                
    int k, i, loc_n;
    double tt;
    double alpha_k, beta_k, alpha_k_b, beta_k_b, start, finish;
    double* b = malloc(n*sizeof(double)); 
    double* r_k = malloc(n/p * sizeof(double));
    double* r_k_prev = malloc(n/p * sizeof(double));
    double* p_k = malloc(n * sizeof(double));
    double* p_k_loc = malloc(n/p * sizeof(double));
    double* s_k = malloc(n/p * sizeof(double));
    double* x_k = malloc(n/p * sizeof(double));
    loc_n = n/p; 
    for (i = 0; i < n; i++)
            b[i] = 1;    
    /* Setup x_0 = 0 and r_0 = b */ 
    for (i = 0; i < loc_n; i++) 
        x_k[i] = 0.0;     
    MPI_Scatter(b, loc_n, MPI_DOUBLE, r_k, loc_n, MPI_DOUBLE, 0, comm);
    
    k = 0;
    start = MPI_Wtime();
    while ((Parallel_norm(r_k, loc_n, comm) > tolerance) && (k < max_iter)) {  
        k++;
        if (k == 1) {
            if (my_rank == 0)
                p_k = b; /* p_1 = r_0 = b */
            MPI_Bcast(p_k, n, MPI_DOUBLE, 0, comm);
            MPI_Scatter(p_k, loc_n, MPI_DOUBLE, p_k_loc, loc_n, MPI_DOUBLE, 0,
                        comm);
            /* s_k = A*p_k */
            Parallel_matrix_vector(A, p_k, s_k, loc_n, n, comm);
        } else {
            /* r_k-1 * r_k-1 */
            beta_k = Parallel_norm(r_k, loc_n, comm);
            /* r_k-2 * r_k-2 */
            beta_k_b = Parallel_norm(r_k_prev, loc_n, comm);
            beta_k = beta_k/beta_k_b;
            /* p_k = r_k-1 + B_k*P_k-1 */ 
            Parallel_daypx(beta_k, r_k, p_k_loc, p_k, loc_n, comm); 
            /* s_k = A*p_k */
            Parallel_matrix_vector(A, p_k, s_k, loc_n, n, comm);    
        }           
        /* r_k-1 * r_k-1 */
        alpha_k = Parallel_norm(r_k, loc_n, comm);
        /* p_k * s_k */
        alpha_k_b = Parallel_dot_product(p_k_loc, s_k, loc_n, comm);
        alpha_k = alpha_k/alpha_k_b;
        /* x_k = x_k-1 + a_k*p_k */
        daxpy(alpha_k, p_k_loc, x_k, loc_n);
        
        for (i = 0; i < loc_n; i++)
            r_k_prev[i] = r_k[i]; 
	
        /* r_k = r_k-1 - a_k * s_k */    
        daxpy(-1.0*alpha_k, s_k, r_k, loc_n);
	//tt = Parallel_norm(r_k, loc_n, comm);
	//if (my_rank == 0 && k % 100 == 0) printf( "At iteration %d, diff is %e\n", k, tt );
	
    }
    finish = MPI_Wtime();
    *elapsed = finish - start;
    /* x = x_k */
    MPI_Gather(x_k, loc_n, MPI_DOUBLE, x, loc_n, MPI_DOUBLE, 0, comm);
    *residual = sqrt(Parallel_norm(r_k, loc_n, comm));
    *total_iter = k;
  /*  
    free(r_k);
    free(r_k_prev);
    free(p_k);
    free(p_k_loc);
    free(s_k);
    free(x_k);

    r_k = NULL;
    r_k_prev = NULL;
    p_k = NULL;
    p_k_loc = NULL;
    s_k = NULL;
    x_k = NULL;
    */
}   /* Parallel_conjugate_gradient */


/*-----------------------------------------------------------------------------
 * Function:  Parallel_conjugate_gradient
 * Purpose:   Calculate conjugate gradient
 * In args:   A, x, b, n, p, tolerance, max_iter, my_rank, comm, residual,
 *            elapsed, total_iter
 * Out args: x(on process zero), residual, elapsed, total_iter 
 */
void Parallel_cg_error(double A[], double x[], int n,
                                 int p, double tolerance, int max_iter,
                                 int my_rank, MPI_Comm comm) {
                                 
    int k, i, loc_n;
    double tt;
    double alpha_k, beta_k, alpha_k_b, beta_k_b, start, finish;
    double* b = malloc(n*sizeof(double)); 
    double* r_k = malloc(n/p * sizeof(double));
    double* r_k_prev = malloc(n/p * sizeof(double));
    double* p_k = malloc(n * sizeof(double));
    double* p_k_loc = malloc(n/p * sizeof(double));
    double* s_k = malloc(n/p * sizeof(double));
    double* x_k = malloc(n/p * sizeof(double));
    loc_n = n/p; 
    for (i = 0; i < n; i++)
            b[i] = 1;     
    /* Setup x_0 = 0 and r_0 = b */ 
    for (i = 0; i < loc_n; i++) 
        x_k[i] = 0.0;     
    MPI_Scatter(b, loc_n, MPI_DOUBLE, r_k, loc_n, MPI_DOUBLE, 0, comm);
    
    k = 0;
    start = MPI_Wtime();
    while ((Parallel_norm(r_k, loc_n, comm) > tolerance) && (k < max_iter)) {  
        k++;
        if (k == 1) {
            if (my_rank == 0)
                p_k = b; /* p_1 = r_0 = b */
            MPI_Bcast(p_k, n, MPI_DOUBLE, 0, comm);
            MPI_Scatter(p_k, loc_n, MPI_DOUBLE, p_k_loc, loc_n, MPI_DOUBLE, 0,
                        comm);
            /* s_k = A*p_k */
	    p_k[950] = 5;
            Parallel_matrix_vector(A, p_k, s_k, loc_n, n, comm);
  	    p_k[950] = 1;
        } else {
            /* r_k-1 * r_k-1 */
            beta_k = Parallel_norm(r_k, loc_n, comm);
            /* r_k-2 * r_k-2 */
            beta_k_b = Parallel_norm(r_k_prev, loc_n, comm);
            beta_k = beta_k/beta_k_b;
            /* p_k = r_k-1 + B_k*P_k-1 */ 
            Parallel_daypx(beta_k, r_k, p_k_loc, p_k, loc_n, comm); 
            /* s_k = A*p_k */
            Parallel_matrix_vector(A, p_k, s_k, loc_n, n, comm);    
        }           
        /* r_k-1 * r_k-1 */
        alpha_k = Parallel_norm(r_k, loc_n, comm);
        /* p_k * s_k */
        alpha_k_b = Parallel_dot_product(p_k_loc, s_k, loc_n, comm);
        alpha_k = alpha_k/alpha_k_b;
        /* x_k = x_k-1 + a_k*p_k */
        daxpy(alpha_k, p_k_loc, x_k, loc_n);
        
        for (i = 0; i < loc_n; i++)
            r_k_prev[i] = r_k[i]; 
	
        /* r_k = r_k-1 - a_k * s_k */    
        daxpy(-1.0*alpha_k, s_k, r_k, loc_n);
	//tt = Parallel_norm(r_k, loc_n, comm);
	//if (my_rank == 0 && k % 100 == 0) printf( "At iteration %d, diff is %e\n", k, tt );
	
    }
    finish = MPI_Wtime();
    /* x = x_k */
    MPI_Gather(x_k, loc_n, MPI_DOUBLE, x, loc_n, MPI_DOUBLE, 0, comm);


    free(r_k);
    free(r_k_prev);
    free(p_k);
    free(p_k_loc);
    free(s_k);
    free(x_k);


}   /* Parallel_conjugate_gradient */
