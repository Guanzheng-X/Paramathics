This is Chaoran's folder.

Propose: Implement Conjugate Gradient Method and calculate the Relative A Norm. 

File:     mpi_conjugate_gradient.c

Compile:  mpicc -g -Wall -o mpi_conjugate_gradient mpi_conjugate_gradient.c -lm

Run:      mpiexec -n <number of processes> ./mpi_conjugate_gradient \
           <the row number of matrix> <tolerance> <max number of iterations> \
