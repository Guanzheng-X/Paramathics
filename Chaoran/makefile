
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules



ex4: ex4.o  chkopts
	-${CLINKER}  -o ex4 ex4.o ${PETSC_KSP_LIB}
	${RM} ex4.o

ex5: ex5.o  chkopts
	-${CLINKER}  -o ex5 ex5.o ${PETSC_KSP_LIB}
	${RM} ex5.o


#------------------------------------------------------------------------------------

runex4:
	-@${MPIEXEC} -n 1 ./ex4  -fin test.mtx -ksp_max_it  100
	-@${RM} -f ex4.tmp
runex5:
	-@${MPIEXEC} -n 1 ./ex5  -fin test.mtx 
	-@${RM} -f ex5.tmp


include ${PETSC_DIR}/lib/petsc/conf/test

