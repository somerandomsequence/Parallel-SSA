build:
	mpicc master_slave.c -o master_slave

test:
	mpirun -np 4 master_slave

prep:
	R --vanilla < preprocess.R
