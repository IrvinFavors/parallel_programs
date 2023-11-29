# parallel_programs

**Jacobi Method** <br />
Compile: mpicc mpi_jacobi.c -lm -o jacobi
Run: srun -n 4 ./jacobi

**(Temp)** <br />
Compile: gcc -g -o trap trap.c -lpthread
Run: ./trap 4 2 4 6

