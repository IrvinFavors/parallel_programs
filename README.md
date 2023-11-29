# parallel_programs

**Jacobi Method** <br />
Compile: mpicc mpi_jacobi.c -lm -o jacobi <br />
Run: srun -n 4 ./jacobi

**Trapezoidal Rule** <br />
Compile: gcc -g -o trap pth_trap.c -lpthread <br />
Run: ./trap 4 2 4 6

