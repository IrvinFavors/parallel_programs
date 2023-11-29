#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "timer.h"

/*** TODO 1: include the required MPI header file ***/
#include <mpi.h>

#define REAL float /* define for single or double precision */
#define DEFAULT_DIMSIZE 256

/************************************************************
 * program to solve a finite difference
 * discretization of the Helmholtz equation
 * (d2/dx2)u + (d2/dy2)u - alpha u = f
 * using Jacobi's iterative method.
 *
 * Input:   n - grid dimension in x direction
 *          m - grid dimension in y direction
 *          alpha - Helmholtz constant (always greater than 0.0)
 *          tol   - error tolerance for iterative solver
 *          relax - Successice over relaxation parameter
 *          mits  - Maximum iterations for iterative solver
 *
 * On output
 *       : u(n,m) - Dependent variable (solutions)
 *       : f(n,m) - Right hand side function
 *************************************************************/

/* subroutine initialize (n,m,alpha,dx,dy,u,f)
 ******************************************************
 * Initializes u and f
 * Assumes exact solution is u(x,y) = (1-x^2)*(1-y^2)
 ******************************************************/
void initialize(long n, long m, REAL alpha, REAL dx, REAL dy, REAL *u_p, REAL *f_p) {
    REAL (*u)[m] = (REAL(*)[m])u_p;
    REAL (*f)[m] = (REAL(*)[m])f_p;
	REAL x, y;

	/* Initialize initial condition and RHS */
	for (long i = 0; i < n; i++) {
		for (long j = 0; j < m; j++) {
			x = -1.0 + (i-1)*dx;
			y = -1.0 + (j-1)*dy;
			u[i][j] = 0.0;
			f[i][j] = -2.0*(1.0 - x*x) - 2.0*(1.0 - y*y) - alpha * (1.0 - x*x) * (1.0 - y*y);
		}
	}

}

/* subroutine error_check (n,m,alpha,dx,dy,u,f)
 ******************************************************
 * Checks error between numerical and exact solution
 ******************************************************/
double error_check(long n, long m, REAL alpha, REAL dx, REAL dy, REAL * u_p, REAL * f_p) {
	REAL (*u)[m] = (REAL(*)[m])u_p;
	REAL (*f)[m] = (REAL(*)[m])f_p;
	REAL x, y;
	REAL resid;
	double error = 0.0;

	for (long i = 0; i < n; i++) {
		for (long j = 0; j < m; j++) {
			x = -1.0 + (i-1)*dx;
			y = -1.0 + (j-1)*dy;
			resid = u[i][j] - (1.0 - x*x) * (1.0 - y*y);
			error = error + resid * resid;
		}
	}
	error = sqrt(error) / (n*m);

	return error;
}

/* subroutine jacobi (n,m,dx,dy,alpha,omega,u,f,tol,mits)
 ******************************************************************
 * Solves poisson equation on rectangular grid assuming:
 * (1) Uniform discretization in each direction, and
 * (2) Dirichlect boundary conditions
 *
 * Jacobi method is used in this routine
 *
 * Input:  n,m   Number of grid points in the X/Y directions
 *         dx,dy Grid spacing in the X/Y directions
 *         alpha Helmholtz eqn. coefficient
 *         omega Relaxation factor
 *         f(n,m) Right hand side function
 *         u(n,m) Dependent variable/Solution
 *         tol    Tolerance for iterative solver
 *         mits  Maximum number of iterations
 *
 * Output: u(n,m) - Solution
 *****************************************************************/
void jacobi_seq(long n, long m, REAL dx, REAL dy, REAL alpha, REAL omega, REAL tol, long mits, REAL *u_p, REAL *f_p) {
    REAL (*u)[m] = (REAL(*)[m])u_p;
    REAL (*f)[m] = (REAL(*)[m])f_p;
	REAL uold[n][m];

	REAL resid, error;

	/* Initialize coefficients */
	REAL ax = 1.0 / (dx*dx); // x-direction coef
	REAL ay = 1.0 / (dy*dy); // y-direction coeff
	REAL b = -2.0 / (dx*dx) - 2.0 / (dy*dy) - alpha; // central coeff

	/* Jacobi iteration */
	error = (10.0 * tol);
	long k = 1;
	while ((k <= mits) && (error > tol)) {
		error = 0.0;

		/* Copy new solution into old */
		for (long i = 0; i < n; i++) {
			for (long j = 0; j < m; j++) {
				uold[i][j] = u[i][j];
			}
		}

		for (long i = 1; i < n-1; i++) {
			for (long j = 1; j < m-1; j++) {
				resid = (ax * (uold[i-1][j] + uold[i+1][j]) + ay * (uold[i][j-1] + uold[i][j+1]) + b * uold[i][j] - f[i][j]) / b;
				u[i][j] = uold[i][j] - omega * resid;
				error = error + resid * resid;
			}
		}
		error = sqrt(error) / (n*m);
	
		/* Error check */
		if (k % 500 == 0)
		printf("Finished %ld iterations with error: %g\n", k, error);
	
		k = k + 1;
	} /*  End iteration loop */
	
	printf("Total Number of Iterations: %ld\n", k);
	printf("Residual: %.15g\n", error);

}

void jacobi_mpi(long n, long m, REAL dx, REAL dy, REAL alpha, REAL omega, REAL tol, long mits, REAL *u_p, REAL *f_p) {
	REAL (*u_global)[m] = (REAL (*)[])u_p;
	REAL (*f_global)[m] = (REAL (*)[])f_p;

	int myrank=0, numprocs=1;
	
	/*** TODO 2: add code to set numprocs and myrank to the appropriate values ***/
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	/* number of rows on each process */
	/*** TODO 3: add code to determine the number of rows that are handled by each process ***/
	long my_n = n / numprocs;

	/* memory management */
	REAL (*u)[m] = malloc((my_n+2)*m*sizeof(REAL));
	REAL (*f)[m] = malloc((my_n+2)*m*sizeof(REAL));
	REAL uold[my_n+2][m];

	/*** TODO 4: replace the following two lines by code that distributes the rows of u_global and f_global to the processes ***/
	MPI_Scatter(u_global, my_n*m, MPI_DOUBLE, &u[1][0], my_n*m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(f_global, my_n*m, MPI_DOUBLE, &f[1][0], my_n*m, MPI_DOUBLE, 0, MPI_COMM_WORLD);


	int li = /*** TODO 5a: add code to determine the index of the first row that this process has to update ***/ 2;
	int ui = /*** TODO 5b: add code to determine the index of the last row that this process has to update ***/ my_n;

	/* Initialize coefficients */
	REAL ax = 1.0 / (dx*dx); // x-direction coef
	REAL ay = 1.0 / (dy*dy); // y-direction coeff
	REAL b = -2.0 / (dx*dx) - 2.0 / (dy*dy) - alpha; // central coeff

	/* Jacobi iteration */
	REAL resid, my_error, error=10.0*tol;
	long k = 1;
	while ((k <= mits) && (error > tol)) {
		my_error = 0.0;

		/* Copy new solution into old */
		memcpy(uold[1], u[1], my_n*m*sizeof(REAL));

		/*** TODO 6: add code for the communication between processes
		if (myrank > 0) {
		 * all processes except rank 0 send their first row to the previous rank
		 * all processes except rank 0 receive a copy of the previous rank's last row
		}
		if (myrank < numprocs-1) {
		 * all processes ecept the last rank send their last row to the subsequent rank
		 * all processes ecept the last rank receive a copy of the subsequent rank's first row
		}
		***/
		MPI_Request send_req, recv_req;

		if (myrank > 0) {
			MPI_Isend(&u[1][0], m, MPI_REAL, myrank-1, 0, MPI_COMM_WORLD, &send_req);
			MPI_Irecv(&u[0][0], m, MPI_REAL, myrank-1, 0, MPI_COMM_WORLD, &recv_req);
			MPI_Wait(&send_req, MPI_STATUS_IGNORE);
			MPI_Wait(&recv_req, MPI_STATUS_IGNORE);
		}

		if (myrank < numprocs-1) {
			MPI_Isend(&u[my_n][0], m, MPI_REAL, myrank+1, 0, MPI_COMM_WORLD, &send_req);
			MPI_Irecv(&u[my_n+1][0], m, MPI_REAL, myrank+1, 0, MPI_COMM_WORLD, &recv_req);
			MPI_Wait(&send_req, MPI_STATUS_IGNORE);
			MPI_Wait(&recv_req, MPI_STATUS_IGNORE);
		}

		/* perform Jacobi iteration */		
		for (long i = li; i < ui; i++) {
			for (long j = 1; j < m-1; j++) {
				resid = (ax * (uold[i-1][j] + uold[i+1][j]) + ay * (uold[i][j-1] + uold[i][j+1]) + b * uold[i][j] - f[i][j]) / b;
				u[i][j] = uold[i][j] - omega * resid;
				my_error = my_error + resid * resid;
			}
		}
		
		/*** TODO 7: replace the following line by summing my_error across all processes */
		//error = my_error;
		MPI_Allreduce(&my_error, &error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		
		error = sqrt(error) / (n*m);

		if (myrank == 0) {
			if (k % 500 == 0)	
			printf("Finished %ld iterations with error: %g\n", k, error);
		}

		k = k + 1;
	} /*  End iteration loop */

	if (myrank ==0) {
		printf("Total Number of Iterations: %ld\n", k);
		printf("Residual: %.15g\n", error);
	}

	/*** TODO 8: replace the following line by code that collects the rows from the processes into u_global  ***/
	//memcpy(u_global, u[1], n*m*sizeof(REAL));
	MPI_Gather(&u[1][0], my_n*m, MPI_REAL, u_global, my_n*m, MPI_REAL, 0, MPI_COMM_WORLD);

	free(u);
	free(f);

}

#ifdef DEBUG
void print_array(char *title, char *name, REAL *A, long n, long m) {
	printf("%s:\n", title);
    for (long i = 0; i<n; i++) {
        for (long j = 0; j<m; j++) {
            printf("%s[%ld][%ld]:%f  ", name, i, j, A[i * m + j]);
        }
        printf("\n");
    }
    printf("\n");
}
#endif

int main(int argc, char * argv[]) {
	int numprocs=1, myrank=0;
	double start, finish, elapsed_seq, elapsed_mpi;
	double error_seq, error_mpi;

	/* default values */
    long n = DEFAULT_DIMSIZE;
    long m = DEFAULT_DIMSIZE;
    REAL alpha = 0.0543;
    REAL tol = 0.0000000001;
    REAL relax = 1.0;
    int mits = 5000;

	/*** TODO 9a: add code to set up the MPI execution environment ***/
	MPI_Init(&argc, &argv);

	/*** TODO 9b: add code to set numprocs and myrank to the appropriate values ***/
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	int err_exit = 0;
	if (myrank == 0) {
		fprintf(stderr,"Usage: jacobi [<n> <m> <alpha> <tol> <relax> <mits>]\n");
		fprintf(stderr, "\tn - grid dimension in x direction, default: %ld\n", n);
		fprintf(stderr, "\tm - grid dimension in y direction, default: n if provided or %ld\n", m);
		fprintf(stderr, "\talpha - Helmholtz constant (always greater than 0.0), default: %g\n", alpha);
		fprintf(stderr, "\ttol   - error tolerance for iterative solver, default: %g\n", tol);
		fprintf(stderr, "\trelax - Successice over relaxation parameter, default: %g\n", relax);
		fprintf(stderr, "\tmits  - Maximum iterations for iterative solver, default: %d\n", mits);

		if (argc == 2)      { sscanf(argv[1], "%ld", &n); m = n; }
		else if (argc == 3) { sscanf(argv[1], "%ld", &n); sscanf(argv[2], "%ld", &m); }
		else if (argc == 4) { sscanf(argv[1], "%ld", &n); sscanf(argv[2], "%ld", &m); sscanf(argv[3], "%g", &alpha); }
		else if (argc == 5) { sscanf(argv[1], "%ld", &n); sscanf(argv[2], "%ld", &m); sscanf(argv[3], "%g", &alpha); sscanf(argv[4], "%g", &tol); }
		else if (argc == 6) { sscanf(argv[1], "%ld", &n); sscanf(argv[2], "%ld", &m); sscanf(argv[3], "%g", &alpha); sscanf(argv[4], "%g", &tol); sscanf(argv[5], "%g", &relax); }
		else if (argc == 7) { sscanf(argv[1], "%ld", &n); sscanf(argv[2], "%ld", &m); sscanf(argv[3], "%g", &alpha); sscanf(argv[4], "%g", &tol); sscanf(argv[5], "%g", &relax); sscanf(argv[6], "%d", &mits); }
		else {
			/* any remaining command line arguments are ignored */
		}
		
		printf("\njacobi %ld %ld %g %g %g %d\n\n", n, m, alpha, tol, relax, mits);

		if (n%numprocs) {
			fprintf(stderr, "Number of rows %ld is not evenly divisible by number of processes %d.\n", n, numprocs);
			err_exit = 1;
		}
	}

	/*** TODO 10: add code to distribute err_exit to all processes ***/
	MPI_Bcast(&err_exit, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (err_exit) {
		/*** TODO 11: add code to terminate the MPI execution environment ***/
		MPI_Finalize();
		exit(-1);
	}

	REAL dx = 2.0 / (n-1); /* grid spacing in x direction */
    REAL dy = 2.0 / (m-1); /* grid spacing in y direction */

	REAL *u;
	REAL *f;

	if (myrank == 0) {
		printf("================================= Sequential Execution ======================================\n");
		
		u = malloc(n*m*sizeof(REAL));
		f = malloc(n*m*sizeof(REAL));
		
		initialize(n, m, alpha, dx, dy, u, f);
		GET_TIME(start);
		jacobi_seq(n, m, dx, dy, alpha, relax, tol, mits, u, f);
		GET_TIME(finish);
		elapsed_seq = finish - start;
		error_seq = error_check(n, m, alpha, dx, dy, u, f);

#if DEBUG
		print_array("Sequential Run", "u", u, n, m);
#endif
	}

	/*** TODO 12: add code to synchronize the processes ***/
	MPI_Barrier(MPI_COMM_WORLD);

    if (myrank == 0) {
     	printf("\n========================== Parallel MPI Execution (%d processes) =============================\n", numprocs);
		initialize(n, m, alpha, dx, dy, u, f); // re-initialize for parallel run
	}	

	start = /*** TODO 13a: add code to time the parallel execution ***/ MPI_Wtime();
	jacobi_mpi(n, m, dx, dy, alpha, relax, tol, mits, u, f);
	finish = /*** TODO 13b: add code to time the parallel execution ***/ MPI_Wtime();
	double elapsed_local = /*** TODO 13c: add code to calculate the elapsed time ***/ finish - start;

	/*** TODO 14: replace the following line to set elapsed_mpi to the maximum time across all processes ***/
	MPI_Reduce(&elapsed_local, &elapsed_mpi, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	if (myrank == 0) {
#if DEBUG
		print_array("Parallel Run", "u", u, n, m);
#endif

		double flops = mits*(n-2)*(m-2)*13;

		printf("------------------------------------------------------------------------------------------------------\n");
    	printf("Performance:\tRuntime(s)\tGFLOPS\t\tError\n");
    	printf("------------------------------------------------------------------------------------------------------\n");
    	printf("sequential:\t%.2f\t\t%.3f\t\t%g\n", elapsed_seq, flops / (1.0e9 * elapsed_seq), error_seq);
    	printf("mpi (%d procs):\t%.2f\t\t%.3f\t\t%g\n", numprocs, elapsed_mpi, flops / (1.0e9 * elapsed_mpi), error_check(n, m, alpha, dx, dy, u, f));
    	printf("------------------------------------------------------------------------------------------------------\n");
		printf("speedup:\t%.2f\n", elapsed_seq/elapsed_mpi);

		free(u);
		free(f);
	}

	/*** TODO 15: don't forget to terminate the MPI execution environment */
	MPI_Finalize();

    return 0;
}

/* you can compile and run the program on cisc372.cis.udel.edu using the following commands:
 *   mpicc mpi_jacobi.c -lm -o jacobi
 *   srun -n <numprocs> ./jacobi
 */
