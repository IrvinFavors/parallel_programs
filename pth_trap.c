/* File:    trap.c
 * Purpose: Calculate definite integral using trapezoidal 
 *          rule.
 *
 * Input:   a, b, n
 * Output:  Estimate of integral from a to b of f(x)
 *          using n trapezoids.
 *
 * Compile: gcc -g -Wall -o trap trap.c
 * Usage:   ./trap
 *
 * Note:    The function f(x) is hardwired.
 *
 * IPP2:    Section 3.2.1 (pp. 101 and ff.) and 5.2 (p. 228)
 */

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

int thread_count;

double f(double x);    /* Function we're integrating */
double  integral;   /* Store result in integral   */
double  a, b;       /* Left and right endpoints   */
int     n;          /* Number of trapezoids       */
double  h;          /* Height of trapezoids       */
pthread_mutex_t mutex;

void* Trap(void* rank);

int main(int argc, char* argv[]) {
   long thread;
   pthread_t* thread_handles;
   thread_count = strtol(argv[1], NULL, 10);
   thread_handles = malloc(thread_count*sizeof(pthread_t));


   printf("Enter a, b, and n\n");
   scanf("%lf", &a);
   scanf("%lf", &b);
   scanf("%d", &n);

   h = (b-a)/n;
   integral = 0.0;
   pthread_mutex_init(&mutex, NULL);

   for (thread = 0; thread < thread_count; thread++) {
      pthread_create(&thread_handles[thread], NULL, &Trap, (void*)thread);
   }
   for (thread = 0; thread < thread_count; thread++) {
      pthread_join(thread_handles[thread], NULL);
   }
   
   integral *= h;
   
   printf("With n = %d trapezoids, our estimate\n", n);
   printf("of the integral from %f to %f = %.15f\n",
      a, b, integral);

   free(thread_handles);
   pthread_mutex_destroy(&mutex);

   return 0;
}  /* main */

/*------------------------------------------------------------------
 * Function:    Trap
 * Purpose:     Estimate integral from a to b of f using trap rule and
 *              n trapezoids
 * Input args:  a, b, n, h
 * Return val:  Estimate of the integral 
 */
void* Trap(void* rank) {
   long my_rank = (long)rank;
   double local_a = a + my_rank * h * n / thread_count;
   double local_b = a +  (my_rank + 1) * h * n / thread_count;
   int local_n = n / thread_count;
   double my_integral = (f(local_a) + f(local_b)) / 2.0;
   int k;

   for (k = 1; k <= local_n-1; k++) {
     my_integral += f(local_a+k*h);
   }

   pthread_mutex_lock(&mutex);
   integral += my_integral;
   pthread_mutex_unlock(&mutex);

   return NULL;
}  /* Trap */

/*------------------------------------------------------------------
 * Function:    f
 * Purpose:     Compute value of function to be integrated
 * Input args:  x
 */
double f(double x) {
   double return_val;

   return_val = x*x;
   return return_val;
}  /* f */
