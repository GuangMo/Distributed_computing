#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "julia.h"

int main(int argc, char *argv[])
{
	int thread_count; // Change this to max threads before submitting

  int width, height, maxiter, flag;
  double x[2], y[2], c[2];
  char *image, *stats;
	double time;
  getParams(argv, &flag, c, x, y, &width, &height, &maxiter, &image, &stats);
	thread_count = omp_get_max_threads();
	double * end_times = (double *)calloc(thread_count,sizeof( double));
  double * start_times = (double *)calloc(thread_count,sizeof( double));
  int * max_iterations = (int *)calloc(thread_count,sizeof(int));


  int *iterations = (int*)malloc( sizeof(int) * width * height );
  assert(iterations);

  int n,i,tmp, phase;
  int * a;

  	#ifdef _OPENMP
  time = omp_get_wtime();
	int maxCount = julia(x, width, y, height, c, flag, maxiter, iterations,end_times, start_times, max_iterations, thread_count);
  time = omp_get_wtime() - time;
	#endif


	/* save our picture for the viewer */
  saveBMP(image, iterations, width, height);
 //Write stats file.

 	FILE *f;
	f = fopen(stats, "w");
	for (int i =0 ; i < thread_count; i ++){
		fprintf(f, "%d %f %d\n",i, end_times[i]-start_times[i], max_iterations[i]);
	}

	fclose(f);

	printf("Number of threads: %d\n", thread_count);
	printf("Execution time: %f\n", time);

  free(iterations);
  return 0;
}




