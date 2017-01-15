#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#ifdef _OPENACC
#include <openacc.h>
#endif


#include "julia.h"

int main(int argc, char *argv[])
{
  int width, height, maxiter, flag;
  double x[2], y[2], c[2];
  char *image, *stats;
	double time;
  getParams(argv, &flag, c, x, y, &width, &height, &maxiter, &image, &stats);

  int *iterations = (int*)malloc( sizeof(int) * width * height );
  assert(iterations);

  int n,i,tmp, phase;
  int * a;

  /* compute set */
	#ifdef _OPENACC
	struct timeval start, end;
	gettimeofday(&start, NULL);
  int maxCount = julia(x, width, y, height, c, flag, maxiter, iterations);
	gettimeofday(&end, NULL);
  time = end.tv_sec-start.tv_sec+(end.tv_usec - start.tv_usec)*1.e-6;
	#endif


  FILE *f;
  f = fopen(stats, "w");

  fprintf(f, "%d %f %d\n",1, time,maxiter);


  fclose(f);
	/* save our picture for the viewer */
  saveBMP(image, iterations, width, height);

	printf("Execution time: %f\n", time);

  free(iterations);
  return 0;
}




