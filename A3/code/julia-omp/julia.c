#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "julia.h"

int julia(const double *x, int xres, const double *y, int yres, const double *c,
	  int flag, int maxIterations, int *iterations, double * end_times, \
	  double * start_times, int *max_iterations, int thread_count)
	{
	int i,j;
	double xi, yi, xgap, ygap, savex, savey, radius;
	int count;
	xgap = (x[1] - x[0]) / xres;
	ygap = (y[1] - y[0]) / yres;

	/* Open MP */
	#pragma omp parallel  num_threads(thread_count) \
	private(i,j,xi,yi,radius, count,savex,savey)
	{
	int maxIterationCount = 0;
	start_times[omp_get_thread_num()] = omp_get_wtime();

	for (j = 0; j < yres; j++)
    {
		#pragma omp for schedule(static,1)

		for (i = 0; i < xres; i++)
		{
		/* pixel to coordinates */
		xi = x[0] + i * xgap;
		yi = y[0] + j * ygap;

		/* initial value for the iteration */
		savex = flag*c[0] + (1-flag)*xi;
		savey = flag*c[1] + (1-flag)*yi;

		radius = 0.0;
		count = 0;
		while ( radius <= 4.0 && count < maxIterations )
		{
			double savex2 = xi;
			xi = xi * xi - yi * yi + savex;
			yi = 2.0f * savex2 * yi + savey;
			radius = xi * xi + yi * yi;
			count++;
		}

		if(count > maxIterationCount ){
		    maxIterationCount = count;
		}
		int *p = iterations + j*xres+i;

		/* If radius <= 4.0, we have hit maxIterations. The point is
		 likely in the set. */
		if (radius <= 4.0)
		{
			//assert(count==maxIterations);
			*p = 0;
		}
		else
		/* radius > 4.0. The point tends to infinity. We recognize
		   this at iteration number count */
		{
			*p = count;
		}
	}

}
		end_times[omp_get_thread_num()] = omp_get_wtime();
		max_iterations[omp_get_thread_num()] = maxIterationCount;
}
  return 1;
}
