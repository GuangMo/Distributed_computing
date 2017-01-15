#include <stdio.h>
#include <stdlib.h>

void getParams(char **argv, int *flag, double *c, double *x, double *y, int *width, int *height,
	       int *maxiter, char **image, char **stats)
{
  *flag = strtol(argv[1], NULL, 0);

  c[0]  = strtod(argv[2], NULL);
  c[1]  = strtod(argv[3], NULL);

  x[0] = strtod(argv[4], NULL);
  x[1] = strtod(argv[5], NULL);
  y[0] = strtod(argv[6], NULL);
  y[1] = strtod(argv[7], NULL);

  *width  = strtol(argv[8], NULL, 0);
  *height = strtol(argv[9], NULL, 0);
  *maxiter= strtol(argv[10], NULL, 0);

  *image = argv[11];
  *stats = argv[12];
}
