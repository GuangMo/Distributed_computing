
int julia(const float *x, int xres, const float *y, int yres, const float *c, int flag, int maxIterations,
	  int *iterations);

void getParams(char **argv,  int *flag, float *c, float *x, float *y, int *width, int *height,
	       int *maxiter, char **image, char **stats);


void iterations2color(int width, int height, const int *iterations, int max_iterations, int *image);

void saveBMP(char* filename, int* result, int width, int height);
