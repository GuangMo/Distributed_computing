
void iterations2color(int width, int height, int *iterations, int max_iterations, int *image)
{
  /* map iteratin number to color value */
  int num_colors = 255;
  int i,j;
  
  for (i=0; i<width; i++)
    for (j=0; j<height; j++)
      {
	int offset = j*width+i;
	int iter = iterations[offset];
	if (iter > max_iterations)
	  image[offset] = 0;
	else
	  image[offset] = iter % num_colors;
      }
}
