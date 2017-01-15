#include<stdio.h>
#include<stdlib.h>

typedef struct {
    unsigned char r;
    unsigned char g;
    unsigned char b;
} RGB;

static RGB table[256] = {{0,0,0}}; /* Initialize to all black */

void initColours()
{
  
  int i;
  int red, green, blue;
  
  RGB black = {0,0,0};
  RGB white = {255,255,255};
  
  table[0]   = white; /* Windows reserves first color for white */
  table[255] = black; /* Windows reserves last color as black   */
  
  i = 20; /* first 20 and last 20 are reserved */
  for (red = 0; red <= 255; red+= 51) 
    {/* the six values of red */
      for (green = 0; green <= 255; green += 51)
	{
	  for (blue = 0; blue <= 255; blue+= 51) 
	    {
	      table[i].r = red;
	      table[i].g = green;
	      table[i].b = blue;
	      ++i;
	    }
	}
    }
}

void saveBMP(char* filename, int* result, int w, int h, int max_colors){
        initColours();
	FILE *f;
	unsigned char *img = NULL;
	int filesize = 54 + 3*w*h;  //w is your image width, h is image height, both int

	unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
	unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
	unsigned char bmppad[3] = {0,0,0};

	bmpfileheader[ 2] = (unsigned char)(filesize    );
	bmpfileheader[ 3] = (unsigned char)(filesize>> 8);
	bmpfileheader[ 4] = (unsigned char)(filesize>>16);
	bmpfileheader[ 5] = (unsigned char)(filesize>>24);

	bmpinfoheader[ 4] = (unsigned char)(       w    );
	bmpinfoheader[ 5] = (unsigned char)(       w>> 8);
	bmpinfoheader[ 6] = (unsigned char)(       w>>16);
	bmpinfoheader[ 7] = (unsigned char)(       w>>24);
	bmpinfoheader[ 8] = (unsigned char)(       h    );
	bmpinfoheader[ 9] = (unsigned char)(       h>> 8);
	bmpinfoheader[10] = (unsigned char)(       h>>16);
	bmpinfoheader[11] = (unsigned char)(       h>>24);

	f = fopen(filename,"wb");
	fwrite(bmpfileheader,1,14,f);
	fwrite(bmpinfoheader,1,40,f);

	
	img = (unsigned char *)malloc(3*w);

	int i,j;
	for(j=0; j<h; j++)
	{
	    for(i=0; i<w; i++)
		{
		  int index = result[j*w+i]%255;
		  int r = table[index].r;
		  int g = table[index].g;
		  int b = table[index].b;
		  if (r > 255) r=255;
		    if (g > 255) g=255;
		    if (b > 255) b=255;
		    img[i*3+2] = (unsigned char)(r);
		    img[i*3+1] = (unsigned char)(g);
		    img[i*3+0] = (unsigned char)(b);
		}
		fwrite(img,3,w,f);
	    fwrite(bmppad,1,(4-(w*3)%4)%4,f);
	}
	fclose(f);
}
