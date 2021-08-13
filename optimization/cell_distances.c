#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "omp.h"
#include <string.h>

char* filename = "cells";
int num_threads;

int main(int argc, char const *argv[])
{
  if(argc != 2)
   {
    printf("\n!!! Argument count not equal to 2. \n");

    exit(1);
   }
  else
    {
     num_threads = atoi(argv[1]+2);
    }

  // read the files
  FILE* fp = fopen(filename, "r");
  int count = 0;
  char c = 0;
  
  // get the characters and store in char c, increase the count if the character is new line
  while(!feof(fp))
  {
    c = fgetc(fp);
    if(c == '\n')
     {
       ++count;
     }
  }
  fclose(fp);



  int lines = count + 1;
  float coord[lines][3];

//get the x,y,z from file
  FILE* fpr = fopen(filename, "r");

  float xt;
  float yt;
  float zt;

  for(size_t ix = 0; ix < lines ; ++ix)
    {
     fscanf(fp, "%f %f %f", &xt, &yt, &zt);
     coord[ix][0] = (xt);
     coord[ix][1] = (yt);
     coord[ix][2] = (zt);
    }
   fclose(fpr);

  
  int MAX_DIST = 3465;//as the coordinates are between -10 and 10, the maximum distance between two points will be 20*sqrt(3) = 34.64, and there will be 3465 possible distancs to be counted
  unsigned long dis_count[3465] = {0};

  size_t cell_1;
  size_t cell_2;

  omp_set_num_threads(num_threads);
#pragma omp parallel shared(coord, lines)
  {
    unsigned int local_distCounter[3465] = {0};
#pragma omp for private(cell_1, cell_2) schedule (dynamic, 30)
    for(cell_1 = 0; cell_1 < (lines - 1); ++cell_1)
      {
       float x1 = coord[cell_1][0];
       float y1 = coord[cell_1][1];
       float z1 = coord[cell_1][2];

      // compute the distances
      for(cell_2 = (cell_1 + 1); cell_2 < lines; ++cell_2)
	{
 	 float dis_x = x1 - coord[cell_2][0];
	 float dis_y = y1 - coord[cell_2][1];
	 float dis_z = z1 - coord[cell_2][2];
	 ++local_distCounter[(unsigned int) (sqrtf((dis_x*dis_x + dis_y*dis_y + dis_z*dis_z)) * 100)];
        }
     }

#pragma omp critical
    {
      for(size_t ix = 0; ix < MAX_DIST ; ++ix)
	{
	 dis_count[ix] += local_distCounter[ix];
       }
    }
  }

  for(size_t ix = 0; ix < MAX_DIST; ++ix)
    {
     if(dis_count[ix] >= 1)
       {
	 printf("%05.2f %d\n", ((float) ix ) / 100, dis_count[ix]);
       }
    }
  return 0;
}
