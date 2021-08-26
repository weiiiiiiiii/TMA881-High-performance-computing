#include <stdlib.h>
#include <stdio.h>
#include <threads.h>
#include <string.h>
#include <complex.h>
#include <math.h>

// number of threads, picture size and exponent degree
int nthrds, img_size, degree;

// epsilon for breaking the computation loop
const double eps = 0.001;
const long upper_bound = 10000000000;
const int max_conv = 50;

// roots for the exponent expression
double complex **roots;

// color map for drawing attractor image
char *colormap[10] = {
		      "180 000 030", "000 180 030", "000 030 180", "000 190 180", "180 000 175",
		      "180 255 000", "155 170 180", "070 050 000", "150 060 000", "000 150 060"
};
char *colormap_conv[50] = {
			   "005 005 005 ", "010 010 010 ", "015 015 015 ", "020 020 020 ", "025 025 025 ",
			   "030 030 030 ", "035 035 035 ", "040 040 040 ", "045 045 045 ", "050 050 050 ",
			   "056 056 056 ", "061 061 061 ", "066 066 066 ", "071 071 071 ", "076 076 076 ",
			   "081 081 081 ", "086 086 086 ", "091 091 091 ", "096 096 096 ", "100 100 100 ",
			   "107 107 107 ", "112 112 112 ", "117 117 117 ", "122 122 122 ", "127 127 127 ",
			   "132 132 132 ", "137 137 137 ", "142 142 142 ", "147 147 147 ", "153 153 153 ",
			   "158 158 158 ", "163 163 163 ", "168 168 168 ", "173 173 173 ", "178 178 178 ",
			   "183 183 183 ", "188 188 188 ", "193 193 193 ", "198 198 198 ", "204 204 204 ",
			   "209 209 209 ", "214 214 214 ", "219 219 219 ", "224 224 224 ", "229 229 229 ",
			   "234 234 234 ", "239 239 239 ", "244 244 244 ", "249 249 249 ", "255 255 255 "
};

// argument type for computation thread
typedef struct {
  int thrd_idx;
  char *attr;
  char *conv;
  char *row_done;
} comp_thrd_info_t;

// argument type for writing thread
typedef struct {
  FILE *attrfile;
  FILE *convfile;
  char *attr;
  char *conv;
  char *row_done;
} write_thrd_info_t;

// given a complex number, calculate which root it will be converged to
// and set the attr and conv with the specific value
void compute(double complex c, char *attr, char *conv) {
  int iter;
  // start iteration
  for (iter = 0, *attr = -1;;iter++) {
    // break when x is closer than 10^-3 to the origin
    if (creal(c)*creal(c) + cimag(c)*cimag(c) <= eps * eps) {
      *attr = 9;
      break;
    }
    // break when the absolute value of the real and imaginary part of
    // x reaches an upper bound
    if (fabs(creal(c)) > upper_bound || fabs(cimag(c)) > upper_bound) {
      *attr = 9;
      break;
    }
    // If x_i is closer than 10^-3 to one of the roots of f(x),
    // then abort the iteration.
    double complex *rs = roots[degree-1];
    double complex diff;
    for (int ix = 0; ix < degree; ix++) {
      diff = c - rs[ix];
      if (creal(diff)*creal(diff) + cimag(diff)*cimag(diff) <= eps * eps) {
	*attr = ix;
	break;
      }
    }
    if (*attr != -1)
      break;

    double complex fx, dfx;
    switch ( degree ) {
    case 1: // hardcoded computation
      c -= c-1;
      break;
    case 2: // hardcoded computation
      c -= (c * c - 1) / (2 * c);
      break;
    case 3:
      c -= (c * c * c - 1) /(3 * c * c);
      break;
    case 4:
      c -=(c * c * c * c - 1) / (4 * c * c * c);
      break;
    case 5:
      c -=(c * c * c * c * c - 1) / (5 * c * c * c * c);
      break;
    case 6:
      c -=(c * c * c * c * c * c - 1) / ( 6 * c * c * c * c * c);
      break;
    case 7:
      c -= (c * c * c * c * c * c * c -1) / (7 * c * c * c * c * c * c);
      break;
    case 8:
      c -=(c * c * c * c * c * c * c * c -1) / (8 * c * c * c * c * c * c * c);
      break;
    case 9:
      c -=(c * c * c * c * c * c * c * c * c - 1) / (9 * c * c * c * c * c *c * c *c);
      break; 
    default:
      fprintf(stderr, "unexpected degree\n");
      exit(1);
    }

  }
  // set convergence value to the number of iterations
  *conv = iter < max_conv ? iter : max_conv-1;

}

// compute thread
int comp_thrd(void *args) {
  // parse arguments
  const comp_thrd_info_t *thrd_info  = (comp_thrd_info_t*) args;
  int thrd_idx = thrd_info->thrd_idx;
  char *attr = thrd_info->attr;
  char *conv = thrd_info->conv;

  // process rows
  for(int ix = thrd_idx; ix < img_size; ix += nthrds) {
    // process each point
    // the value for the real part of the complex number, reduced to [-2, 2]
    double re = ix * (2.0 - (-2.0)) / (img_size - 1.0) - 2.0;
    for (int jx = 0; jx < img_size; jx++) {
      // the value for the imaginary part of the complex number, reduced to [-2, 2]
      double im = jx * (2.0 - (-2.0)) / (img_size - 1.0) - 2.0;
      // construct the complex number
      double complex c = re + im * I;
      compute(c, attr + ix * img_size + jx, conv + ix * img_size + jx);
    }
    thrd_info->row_done[ix] = 1;
  }
}

// writing thread
int writefile(void *args) {
  // parse arguments
  const write_thrd_info_t *thrd_info  = (write_thrd_info_t*) args;
  FILE *attrfile = thrd_info->attrfile;
  FILE *convfile = thrd_info->convfile;

  struct timespec sleep_timespec;
  sleep_timespec.tv_sec = 0;
  sleep_timespec.tv_nsec = 100;

  // representing the row string of attractor image and convergence image
  char attr_color_row[12 * img_size];
  char con_color_row[12 * img_size];
  // write row by row
  for (int ix = 0; ix < img_size; ix++) {
    // only write a row when a compute thread finish this row, otherwise wait for timespec
    while (!thrd_info->row_done[ix]) {
      nanosleep(&sleep_timespec, NULL);
      continue;
    }
    // get the row from the attractor array and convergence array
    char *attr_row = thrd_info->attr + ix * img_size;
    char *conv_row = thrd_info->conv + ix * img_size;
    // write pixel by pixel
    for (int jx = 0; jx < img_size; jx++) {
      // choose the color according to the attr and conv value
      char *attr_color = colormap[attr_row[jx]];
      memcpy(attr_color_row+12*jx, attr_color, 12);

      char *conv_color = colormap_conv[conv_row[jx]];
      memcpy(con_color_row+12*jx, conv_color, 12);
    }
    attr_color_row[12 * img_size - 1] = '\n';
    con_color_row[12 * img_size - 1] = '\n';

    // write to the file
    fwrite(attr_color_row, sizeof(char), 12 * img_size, attrfile);
    fwrite(con_color_row, sizeof(char), 12 * img_size, convfile);
  }
}

// roots for the exponent expression x^degree - 1
void init_roots() {
  roots = (double complex**) malloc(sizeof(double complex*) * 9);
  for (int ix = 0; ix < 9; ix++) {
    roots[ix] = (double complex*) malloc(sizeof(double complex) * (ix + 1));
  }
  // roots for x - 1
  roots[0][0] = 1 + 0 * I;
  // roots for x^2 - 1
  roots[1][0] = 1 + 0 * I;
  roots[1][1] = -1 + 0 * I;
  // roots for x^3 - 1
  roots[2][0] = 1 + 0 * I;
  roots[2][1] = -0.5 + 0.86603 * I;
  roots[2][2] = -0.5 - 0.86606 * I;
  // roots for x^4 - 1
  roots[3][0] = 1 + 0 * I;
  roots[3][1] = 0 + 1 * I;
  roots[3][2] = -1 + 0 * I;
  roots[3][3] = 0 - 1 * I;
  // roots for x^5 - 1
  roots[4][0] = 1 + 0 * I;
  roots[4][1] = 0.309017 + 0.951057 * I;
  roots[4][2] = -0.809017 + 0.587785 * I;
  roots[4][3] = -0.809017 - 0.587785 * I;
  roots[4][4] = 0.309017 - 0.951057 * I;
  // roots for x^6 - 1
  roots[5][0] = 1 + 0 * I;
  roots[5][1] = 0.5 + 0.866025 * I;
  roots[5][2] = -0.5 + 0.866025 * I;
  roots[5][3] = -1 - 0 * I;
  roots[5][4] = -0.5 - 0.866025 * I;
  roots[5][5] = 0.5 - 0.866025 * I;
  // roots for x^7 - 1
  roots[6][0] = 1 + 0 * I;
  roots[6][1] = 0.62349 + 0.781831 * I;
  roots[6][2] = -0.222521 + 0.974928 * I;
  roots[6][3] = -0.900969 + 0.433884 * I;
  roots[6][4] = -0.900969 - 0.433884 * I;
  roots[6][5] = -0.222521 - 0.974928 * I;
  roots[6][6] = 0.62349 - 0.781831 * I;
  // roots for x^8 - 1
  roots[7][0] = 1 + 0 * I;
  roots[7][1] = 0.707107 + 0.707107 * I;
  roots[7][2] = 0 + 1 * I;
  roots[7][3] = -0.707107 + 0.707107 * I;
  roots[7][4] = -1 + 0 * I;
  roots[7][5] = -0.707107 - 0.707107 * I;
  roots[7][6] = 0 - 1 * I;
  roots[7][7] = 0.707107 - 0.707107 * I;
  // roots for x^9 - 1
  roots[8][0] = 1 + 0 * I;
  roots[8][1] = 0.766044 + 0.642788 * I;
  roots[8][2] = 0.173648 + 0.984808 * I;
  roots[8][3] = -0.5 + 0.866025 * I;
  roots[8][4] = -0.939693 + 0.34202 * I;
  roots[8][5] = -0.939693 - 0.34202 * I;
  roots[8][6] = -0.5 - 0.866025 * I;
  roots[8][7] = 0.173648 - 0.984808 * I;
  roots[8][8] = 0.766044 - 0.642788 * I;
}

int main(int argc, char *argv[]) {
  if (argc < 4) {
    printf("Usage: newton -t[NumberOfThreads] -l[ImageSize] degreeonent\n");
    exit(1);
  }

  // Parsing command line arguments
  for (int ix = 1; ix < argc - 1; ix++) {
    if (strncmp(argv[ix], "-t", 2) == 0)
      nthrds = atoi(argv[ix]+2);
    if (strncmp(argv[ix], "-l", 2) == 0)
      img_size = atoi(argv[ix]+2);
  }
  degree = atoi(argv[argc-1]);

  // create attractor file and convergence file
  // and write the required file header
  FILE *attrfile, *convfile;
  char filename[26];

  sprintf(filename, "newton_attractors_x%d.ppm", degree);
  attrfile = fopen(filename, "w");
  fprintf(attrfile, "P3\n");
  fprintf(attrfile, "%d %d \n", img_size, img_size);
  fprintf(attrfile,"255\n");

  sprintf(filename, "newton_convergence_x%d.ppm", degree);
  convfile = fopen(filename, "w");
  fprintf(convfile, "P3\n");
  fprintf(convfile, "%d %d \n", img_size, img_size);
  fprintf(convfile,"255\n");

  // initialized roots for degreeonent degreeression
  init_roots();

  // allocate attr and conv array, the size equals img_size * img_size
  char* attr = (char*) malloc(sizeof(char) * img_size * img_size);
  char* conv = (char*) malloc(sizeof(char) * img_size * img_size);
  char* row_done = (char*) calloc(img_size, sizeof(char));
  // Synchronization of compute and write threads.
  thrd_t comp_thrds[nthrds];
  thrd_t write_thrd;

  // computation thread arguments array
  comp_thrd_info_t comp_thrds_info[nthrds];
  int r;
  for (int tx = 0; tx < nthrds; ++tx) {
    comp_thrds_info[tx].thrd_idx = tx;
    comp_thrds_info[tx].attr = attr;
    comp_thrds_info[tx].conv = conv;
    comp_thrds_info[tx].row_done = row_done;
    // start computation thread
    r = thrd_create(comp_thrds+tx, comp_thrd, (void *)(&comp_thrds_info[tx]));
    if (r != thrd_success) {
      fprintf(stderr, "failed to create thread\n");
      exit(1);
    }
  }

  // writing thread argument
  write_thrd_info_t write_thrd_info;
  write_thrd_info.attrfile = attrfile;
  write_thrd_info.convfile = convfile;
  write_thrd_info.attr = attr;
  write_thrd_info.conv = conv;
  write_thrd_info.row_done = row_done;
  // start writing thread
  r = thrd_create(&write_thrd, writefile, (void *)(&write_thrd_info));
  if (r) {
    fprintf(stderr, "failed to create thread\n");
    exit(1);
  }

  for (int tx = 0; tx < nthrds; tx++) {
    thrd_join(comp_thrds[tx], NULL);
  }
  thrd_join(write_thrd, NULL);

  // close file
  fclose(attrfile);
  fclose(convfile);

  // release allocated memory
  free(attr);
  free(conv);
  free(row_done);
  for (int ix = 0; ix < 9; ix++)
    free(roots[ix]);
  free(roots);
}



