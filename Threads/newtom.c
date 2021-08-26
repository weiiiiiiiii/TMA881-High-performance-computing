#include <stdlib.h>
#include <stdio.h>
#include <threads.h>
#include <string.h>
#include <complex.h>
#include <math.h>

int nthrds, img_size, degree;

const double eps = 0.001;
const long upper_bound = 10000000000;
const int max_conv = 100;

double complex **roots;
int colormap[10][3] = {
    {180, 0, 30}, {0, 180, 30}, {0, 30, 180}, {0, 190, 180}, {180, 0, 175},
    {180, 255, 0}, {155, 170, 180}, {70, 50, 0}, {150, 60, 0}, {0, 150, 60}
};

typedef struct {
    int thrd_idx;
    char **attr;
    char **conv;
    char *row_done;
} comp_thrd_info_t;

typedef struct {
    FILE *attrfile;
    FILE *convfile;
    char **attr;
    char **conv;
    char *row_done;
} write_thrd_info_t;


void compute(double complex c, char *attr, char *conv) {
    int iter;
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
            case 1:
                c -= c-1;
                break;
            case 2:
                c -= (c * c - 1) / (2 * c);
                break;
            case 3:
            case 4:
            case 5:
            case 6:
            case 7:
            case 8:
            case 9:
                // computation
                fx = cpow(c, degree) - 1;
                dfx = degree * cpow(c, degree-1);
                diff = fx / dfx;
                c = c - diff;
                break;
            default:
                fprintf(stderr, "unexpected degree\n");
                exit(1);
        }
        
        
    }
    *conv = iter <= max_conv ? iter : max_conv;
    
}

int comp_thrd(void *args) {
    const comp_thrd_info_t *thrd_info  = (comp_thrd_info_t*) args;
    int thrd_idx = thrd_info->thrd_idx;
    char **attr = thrd_info->attr;
    char **conv = thrd_info->conv;
    
    // process rows
    for(int ix = thrd_idx; ix < img_size; ix += nthrds) {
        // process each point
        double re = ix * (2.0 - (-2.0)) / (img_size - 1.0) - 2.0;
        for (int jx = 0; jx < img_size; jx++) {
            double im = jx * (2.0 - (-2.0)) / (img_size - 1.0) - 2.0;
            double complex c = re + im * I;
            compute(c, &attr[ix][jx], &conv[ix][jx]);
        }
        thrd_info->row_done[ix] = 1;
    }
}

int writefile(void *args) {
    const write_thrd_info_t *thrd_info  = (write_thrd_info_t*) args;
    FILE *attrfile = thrd_info->attrfile;
    FILE *convfile = thrd_info->convfile;
    
    struct timespec sleep_timespec;
    sleep_timespec.tv_sec = 0;
    sleep_timespec.tv_nsec = 100000;
    
    
    char attr_color_row[12 * img_size];
    char con_color_row[12 * img_size];
    for (int ix = 0; ix < img_size; ix++) {
        if (!thrd_info->row_done[ix]) {
            nanosleep(&sleep_timespec, NULL);
            continue;
        }
        char *attr_row = thrd_info->attr[ix];
        char *conv_row = thrd_info->conv[ix];
        for (int jx = 0; jx < img_size; jx++) {
            int *attr_color = colormap[attr_row[jx]];
            sprintf(attr_color_row+12*jx, "%03d ", attr_color[0]);
            sprintf(attr_color_row+12*jx+4, "%03d ", attr_color[1]);
            sprintf(attr_color_row+12*jx+8, "%03d ", attr_color[2]);
            
            int conv_color = 255 * (conv_row[jx] + 1) / 100;
            sprintf(con_color_row+12*jx, "%03d ", conv_color);
            sprintf(con_color_row+12*jx+4, "%03d ", conv_color);
            sprintf(con_color_row+12*jx+8, "%03d ", conv_color);
        }
        attr_color_row[12 * img_size - 1] = '\n';
        con_color_row[12 * img_size - 1] = '\n';
        
        fwrite(attr_color_row, sizeof(char), 12 * img_size, attrfile);
        fwrite(con_color_row, sizeof(char), 12 * img_size, convfile);
    }
}

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
    printf(attrfile, "P3\n");
    printf(attrfile, "%d %d \n", img_size, img_size);
    printf(attrfile,"255\n");
    
    sprintf(filename, "newton_convergence_x%d.ppm", degree);
    convfile = fopen(filename, "w");
    printf(convfile, "P3\n");
    printf(convfile, "%d %d \n", img_size, img_size);
    printf(convfile,"255\n");
    
    // initialized roots for degreeonent degreeression
    init_roots();
    
    char** attr = (char**) malloc(sizeof(char*) * img_size);
    for (int ix = 0; ix < img_size; ix++)
        attr[ix] = (char *) malloc(sizeof(char) *img_size);
    
    char** conv = (char**) malloc(sizeof(char*) * img_size);
    for (int ix = 0; ix < img_size; ix++)
        conv[ix] = (char *) malloc(sizeof(char) *img_size);
    
    char* row_done = (char*) calloc(img_size, sizeof(char));
    
    // Synchronization of compute and write threads.
    thrd_t comp_thrds[nthrds];
    thrd_t write_thrd;
    
    comp_thrd_info_t comp_thrds_info[nthrds];
    int r;
    for (int tx = 0; tx < nthrds; ++tx) {
        comp_thrds_info[tx].thrd_idx = tx;
        comp_thrds_info[tx].attr = attr;
        comp_thrds_info[tx].conv = conv;
        comp_thrds_info[tx].row_done = row_done;
        r = thrd_create(comp_thrds+tx, comp_thrd, (void *)(&comp_thrds_info[tx]));
        if (r != thrd_success) {
            fprintf(stderr, "failed to create thread\n");
            exit(1);
        }
    }
    
    write_thrd_info_t write_thrd_info;
    write_thrd_info.attrfile = attrfile;
    write_thrd_info.convfile = convfile;
    write_thrd_info.attr = attr;
    write_thrd_info.conv = conv;
    write_thrd_info.row_done = row_done;
    r = thrd_create(&write_thrd, writefile, (void *)(&write_thrd_info));
    if (r) {
        fprintf(stderr, "failed to create thread\n");
        exit(1);
    }
    
    
    for (int tx = 0; tx < nthrds; tx++) {
        thrd_join(comp_thrds[tx], NULL);
    }
    thrd_join(write_thrd, NULL);
    
    fclose(attrfile);
    fclose(convfile);

    // release memory
    for (int ix = 0; ix < img_size; ix++) {
        free(attr[ix]);
        free(conv[ix]);
    }
    free(attr);
    free(conv);
    free(row_done);
    for (int ix = 0; ix < 9; ix++) 
        free(roots[ix]);
    free(roots);
}

