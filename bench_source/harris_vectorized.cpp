#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <cmath>
#include <string.h>
#include <immintrin.h>
#include "harris.h"

#ifdef __AVX__
#define VF 4
#define VTYPE __m256d
#define VSETZERO _mm256_setzero_pd
#define VSET1 _mm256_set1_pd
#define VLOAD _mm256_loadu_pd
#define VSTORE _mm256_storeu_pd
#define VLOADU _mm256_loadu_pd
#define VSTOREU _mm256_storeu_pd
#define VMUL _mm256_mul_pd
#define VADD _mm256_add_pd
#ifdef __AVX2__
#define VFMADD _mm256_fmadd_pd
#else
#define VFMADD(a, b, c) VADD(VMUL(a, b), c)
#endif
#else
#define VF 2
#define VTYPE __m128d
#define VSETZERO _mm_setzero_pd
#define VSET1 _mm_set1_pd
#define VLOAD _mm_loadu_pd
#define VSTORE _mm_storeu_pd
#define VLOADU _mm_loadu_pd
#define VSTOREU _mm_storeu_pd
#define VMUL _mm_mul_pd
#define VADD _mm_add_pd
#endif


#define _mm_sobelY(A,i,j) (tab_cell(A, i-1, j-1) *( -0.0833333333333f) + \
                        (tab_cell(A, i-1, j+1) * 0.0833333333333f) + \
                      (tab_cell(A, i, j-1) * -0.166666666667f) + \
                       (tab_cell(A, i, j+1) * 0.166666666667f) + \
                       (tab_cell(A, i+1, j-1) * -0.0833333333333f) +\
                       (tab_cell(A, i+1, j+1) * 0.0833333333333f))

#define _mm_sobelX(A,i,j) 
tab_cell(A, i-1, j-1) *( -0.0833333333333f) +\
                         (tab_cell(A, i+1, j-1) * 0.0833333333333f) + \
                      (tab_cell(A, i+1, j) * 0.166666666667f) + \
                         (tab_cell(A, i-1, j) * (-0.166666666667f)) + \
                       (tab_cell(A, i-1, j+1) * (-0.0833333333333f)) + \
                       (tab_cell(A, i+1, j+1) * 0.0833333333333f))



extern "C" void  pipeline_harris(int  C, int  R, void * img_void_arg, void * harris_void_arg)
{
    float * img;
    img = (float *) (img_void_arg);
    float * harris;
    harris = (float *) (harris_void_arg);

    float constants[] = {-0.0833333333333f,-0.166666666667f,0.0833333333333f,
        -0-0.0833333333333f};
    VTYPE sobel1 = VLOADU(&constants[0]);
    VTYPE sobel2 = VLOADU(&constants[1]);
    VTYPE sobel3 = VLOADU(&constants[2]);
    VTYPE sobel4 = VLOADU(&constants[3]);

    // Tile size
    static int TSIZEX = 32;
    static int TSIZEY = 256;
    // Filter size
    static int ft_size = 1;

    #pragma omp parallel for schedule(static)
    for (int  Ti = 0; (Ti <= (R / TSIZEX)); Ti ++)
    {
        float Ix = (float*) malloc();
        float Iy[TSIZEX+2*ft_size][TSIZEY+2*ft_size];

        for (int  Tj = 0; (Tj <= (C / TSIZEY)); Tj ++)
        {

            int bot, top, right, left;
            int bot0, top0, right0, left0;
            int height,width;
            // Tile bounds
            bot0 = isl_min(isl_max(Ti * TSIZEX, ft_size), R-ft_size);
            top0 = isl_min( (Ti + 1) * TSIZEX , R-ft_size);
            left0 = isl_min(isl_max(Tj * TSIZEY, ft_size), C-ft_size);
            right0 = isl_min( (Tj + 1) * TSIZEY, C-ft_size);

            width = right0 - left0;
            height = top0 - bot0;

            for (int  i = bot0; i <= top0 ; i++)
            {
                #pragma ivdep
                for (int  j = left0 ; j <= right0 ; j ++)
                {
                    tile_cell(Iy,i,j) = sobelY(img,i,j);
                    VMUL
                }
            }
            for (int  i = bot0; i <= top0 ; i++)
            {
                #pragma ivdep
                for (int  j = left0 ; j <= right0 ; j ++)
                {
                  tile_cell(Ix,i,j) = sobelX(img, i, j);
                }
            }

            for (int  i = bot0; i < top0 ; i++)
            {
                #pragma ivdep
                for (int  j = left0 ; j < right0 ; j ++)
                {
                    tab_cell(harris,i,j) = (filter2sq(Ix,Ix,i,j)*filter2sq(Iy,Iy,i,j) - filter2sq(Ix,Iy,i,j)*filter2sq(Ix,Iy,i,j))- (0.04f * ( filter2sq(Ix,Ix,i,j) + filter2sq(Iy,Iy,i,j)) * ( filter2sq(Ix,Ix,i,j) + filter2sq(Iy,Iy,i,j)));
                }

            }

        }
    }
}
