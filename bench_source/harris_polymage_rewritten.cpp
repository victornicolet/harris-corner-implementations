#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <cmath>
#include <string.h>
#include "harris.h"

/* This is only the polymage output alogirthm rewritten in a more
* human-readable fashion. It also seems there was segmentation faults in
* the output.
* Here we have only one omp parallel for construct on the outermost loop
* enclosing the tiles. The different steps of the pipeline are strongly grouped
* and executed sequentially in each tile. This version also benefits from
* compiler optimizations since it relies on few but costly calculatory steps.
*/

extern "C" void  pipeline_harris(int  C, int  R, void * img_void_arg, void * harris_void_arg)
{
    float * img;
    img = (float *) (img_void_arg);
    float * harris;
    harris = (float *) (harris_void_arg);



    // Tile size
    static int TSIZEX = 32;
    static int TSIZEY = 256;
    // Filter size
    static int ft_size = 1;

    int i, j;

    #pragma omp parallel for schedule(static) private(i, j) \
        shared(img, harris)
    for (int  Ti = 0; (Ti <= (R / TSIZEX)); Ti ++)
    {
        float Ix[TSIZEX+2*ft_size][TSIZEY+2*ft_size];
        float Iy[TSIZEX+2*ft_size][TSIZEY+2*ft_size];

        #pragma omp private(Ix[TSIZEX+2*ft_size], Iy[TSIZEX+2*ft_size]))
        for (int  Tj = 0; (Tj <= (C / TSIZEY)); Tj ++)
        {
            int bot0, top0, right0, left0;
            int height,width;
            // Tile bounds
            bot0 = isl_min(isl_max(Ti * TSIZEX, ft_size), R-ft_size);
            top0 = isl_min( (Ti + 1) * TSIZEX , R-ft_size);
            left0 = isl_min(isl_max(Tj * TSIZEY, ft_size), C-ft_size);
            right0 = isl_min( (Tj + 1) * TSIZEY, C-ft_size);

            width = right0 - left0;
            height = top0 - bot0;

            for (i = bot0; i <= top0 ; i++)
            {
                #pragma ivdep
                for (int  j = left0 ; j <= right0 ; j ++)
                {
                    tile_cell(Iy,i,j) = sobelY(img,i,j);
                }
            }
            for (i = bot0; i <= top0 ; i++)
            {
                #pragma ivdep
                for (int  j = left0 ; j <= right0 ; j ++)
                {
                  tile_cell(Ix,i,j) = sobelX(img, i, j);
                }
            }
            for (i = bot0; i < top0 ; i++)
            {
                #pragma ivdep
                for (j = left0 ; j < right0 ; j ++)
                {
                    tab_cell(harris,i,j) =
                    (filter2sq_t(Ix,Ix,i,j)*filter2sq_t(Iy,Iy,i,j) -
                        filter2sq_t(Ix,Iy,i,j)*filter2sq_t(Ix,Iy,i,j))-
                    (0.04f * ( filter2sq_t(Ix,Ix,i,j) + filter2sq_t(Iy,Iy,i,j)) *
                        ( filter2sq_t(Ix,Ix,i,j) + filter2sq_t(Iy,Iy,i,j)));
                }

            }

        }
    }
}
