#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <cmath>
#include <string.h>
#include "harris_tiles.h"

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

    #pragma omp parallel for schedule(static)
    for (int  Ti = 0; (Ti <= (R / TSIZEX)); Ti ++)
    {
        float Ix[TSIZEX+2*ft_size][TSIZEY+2*ft_size];
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
