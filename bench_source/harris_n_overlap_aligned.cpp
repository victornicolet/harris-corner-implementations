/* Victor Nicolet */

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <cmath>
#include <string.h>
#include <omp.h>
#include "harris.h"

extern "C" void  pipeline_harris_aligned(int  C, int  R, float * img, float * harris)
{
  // Tile size
  static int TSIZEX = 32;
  static int TSIZEY = 256;

  float * Ix = falloc_aligned_padded(R,C,CACHE_LINE_SIZE);
  float * Iy = falloc_aligned_padded(R,C,CACHE_LINE_SIZE);
  float * Sxx = falloc_aligned_padded(R,C,CACHE_LINE_SIZE);
  float * Sxy = falloc_aligned_padded(R,C,CACHE_LINE_SIZE);
  float * Syy = falloc_aligned_padded(R,C,CACHE_LINE_SIZE);

  // Filter size
  // filter2 -> ft_size =1 or filter3 -> ft_size = 2
  static int ft_size = 1;

  #pragma omp parallel
  {
    #pragma omp single
    {

      for (int  Ti = 0; (Ti <= (R / TSIZEX)); Ti ++)
      {

        for (int  Tj = 0; (Tj <= (C / TSIZEY)); Tj ++)
        {
          int bot, top, right, left;
          int bot0, top0, right0, left0;
          int height,width;

          // Tile bounds :
          // Level 0 : margins before smoothing filter
          bot0 = isl_max(Ti * TSIZEX, 1);
          top0 = isl_min( (Ti + 1) * TSIZEX , R-1);
          left0 = isl_max(Tj * TSIZEY, 1);
          right0 = isl_min( (Tj + 1) * TSIZEY, C-1);
          // Margins after smoothing filter
          bot = isl_max(Ti * TSIZEX, ft_size);
          top = isl_min( (Ti + 1) * TSIZEX+1 , R-ft_size-1);
          left = isl_max(Tj * TSIZEY, ft_size+1);
          right = isl_min( (Tj + 1) * TSIZEY, C-ft_size-1);

          // Tile size after smoothing filter
          width = right - left;
          height = top - bot;

          // Loading task + sobel filter
          #pragma omp task depend(out : Iy,\
            Ix)
          {
            for (int  i = bot0; i < top0 ; i++)
            {

              #pragma ivdep
              for (int  j = left0 ; j < right0 ; j ++)
              {
                tab_cell(Iy,i,j) = sobelY(img, i, j);
                tab_cell(Ix,i,j) = sobelX(img, i, j);
              }
            }
          }

          // Computational tasks : Ixx & Sxx fused

          #pragma omp task depend(out : Syy) \
              depend(in : Iy)
          {
            for (int  i = bot; i < top ; i++)
            {
              #pragma ivdep
              for (int  j = left ; j < right ; j ++)
              {
                tab_cell(Syy,i,j) = t_filter2sq( Iy, Iy, i,j) ;
              }
            }
          }

          #pragma omp task depend(in : Ix,\
               Iy) \
            depend(out : Sxy)
          {
            for (int  i = bot; i < top ; i++)
            {
              #pragma ivdep
              for (int  j = left ; j < right ; j ++)
              {
                tab_cell(Sxy, i, j) = t_filter2sq( Ix, Iy, i,j) ;
              }
            }
          }

          #pragma omp task depend(in : Ix) \
            depend( out : Sxx)
          {
            for (int  i = bot; i < top ; i++)
            {
              #pragma ivdep
              for (int  j = left ; j < right ; j ++)
              {
                tab_cell(Sxx, i, j) = t_filter2sq( Ix,Ix, i,j) ;
              }
            }

          }

          // det + trace
          #pragma omp task depend( in : Sxx, \
             Sxy, \
             Syy)
          for (int  i = bot; i < top ; i++)
          {
            #pragma ivdep
            for (int  j = left ; j < right ; j ++)
            {
                tab_cell(harris,i,j) =  t_det(i, j) - (0.04f * ((tab_cell(Sxx, i, j) + tab_cell(Syy, i, j)) * (tab_cell(Sxx, i, j) + tab_cell(Syy, i, j))));
            }
          }

        }
      }
    }
  }

  free(Ix);
  free(Iy);
  free(Sxx);
  free(Sxy);
  free(Syy);
}
