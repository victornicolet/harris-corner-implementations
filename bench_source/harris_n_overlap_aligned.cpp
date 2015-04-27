/* Victor Nicolet */

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <cmath>
#include <string.h>
#include <omp.h>
#include "harris.h"

extern "C" void  pipeline_harris(int  C, int  R, void * img_void_arg, void * harris_void_arg)
{
  float * img;
  img = (float *) (img_void_arg);
  float * harris;
  harris = (float *) (harris_void_arg);


  float ** Ix = aligned_allocmatrix(R,C);
  float ** Iy = aligned_allocmatrix(R,C);
  float ** Sxx = aligned_allocmatrix(R,C);
  float ** Sxy = aligned_allocmatrix(R,C);
  float ** Syy = aligned_allocmatrix(R,C);



  // Tile size
  static int TSIZEX = 32;
  static int TSIZEY = 256;
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
          #pragma omp task depend(out : Iy[bot0:height][ left0:width],\
            Ix[bot0:height][left0:width])
          {
            for (int  i = bot0; i < top0 ; i++)
            {

              #pragma ivdep
              for (int  j = left0 ; j < right0 ; j ++)
              {
                mat_cell(Iy,i,j) = sobelY(img, i, j);
                mat_cell(Ix,i,j) = sobelX(img, i, j);
              }
            }
          }

          // Computational tasks : Ixx & Sxx fused

          #pragma omp task depend(out : Syy[bot:height][left:width]) \
              depend(in : Iy[bot-ft_size:height+ft_size][left-ft_size:width+ft_size])
          {
            for (int  i = bot; i < top ; i++)
            {
              #pragma ivdep
              for (int  j = left ; j < right ; j ++)
              {
                mat_cell(Syy,i,j) = filter2sq( Iy, Iy, i,j) ;
              }
            }
          }

          #pragma omp task depend(in : Ix[bot-ft_size:height+ft_size][left-ft_size:width+ft_size],\
               Iy[bot-ft_size:height+ft_size][left-ft_size:width+ft_size]) \
            depend(out : Sxy[bot:height][left:width])
          {
            for (int  i = bot; i < top ; i++)
            {
              #pragma ivdep
              for (int  j = left ; j < right ; j ++)
              {
                mat_cell(Sxy, i, j) = filter2sq( Ix, Iy, i,j) ;
              }
            }
          }

          #pragma omp task depend(in : Ix[bot-ft_size:height+ft_size][left-ft_size:width+ft_size]) \
            depend( out : Sxx[bot:height][left:width])
          {
            for (int  i = bot; i < top ; i++)
            {
              #pragma ivdep
              for (int  j = left ; j < right ; j ++)
              {
                mat_cell(Sxx, i, j) = filter2sq( Ix,Ix, i,j) ;
              }
            }

          }

          // det + trace
          #pragma omp task depend( in : Sxx[bot:height][left:width], \
             Sxy[bot:height][left:width], \
             Syy[bot:height][left:width])
          for (int  i = bot; i < top ; i++)
          {
            #pragma ivdep
            for (int  j = left ; j < right ; j ++)
            {
              tab_cell(harris,i,j) =  det(i, j) - (0.04f * ((mat_cell(Sxx, i, j) + mat_cell(Syy, i, j)) * (mat_cell(Sxx, i, j) + mat_cell(Syy, i, j))));
            }
          }

        }
      }
    }
  }

  freematrix(Ix, R);
  freematrix(Iy, R);
  freematrix(Sxx, R);
  freematrix(Sxy, R);
  freematrix(Syy, R);
}
