/* Victor Nicolet */

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <cmath>
#include <string.h>
#include <omp.h>
#define isl_min(x,y) ((x) < (y) ? (x) : (y))
#define isl_max(x,y) ((x) > (y) ? (x) : (y))
#define isl_floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))

#define mat_cell(A,i,j) A[i][j]
#define tab_cell(A,i,j) A[((i)*C + (j))]
#define sobelY(A,i,j) (tab_cell(A, i-1, j-1) *( -0.0833333333333f) + (tab_cell(A, i-1, j+1) * 0.0833333333333f) + \
                      (tab_cell(A, i, j-1) * -0.166666666667f) + (tab_cell(A, i, j+1) * 0.166666666667f) + \
                       (tab_cell(A, i+1, j-1) * -0.0833333333333f) +(tab_cell(A, i+1, j+1) * 0.0833333333333f))

#define sobelX(A,i,j) (tab_cell(A, i-1, j-1) *( -0.0833333333333f) + (tab_cell(A, i+1, j-1) * 0.0833333333333f) + \
                      (tab_cell(A, i+1, j) * 0.166666666667f) + (tab_cell(A, i-1, j) * (-0.166666666667f)) + \
                       (tab_cell(A, i-1, j+1) * (-0.0833333333333f)) + (tab_cell(A, i+1, j+1) * 0.0833333333333f))

// Macros implementing various smoothing filters. Use regexes to edit !

#define filter3(A,i,j) mat_cell(A, i-2, j-2) + mat_cell(A, i-2, j-1) + mat_cell(A, i-2, j) + \
                        mat_cell(A, i-2, j+1) + mat_cell(A, i-2, j+2) + mat_cell(A, i-1, j-2) + \
                        mat_cell(A, i-1, j-1) + mat_cell(A, i-1, j) + mat_cell(A, i-1, j+1) + \
                        mat_cell(A, i-1, j+2) + mat_cell(A, i, j-2)   + mat_cell(A, i, j-1)   + \
                        mat_cell(A, i, j)   + mat_cell(A, i, j+1)   + mat_cell(A, i, j+2) + \
                        mat_cell(A, i+1, j-2) + mat_cell(A, i+1, j-1) + mat_cell(A, i+1, j) + \
                        mat_cell(A, i+1, j+1) + mat_cell(A, i+1, j+2)+ mat_cell(A, i+2, j-2) + \
                        mat_cell(A, i+1, j-1) + mat_cell(A, i+2, j) + mat_cell(A, i+2, j+1) + \
                         mat_cell(A, i+2, j+2)

//([mat_cell\(]+)[A,]+ ([i]+[-|+]*[1-2]*[, j]+[-|+]*[1-2]*[\)]+) -> $1 A, $2 * $1 B, $2
#define filter3sq(A,B,i,j) mat_cell( A, i-2, j-2) * mat_cell( B, i-2, j-2) + \
                      mat_cell( A, i-2, j-1) * mat_cell( B, i-2, j-1) +\
                      mat_cell( A, i-2, j) * mat_cell( B, i-2, j) + \
                      mat_cell( A, i-2, j+1) * mat_cell( B, i-2, j+1) + \
                      mat_cell( A, i-2, j+2) * mat_cell( B, i-2, j+2) + \
                      mat_cell( A, i-1, j-2) * mat_cell( B, i-1, j-2) + \
                      mat_cell( A, i-1, j-1) * mat_cell( B, i-1, j-1) + \
                      mat_cell( A, i-1, j) * mat_cell( B, i-1, j) +  \
                      mat_cell( A, i-1, j+1) * mat_cell( B, i-1, j+1) + \
                      mat_cell( A, i-1, j+2) * mat_cell( B, i-1, j+2) + \
                      mat_cell( A, i, j-2) * mat_cell( B, i, j-2)   + \
                      mat_cell( A, i, j-1) * mat_cell( B, i, j-1)   + \
                      mat_cell( A, i, j) * mat_cell( B, i, j)   + \
                      mat_cell( A, i, j+1) * mat_cell( B, i, j+1)   + \
                      mat_cell( A, i, j+2) * mat_cell( B, i, j+2) + \
                      mat_cell( A, i+1, j-2) * mat_cell( B, i+1, j-2) + \
                      mat_cell( A, i+1, j-1) * mat_cell( B, i+1, j-1) + \
                      mat_cell( A, i+1, j) * mat_cell( B, i+1, j) + \
                      mat_cell( A, i+1, j+1) * mat_cell( B, i+1, j+1) + \
                      mat_cell( A, i+1, j+2) * mat_cell( B, i+1, j+2)+ \
                      mat_cell( A, i+2, j-2) * mat_cell( B, i+2, j-2) + \
                      mat_cell( A, i+1, j-1) * mat_cell( B, i+1, j-1) + \
                      mat_cell( A, i+2, j) * mat_cell( B, i+2, j) + \
                      mat_cell( A, i+2, j+1) * mat_cell( B, i+2, j+1) + \
                      mat_cell( A, i+2, j+2) * mat_cell( B, i+2, j+2)

#define filter2(A,i,j)  mat_cell(A, i-1, j-1) + mat_cell(A, i-1, j) + mat_cell(A, i-1, j+1) + \
                        mat_cell(A, i, j-1)   + mat_cell(A, i, j)   + mat_cell(A, i, j+1)  + \
                        mat_cell(A, i+1, j-1) + mat_cell(A, i+1, j) + mat_cell(A, i+1, j+1)

#define filter2sq(A,B,i,j) mat_cell(A, i-1, j-1) * mat_cell(B, i-1, j-1) +\
                        mat_cell(A, i-1, j) * mat_cell(B, i-1, j) +\
                        mat_cell(A, i-1, j+1) * mat_cell(B, i-1, j+1) + \
                        mat_cell(A, i, j-1) * mat_cell(B, i, j-1) + \
                        mat_cell(A, i, j) * mat_cell(B, i, j)  + \
                        mat_cell(A, i, j+1) * mat_cell(B, i, j+1)  + \
                        mat_cell(A, i+1, j-1) * mat_cell(B, i+1, j-1) + \
                        mat_cell(A, i+1, j) * mat_cell(B, i+1, j) +\
                        mat_cell(A, i+1, j+1) * mat_cell(B, i+1, j+1)

#define det(i,j)        mat_cell(Sxx, i, j) * mat_cell(Syy, i, j) - \
                        mat_cell(Sxy, i, j) * mat_cell(Sxy, i, j)



float ** allocmatrix(int rows, int cols){
  float ** Res = (float **) malloc(sizeof(float *) * rows);
  for(int i = 0; i < rows ; i++){
    Res[i] = (float *) malloc(sizeof(float) * cols);

    if( Res[i] == NULL){
      printf("Error while allocating two dimensional matrix at line %i \n",i);
      return NULL;
    }
  }
  if( Res == NULL ){
    printf("Error while allocating two dimensionnal matrix \n");
    return NULL;
  }
  return Res;
}

int freematrix(float ** Mat, int rows){
  for(int i = 0; i < rows; i++){
    free(Mat[i]);
  }
  free(Mat);
  return 0;
}

extern "C" void  pipeline_harris(int  C, int  R, void * img_void_arg, void * harris_void_arg)
{
  float * img;
  img = (float *) (img_void_arg);
  float * harris;
  harris = (float *) (harris_void_arg);


  float ** Ix = allocmatrix(R,C);
  float ** Iy = allocmatrix(R,C);
  float ** Sxx = allocmatrix(R,C);
  float ** Sxy = allocmatrix(R,C);
  float ** Syy = allocmatrix(R,C);



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
