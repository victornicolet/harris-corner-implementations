#ifndef HARRIS_H
#define HARRIS_H

#include <errno.h>
#include <unistd.h>

#define isl_min(x,y) ((x) < (y) ? (x) : (y))
#define isl_max(x,y) ((x) > (y) ? (x) : (y))
#define isl_floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))

#define mat_cell(A,i,j) A[i][j]
#define tab_cell(A,i,j) A[((i)*C + (j))]
#define tile_cell(A,i,j) mat_cell(A,i-bot0,j-left0)

static inline int get_cache_line_size(){
  int cache_line_size;
  #if defined(linux)
    cache_line_size = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);
  #else
    cache_line_size = 64;
    printf("Warning using default cache line size : %i", cache_line_size);
  #endif
  return cache_line_size;
}

static inline float ** allocmatrix(int rows, int cols){
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

// Allocate an array where each line is aligned in memory.
static inline float ** alloc_array_lines(int R, int C, int cache_line_size){
  int memalign;
  float ** mat;
  mat = (float **) malloc(R * sizeof(float *));

  for(int i = 0; i < R; i++){

      memalign = posix_memalign((void **)&mat[i], cache_line_size,
        sizeof(float) * C);

      if(memalign != 0){
        printf("Error while allocating tile at %i\n" , i);
        return NULL;
      } else if (memalign == EINVAL) {
        printf("alloc_array_lines error : alignment parameter is not a power"
          " of two \n");
        return NULL;
      } else if (memalign == ENOMEM) {
        printf("alloc_array_lines error : insufficient memory\n");
        return NULL;
      }
  }
  if( mat == NULL){
    printf("Error while allocating matrix !\n");
    return NULL;
  }

  return mat;
}

// Array [i][j] is j'th element of i'th tile, each tile is aligned in memory
static inline float ** alloc_array_tiles(int R, int C, int TSIZEX, int TSIZEY){
  int cache_line_size;
  int memalign;
  int x,y;
  float ** mat;

  #if defined(linux)
    cache_line_size = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);
  #else
    cache_line_size = 64;
    printf("Warning using default cache line size : %i", cache_line_size);
  #endif

  int n_tiles_x = R / TSIZEX;
  int n_tiles_y = R / TSIZEY;
  mat = (float **) malloc(n_tiles_x * n_tiles_y * sizeof(float *));

  for(int i = 0; i < n_tiles_x; i++){
    for(int j = 0; i < n_tiles_y; j++){
      x = isl_min( (i + 1) * TSIZEX , R) - isl_max(i * TSIZEX, 1);
      y = isl_min( (j + 1) * TSIZEY, C-1) - isl_max(j * TSIZEY, 1);
      memalign = posix_memalign((void **)&mat[i], cache_line_size,
        sizeof(float) * x * y);

      if(memalign != 0){
        printf("Error while allocating tile at %i : %i\n" , i, j);
        return NULL;
      }
    }
  }

  if( mat == NULL){
    printf("Error while allocating matrix !\n");
    return NULL;
  }

  return mat;
}

static inline int freematrix(float ** Mat, int rows){
  for(int i = 0; i < rows; i++){
    free(Mat[i]);
  }
  free(Mat);
  return 0;
}

#ifndef VERSION_ALIGNED
  extern "C" void  pipeline_harris(int  C, int  R, void * img_void_arg, void * harris_void_arg);
#else
  extern "C" void  pipeline_harris(int  C, int  R, float ** img_void_arg, float ** harris_void_arg);
#endif

// Macros

#define sobelY(A,i,j) (tab_cell(A, i-1, j-1) *( -0.0833333333333f) + (tab_cell(A, i-1, j+1) * 0.0833333333333f) + \
                      (tab_cell(A, i, j-1) * -0.166666666667f) + (tab_cell(A, i, j+1) * 0.166666666667f) + \
                       (tab_cell(A, i+1, j-1) * -0.0833333333333f) +(tab_cell(A, i+1, j+1) * 0.0833333333333f))

#define sobelX(A,i,j) (tab_cell(A, i-1, j-1) *( -0.0833333333333f) + (tab_cell(A, i+1, j-1) * 0.0833333333333f) + \
                      (tab_cell(A, i+1, j) * 0.166666666667f) + (tab_cell(A, i-1, j) * (-0.166666666667f)) + \
                       (tab_cell(A, i-1, j+1) * (-0.0833333333333f)) + (tab_cell(A, i+1, j+1) * 0.0833333333333f))

#define m_sobelY(A,i,j) (mat_cell(A, i-1, j-1) *( -0.0833333333333f) + (mat_cell(A, i-1, j+1) * 0.0833333333333f) + \
                     (mat_cell(A, i, j-1) * -0.166666666667f) + (mat_cell(A, i, j+1) * 0.166666666667f) + \
                      (mat_cell(A, i+1, j-1) * -0.0833333333333f) +(mat_cell(A, i+1, j+1) * 0.0833333333333f))

#define m_sobelX(A,i,j) (mat_cell(A, i-1, j-1) *( -0.0833333333333f) + (mat_cell(A, i+1, j-1) * 0.0833333333333f) + \
                   (mat_cell(A, i+1, j) * 0.166666666667f) + (mat_cell(A, i-1, j) * (-0.166666666667f)) + \
                    (mat_cell(A, i-1, j+1) * (-0.0833333333333f)) + (mat_cell(A, i+1, j+1) * 0.0833333333333f))

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

#endif /*HARRIS_H*/
