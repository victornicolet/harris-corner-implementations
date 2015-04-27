#ifndef HARRIS
#define HARRIS

#define isl_min(x,y) ((x) < (y) ? (x) : (y))
#define isl_max(x,y) ((x) > (y) ? (x) : (y))
#define isl_floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))

#define mat_cell(A,i,j) A[i][j]
#define tab_cell(A,i,j) A[((i)*C + (j))]
#define tile_cell(A,i,j) mat_cell(A,i-bot0,j-left0)

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

inline float ** allocmatrix(int rows, int cols){
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

inline int freematrix(float ** Mat, int rows){
  for(int i = 0; i < rows; i++){
    free(Mat[i]);
  }
  free(Mat);
  return 0;
}

extern "C" void  pipeline_harris(int  C, int  R, void * img_void_arg, void * harris_void_arg);

#endif /*HARRIS*/
