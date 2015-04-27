#ifndef HARRIS_POLYMAGE
#define HARRIS_POLYMAGE

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
#define filter3(A,i,j) tile_cell(A, i-2, j-2) + tile_cell(A, i-2, j-1) + tile_cell(A, i-2, j) + \
                      tile_cell(A, i-2, j+1) + tile_cell(A, i-2, j+2) + tile_cell(A, i-1, j-2) + \
                      tile_cell(A, i-1, j-1) + tile_cell(A, i-1, j) + tile_cell(A, i-1, j+1) + \
                      tile_cell(A, i-1, j+2) + tile_cell(A, i, j-2)   + tile_cell(A, i, j-1)   + \
                      tile_cell(A, i, j)   + tile_cell(A, i, j+1)   + tile_cell(A, i, j+2) + \
                      tile_cell(A, i+1, j-2) + tile_cell(A, i+1, j-1) + tile_cell(A, i+1, j) + \
                      tile_cell(A, i+1, j+1) + tile_cell(A, i+1, j+2)+ tile_cell(A, i+2, j-2) + \
                      tile_cell(A, i+1, j-1) + tile_cell(A, i+2, j) + tile_cell(A, i+2, j+1) + \
                      tile_cell(A, i+2, j+2)

//([tile_cell\(]+)[A,]+ ([i]+[-|+]*[1-2]*[, j]+[-|+]*[1-2]*[\)]+) -> $1 A, $2 * $1 B, $2
#define filter3sq(A,B,i,j)  tile_cell( A, i-2, j-2) * tile_cell( B, i-2, j-2) + \
                            tile_cell( A, i-2, j-1) * tile_cell( B, i-2, j-1) +\
                            tile_cell( A, i-2, j) * tile_cell( B, i-2, j) + \
                            tile_cell( A, i-2, j+1) * tile_cell( B, i-2, j+1) + \
                            tile_cell( A, i-2, j+2) * tile_cell( B, i-2, j+2) + \
                            tile_cell( A, i-1, j-2) * tile_cell( B, i-1, j-2) + \
                            tile_cell( A, i-1, j-1) * tile_cell( B, i-1, j-1) + \
                            tile_cell( A, i-1, j) * tile_cell( B, i-1, j) +  \
                            tile_cell( A, i-1, j+1) * tile_cell( B, i-1, j+1) + \
                            tile_cell( A, i-1, j+2) * tile_cell( B, i-1, j+2) + \
                            tile_cell( A, i, j-2) * tile_cell( B, i, j-2)   + \
                            tile_cell( A, i, j-1) * tile_cell( B, i, j-1)   + \
                            tile_cell( A, i, j) * tile_cell( B, i, j)   + \
                            tile_cell( A, i, j+1) * tile_cell( B, i, j+1)   + \
                            tile_cell( A, i, j+2) * tile_cell( B, i, j+2) + \
                            tile_cell( A, i+1, j-2) * tile_cell( B, i+1, j-2) + \
                            tile_cell( A, i+1, j-1) * tile_cell( B, i+1, j-1) + \
                            tile_cell( A, i+1, j) * tile_cell( B, i+1, j) + \
                            tile_cell( A, i+1, j+1) * tile_cell( B, i+1, j+1) + \
                            tile_cell( A, i+1, j+2) * tile_cell( B, i+1, j+2)+ \
                            tile_cell( A, i+2, j-2) * tile_cell( B, i+2, j-2) + \
                            tile_cell( A, i+1, j-1) * tile_cell( B, i+1, j-1) + \
                            tile_cell( A, i+2, j) * tile_cell( B, i+2, j) + \
                            tile_cell( A, i+2, j+1) * tile_cell( B, i+2, j+1) + \
                            tile_cell( A, i+2, j+2) * tile_cell( B, i+2, j+2)

#define filter2(A,i,j)  tile_cell(A, i-1, j-1) + tile_cell(A, i-1, j) + tile_cell(A, i-1, j+1) + \
                        tile_cell(A, i, j-1)   + tile_cell(A, i, j)   + tile_cell(A, i, j+1)  + \
                        tile_cell(A, i+1, j-1) + tile_cell(A, i+1, j) + tile_cell(A, i+1, j+1)

#define filter2sq(A,B,i,j) ((tile_cell(A, i-1, j-1) * tile_cell(B, i-1, j-1)) + tile_cell(A, i-1, j) * tile_cell(B, i-1, j) +(tile_cell(A, i-1, j+1) * tile_cell(B, i-1, j+1)) + (tile_cell(A, i, j-1) * tile_cell(B, i, j-1)) + (tile_cell(A, i, j) * tile_cell(B, i, j))  + (tile_cell(A, i, j+1) * tile_cell(B, i, j+1))  + (tile_cell(A, i+1, j-1) * tile_cell(B, i+1, j-1)) + (tile_cell(A, i+1, j) * tile_cell(B, i+1, j)) + (tile_cell(A, i+1, j+1) * tile_cell(B, i+1, j+1)))

#endif /*HARRIS_POLYMAGE*/
