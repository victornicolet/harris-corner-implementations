#include <cv.h>
#include <highgui.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "bench_source/harris.h"

// CHECK_FINAL_RESULT toggles the display of the input and output images
#define CHECK_FINAL_RESULT
// Experiment : try to run multiples pipelines in parallel
//#define RUN_PARALLEL = true
// Exp. : align lines of matrixes in memory ( assuming tile size &- cache line size)
#define VERSION_ALIGNED
// Check if image to matrix translation produces the correst output
//#define CHECK_LOADING_DATA
using namespace std;

static int minruns = 1;

int main(int argc, char ** argv)
{
  #ifdef VERSION_ALIGNED
    printf("VERSION_ALIGNED\n");
    int cache_line_size = get_cache_line_size();
  #endif

  double begin, end;
  double stime, avgt;
  int R, C, nruns;
  #ifdef VERSION_ALIGNED
    float ** res;
    float ** data;
  #else
    float * res;
    float * data;
  #endif

  if ( argc != 3 )
  {
      printf("usage: harris <Image_Path> <Nruns>\n");
     return -1;
  }

  cv::Mat image;
  printf("Loading image ....\n");
  image = cv::imread( argv[1], 1 );
  printf("%s successfully loaded !\n\n", argv[1]);

  if ( !image.data )
  {
      printf("No image data ! Are you sure %s is an image ?\n", argv[1]);
      return -1;
  }

  // Convert image input to grayscale floating point
  cv::cvtColor(image, image, cv::COLOR_BGR2GRAY);
  cv::Size size = image.size();
  C = size.width;
  R = size.height;
  nruns = max(minruns, atoi(argv[2]));

  printf("_________________________________________\n");
  printf("Values settings :\n");
  printf("Nruns : %i || %s [%i, %i]\n", nruns, argv[1], R, C);
  printf("_________________________________________\n");

  #ifdef VERSION_ALIGNED
    res = alloc_array_lines(R, C, cache_line_size);
  #else
    res = (float *) calloc(R*C, sizeof(float));
  #endif

  if(res == NULL)
  {
    printf("Error while allocating result table of size %ld B\n", (sizeof(float))*C*R );
    return -1;
  }

  cv::Scalar sc;

  #ifdef VERSION_ALIGNED
    data = (float **) alloc_array_lines(R, C, cache_line_size);
    for(int i= 0; i < R;i++){
      for(int j = 0; j < C;j++){
        sc = image.at<uchar>(i, j) ;
        data[i][j] = (float) sc.val[0]/255;
      }
    }
  #else
    data = (float *) malloc(R*C*sizeof(float));
    for(int i= 0; i < R;i++){
      for(int j = 0; j < C;j++){
        sc = image.at<uchar>(i, j) ;
        data[i*C+j] = (float) sc.val[0]/255;
      }
    }
  #endif

  // Running tests
  avgt = 0.0f;
  double init,finish;
  /*
  Do not use clock here we need elapsed "wall clock time", not total CPU time.
  */
  init =  omp_get_wtime();
  #ifdef RUN_PARALLEL
    #pragma omp parallel for shared(avgt)
  #endif
  for(int run = 0; run <= nruns; run++)
  {
    begin = omp_get_wtime();
    #ifdef VERSION_ALIGNED
      pipeline_harris_aligned(C, R, data, res);
    #else
      pipeline_harris(C, R, data, res);
    #endif
    end = omp_get_wtime();
    stime = end - begin;
    if(run !=0){
      printf("Run %i : \t\t %f ms\n", run, (double) stime * 1000.0 );
      #ifdef RUN_PARALLEL
        #pragma omp atomic
      #endif
      avgt += stime;
    }
  }
  finish =  omp_get_wtime();
  if(avgt == 0)
  {
    printf("Error : running didn't take time !");
    return -1;
  }
  printf("Average time : %f ms\n", (double) (1000.0*avgt / (nruns)));
  printf("Total time : %f ms\n", (double) (finish-init) * 1000.0);

  #ifdef RUN_PARALLEL
    printf("Gain total times to run %i instances in parallel / serial time :\n ", nruns);
    printf("\t %f\n",(double) (finish-init)/(avgt));
  #endif

  // Checking images using OpenCV

  #ifdef VERSION_ALIGNED
    float * t_res = (float *) malloc(sizeof(float)*R*C);
    for(int i = 0; i < R; i++){
      for(int j = 0; j < C; j++){
        t_res[i * C + j] = res[i][j];
      }
    }
  #else
   float * t_res = res;
  #endif

  #ifdef CHECK_LOADING_DATA
    #ifdef VERSION_ALIGNED
    float * t_data = (float *) malloc(sizeof(float)*R*C);
    for(int i = 0; i < R; i++){
      for(int j = 0; j < C; j++){
        t_data[i * C + j] = data[i][j];
      }
    }
    cv::Mat loaded_data = cv::Mat(R,C,CV_32F, t_data);
    #else
    cv::Mat loaded_data = cv::Mat(R,C,CV_32F, data);
    #endif

    cv::namedWindow( "Check data", cv::WINDOW_NORMAL);
    cv::imshow( "Check data", loaded_data);
    cv::waitKey(0);
    cv::destroyAllWindows();
    loaded_data.release();

    #ifdef VERSION_ALIGNED
      free(t_data);
    #endif

  #endif

  #ifdef CHECK_FINAL_RESULT
    cv::Mat imres = cv::Mat(R, C, CV_32F, t_res);
    cv::namedWindow( "Input", cv::WINDOW_NORMAL );
    cv::imshow( "Input", image );
    cv::namedWindow( "Output", cv::WINDOW_NORMAL );
    cv::imshow( "Output", imres * 65535.0 );
    cv::waitKey(0);
    cv::destroyAllWindows();
    imres.release();
    free(t_res);
  #endif // CHECK_FINAL_RESULT

  #ifdef VERSION_ALIGNED
    free(res);
  #endif

  image.release();
  free(data);
  return 0;

}
