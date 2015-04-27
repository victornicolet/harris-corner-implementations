#include <cv.h>
#include <highgui.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "bench_source/harris.h"

#define CHECK_FINAL_RESULT = true;
//#define RUN_PARALLEL = true
using namespace std;

static int minruns = 1;

int main(int argc, char ** argv)
{
  double begin, end;
  double stime, avgt;
  int R, C, nruns;
  float * res;
  float * data;

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

  //res = (float *) calloc(R*C, sizeof(float));
  posix_memalign((void **)&res, 64, sizeof(float)*R*C);

  if(res == NULL)
  {
    printf("Error while allocating result table of size %ld B\n", (sizeof(float))*C*R );
    return -1;
  }
  cv::Scalar sc;
  //data = (float *) image.data ;
  //data = (float *) malloc(sizeof(float) * R * C);
  posix_memalign((void **)&data, 64, sizeof(float)*R*C);
  for(int i= 0; i < R;i++){
    for(int j = 0; j < C;j++){
      sc = image.at<uchar>(i, j) ;
      data[i * C + j] = (float) sc.val[0]/255;
    }
  }

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
  for(int run = 1; run <= nruns; run++)
  {
    begin = omp_get_wtime();
    pipeline_harris(C, R, data, res);
    end = omp_get_wtime();
    stime = end - begin;
    printf("Run %i : \t\t %f ms\n", run, (double) stime * 1000.0 );

    #ifdef RUN_PARALLEL
      #pragma omp atomic
    #endif
    avgt += stime;
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

  #ifdef CHECK_FINAL_RESULT
    cv::Mat imres = cv::Mat(R, C, CV_32F, res);
    cv::namedWindow( "Input", cv::WINDOW_NORMAL );
    cv::imshow( "Input", image );
    cv::namedWindow( "Output", cv::WINDOW_NORMAL );
    cv::imshow( "Output", imres * 65535.0 );
    cv::waitKey(0);
    cv::destroyAllWindows();
    imres.release();
  #endif // CHECK_FINAl_RESULT

  image.release();
  free(res);
  free(data);
  return 0;

}
