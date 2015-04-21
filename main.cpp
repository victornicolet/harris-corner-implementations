#include <opencv2/opencv.hpp>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "bench_source/harris.h"

//#define CHECK_FINAL_RESULT = true;

using namespace std;

static int minruns = 5;

int main(int argc, char ** argv)
{
  clock_t begin, end;
  int stime, avgt;
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
      printf("No image data \n");
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

  res = (float *) malloc((sizeof(float))*C*R);

  if(res == NULL)
  {
    printf("Error while allocating result table of size %ld B\n", (sizeof(float))*C*R );
    return -1;
  }
  cv::Scalar sc;
  //data = (float *) image.data ;
  data = (float *) malloc(sizeof(float) * R * C);
  for(int i= 0; i < R;i++){
    for(int j = 0; j < C;j++){
      sc = image.at<uchar>(i, j) ;
      data[i * C + j] = sc.val[0]/255;
    }
  }

  // Running tests
  avgt = 0.0f;
  int init,finish;
  init = clock();
  for(int run = 1; run <= nruns; run++)
  {
    printf("Run %i : ",run);
    begin = clock();
    pipeline_harris(C, R, data, res);
    end = clock();
    stime = end - begin;
    printf("\t(%i)\t %f\n", omp_get_thread_num(),(float) stime  / CLOCKS_PER_SEC );
    avgt += stime;
  }
  finish = clock();
  if(avgt == 0)
  {
    printf("Error : running took no time !");
    return -1;
  }
  printf("Average time : %f ms\n", (float) (avgt / (nruns * CLOCKS_PER_SEC)));
  printf("Diff beteween total loop time and sum of running times :\n ");
  printf("\t %f\n",(float) (finish-init-avgt)/ CLOCKS_PER_SEC);

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
