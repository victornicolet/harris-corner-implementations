/*Copyright (c) 2015-, Ravi Teja Mullapudi, Vinay Vasista, CSA Indian
  Institute of Science
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:
  1. Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
  2. Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.
  3. Neither the name of the Indian Institute of Science nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY Ravi Teja Mullapudi, Vinay Vasista, CSA Indian
  Institute of Science ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
  BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL Ravi Teja Mullapudi, Vinay
  Vasista, CSA Indian Institute of Science BE LIABLE FOR ANY DIRECT, INDIRECT,
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
  OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
  ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <cmath>
#include <string.h>
#define isl_min(x,y) ((x) < (y) ? (x) : (y))
#define isl_max(x,y) ((x) > (y) ? (x) : (y))
#define isl_floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
extern "C" void  pipeline_harris(int  C, int  R, void * img_void_arg, void * harris_void_arg)
{
  float * img;
  img = (float *) (img_void_arg);
  float * harris;
  harris = (float *) (harris_void_arg);
  #pragma omp parallel for schedule(static)
  for (int  Ti = -1; (Ti <= (R / 32)); Ti = (Ti + 1))
  {
    float  Iy[34][258];
    float  Ix[34][258];
    for (int  Tj = -1; (Tj <= (C / 256)); Tj = (Tj + 1))
    {

      int  _ct0 = ((R < ((32 * Ti) + 33))? R: ((32 * Ti) + 33));
      int  _ct1 = ((1 > (32 * Ti))? 1: (32 * Ti));
      int  _ct2 = ((C < ((256 * Tj) + 257))? C: ((256 * Tj) + 257));
      int  _ct3 = ((1 > (256 * Tj))? 1: (256 * Tj));

      for (int  _i0 = _ct1; (_i0 <= _ct0); _i0 = (_i0 + 1))
      {
        #pragma ivdep
        for (int  _i1 = _ct3; (_i1 <= _ct2); _i1 = (_i1 + 1))
        {
          Iy[((-32 * Ti) + _i0)][((-256 * Tj) + _i1)] = ((((((img[(((-1 + _i0) * (C + 2)) + (-1 + _i1))] * -0.0833333333333f) + (img[(((-1 + _i0) * (C + 2)) + (1 + _i1))] * 0.0833333333333f)) + (img[((_i0 * (C + 2)) + (-1 + _i1))] * -0.166666666667f)) + (img[((_i0 * (C + 2)) + (1 + _i1))] * 0.166666666667f)) + (img[(((1 + _i0) * (C + 2)) + (-1 + _i1))] * -0.0833333333333f)) + (img[(((1 + _i0) * (C + 2)) + (1 + _i1))] * 0.0833333333333f));
        }
      }
      for (int  _i0 = _ct1; (_i0 <= _ct0); _i0 = (_i0 + 1))
      {
        #pragma ivdep
        for (int  _i1 = _ct3; (_i1 <= _ct2); _i1 = (_i1 + 1))
        {
          Ix[((-32 * Ti) + _i0)][((-256 * Tj) + _i1)] = ((((((img[(((-1 + _i0) * (C + 2)) + (-1 + _i1))] * -0.0833333333333f) + (img[(((1 + _i0) * (C + 2)) + (-1 + _i1))] * 0.0833333333333f)) + (img[(((-1 + _i0) * (C + 2)) + _i1)] * -0.166666666667f)) + (img[(((1 + _i0) * (C + 2)) + _i1)] * 0.166666666667f)) + (img[(((-1 + _i0) * (C + 2)) + (1 + _i1))] * -0.0833333333333f)) + (img[(((1 + _i0) * (C + 2)) + (1 + _i1))] * 0.0833333333333f));
        }
      }

      int  _ct8 = (((R - 1) < ((32 * Ti) + 32))? (R - 1): ((32 * Ti) + 32));
      int  _ct9 = ((2 > ((32 * Ti) + 1))? 2: ((32 * Ti) + 1));
      for (int  _i0 = _ct9; (_i0 <= _ct8); _i0 = (_i0 + 1))
      {
        int  _ct10 = (((C - 1) < ((256 * Tj) + 256))? (C - 1): ((256 * Tj) + 256));
        int  _ct11 = ((2 > ((256 * Tj) + 1))? 2: ((256 * Tj) + 1));
        #pragma ivdep
        for (int  _i1 = _ct11; (_i1 <= _ct10); _i1 = (_i1 + 1))
        {
          harris[((_i0 * (2 + C)) + _i1)] = ((((((((((((Ix[(-1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))] * Ix[(-1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))]) + (Ix[(-1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)] * Ix[(-1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)])) + (Ix[(-1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))] * Ix[(-1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))])) + (Ix[((-32 * Ti) + _i0)][(-1 + ((-256 * Tj) + _i1))] * Ix[((-32 * Ti) + _i0)][(-1 + ((-256 * Tj) + _i1))])) + (Ix[((-32 * Ti) + _i0)][((-256 * Tj) + _i1)] * Ix[((-32 * Ti) + _i0)][((-256 * Tj) + _i1)])) + (Ix[((-32 * Ti) + _i0)][(1 + ((-256 * Tj) + _i1))] * Ix[((-32 * Ti) + _i0)][(1 + ((-256 * Tj) + _i1))])) + (Ix[(1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))] * Ix[(1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))])) + (Ix[(1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)] * Ix[(1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)])) + (Ix[(1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))] * Ix[(1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))])) * (((((((((Iy[(-1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))] * Iy[(-1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))]) + (Iy[(-1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)] * Iy[(-1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)])) + (Iy[(-1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))] * Iy[(-1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))])) + (Iy[((-32 * Ti) + _i0)][(-1 + ((-256 * Tj) + _i1))] * Iy[((-32 * Ti) + _i0)][(-1 + ((-256 * Tj) + _i1))])) + (Iy[((-32 * Ti) + _i0)][((-256 * Tj) + _i1)] * Iy[((-32 * Ti) + _i0)][((-256 * Tj) + _i1)])) + (Iy[((-32 * Ti) + _i0)][(1 + ((-256 * Tj) + _i1))] * Iy[((-32 * Ti) + _i0)][(1 + ((-256 * Tj) + _i1))])) + (Iy[(1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))] * Iy[(1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))])) + (Iy[(1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)] * Iy[(1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)])) + (Iy[(1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))] * Iy[(1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))]))) - ((((((((((Ix[(-1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))] * Iy[(-1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))]) + (Ix[(-1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)] * Iy[(-1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)])) + (Ix[(-1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))] * Iy[(-1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))])) + (Ix[((-32 * Ti) + _i0)][(-1 + ((-256 * Tj) + _i1))] * Iy[((-32 * Ti) + _i0)][(-1 + ((-256 * Tj) + _i1))])) + (Ix[((-32 * Ti) + _i0)][((-256 * Tj) + _i1)] * Iy[((-32 * Ti) + _i0)][((-256 * Tj) + _i1)])) + (Ix[((-32 * Ti) + _i0)][(1 + ((-256 * Tj) + _i1))] * Iy[((-32 * Ti) + _i0)][(1 + ((-256 * Tj) + _i1))])) + (Ix[(1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))] * Iy[(1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))])) + (Ix[(1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)] * Iy[(1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)])) + (Ix[(1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))] * Iy[(1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))])) * (((((((((Ix[(-1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))] * Iy[(-1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))]) + (Ix[(-1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)] * Iy[(-1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)])) + (Ix[(-1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))] * Iy[(-1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))])) + (Ix[((-32 * Ti) + _i0)][(-1 + ((-256 * Tj) + _i1))] * Iy[((-32 * Ti) + _i0)][(-1 + ((-256 * Tj) + _i1))])) + (Ix[((-32 * Ti) + _i0)][((-256 * Tj) + _i1)] * Iy[((-32 * Ti) + _i0)][((-256 * Tj) + _i1)])) + (Ix[((-32 * Ti) + _i0)][(1 + ((-256 * Tj) + _i1))] * Iy[((-32 * Ti) + _i0)][(1 + ((-256 * Tj) + _i1))])) + (Ix[(1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))] * Iy[(1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))])) + (Ix[(1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)] * Iy[(1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)])) + (Ix[(1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))] * Iy[(1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))])))) - ((0.04f * ((((((((((Ix[(-1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))] * Ix[(-1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))]) + (Ix[(-1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)] * Ix[(-1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)])) + (Ix[(-1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))] * Ix[(-1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))])) + (Ix[((-32 * Ti) + _i0)][(-1 + ((-256 * Tj) + _i1))] * Ix[((-32 * Ti) + _i0)][(-1 + ((-256 * Tj) + _i1))])) + (Ix[((-32 * Ti) + _i0)][((-256 * Tj) + _i1)] * Ix[((-32 * Ti) + _i0)][((-256 * Tj) + _i1)])) + (Ix[((-32 * Ti) + _i0)][(1 + ((-256 * Tj) + _i1))] * Ix[((-32 * Ti) + _i0)][(1 + ((-256 * Tj) + _i1))])) + (Ix[(1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))] * Ix[(1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))])) + (Ix[(1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)] * Ix[(1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)])) + (Ix[(1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))] * Ix[(1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))])) + (((((((((Iy[(-1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))] * Iy[(-1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))]) + (Iy[(-1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)] * Iy[(-1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)])) + (Iy[(-1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))] * Iy[(-1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))])) + (Iy[((-32 * Ti) + _i0)][(-1 + ((-256 * Tj) + _i1))] * Iy[((-32 * Ti) + _i0)][(-1 + ((-256 * Tj) + _i1))])) + (Iy[((-32 * Ti) + _i0)][((-256 * Tj) + _i1)] * Iy[((-32 * Ti) + _i0)][((-256 * Tj) + _i1)])) + (Iy[((-32 * Ti) + _i0)][(1 + ((-256 * Tj) + _i1))] * Iy[((-32 * Ti) + _i0)][(1 + ((-256 * Tj) + _i1))])) + (Iy[(1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))] * Iy[(1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))])) + (Iy[(1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)] * Iy[(1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)])) + (Iy[(1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))] * Iy[(1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))])))) * ((((((((((Ix[(-1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))] * Ix[(-1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))]) + (Ix[(-1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)] * Ix[(-1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)])) + (Ix[(-1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))] * Ix[(-1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))])) + (Ix[((-32 * Ti) + _i0)][(-1 + ((-256 * Tj) + _i1))] * Ix[((-32 * Ti) + _i0)][(-1 + ((-256 * Tj) + _i1))])) + (Ix[((-32 * Ti) + _i0)][((-256 * Tj) + _i1)] * Ix[((-32 * Ti) + _i0)][((-256 * Tj) + _i1)])) + (Ix[((-32 * Ti) + _i0)][(1 + ((-256 * Tj) + _i1))] * Ix[((-32 * Ti) + _i0)][(1 + ((-256 * Tj) + _i1))])) + (Ix[(1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))] * Ix[(1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))])) + (Ix[(1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)] * Ix[(1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)])) + (Ix[(1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))] * Ix[(1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))])) + (((((((((Iy[(-1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))] * Iy[(-1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))]) + (Iy[(-1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)] * Iy[(-1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)])) + (Iy[(-1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))] * Iy[(-1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))])) + (Iy[((-32 * Ti) + _i0)][(-1 + ((-256 * Tj) + _i1))] * Iy[((-32 * Ti) + _i0)][(-1 + ((-256 * Tj) + _i1))])) + (Iy[((-32 * Ti) + _i0)][((-256 * Tj) + _i1)] * Iy[((-32 * Ti) + _i0)][((-256 * Tj) + _i1)])) + (Iy[((-32 * Ti) + _i0)][(1 + ((-256 * Tj) + _i1))] * Iy[((-32 * Ti) + _i0)][(1 + ((-256 * Tj) + _i1))])) + (Iy[(1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))] * Iy[(1 + ((-32 * Ti) + _i0))][(-1 + ((-256 * Tj) + _i1))])) + (Iy[(1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)] * Iy[(1 + ((-32 * Ti) + _i0))][((-256 * Tj) + _i1)])) + (Iy[(1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))] * Iy[(1 + ((-32 * Ti) + _i0))][(1 + ((-256 * Tj) + _i1))])))));
        }
      }
    }
  }
}
