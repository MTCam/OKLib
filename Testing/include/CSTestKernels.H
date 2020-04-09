#ifndef __CS_TEST_KERNELS_H__
#define __CS_TEST_KERNELS_H__

#include "FC.h"

extern "C" {
  void FC_MODULE(testkernels,testkernel,TESTKERNELS,TESTKERNEL)(int *result);
};

extern "C" {
  void FC_GLOBAL(iceunrollfortrantestkernel,ICEUNROLLFORTRANTESTKERNEL)(int *result);
};
extern "C" {
  void FC_GLOBAL(iceinterchangefortrantestkernel,ICEINTERCHANGEFORTRANTESTKERNEL)(int *result);
};
extern "C" {
  void FC_GLOBAL(icetilefortrantestkernel,ICETILEFORTRANTESTKERNEL)(int *result);
};
extern "C" {
  void FC_GLOBAL(icestripminefortrantestkernel,ICESTRIPMINEFORTRANTESTKERNEL)(int *result);
};

int ICEUnrollCxxTestKernel(int &result);
int ICEInterchangeCxxTestKernel(int &result);
int ICETileCxxTestKernel(int &result);
int ICEStripMineCxxTestKernel(int &result);

#endif
