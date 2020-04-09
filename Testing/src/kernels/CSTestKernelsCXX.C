#include <CSTestKernels.H>
#include <iostream>
//
// matrix multiplication
//
//       [ 1  2  3  4  5]
//       [ 6  7  8  9 10]
// [A] = [11 12 13 14 15]
//       [16 17 18 19 20]
//       [21 22 23 24 25]
//
//             [ 215  230  245  260  275]
//             [ 490  530  570  610  650]
// [A]*[A] =   [ 765  830  895  960 1025]
//             [1040 1130 1220 1310 1400]
//             [1315 1430 1545 1660 1775]
//
int ICEUnrollCxxTestKernel(int &result) {

  int m = 5, p = 5, q = 5, c, d, k, sum = 0;
  int first[m][p], second[p][q], multiply[m][q];

  sum = 0;
  for (c = 0; c < m; c++) {
    for (d = 0; d < p; d++) {
      sum++;
      first[c][d] = sum;
    }
  }

  sum = 0;
  for (c = 0; c < p; c++) {
    for (d = 0; d < q; d++) {
      sum++;
      second[c][d] = sum;
    }
  }

  //sum = 0;
#pragma @ICE loop=unrollTestCXX
  for (c = 0; c < m; c++) {
    for (d = 0; d < q; d++) {
      for (k = 0; k < p; k++) {
        //sum = sum + first[c][k]*second[k][d];
        multiply[c][d] = multiply[c][d] + first[c][k]*second[k][d];
      }
      //multiply[c][d] = sum;
      //sum = 0;
    }
  }
#pragma @ICE endloop

  result=multiply[2][2];
  return 0;
}

//
// matrix multiplication
//
//       [ 1  2  3  4  5]
//       [ 6  7  8  9 10]
// [A] = [11 12 13 14 15]
//       [16 17 18 19 20]
//       [21 22 23 24 25]
//
//             [ 215  230  245  260  275]
//             [ 490  530  570  610  650]
// [A]*[A] =   [ 765  830  895  960 1025]
//             [1040 1130 1220 1310 1400]
//             [1315 1430 1545 1660 1775]
//
int ICEInterchangeCxxTestKernel(int &result) {

  int m = 5, p = 5, q = 5, c, d, k, sum = 0;
  int first[m][p], second[p][q], multiply[m][q];

  sum = 0;
  for (c = 0; c < m; c++) {
    for (d = 0; d < p; d++) {
      sum++;
      first[c][d] = sum;
    }
  }

  sum = 0;
  for (c = 0; c < p; c++) {
    for (d = 0; d < q; d++) {
      sum++;
      second[c][d] = sum;
    }
  }

  for (c = 0; c < m; c++) {
    for (d = 0; d < q; d++) {
      multiply[c][d] = 0;
    }
  }

#pragma @ICE loop=interchangeTestCXX
  for (c = 0; c < m; c++) {
    for (d = 0; d < q; d++) {
      for (k = 0; k < p; k++) {
        //sum = sum + first[c][k]*second[k][d];
        multiply[c][d] = multiply[c][d] + first[c][k]*second[k][d];
      }
    }
  }
#pragma @ICE endloop

  result=multiply[2][2];
  return 0;
}

//
// matrix multiplication
//
//       [ 1  2  3  4  5]
//       [ 6  7  8  9 10]
// [A] = [11 12 13 14 15]
//       [16 17 18 19 20]
//       [21 22 23 24 25]
//
//             [ 215  230  245  260  275]
//             [ 490  530  570  610  650]
// [A]*[A] =   [ 765  830  895  960 1025]
//             [1040 1130 1220 1310 1400]
//             [1315 1430 1545 1660 1775]
//
int ICETileCxxTestKernel(int &result) {

  int m = 5, p = 5, q = 5, c, d, k, sum = 0;
  int first[m][p], second[p][q], multiply[m][q];

  sum = 0;
  for (c = 0; c < m; c++) {
    for (d = 0; d < p; d++) {
      sum++;
      first[c][d] = sum;
    }
  }

  sum = 0;
  for (c = 0; c < p; c++) {
    for (d = 0; d < q; d++) {
      sum++;
      second[c][d] = sum;
    }
  }

  sum = 0;
#pragma @ICE loop=tileTestCXX
  for (c = 0; c < m; c++) {
    for (d = 0; d < q; d++) {
      for (k = 0; k < p; k++) {
        sum = sum + first[c][k]*second[k][d];
      }
      multiply[c][d] = sum;
      sum = 0;
    }
  }
#pragma @ICE endloop

  result=multiply[2][2];
  return 0;
}
//
// matrix multiplication
//
//       [ 1  2  3  4  5]
//       [ 6  7  8  9 10]
// [A] = [11 12 13 14 15]
//       [16 17 18 19 20]
//       [21 22 23 24 25]
//
//             [ 215  230  245  260  275]
//             [ 490  530  570  610  650]
// [A]*[A] =   [ 765  830  895  960 1025]
//             [1040 1130 1220 1310 1400]
//             [1315 1430 1545 1660 1775]
//
int ICEStripMineCxxTestKernel(int &result) {

  int m = 5, p = 5, q = 5, c, d, k, sum = 0;
  int first[m][p], second[p][q], multiply[m][q];

  sum = 0;
  for (c = 0; c < m; c++) {
    for (d = 0; d < p; d++) {
      sum++;
      first[c][d] = sum;
    }
  }

  sum = 0;
  for (c = 0; c < p; c++) {
    for (d = 0; d < q; d++) {
      sum++;
      second[c][d] = sum;
    }
  }

  sum = 0;
#pragma @ICE loop=stripMineTestCXX
  for (c = 0; c < m; c++) {
    for (d = 0; d < q; d++) {
      for (k = 0; k < p; k++) {
        sum = sum + first[c][k]*second[k][d];
      }
      multiply[c][d] = sum;
      sum = 0;
    }
  }
#pragma @ICE endloop

  result=multiply[2][2];
  return 0;
}
