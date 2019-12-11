#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cblas.h"

/* A^(T)b = xを求める */
int main(int argc, const char *argv[]){
  double A[2*3] = {1,4,2,5,3,6};
  double b[2*1] = {1,2};
  double x[3*1] = {0,0,0};
  int m = 2;
  int n = 3;
  int lda = m;
  int incb = 1;
  int incx = 1;
  int i;
  /*
  CblasNoTransの場合、 x := alpha*A*b + beta*x
  y := alpha*A^T*x + beta*y
  cblas_dgemv(CblasColMajor, CblasNoTrans, Aの行数, Aの列数,
     Aの定数倍(alpha), 行列Aへのポインタ, LDA（Aの行数）,
     ベクトルbへのポインタ, bの次の要素への増分(incb),
     出力xの定数倍(beta), 計算結果を書き込むベクトルxへのポインタ,
    xの次の要素への増分(incx) );
    */
  cblas_dgemv(CblasColMajor,CblasTrans,m,n,
      1.0, A, lda, b, incb,
      0.0, x, incx);
  /* output */
  for(i=0;i<n;i++){
    printf("%f  ",x[i]);
    printf("\n");
  }

}
