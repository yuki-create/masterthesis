/*
gcc -I/usr/local/include -c least_squares.c -o least_squares.o
gcc -L/Users/saki/lapack-3.8.0 least_squares.o  -llapacke -llapack -lcblas -lblas -lgfortran -lm -o least_squares
./least_squares
*/
/* get x minimizing ||Ax-B||_2 */
#include <stdio.h>
#include "lapacke.h"

int main(int argc, const char *argv[]){
  lapack_int m = 5;  // 係数行列Aの実際の行数（実装ではなく計算上考えている行列）
  lapack_int n = 3;  // 係数行列Aの実際の列数
  lapack_int nrhs = 2;  // 右辺の個数（行列Bとxの列数）
  lapack_int ldA = m;   // 配列Aの第一次元
  // ldA = m （Aの実際の行数）と置いた方が良い
  lapack_int ldB = m;   // 配列Bの第一次元。ldB = mで問題なさそう。
  lapack_int info;    // LAPACKEルーチンの出力用
  double A[5*3] = {1,2,3,4,5,1,3,5,2,4,1,4,2,5,3};
  // COL_MAJORの場合、計算上考えている行列Aの要素を列に沿って順に代入する
  double B[5*2] = {-10,12,14,16,18,-3,14,12,16,16};
 // 出力時、Bは解ベクトルで上書きされる
  int i,j;

  info = LAPACKE_dgels(LAPACK_COL_MAJOR,'N',m,n,nrhs,A,ldA,B,ldB);
  for(i=0;i<n;i++) // Bのn行目までが解ベクトルで更新されている
  {
     for(j=0;j<nrhs;j++) // 各行ごとに出力
     {
        printf("%lf ",B[i+ldB*j]); //
     }
     printf("\n");
  }
  /*for(i=0;i<10;i++){
      printf("%lf ",B[i]);
  }*/
  return(info);
}
