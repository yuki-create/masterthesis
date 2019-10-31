#include <stdio.h>
#include <stdlib.h>

/* #include <cblas.h> このプログラムでは要らない */
#include <lapacke.h>

/* 行列の中身を表示する */
static void print_mat(const int m, const int n, const double *a, const int lda)
{
  int i, j;

  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j) {
      printf("% 10.3g ", a[j*lda + i]);
    }
    printf("\n");
  }
  printf("\n");
}

/* メイン関数 */
int main()
{
  const int n = 3;
  const int lda = n;

  double a[lda*n], ev[n];
  lapack_int info;
  int i, j;

  /* 行列に値を入れる */
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      a[lda*i + j] = (i == j || i*j != 0) ? 1 : 0;
    }
  }

  /* ちゃんと入ったか確認 */
  printf("a:\n");
  print_mat(n, n, a, lda);

  /* lapacke経由でdsyev(対称行列の固有値・固有ベクトルを計算する関数)を呼び出す */
  info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', n, a, lda, ev);
  if (info != 0) {
    printf("dsyev failed. info = %d\n", info);
    return EXIT_FAILURE;
  }

  /* 固有ベクトル表示 */
  printf("eigen vectors:\n");
  print_mat(n, n, a, lda);

  /* 固有値表示 */
  printf("eigen values:\n");
  for (i = 0; i < n; ++i) {
    printf("% 10.3g ", ev[i]);
  }
  printf("\n");

  return EXIT_SUCCESS;
}
