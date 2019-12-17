// 実行は ./ms N NSTEP
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "lapacke.h"
#include "cblas.h"
// decralation of functions
void genGraph(int size, int graph[size][size]);
void printGraph(int size, int graph[size][size], int table[size][size]);
void initArrayZero(int size, double array[size]);
void initCoordinates(int size, double array_x[size], double array_y[size], double array_x_d[size], double array_y_d[size]);
void initSpringLength(int size, int natu_l, double array_l[size], double array_l_d[size]);
void initSpringParameters(int size, double mu_k, double sigma_k,double mu_g, double sigma_g, double array_k[size], double array_g[size]);
double rand_normal(double mu, double sigma);
double Uniform( void );
void initFiles();
void initForLapack();
void initMinuteInitialStates();


int main(int argc, char *argv[]){
  //////// parameters to be costomized ////////
  const int N = atoi(argv[1]); // number of mass points
  const int M = (sqrt(N)-1)*sqrt(N)*2; // number of springs (root_N-1)*root_N*2
  const char *dirname = ".";
// washout, learning, evaluating term (time steps)
  const int WASHOUT = 5000;
  const int LEARNING = 5000;
  const int EVAL = 1000;
  // simulating time steps
  const int NSTEP = WASHOUT+LEARNING+EVAL;

  const double dt = 0.0025;
  const int T_input = 1; // adjust frequency of input signal
  const double natu_l = 1.0;
  const double w_in[] = {1.0};
  int fixed_p[] = {4,20}; // index array of fixed points
  int in_p[] = {N-1}; // index array of input points
  // initial purtubation to index 0
  int seed_flag = 0;

  //////// variables dosen't need to be costomized ////////
  double k[M];
  double gamma1[M];
  double mu_k = strtod(argv[2],NULL);
  double sigma_k = strtod(argv[3],NULL);
  double mu_g = strtod(argv[4],NULL);
  double sigma_g = strtod(argv[5],NULL);
  int fixed_num = 0; // number of fixed points (elements of fixed_p)
  int in_num = 0; // number of inputted points (elements of in_p)
  double input[20]={0.0}; // input timesiries
  double x[N];
  double u[N]; // dx/dt
  double y[N];
  double v[N]; // dy/dt
  double m[N];
  double l[M]; // length of spring at time n*dt as outpus.
  int G[N][N]; // adjacency matrix (symmetric matrix)
  int p2l_mat[N][N]; // point index x[idx1],x[idx2] -> spring index p2l_mat[idx1][idx2] (upper triangular matrix)
  double force_x[N];

  // variables and parameters for RK4
  double ku1[N], ku2[N], ku3[N], ku4[N], kx1[N], kx2[N], kx3[N], kx4[N];
  double kv1[N], kv2[N], kv3[N], kv4[N], ky1[N], ky2[N], ky3[N], ky4[N];

  // variables for mass-spring system with minute initial state
  double x_d[N];
  double u_d[N]; // dx/dt
  double y_d[N];
  double v_d[N]; // dy/dt
  double ku1_d[N], ku2_d[N], ku3_d[N], ku4_d[N], kx1_d[N], kx2_d[N], kx3_d[N], kx4_d[N];
  double kv1_d[N], kv2_d[N], kv3_d[N], kv4_d[N], ky1_d[N], ky2_d[N], ky3_d[N], ky4_d[N];
  double l_d[N];

  // variables for NARMA models
  double o_nrm2[3]={0.0}; // NARMA2 output
  double o_nrm10[11]={0.0}; // NARMA10 output
  double o_nrm20[21]={0.0}; // NARMA20 output

  // variables for learning
  double *T, *W_out, *L;

  // variables for EVAL phase
  double o_ms[3]={0.0}; // outputs of MS [O_ms(t) for narma2, narma10, narma20]
  double squared_err[3]={0.0}; // squared errors [narma2, narma10, narma20]
  double normalizer[3]={0.0}; // normalizers [narma2, narma10, narma20]
  double err[3]={0.0}; // normilized mean squaered errors

  // variables for lyapunov exponent
  double lyapunov = 0.0;
  double sum_log = 0.0;
  double initial_d = 0.0; // initial distance of two systems
  double norm2_pre = 0.0;

  // parameters for file I/O
  FILE *fp1; // for export coodinates
  FILE *fp2; // for export length of springs
  FILE *fp3; // for outputs of MS and NARMA after learning
  FILE *fp4; // for parameters and results file
  FILE *fp5; // for export lyapunov exponent
  FILE *fp6; // for export many results
  char filename1[40],filename2[40],filename3[40],filename4[40],filename5[40],filename6[40];
  int i,j;

  // 初期化
  initArrayZero(N,u);
  initArrayZero(N,v);
  initArrayZero(N,u_d);
  initArrayZero(N,v_d);
  initArrayZero(N,m);
  initArrayZero(N,force_x);
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      p2l_mat[i][j] = -1;
    }
  }
  initCoordinates(N,x,y,x_d,y_d);
  initSpringLength(M,natu_l,l,l_d);
  initSpringParameters(M,mu_k,sigma_k,mu_g,sigma_g,k,gamma1);

  genGraph(N,G);
  //printGraph(N,G,p2l_mat);

  srand((unsigned)time(NULL));
  double tmp = 0;
  int tmp2=0;
  for(i=0;i<10000;i++){
    tmp = round(rand_normal(100,3)*10);
    tmp2 = (int)tmp;
    printf("%3d\n",tmp2);
  }
}

//隣接行列を生成
void genGraph(int size, int graph[size][size]){
  int root_size = sqrt(size);
  int i,j;
  for(i=0;i<size;i++){
    for(j=0;j<size;j++){
      graph[i][j] = 0;
    }
  }
  for(i=0;i<size;i++){
    if(i-root_size >= 0)
    graph[i][i-root_size] = 1;
    if(i % root_size != root_size-1)
    graph[i][i+1] = 1;
    if(i % root_size != 0)
    graph[i][i-1] = 1;
    if(i+root_size <= size-1)
    graph[i][i+root_size] = 1;
  }
}
//隣接行列の要素出力
void printGraph(int size, int graph[size][size], int table[size][size]){
  int i,j;
  printf("adjacency matrix G =\n");
  //列のindexを表示
  printf("   ");
  for(j=0;j<size;j++){
    printf("%d ",j);
  }
  printf("\n");

  printf("   ");
  for(j=0;j<size;j++){
    printf("－");
  }
  printf("\n");

  for(i=0;i<size;i++){
    //行のindexを表示
    printf("%d | ",i);
    //行毎に要素を出力
    for(j=0;j<size;j++){
      printf("%d ",graph[i][j]);
    }
    printf("\n");
  }
  //indexなし
  printf("with no index\n");
  for(i=0;i<size;i++){
    //行毎に要素を出力
    for(j=0;j<size;j++){
      printf("%d, ",graph[i][j]);
    }
    printf("\n");
  }
  // p2l_matを出力
  printf("p2l_mat =\n");
  //列のindexを表示
  printf("   ");
  for(j=0;j<size;j++){
    printf("%d ",j);
  }
  printf("\n");

  printf("   ");
  for(j=0;j<size;j++){
    printf("－");
  }
  printf("\n");

  for(i=0;i<size;i++){
    //行のindexを表示
    printf("%d | ",i);
    //行毎に要素を出力
    for(j=0;j<size;j++){
      printf("%d ",table[i][j]);
    }
    printf("\n");
  }
}

void initArrayZero(int size, double array[size]){
  int i;
  for(i=0;i<size;i++){
    array[i] = 0.0;
  }
}

void initCoordinates(int size, double array_x[size], double array_y[size], double array_x_d[size], double array_y_d[size])
{
  // init coordinates of mass points
  int root_size = sqrt(size);
  int i,j;
  for(i=0;i<root_size;i++){
    for(j=0;j<root_size;j++){
      array_x[root_size*i+j] = j;
      array_y[root_size*i+j] = i;
      array_x_d[root_size*i+j] = j;
      array_y_d[root_size*i+j] = i;
    }
  }
    array_x_d[0] = array_x_d[0] + 0.000000000001;
    array_y_d[0] =array_y_d[0] + 0.000000000001;
}

void initSpringLength(int size, int natu_l, double array_l[size], double array_l_d[size]){
  int i;
  for(i=0;i<size;i++){
    array_l[i] = natu_l;
    array_l_d[i]  = natu_l;
  }
}
void initSpringParameters(int size, double mu_k, double sigma_k,double mu_g, double sigma_g, double array_k[size], double array_g[size]){
  int i;
  for(i=0;i<size;i++){
    array_k[i] = rand_normal(mu_k,sigma_k);
    array_g[i] = rand_normal(mu_g,sigma_g);
  }
}

double rand_normal(double mu, double sigma){
      double z=sqrt( -2.0*log(Uniform()) ) * sin( 2.0*M_PI*Uniform() );
      return mu + sigma*z;
}
double Uniform( void ){
    return ((double)rand()+1.0)/((double)RAND_MAX+2.0);
}
