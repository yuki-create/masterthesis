// 実行は ./ms N NSTEP
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "lapacke.h"
#include "cblas.h"
//////// parameters to be costomized ////////
#define N 36 // number of mass points
#define M 60 // number of springs (root_N-1)*root_N*2
const char *dirname = ".";
/* washout, learning, evaluating term (time steps) */
const int WASHOUT = 1000;
const int LEARNING = 0;
const int EVAL = 0;
/* simulating time steps */
const int NSTEP = WASHOUT+LEARNING+EVAL;

const double dt = 0.0025;
const int T_input = 1; // adjust frequency of input signal
// const double gamma1 = 0.1;
// const double k = 100.0;
const double natu_l = 1.0;
const double w_in[] = {1.0};
const double k_bottom = 1.0;
const double k_top = 2000.0;
const double gamma1_bottom = 0.001;
const double gamma1_top = 0.01;
int fixed_p[] = {}; // index array of fixed points
int in_p[] = {0}; // index array of input points
// initial purtubation to index 0
double k[M];
double gamma1[M];
int seed_flag = 0;

//////// variables dosen't need to be costomized ////////
int fixed_num = 0; // number of fixed points (elements of fixed_p)
int in_num = 0; // number of inputted points (elements of in_p)
double input[20]; // input timesiries
double x[N];
double u[N]; // dx/dt
double y[N];
double v[N]; // dy/dt
double m[N];
double l[M]; // length of spring at time n*dt as outpus.
int G[N][N]; // adjacency matrix (symmetric matrix)
int p2l_mat[N][N]; // point index x[idx1],x[idx2] -> spring index p2l_mat[idx1][idx2] (upper triangular matrix)
double force_x[N];

/* variables and parameters for RK4 */
double ku1[N], ku2[N], ku3[N], ku4[N], kx1[N], kx2[N], kx3[N], kx4[N];
double kv1[N], kv2[N], kv3[N], kv4[N], ky1[N], ky2[N], ky3[N], ky4[N];

/* variables for mass-spring system with minute initial state */
double x_d[N];
double u_d[N]; // dx/dt
double y_d[N];
double v_d[N]; // dy/dt
double ku1_d[N], ku2_d[N], ku3_d[N], ku4_d[N], kx1_d[N], kx2_d[N], kx3_d[N], kx4_d[N];
double kv1_d[N], kv2_d[N], kv3_d[N], kv4_d[N], ky1_d[N], ky2_d[N], ky3_d[N], ky4_d[N];
double l_d[M];

/* variables for NARMA models */
double o_nrm2[3]; // NARMA2 output
double o_nrm10[11]; // NARMA10 output
double o_nrm20[21]; // NARMA20 output

/* variables for learning */
double *T, *W_out, *L;

/* variables for EVAL phase */
double o_ms[3]; // outputs of MS [O_ms(t) for narma2, narma10, narma20]
double squared_err[3]; // squared errors [narma2, narma10, narma20]
double normalizer[3]; // normalizers [narma2, narma10, narma20]
double err[3]; // normilized mean squaered errors

/* variables for lyapunov exponent */
double lyapunov = 0.0;
double sum_log = 0.0;
double initial_d = 0.0; // initial distance of two systems
double norm2_pre = 0.0;
double norm2=0.0;

/* parameters for file I/O */
FILE *fp1; // for export coodinates
FILE *fp2; // for export length of springs
FILE *fp3; // for outputs of MS and NARMA after learning
FILE *fp4; // for parameters and results file
FILE *fp5; // for export lyapunov exponent
FILE *fp6; // for export many results
char filename1[40]; // for export coodinates
char filename2[40];
char filename3[40];
char filename4[40];
char filename5[40];
char filename6[40];
/* decralation of functions */
double rand_normal(double mu, double sigma);
double Uniform(void);
void genGraph();
void printGraph();
void init();
void initFiles();
void initForLapack();
void initMinuteInitialStates();
void updateInput(int time_steps);
void updateNarma2();
void updateNarma10();
void updateNarma20();
void updateLearnigData(int time_steps);
void getWeights();
void printWeights();
void updateOutputsMS();
void updateErr();
void getErr();
void rk4(int n);
void rk4MinuteInitialStates(int n);
void updateLyapunovExponent(int time_steps,double *array1, double *array2);
void getSpringLength(double *array_l, double *array_x, double *array_y);
void exportCoordinates(int time_steps);
void exportLength(int time_steps);
void exportOutputs(int time_steps);
void exportResults(); // results and parameters
void exportLyapunovExponent(int time_steps);
double Fx(double *array1, double *array2, double *array3, int idx1);
double Fy(double *array1, double *array2, double *array3, int idx1);
double f(double *array1, double *array2, int idx1, int idx2);
void eular();
void dammy_genGraph();
void test_getSpringLength();
void test_updateNarma(int time_steps);
void test_updateLearningData();

//main(コマンドライン引数の個数，引数を文字列として保存する配列)
// argv[] = {"./ms", "N", "NSTEP"}
int main(int argc, char *argv[]){
  int i,n;
  fixed_num = sizeof fixed_p / sizeof fixed_p[0];
  in_num = sizeof in_p / sizeof in_p[0];
  if(seed_flag == 1) srand((unsigned)time(NULL));
  init();
  getSpringLength(l,x,y); //系の出力となるばねの長さを求め、配列l[]を更新
  initFiles();
  initForLapack();
  initMinuteInitialStates(); //初期値決定とinitial_dの計算
  // test_getSpringLength();
   // printGraph();
  // printf("k[-100]=%f\n",k[-100]);  //アクセスできてしまう
  // printf("%14.12lf %.12lf %.12lf\n",initial_d,l[0],l_d[0]);
  //for(i=0;i<M;i++){
  //  printf("%.12lf\n",k[i]);
  //printf("%.12lf  %.12lf  %.12lf  %.12lf  %.12lf  %.12lf  %.12lf  %.12lf  \n",x[i],y[i],u[i],v[i],kx2[i],ky2[i],kv2[i],ku2[i]);
//}
  /* washout phase */
  for(n=0;n<WASHOUT;n++){
    updateInput(n); //inputノードに外力を入力
  //  updateNarma2();
  //  updateNarma10();
  //  updateNarma20();
  //printf("%d  ",n);
    rk4(n); //全ての質点の座標の更新
    rk4MinuteInitialStates(n); // 初期値をずらした質点ばね系の更新
    getSpringLength(l,x,y); //系の出力となるばねの長さを求め、配列l[]を更新
   getSpringLength(l_d,x_d,y_d);
    updateLyapunovExponent(n,l,l_d);

//    printf("%.12lf\n",l_d[0]);
  //  exportLyapunovExponent(n);
  //  exportLength(n);
   printf("%d %.12lf\n",n,l[0]);
  //  printf("%d %.12lf\n",n,l[0]);
  //  printf("%d    %14.12lf   %14.12lf\n", time_steps, l[0],l_d[0]);
  }
  printf("N=%d, M=%d\n",N,M);
  printf("initial_d=%.12lf\n",initial_d);
  printf("le=%.12lf\n",lyapunov);
  printf("gamma1[-1]=%.12lf\n",gamma1[-1]);

  /* learning phase */
  for(n=WASHOUT;n<WASHOUT+LEARNING;n++){
    //  test_updateNarma(n);
    updateInput(n); //inputノードに外力を入力
    updateNarma2();
    updateNarma10();
    updateNarma20();
    rk4(n); //全ての質点の座標の更新
    rk4MinuteInitialStates(n); // 初期値をずらした質点ばね系の更新
    getSpringLength(l,x,y); //系の出力となるばねの長さを求め、配列l[]を更新
     getSpringLength(l_d,x_d,y_d);
    updateLyapunovExponent(n,l,l_d);
    //    exportCoordinates(n); //座標のデータをファイル出力
    //    exportLength(n); //ばねの長さをファイル出力
    updateLearnigData(n);
    exportLyapunovExponent(n);
  }

  /* determin W_out[M] */
  getWeights();
  //printWeights();

  /* evaluation phase */
  for(n=WASHOUT+LEARNING;n<WASHOUT+LEARNING+EVAL;n++){
    //  test_updateNarma(n);
    updateInput(n); //inputノードに外力を入力
    updateNarma2();
    updateNarma10();
    updateNarma20();
    rk4(n); //全ての質点の座標の更新
    rk4MinuteInitialStates(n); // 初期値をずらした質点ばね系の更新
    getSpringLength(l,x,y); //系の出力となるばねの長さを求め、配列l[]を更新
    getSpringLength(l_d,x_d,y_d);
    //    exportCoordinates(n); //座標のデータをファイル出力
    //    exportLength(n); //ばねの長さをファイル出力
    updateOutputsMS();
    updateLyapunovExponent(n,l,l_d);
    updateErr();
    exportOutputs(n); //近似結果を出力
    exportLyapunovExponent(n);
  }
  getErr();
  exportResults(); // 近似誤差と設定パラメータの出力
  printf("Lyapunov exponent: %f\n",lyapunov);
  //  test_updateLearningData();
  free(T);
  free(L);
  free(W_out);
  return (0);
}

double rand_normal(double mu, double sigma){
      double z=sqrt( -2.0*log(Uniform()) ) * sin( 2.0*M_PI*Uniform() );
      return mu + sigma*z;
}
double Uniform( void ){
    return ((double)rand()+1.0)/((double)RAND_MAX+2.0);
}

void initForLapack(){
  int i,j=0;
  T = malloc(LEARNING*3*sizeof(double));
  L = malloc(LEARNING*M*sizeof(double));
  W_out = malloc(M*3*sizeof(double));
  /* 考えている行列の、列ごとに要素を代入する */
  for(j=0;j<3;j++){
    for(i=0;i<LEARNING;i++){
      T[j*LEARNING+i] = 0.0;
    }
  }
  for(j=0;j<M;j++){
    for(i=0;i<LEARNING;i++){
      L[j*LEARNING+i] = 0.0;
    }
  }
  for(j=0;j<3;j++){
    for(i=0;i<M;i++){
      W_out[j*M+i] = 0.0;
    }
  }
}

void init(){
  /* init coordinates of mass points */
  int root_N = sqrt(N);
  int i=0,j=0;

  genGraph();
  // dammy_genGraph();
  for(i=0;i<root_N;i++){
    for(j=0;j<root_N;j++){
      x[root_N*i+j] = j;
      y[root_N*i+j] = i;
    }
  }
  for(i=0;i<M;i++){
  //  l[i] = natu_l;
  //  l_d[i]  = natu_l;
    // k_bottomからk_topまでのランダムな実数
    //k[i] = (double)rand()/RAND_MAX*(k_top-k_bottom)+k_bottom;
    // gamma1_bottomからgamma1_topまでのランダムな実数
    //gamma1[i] = (double)rand()/RAND_MAX*(gamma1_top-gamma1_bottom)+gamma1_bottom;
    /* k[i] = 1000.0;
    gamma1[i] = 0.005; */
  //  printf("gamma=%f, k=%f\n",gamma1[i],k[i]);
  do { k[i] = rand_normal(1000,300); } while( k[i]<0 );
  do { gamma1[i] = rand_normal(0.001,0.0003); } while( gamma1[i]<0 );

  }
  for(i=0;i<N;i++){
    m[i] = 1.0;
    u[i] = 0.0;
    v[i] = 0.0;
    u_d[i] = 0.0;
    v_d[i] = 0.0;
    force_x[i] = 0.0;
    ku1[i] = 0.0;
    ku2[i] = 0.0;
    ku3[i] = 0.0;
    ku4[i] = 0.0;
    kv1[i] = 0.0;
    kv2[i] = 0.0;
    kv3[i] = 0.0;
    kv4[i] = 0.0;
    kx1[i] = 0.0;
    kx2[i] = 0.0;
    kx3[i] = 0.0;
    kx4[i] = 0.0;
    ky1[i] = 0.0;
    ky2[i] = 0.0;
    ky3[i] = 0.0;
    ky4[i] = 0.0;
  }

  for(i=0;i<20;i++){
    input[i] = 0.0;
  }
  for(i=0;i<3;i++){
    o_nrm2[i] = 0.0;
  }
  for(i=0;i<11;i++){
    o_nrm10[i] = 0.0;
  }
  for(i=0;i<21;i++){
    o_nrm20[i] = 0.0;
  }
  for(i=0;i<3;i++){
    o_ms[i] = 0.0;
    squared_err[i] = 0.0;
    normalizer[i] = 0.0;
    err[i] = 0.0;
  }
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      // -1で初期化したところ、l[-1]にアクセスしてエラーが出なかった
      p2l_mat[i][j] = -1;
    }
  }
}

void initMinuteInitialStates(){
  int root_N = sqrt(N);
  int i,j;
  double tmp = 0.0;
  double delta1,delta2;
  for(i=0;i<root_N;i++){
    for(j=0;j<root_N;j++){
    /*  delta1 = (double)rand()/RAND_MAX*0.0000000002 - 0.0000000001;
      delta2 = (double)rand()/RAND_MAX*0.0000000002 - 0.0000000001;
      x_d[root_N*i+j] = j + delta1;
      y_d[root_N*i+j] = i + delta2; */
      x_d[root_N*i+j] = j;
      y_d[root_N*i+j] = i;
      if(i==0&&j==0){
        x_d[root_N*i+j] = j + 0.000000000001;
        y_d[root_N*i+j] = i + 0.000000000001;
      }
    }
  }
  getSpringLength(l_d,x_d,y_d);
    for(i=0;i<M;i++){
        tmp += pow( (l[i]-l_d[i]), 2.0 );
    }
    initial_d = sqrt(tmp);
    norm2_pre = sqrt(tmp);
}

void initFiles(){
  int looked_idx = 0;
  int arrayl_idx = 0;
  int i,j;
  /* make coodinates data files */
  for(i=0;i<N;i++){
    sprintf(filename1,"%s/points/point%d.dat",dirname,i);
    fp1 = fopen(filename1,"w");
    /* error handling */
    if( fp1 == NULL ){
      printf("cannot open file %s\n",filename1);
    }
    else{
      //    printf("open file %s\n",filename1);
    }
    fprintf(fp1,"x[%d] y[%d]\n",i,i);
    fprintf(fp1,"%f %f\n",x[i],y[i]); //coordinates at n=0 (t=0)初期値を書き込み
    fclose(fp1);
  }

  /* make length of springs data files and set p2l_mat[i][j] */
  for(i=0;i<N;i++){
    for(j=looked_idx;j<N;j++){
      if(G[i][j] == 1){

        sprintf(filename2,"%s/springs/length%d.dat",dirname,arrayl_idx);
        fp2 = fopen(filename2,"w");
        /* error handling */
        if( fp2 == NULL ){
          printf("cannot open file %s\n",filename2);
        }
        else{
          //    printf("open file %s\n",filename2);
        }
        fprintf(fp2,"# l[%d]( l_%d(t) ) connection: point%d-point%d\n",arrayl_idx,arrayl_idx,i,j);
        fprintf(fp2,"n t l_%d(t)\n",arrayl_idx);
        fprintf(fp2,"0 0 %f\n",l[arrayl_idx]);
        fclose(fp2);
        /* 質点idx1とidx2と、バネのインデックスarrayl_idxを対応づける行列 */
        /* p2l_matを対称行列にする */
        p2l_mat[i][j] = arrayl_idx;
        p2l_mat[j][i] = arrayl_idx;
        arrayl_idx++;
      }
    }
    looked_idx++;
  }

  /* make outputs.dat file */
  sprintf(filename3,"%s/results/outputs.dat",dirname);
  fp3 = fopen(filename3,"w");
  /* error handling */
  if( fp3 == NULL ){
    printf("cannot open file %s\n",filename3);
  }
  else{
    //    printf("open file %s\n",filename3);
  }
  fprintf(fp3,"#n  ms_narma2 ms_narma10  ms_narma20  narma2  narma10  narma20\n");
  fclose(fp3);

  /* make le.dat file */
  sprintf(filename5,"%s/results/le.dat",dirname);
  fp5 = fopen(filename5,"w");
  /* error handling */
  if( fp5 == NULL ){
    printf("cannot open file %s\n",filename5);
  }
  else{
    //    printf("open file %s\n",filename3);
  }
  fprintf(fp5,"#n  lyapunov_exponential\n");
  fclose(fp5);
}

//隣接行列を生成
void genGraph(){
  int root_N = sqrt(N);
  int i=0;
  int j=0;
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      G[i][j] = 0;
    }
  }
  for(i=0;i<N;i++){
    if(i-root_N >= 0)
    G[i][i-root_N] = 1;
    if(i % root_N != root_N-1)
    G[i][i+1] = 1;
    if(i % root_N != 0)
    G[i][i-1] = 1;
    if(i+root_N <= N-1)
    G[i][i+root_N] = 1;
  }
}

//隣接行列の要素出力
void printGraph(){
  int i=0, j=0;
  printf("adjacency matrix G =\n");
  //列のindexを表示
  printf("   ");
  for(j=0;j<N;j++){
    printf("%d ",j);
  }
  printf("\n");

  printf("   ");
  for(j=0;j<N;j++){
    printf("－");
  }
  printf("\n");

  for(i=0;i<N;i++){
    //行のindexを表示
    printf("%d | ",i);
    //行毎に要素を出力
    for(j=0;j<N;j++){
      printf("%d ",G[i][j]);
    }
    printf("\n");
  }
  //indexなし
  printf("with no index\n");
  for(i=0;i<N;i++){
    //行毎に要素を出力
    for(j=0;j<N;j++){
      printf("%d, ",G[i][j]);
    }
    printf("\n");
  }

  // p2l_matを出力
  printf("p2l_mat =\n");
  //列のindexを表示
  printf("   ");
  for(j=0;j<N;j++){
    printf("%d ",j);
  }
  printf("\n");

  printf("   ");
  for(j=0;j<N;j++){
    printf("－");
  }
  printf("\n");

  for(i=0;i<N;i++){
    //行のindexを表示
    printf("%d | ",i);
    //行毎に要素を出力
    for(j=0;j<N;j++){
      printf("%d ",p2l_mat[i][j]);
    }
    printf("\n");
  }
}
// ファイルへの座標書き込み
void exportCoordinates(int time_steps){
  int i=0;
  for(i=0;i<N;i++){
    sprintf(filename1,"%s/points/point%d.dat",dirname,i);
    fp1 = fopen(filename1,"a");
    if( fp1 == NULL ){
      printf("loop count n=%d : cannot open file %s\n",time_steps,filename1);
    }
    //  printf("%f %f\n",x[0],y[0]);
    fprintf(fp1,"%f %f\n",x[i],y[i]);
    fclose(fp1);
  }
}

void exportLength(int time_steps){
  int looked_idx = 0;
  int arrayl_idx = 0;
  int i=0,j=0;
  double real_time=time_steps*dt;
  for(i=0;i<N;i++){
    for(j=looked_idx;j<N;j++){
      if(G[i][j] == 1){
        sprintf(filename2,"%s/springs/length%d.dat",dirname,arrayl_idx);
        fp2 = fopen(filename2,"a");
        if( fp2 == NULL ){
          printf("cannot open file %s\n",filename2);
        }
        fprintf(fp2,"%d %lf %.12lf %.12lf\n",time_steps,real_time,l[arrayl_idx], l_d[arrayl_idx]);
        fclose(fp2);
        arrayl_idx++;
      }
    }
    looked_idx++;
  }
}

void updateInput(int time_steps){
  int i=0;
  int idx=0;
  double t = time_steps*dt/T_input;
  for(i=0;i<19;i++){
    input[19-i] = input[18-i];
  }
  input[0] = 0.2*sin(2.0*M_PI*2.11*t)*sin(2.0*M_PI*3.73*t)*sin(2.0*M_PI*4.33*t);

  for(i=0;i<in_num;i++){
    idx = in_p[i];
    //force_x[idx] = w_in[i]*sin(2.0*M_PI*t); // single sinwave
    force_x[idx] = w_in[i]*input[0];
  }
}

void updateNarma2(){
  // o_nrm2 = {o(t), o(t-1), o(t-2)}
  o_nrm2[2] = o_nrm2[1];
  o_nrm2[1] = o_nrm2[0];
  o_nrm2[0] = 0.4 * o_nrm2[1] + 0.4 * o_nrm2[1] * o_nrm2[2] + 0.6 * pow(input[0], 3.0) + 0.1;
}

void updateNarma10(){
  int i=0;
  double tmp=0;
  for(i=0;i<10;i++){
    o_nrm10[10-i] = o_nrm10[9-i];
    tmp += o_nrm10[10-i];
  }
  o_nrm10[0] = 0.3*o_nrm10[1] + 0.05*o_nrm10[1]*tmp + 1.5*input[9]*input[0] + 0.1;
}

void updateNarma20(){
  int i=0;
  double tmp=0;
  for(i=0;i<20;i++){
    o_nrm20[20-i] = o_nrm20[19-i];
    tmp += o_nrm20[20-i];
  }
  o_nrm20[0] = 0.3*o_nrm20[1] + 0.05*o_nrm20[1]*tmp + 1.5*input[19]*input[0] + 0.1;
}

void getSpringLength(double *array_l, double *array_x, double *array_y){
  int looked_idx = 0;
  int arrayl_idx = 0;
  int i,j=0;
  for(i=0;i<N;i++){
    for(j=looked_idx;j<N;j++){
      if(G[i][j] == 1){
        array_l[arrayl_idx] = sqrt( pow( (array_x[i]-array_x[j]), 2.0 ) + pow( (array_y[i]-array_y[j]), 2.0 ) );
        arrayl_idx++;
      }
    }
    looked_idx++;
  }
}

void updateLearnigData(int time_steps){
  int i=0;
  T[time_steps-WASHOUT] = o_nrm2[0];
  T[time_steps-WASHOUT+LEARNING] = o_nrm10[0];
  T[time_steps-WASHOUT+LEARNING*2] = o_nrm20[0];
  //// Bus error ////
  for(i=0;i<M;i++){
    L[time_steps-WASHOUT+LEARNING*i] = l[i];
  }
}

void updateOutputsMS(){
  int ldW_out = M; //W_outの行数
  int incl = 1;
  int inco_ms = 1;
  // 行列W_outの行数, 列数（転置前）
  cblas_dgemv(CblasColMajor,CblasTrans,M,3,
    1.0, W_out, ldW_out, l, incl,
      0.0, o_ms, inco_ms);
}

void getWeights(){
  lapack_int nrhs, ldL, ldT, info;
  int i,j;
  nrhs = 3;
  ldL = LEARNING;
  ldT = LEARNING;
  info = LAPACKE_dgels(LAPACK_COL_MAJOR,'N',LEARNING,M,nrhs,L,ldL,T,ldT);

  for(j=0;j<3;j++){
    for(i=0;i<M;i++){
      W_out[i+M*j] = T[i+ldT*j];
    }
  }
}

void printWeights(){
  int i,j;
  printf("W_out:\n");
  printf("NARMA2  NARMA10 NARMA20\n");
  for(i=0;i<M;i++){
    for(j=0;j<3;j++){
      printf("%lf  ",W_out[i+M*j]);
    }
    printf("\n");
  }
}

void updateErr(){
  squared_err[0] += pow(o_nrm2[0]-o_ms[0], 2.0);
  squared_err[1] += pow(o_nrm10[0]-o_ms[1], 2.0);
  squared_err[2] += pow(o_nrm20[0]-o_ms[2], 2.0);
  normalizer[0] += pow(o_nrm2[0],2.0);
  normalizer[1] += pow(o_nrm10[0],2.0);
  normalizer[2] += pow(o_nrm20[0],2.0);
}

void getErr(){
  int i=0;
  printf("NMSE:\nNARMA2    10    20\n");
  for(i=0;i<3;i++){
    err[i] = squared_err[i] / normalizer[i];
    printf("%f    ",err[i]);
  }
  printf("\n");
}

/* 2つの出力の差　引数(i,l,l_d) */
void updateLyapunovExponent(int time_steps, double *array1, double *array2){
  int i=0;
  double tmp=0;
  double real_time = dt*time_steps;
  for(i=0;i<M;i++){
      tmp = tmp + pow( (array1[i]-array2[i]), 2.0 );
  }
  norm2 = sqrt(tmp);
  //printf("initial_d=%.12lf    norm2=%.12lf\n",initial_d, norm2);
  // overflowを防ぐ
  //for(i=0;i<M;i++){
  //  array2[i] = array1[i] + (initial_d/norm2)*(array2[i]-array1[i]);
  //}
  sum_log = sum_log + log(norm2/initial_d);
//  printf("%d %.12lf %.12lf\n",time_steps,norm2,initial_d);
  lyapunov = sum_log/time_steps;
  // sum_log = sum_log + log(norm2/norm2_pre) / log(2.0);
//  lyapunov = sum_log/real_time;
//  norm2_pre = norm2;
}

void exportOutputs(int time_steps){
  int i;
/*  printf("outputs of MS for\n");
  printf("NARMA2   NARMA10   NARMA20:\n"); */
  sprintf(filename3,"%s/results/outputs.dat",dirname);
  fp3 = fopen(filename3,"a");
  if( fp3 == NULL ){
    printf("loop count n=%d : cannot open file %s\n",time_steps,filename3);
  }
  // time_step, 質点バネ系による近似, NARMAモデルの出力値を1行に出力
  fprintf(fp3,"%d  ",time_steps);
  for(i=0;i<3;i++){
    fprintf(fp3,"%f  ",o_ms[i]);
  }
  fprintf(fp3,"%f  %f  %f",o_nrm2[0],o_nrm10[0],o_nrm20[0]);
  fprintf(fp3,"\n");
  fclose(fp3);
}

void exportLyapunovExponent(int time_steps){
  sprintf(filename5,"%s/results/le.dat",dirname);
  fp5 = fopen(filename5,"a");
  if( fp5 == NULL ){
    printf("loop count n=%d : cannot open file %s\n",time_steps,filename5);
  }
  fprintf(fp5,"%d  %f\n",time_steps, lyapunov);
  fclose(fp5);
}

void exportResults(){
  int i;
  /* make results.dat file (parameters and results) */
  sprintf(filename4,"%s/results/results.txt",dirname);
  fp4 = fopen(filename4,"w");
  /* error handling */
  if( fp4 == NULL ){
    printf("cannot open file %s\n",filename4);
  }
  else{
    //    printf("open file %s\n",filename4);
  }
  /* erros */
  fprintf(fp4,"NMSE:\nNARMA2    10    20\n");
  for(i=0;i<3;i++){
    fprintf(fp4,"%f    ",err[i]);
  }

  fprintf(fp4,"\n---parameters---\n");
  fprintf(fp4,"WASHOUT = %d, LEARNING = %d, EVAL = %d\n",WASHOUT,LEARNING,EVAL);
//  fclose(fp4);
//  fp4 = fopen(filename4,"a");
  fprintf(fp4,"N = %d, M = %d\n",N,M);
  fprintf(fp4,"dt = %f\nT_input = %d\n",dt,T_input);
  fprintf(fp4,"natural_length = %f\n",natu_l);

  fprintf(fp4,"fixed_points_index: ");
  for(i=0;i<fixed_num;i++){
    fprintf(fp4,"%d ",fixed_p[i]);
  }

  fprintf(fp4,"\ninput_points_index: ");
  for(i=0;i<in_num;i++){
    fprintf(fp4,"%d ",in_p[i]);
  }

  fprintf(fp4,"\ninput_weights: ");
  for(i=0;i<in_num;i++){
    fprintf(fp4,"%f ",w_in[i]);
  }

  fprintf(fp4,"\nrange_of_k: %f~%f",k_bottom,k_top);
  fprintf(fp4,"\nrange_of_gamma: %f~%f",gamma1_bottom,gamma1_top);

  fprintf(fp4,"\n---index_of_l[i]  k  gamma---");
  for(i=0;i<M;i++){
    fprintf(fp4,"\n%d   %f   %f",i,k[i],gamma1[i]);
  }
  fprintf(fp4,"\n");
    fclose(fp4);
}

void rk4(int n){
  int fixed_flag = 0; // 0:false 1:true
  int i=0;
  int j=0;
  double tmp_u[N];
  double tmp_x[N];
  double tmp_v[N];
  double tmp_y[N];
  /* skip calculation of fixed points */
  for(i=0;i<N;i++){
    tmp_u[i] = 0.0;
    tmp_x[i] = 0.0;
    tmp_v[i] = 0.0;
    tmp_y[i] = 0.0;
  }
  /* update k1 vectors */
  for(i=0;i<N;i++){
    for(j=0;j<fixed_num;j++){
      if(i==fixed_p[j]){ fixed_flag = 1; break; }
    }
    if(fixed_flag == 0){
      ku1[i] = dt*Fx(x,y,u,i);
      kx1[i] = dt*u[i];
      kv1[i] = dt*Fy(y,x,v,i);
      ky1[i] = dt*v[i];
    }
    fixed_flag = 0;
  }
  /* update k2 vectors */
  for(i=0;i<N;i++){
    tmp_u[i] = u[i] + ku1[i]/2;
    tmp_x[i] = x[i] + kx1[i]/2;
    tmp_v[i] = v[i] + kv1[i]/2;
    tmp_y[i] = y[i] + ky1[i]/2;
  }
  for(i=0;i<N;i++){
    for(j=0;j<fixed_num;j++){
      if(i==fixed_p[j]){ fixed_flag = 1; break; }
    }
    if(fixed_flag == 0){
      ku2[i] = dt*Fx(tmp_x,tmp_y,tmp_u,i);
      kx2[i] = dt*tmp_u[i];
      kv2[i] = dt*Fy(tmp_y,tmp_x,tmp_v,i);
      ky2[i] = dt*tmp_v[i];
    }
    fixed_flag = 0;
  }
  /* update k3 vectors */
  for(i=0;i<N;i++){
    tmp_u[i] = u[i] + ku2[i]/2;
    tmp_x[i] = x[i] + kx2[i]/2;
    tmp_v[i] = v[i] + kv2[i]/2;
    tmp_y[i] = y[i] + ky2[i]/2;
  }
  for(i=0;i<N;i++){
    for(j=0;j<fixed_num;j++){
      if(i==fixed_p[j]){ fixed_flag = 1; break; }
    }
    if(fixed_flag == 0){
      ku3[i] = dt*Fx(tmp_x,tmp_y,tmp_u,i);
      kx3[i] = dt*tmp_u[i];
      kv3[i] = dt*Fy(tmp_y,tmp_x,tmp_v,i);
      ky3[i] = dt*tmp_v[i];
    }
    fixed_flag = 0;
  }
  /* update k4 vectors */
  for(i=0;i<N;i++){
    tmp_u[i] = u[i] + ku3[i];
    tmp_x[i] = x[i] + kx3[i];
    tmp_v[i] = v[i] + kv3[i];
    tmp_y[i] = y[i] + ky3[i];
  }
  for(i=0;i<N;i++){
    for(j=0;j<fixed_num;j++){
      if(i==fixed_p[j]){ fixed_flag = 1; break; }
    }
    if(fixed_flag == 0){
      ku4[i] = dt*Fx(tmp_x,tmp_y,tmp_u,i);
      kx4[i] = dt*tmp_u[i];
      kv4[i] = dt*Fy(tmp_y,tmp_x,tmp_v,i);
      ky4[i] = dt*tmp_v[i];
    }
    fixed_flag = 0;
  }
  /* update u,x,v,y ( time step n -> (n+1) ) */
  for(i=0;i<N;i++){
    for(j=0;j<fixed_num;j++){
      if(i==fixed_p[j]){ fixed_flag = 1; break; }
    }
    if(fixed_flag == 0){
      u[i] = u[i] + ( ku1[i] + 2*ku2[i] + 2*ku3[i] + ku4[i] )/6;
      x[i] = x[i] + ( kx1[i] + 2*kx2[i] + 2*kx3[i] + kx4[i] )/6;
      v[i] = v[i] + ( kv1[i] + 2*kv2[i] + 2*kv3[i] + kv4[i] )/6;
      y[i] = y[i] + ( ky1[i] + 2*ky2[i] + 2*ky3[i] + ky4[i] )/6;
    }
    fixed_flag = 0;
  }
}

void rk4MinuteInitialStates(int n){
  int fixed_flag = 0; // 0:false 1:true
  int i=0;
  int j=0;
  double tmp_u[N];
  double tmp_x[N];
  double tmp_v[N];
  double tmp_y[N];
  /* skip calculation of fixed points */
  for(i=0;i<N;i++){
    tmp_u[i] = 0.0;
    tmp_x[i] = 0.0;
    tmp_v[i] = 0.0;
    tmp_y[i] = 0.0;
  }
  /* update k1 vectors */
  for(i=0;i<N;i++){
    for(j=0;j<fixed_num;j++){
      if(i==fixed_p[j]){ fixed_flag = 1; break; }
    }
    if(fixed_flag == 0){
      ku1_d[i] = dt*Fx(x_d,y_d,u_d,i);
      kx1_d[i] = dt*u_d[i];
      kv1_d[i] = dt*Fy(y_d,x_d,v_d,i);
      ky1_d[i] = dt*v_d[i];
    }
    fixed_flag = 0;
  }
  /* update k2 vectors */
  for(i=0;i<N;i++){
    tmp_u[i] = u_d[i] + ku1_d[i]/2;
    tmp_x[i] = x_d[i] + kx1_d[i]/2;
    tmp_v[i] = v_d[i] + kv1_d[i]/2;
    tmp_y[i] = y_d[i] + ky1_d[i]/2;
  }
  for(i=0;i<N;i++){
    for(j=0;j<fixed_num;j++){
      if(i==fixed_p[j]){ fixed_flag = 1; break; }
    }
    if(fixed_flag == 0){
      ku2_d[i] = dt*Fx(tmp_x,tmp_y,tmp_u,i);
      kx2_d[i] = dt*tmp_u[i];
      kv2_d[i] = dt*Fy(tmp_y,tmp_x,tmp_v,i);
      ky2_d[i] = dt*tmp_v[i];
    }
    fixed_flag = 0;
  }
  /* update k3 vectors */
  for(i=0;i<N;i++){
    tmp_u[i] = u_d[i] + ku2_d[i]/2;
    tmp_x[i] = x_d[i] + kx2_d[i]/2;
    tmp_v[i] = v_d[i] + kv2_d[i]/2;
    tmp_y[i] = y_d[i] + ky2_d[i]/2;
  }
  for(i=0;i<N;i++){
    for(j=0;j<fixed_num;j++){
      if(i==fixed_p[j]){ fixed_flag = 1; break; }
    }
    if(fixed_flag == 0){
      ku3_d[i] = dt*Fx(tmp_x,tmp_y,tmp_u,i);
      kx3_d[i] = dt*tmp_u[i];
      kv3_d[i] = dt*Fy(tmp_y,tmp_x,tmp_v,i);
      ky3_d[i] = dt*tmp_v[i];
    }
    fixed_flag = 0;
  }
  /* update k4 vectors */
  for(i=0;i<N;i++){
    tmp_u[i] = u_d[i] + ku3_d[i];
    tmp_x[i] = x_d[i] + kx3_d[i];
    tmp_v[i] = v_d[i] + kv3_d[i];
    tmp_y[i] = y_d[i] + ky3_d[i];
  }
  for(i=0;i<N;i++){
    for(j=0;j<fixed_num;j++){
      if(i==fixed_p[j]){ fixed_flag = 1; break; }
    }
    if(fixed_flag == 0){
      ku4_d[i] = dt*Fx(tmp_x,tmp_y,tmp_u,i);
      kx4_d[i] = dt*tmp_u[i];
      kv4_d[i] = dt*Fy(tmp_y,tmp_x,tmp_v,i);
      ky4_d[i] = dt*tmp_v[i];
    }
    fixed_flag = 0;
  }
  /* update u,x,v,y ( time step n -> (n+1) ) */
  for(i=0;i<N;i++){
    for(j=0;j<fixed_num;j++){
      if(i==fixed_p[j]){ fixed_flag = 1; break; }
    }
    if(fixed_flag == 0){
      u_d[i] = u_d[i] + ( ku1_d[i] + 2*ku2_d[i] + 2*ku3_d[i] + ku4_d[i] )/6;
      x_d[i] = x_d[i] + ( kx1_d[i] + 2*kx2_d[i] + 2*kx3_d[i] + kx4_d[i] )/6;
      v_d[i] = v_d[i] + ( kv1_d[i] + 2*kv2_d[i] + 2*kv3_d[i] + kv4_d[i] )/6;
      y_d[i] = y_d[i] + ( ky1_d[i] + 2*ky2_d[i] + 2*ky3_d[i] + ky4_d[i] )/6;
    }
    fixed_flag = 0;
  }
}

double Fx(double *array1, double *array2, double *array3, int idx1){
  int idx2=0;
  double ans = 0;
  /* get sum of elastic forces to array1[idx1] */
  double sum = 0;
  int arrayl_idx=-1;
  int i,j;
//  for(i=0;i<N;i++){
  //  for(j=0;j<N;j++){
  //  printf("%d ",p2l_mat[i][j]);
  //}
//}
//printf("\n");
  for(idx2=0;idx2<N;idx2++){
    if(G[idx1][idx2] == 1){
    // if(arrayl_idx==-1) {printf("idx1=%d, idx2=%d\n",idx1,idx2);}
      arrayl_idx = p2l_mat[idx1][idx2];
      sum += f(array1,array2,idx1,idx2);
      //if(idx1==0&&idx2==1) printf("%.12lf\n",sum);
    }
  }
  //printf("%d  %.12lf\n",arrayl_idx, gamma1[arrayl_idx]);
  // rk4から入力された質点のインデックスidx1に、繋がるばねがあるならarrayl_idx!=-1になっている
  if(arrayl_idx!=-1){
  ans = ( sum - force_x[idx1] - gamma1[arrayl_idx] * array3[idx1] ) / m[idx1];
}
//idx1につながっているばねが1つも存在しない
else {
  ans = 0.0;
}
  return ans;
}
// force_xをx方向にのみ入力するために分ける
double Fy(double *array1, double *array2, double *array3, int idx1){
  int idx2=0;
  double ans = 0;
  /* get sum of elastic forces to array1[idx1] */
  double sum = 0;
  int arrayl_idx=-1;
  for(idx2=0;idx2<N;idx2++){
    if(G[idx1][idx2] == 1){
      arrayl_idx = p2l_mat[idx1][idx2];
      sum += f(array1,array2,idx1,idx2);
    }
  }
  if(arrayl_idx!=-1){
  ans = ( sum - force_x[idx1] - gamma1[arrayl_idx] * array3[idx1] ) / m[idx1];
}
else {
  ans = 0.0;
}
  return ans;
}

double f(double *array1, double *array2, int idx1, int idx2){
  double ans = 0;
  int arrayl_idx = p2l_mat[idx1][idx2];
  // debug
   if(arrayl_idx == -1){
      printf("arrayl_idx = -1: idx1=%d, idx2=%d\n",idx1,idx2);
  }
  double distance = sqrt( pow((array1[idx1]-array1[idx2]), 2.0) + pow((array2[idx1]-array2[idx2]), 2.0) );
//  if(idx1==0&&idx2==1) printf("%.12lf\n",k[arrayl_idx]);
  ans = k[arrayl_idx] * ( natu_l* (array1[idx1]-array1[idx2]) / distance + array1[idx2] - array1[idx1]);
//  if(idx1==1&&idx2==0) printf("%.12lf\n",ans);
//printf("%d  %d\n",idx1,idx2);

  return ans;
}

void eular(){
  int fixed_flag = 0; // 0:false 1:true
  int i=0;
  int j=0;
  double tmp_u[N];
  double tmp_x[N];
  double tmp_v[N];
  double tmp_y[N];
  /* skip calculation of fixed points */
  for(i=0;i<N;i++){
    tmp_u[i] = 0.0;
    tmp_x[i] = 0.0;
    tmp_v[i] = 0.0;
    tmp_y[i] = 0.0;
  }
  /* update u,x,v,y ( time step n -> (n+1) ) */
  for(i=0;i<N;i++){
    for(j=0;j<fixed_num;j++){
      if(i==fixed_p[j]){ fixed_flag = 1; break; }
    }
    if(fixed_flag == 0){
      tmp_u[i] = u[i] + dt*Fx(x,y,u,i);
      tmp_x[i] = x[i] + dt*tmp_u[i];
      tmp_v[i] = v[i] + dt*Fy(y,x,v,i);
      tmp_y[i] = y[i] + dt*tmp_v[i];
    }
    fixed_flag = 0;
  }
  for(i=0;i<N;i++){
    for(j=0;j<fixed_num;j++){
      if(i==fixed_p[j]){ fixed_flag = 1; break; }
    }
    if(fixed_flag == 0){
      u[i] = tmp_u[i];
      x[i] = tmp_x[i];
      v[i] = tmp_v[i];
      y[i] = tmp_y[i];
    }
    fixed_flag = 0;
  }
}

void dammy_genGraph(){
  int root_N = sqrt(N);
  int i=0;
  int j=0;
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      G[i][j] = 0;
    }
  }
  G[0][1]=1;
  G[1][0]=1;
}

void test_updateNarma(int time_steps){
  printf("%d %f %f %f %f\n",time_steps,input[0],o_nrm2[0],o_nrm10[0],o_nrm20[0]);
}

void test_getSpringLength(){
  int looked_idx = 0;
  int arrayl_idx = 0;
  int i,j=0;
  for(i=0;i<N;i++){
    for(j=looked_idx;j<N;j++){
      if(G[i][j] == 1){
        //  l[arrayl_idx] = sqrt( pow( (x[i]-x[j]), 2.0 ) + pow( (y[i]-y[j]), 2.0 ) );
        printf("l[%d] connection x[i]-x[j]: %d-%d\n",arrayl_idx,i,j);
        arrayl_idx++;
      }
    }
    looked_idx++;
  }
}

void test_updateLearningData(){
  /*
  plot "T.dat" using ($1*0.0025):2 w l title"narma2", "T.dat" using ($1*0.0025):3 w l title"narma10", "T.dat" using ($1*0.0025):4 w l title"narma20"
  */
  int i,j;
  // export T
  /*  for(j=0;j<LEARNING;j++){
  printf("%d %f %f %f\n",j+WASHOUT,T[j],T[LEARNING+j],T[LEARNING*2+j]);
} */
// export L
for(j=0;j<LEARNING;j++){
  printf("%d ",j+WASHOUT);
  for(i=0;i<M;i++){
    printf("%f ",L[LEARNING*i+j]);
  }
  printf("\n");
}
}
