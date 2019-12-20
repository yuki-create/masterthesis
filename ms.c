// 実行は ./ms N NSTEP
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "lapacke.h"
#include "cblas.h"
// decralation of functions
double rand_normal(double mu, double sigma);
double Uniform( void );
void genGraph(int size, int graph[size][size]);
void printGraph(int size, int graph[size][size], int table[size][size]);
void initArray1dim(int size, double array[size],double initial_value);
void initCoordinates(double natu_l, int size, double array_x[size], double array_y[size], double array_x_d[size], double array_y_d[size], double delta, int d_idx, int fixed_num, int fixed_p[fixed_num]);
void initSpringParameters(int size, double mu_k, double sigma_k,double mu_g, double sigma_g, double array_k[size], double array_g[size]);
void initP2lMat(int size, int graph[size][size], int table[size][size]);
void updateInput(int time_steps, double dt, int T_input, double input[30], int in_num, int in_p[in_num],int N, double force_x[N],double w_in[in_num], double x[N]);
void updateNarma(double input[30], double o_nrm2[3], double o_nrm10[11], double o_nrm20[21], double o_nrm30[31]);
void rk4(int N, double dt, double x[N], double y[N], double u[N], double v[N], int fixed_num, int fixed_p[fixed_num], int G[N][N], int p2l_mat[N][N],double force_x[N], double m[N], int M, double k[M], double gamma1[M], double natu_l);
double Fx(double *array1, double *array2, double *array3, int idx1, int N, int G[N][N], int p2l_mat[N][N],double force_x[N], double m[N], int M, double k[M], double gamma1[M], double natu_l );
double Fy(double *array1, double *array2, double *array3, int idx1, int N, int G[N][N], int p2l_mat[N][N],double force_x[N], double m[N], int M, double k[M], double gamma1[M], double natu_l );
double f(double *array1, double *array2, int idx1, int idx2, int N, int p2l_mat[N][N], int M, double k[M], double natu_l);
void getSpringLength(double *array_l, double *array_x, double *array_y, int N, int G[N][N]);
double updateLyapunovExponent(int time_steps,int M, double array1[M+1], double array2[M], double initial_d, double *sum_log);
void exportLyapunovExponent(int time_steps, double dt, double lyapunov, char *dirname, FILE *fp5, char filename5[40]);
void updateLearnigData(int time_steps, int M, int WASHOUT, int LEARNING, double l[M+1],double o_nrm2[3], double o_nrm10[11], double o_nrm20[21], double o_nrm30[31], double *T, double *L);
void getWeights(int M, int LEARNING, double *L, double *T, double *W_out);
void printWeights(int M, double *W_out);
void updateOutputsMS(int M, double o_ms[4], double l[M+1], double *W_out);
void exportCoordinates(int time_steps,char *dirname, FILE *fp1,  char *filename1, int N, double x[N], double y[N]);
void exportLength(int time_steps, double dt, char *dirname, FILE *fp2, char *filename2, int N, int G[N][N], int M, double l[M+1], double l_d[M]);
void exportOutputs(int time_steps, double dt, char *dirname, FILE *fp3, char *filename3, double o_ms[4], double o_nrm2[3], double o_nrm10[11], double o_nrm20[21], double o_nrm30[31]);
void exportInput(int time_steps, double dt, char *dirname, FILE *fp7, char *filename7, double input[30]);
void exportGraph(char *dirname, int size, int graph[size][size]);

int main(int argc, char *argv[]){
  //////// parameters to be costomized ////////
  int seed_flag = 1;
  int debug_flag = atoi(argv[1]);
  // const int repeat = atoi(argv[1]);
  const int N = pow( atoi(argv[2]), 2.0); // number of mass points
  const int M = (sqrt(N)-1)*sqrt(N)*2; // number of springs (root_N-1)*root_N*2
// washout, learning, evaluating term (time steps)
  const int WASHOUT = 1000;
  const int LEARNING = 5000;
  const int EVAL = 5000;
  // simulating time steps
  const int NSTEP = WASHOUT+LEARNING+EVAL;

  const double dt = 0.0025;
  const int T_input = 1; // adjust frequency of input signal
  const double natu_l = 1.0;
  double w_in[] = {10.0};
  int fixed_p[] = {}; // index array of fixed points
  int in_p[] = {0}; // index array of input points
  int d_idx = 0; // 初期値をずらす質点のindex
  double delta=0.000000000001; // x座標とy座標のずらす量

  //////// variables dosen't need to be costomized ////////
  double k[M];
  double gamma1[M];
  double mu_k = strtod(argv[3],NULL);
  double sigma_k = strtod(argv[4],NULL);
  double mu_g = strtod(argv[5],NULL);
  double sigma_g = strtod(argv[6],NULL);
  int fixed_num = 0; // number of fixed points (elements of fixed_p)
  int in_num = 0; // number of inputted points (elements of in_p)
  double input[30]={0.0}; // input timesiries
  double x[N];
  double u[N]; // dx/dt
  double y[N];
  double v[N]; // dy/dt
  double m[N];
  double *l; // length of spring at time n*dt as outpus.
  int G[N][N]; // adjacency matrix (symmetric matrix)
  int p2l_mat[N][N]; // point index x[idx1],x[idx2] -> spring index p2l_mat[idx1][idx2] (upper triangular matrix)
  double force_x[N];

  // variables and parameters for RK4
  //double ku1[N], ku2[N], ku3[N], ku4[N], kx1[N], kx2[N], kx3[N], kx4[N];
  //double kv1[N], kv2[N], kv3[N], kv4[N], ky1[N], ky2[N], ky3[N], ky4[N];

  // variables for mass-spring system with minute initial state
  double x_d[N];
  double u_d[N]; // dx/dt
  double y_d[N];
  double v_d[N]; // dy/dt
  //double ku1_d[N], ku2_d[N], ku3_d[N], ku4_d[N], kx1_d[N], kx2_d[N], kx3_d[N], kx4_d[N];
  //double kv1_d[N], kv2_d[N], kv3_d[N], kv4_d[N], ky1_d[N], ky2_d[N], ky3_d[N], ky4_d[N];
  double *l_d;

  // variables for NARMA models
  double o_nrm2[3]={0.0}; // NARMA2 output
  double o_nrm10[11]={0.0}; // NARMA10 output
  double o_nrm20[21]={0.0}; // NARMA20 output
  double o_nrm30[31]={0.0}; // NARMA30 output

  // variables for learning
  double *T, *W_out, *L;

  // variables for EVAL phase
  double o_ms[4]={0.0}; // outputs of MS [O_ms(t) for narma2, narma10, narma20]
  double squared_err[4]={0.0}; // squared errors [narma2, narma10, narma20]
  double normalizer[4]={0.0}; // normalizers [narma2, narma10, narma20]
  double err[4]={0.0}; // normilized mean squaered errors

  // variables for lyapunov exponent
  double lyapunov = 0.0;
  double initial_d = 0.0; // initial distance of two systems
  double sum_log = 0.0;

  // parameters for file I/O
  char *dirname = argv[7];
  FILE *fp1; // for export coodinates
  FILE *fp2; // for export length of springs
  FILE *fp3; // for outputs of MS and NARMA after learning
  FILE *fp4; // for parameters and results file
  FILE *fp5; // for export lyapunov exponent
  FILE *fp6; // for export many results
  FILE *fp7; // export input signal
  char filename1[40],filename2[40],filename3[40],filename4[40],filename5[40],filename6[40],filename7[40];
  int i,j,n;
  double tmp = 0.0;
  double norm2 = 0.0;
  double tmp_array[M];
  int loop=0;
  int skip=sqrt(N);
  double le_avr = 0.0;
  double err_avr[4] = {0.0};

  // 初期値の誤差を導入する点を変える時、一回だけやればいいもの

// 初期化始まり　init start
  if(seed_flag == 1) srand((unsigned)time(NULL));

  l = malloc((M+1)*sizeof(double));
  l_d = malloc(M*sizeof(double));
  // lapackのライブラリを使う配列の初期化。行列の列ごとに一次元配列に格納。
  T = malloc(LEARNING*4*sizeof(double));
  L = malloc(LEARNING*(M+1)*sizeof(double));
  W_out = malloc((M+1)*4*sizeof(double));

  fixed_num = sizeof fixed_p / sizeof fixed_p[0];
  in_num = sizeof in_p / sizeof in_p[0];
  genGraph(N,G);
//  exportGraph(dirname,N,G);
  initP2lMat(N,G,p2l_mat);
  //  printGraph(N,G,p2l_mat);
  initSpringParameters(M,mu_k,sigma_k,mu_g,sigma_g,k,gamma1);
  initArray1dim(N,m,1.0);

  // lyapunovの計算だけする
  for(loop=0;loop<skip;loop++){
      d_idx=skip*loop;

      initArray1dim(N,u,0.0);
      initArray1dim(N,v,0.0);
      initArray1dim(N,u_d,0.0);
      initArray1dim(N,v_d,0.0);
      initArray1dim(N,force_x,0.0);
      initArray1dim(M+1,l,natu_l);
    //  l[M]=1.0; // l[M]はバイアスで常にl[M]=1.0
      initArray1dim(M,l_d,natu_l);
      initArray1dim(30,input,0.0);

      sum_log=0.0;
      lyapunov=0.0;

      initCoordinates(natu_l, N,x,y,x_d,y_d, delta, d_idx, fixed_num, fixed_p);

      getSpringLength(l,x,y,N,G);
      getSpringLength(l_d,x_d,y_d,N,G);
      initial_d=0.0;
      for(i=0;i<M;i++){
          initial_d += pow( (l[i]-l_d[i]), 2.0 );
      }
      initial_d = sqrt(initial_d);
    //  printf("%14.12lf\n",initial_d);
    // printf("l[0]=%.12lf\n",l[0]);

      for(n=0;n<WASHOUT+LEARNING;n++){
      updateInput(n,dt, T_input, input, in_num, in_p ,N, force_x, w_in,x);
       //printf("%.12lf\n",o_nrm30[0]);
       rk4(N, dt, x, y, u, v, fixed_num, fixed_p, G, p2l_mat, force_x, m, M, k, gamma1, natu_l);
      rk4(N, dt, x_d, y_d, u_d, v_d, fixed_num, fixed_p, G, p2l_mat, force_x, m, M, k, gamma1, natu_l);
       getSpringLength(l,x,y,N,G);
       getSpringLength(l_d,x_d,y_d,N,G);
       lyapunov = updateLyapunovExponent(n,M,l,l_d, initial_d, &sum_log);
     }
     // loop回数(初期値入れた点)　lyapunov errorNARMA2 10 20 30
  //   printf("%d    %.12lf\n",d_idx,lyapunov);
     le_avr += lyapunov;
  }

  // 近似誤差を測る　
    d_idx=0;

  initArray1dim(N,u,0.0);
  initArray1dim(N,v,0.0);
//  initArray1dim(N,u_d,0.0);
//  initArray1dim(N,v_d,0.0);
  initArray1dim(N,force_x,0.0);
  initArray1dim(M+1,l,natu_l);
//  l[M]=1.0; // l[M]はバイアスで常にl[M]=1.0
  // initArray1dim(M,l_d,natu_l);
  initArray1dim(30,input,0.0);
  initArray1dim(3,o_nrm2,0.0);
  initArray1dim(11,o_nrm10,0.0);
  initArray1dim(21,o_nrm20,0.0);
  initArray1dim(31,o_nrm30,0.0);
  initArray1dim(4,o_ms,0.0);
  initArray1dim(4,squared_err,0.0);
  initArray1dim(4,normalizer,0.0);
  initArray1dim(4,err,0.0);

  sum_log=0.0;
  lyapunov=0.0;

  initCoordinates(natu_l, N,x,y,x_d,y_d, delta, d_idx, fixed_num, fixed_p);

  getSpringLength(l,x,y,N,G);
  //getSpringLength(l_d,x_d,y_d,N,G);
//  initial_d=0.0;
//  for(i=0;i<M;i++){
//      initial_d += pow( (l[i]-l_d[i]), 2.0 );
//  }
//  initial_d = sqrt(initial_d);
//  printf("%14.12lf\n",initial_d);
// printf("l[0]=%.12lf\n",l[0]);

  // 考えている行列の、列ごとに要素を代入する
  for(j=0;j<4;j++){
    for(i=0;i<LEARNING;i++){
      T[j*LEARNING+i] = 0.0;
    }
  }
  for(j=0;j<(M+1);j++){
    for(i=0;i<LEARNING;i++){
      L[j*LEARNING+i] = 0.0;
    }
  }
  for(j=0;j<4;j++){
    for(i=0;i<(M+1);i++){
      W_out[j*(M+1)+i] = 0.0;
    }
  }

  // 出力波形を見たりする時、ファイルを初期化
  if(debug_flag == 1){
    int looked_idx = 0;
    int arrayl_idx = 0;
    int i,j;
    // make coodinates data files
    for(i=0;i<N;i++){
      sprintf(filename1,"%s/data/points/point%d.dat",dirname,i);
      fp1 = fopen(filename1,"w");
      if( fp1 == NULL ){
        printf("cannot open file %s\n",filename1);
      }
      else{
        //    printf("open file %s\n",filename1);
      }
      fprintf(fp1,"x[%d] y[%d]\n",i,i);
      fprintf(fp1,"%.12lf %.12lf\n",x[i],y[i]); //coordinates at n=0 (t=0)初期値を書き込み
      fclose(fp1);
    }

    // make length of springs data files
    for(i=0;i<N;i++){
      for(j=looked_idx;j<N;j++){
        if(G[i][j] == 1){
          sprintf(filename2,"%s/data/springs/length%d.dat",dirname,arrayl_idx);
          fp2 = fopen(filename2,"w");
          if( fp2 == NULL ){
            printf("cannot open file %s\n",filename2);
          }
          else{
            //    printf("open file %s\n",filename2);
          }
          fprintf(fp2,"# l[%d]( l_%d(t) ) connection: point%d-point%d\n",arrayl_idx,arrayl_idx,i,j);
          fprintf(fp2,"# n t l_%d(t)      l_differential_%d(t)\n",arrayl_idx, arrayl_idx);
          fclose(fp2);
          arrayl_idx++;
        }
      }
      looked_idx++;
    }
    sprintf(filename2,"%s/data/springs/bias.dat",dirname);
    fp2 = fopen(filename2,"w");
    fprintf(fp2,"# n t bias ( l[%d] )\n",M);
    fclose(fp2);

    // make outputs.dat file
    sprintf(filename3,"%s/data/outputs.dat",dirname);
    fp3 = fopen(filename3,"w");
    if( fp3 == NULL ){
      printf("cannot open file %s\n",filename3);
    }
    else{
      //    printf("open file %s\n",filename3);
    }
    fprintf(fp3,"#n  real_time  ms_narma2  ms_narma10  ms_narma20  ms_narma30  narma2  narma10  narma20  narma30\n");
    fclose(fp3);

    // make le.dat file
    sprintf(filename5,"%s/data/le.dat",dirname);
    fp5 = fopen(filename5,"w");
    if( fp5 == NULL ){
      printf("cannot open file %s\n",filename5);
    }
    else{
      //    printf("open file %s\n",filename3);
    }
    fprintf(fp5,"#n  real_time  lyapunov_exponential\n");
    fclose(fp5);

    // make input.dat
    sprintf(filename7,"%s/data/input.dat",dirname);
    fp7 = fopen(filename7,"w");
    if( fp7 == NULL ){
      printf("cannot open file %s\n",filename7);
    }
    fprintf(fp7,"#n  real_time  input\n");
    fclose(fp7);
  }

  // 初期化終わり　init end
  // washout
  for(n=0;n<WASHOUT;n++){
    updateInput(n,dt, T_input, input, in_num, in_p ,N, force_x, w_in,x);
     updateNarma(input, o_nrm2, o_nrm10, o_nrm20,o_nrm30);
     //printf("%.12lf\n",o_nrm30[0]);
     rk4(N, dt, x, y, u, v, fixed_num, fixed_p, G, p2l_mat, force_x, m, M, k, gamma1, natu_l);
  //  rk4(N, dt, x_d, y_d, u_d, v_d, fixed_num, fixed_p, G, p2l_mat, force_x, m, M, k, gamma1, natu_l);
     getSpringLength(l,x,y,N,G);
    // getSpringLength(l_d,x_d,y_d,N,G);
     //lyapunov = updateLyapunovExponent(n,M,l,l_d, initial_d, &sum_log);
     if(debug_flag==1){
       exportCoordinates(n,dirname,fp1,filename1, N, x, y);
       exportLength(n, dt, dirname, fp2, filename2, N, G, M, l, l_d);
       exportLyapunovExponent(n, dt, lyapunov, dirname, fp5, filename5);
       exportInput(n, dt, dirname, fp7, filename7, input);
     }
  }
  // 学習
  for(n=WASHOUT;n<WASHOUT+LEARNING;n++){
    updateInput(n,dt, T_input, input, in_num, in_p ,N, force_x, w_in,x);
     updateNarma(input, o_nrm2, o_nrm10, o_nrm20,o_nrm30);
     rk4(N, dt, x, y, u, v, fixed_num, fixed_p, G, p2l_mat, force_x, m, M, k, gamma1, natu_l);
  //    rk4(N, dt, x_d, y_d, u_d, v_d, fixed_num, fixed_p, G, p2l_mat, force_x, m, M, k, gamma1, natu_l);
      getSpringLength(l,x,y,N,G);
    //  getSpringLength(l_d,x_d,y_d,N,G);
      updateLearnigData(n, M, WASHOUT, LEARNING, l, o_nrm2, o_nrm10, o_nrm20, o_nrm30, T, L);
    //  lyapunov = updateLyapunovExponent(n,M, l,l_d,initial_d, &sum_log);
      if(debug_flag==1){
         exportCoordinates(n,dirname,fp1,filename1, N, x, y);
        exportLength(n, dt, dirname, fp2, filename2, N, G, M, l, l_d);
        exportLyapunovExponent(n, dt, lyapunov, dirname, fp5, filename5);
        exportInput(n, dt, dirname, fp7, filename7, input);
      }
  }
  getWeights(M, LEARNING, L, T, W_out);
  // printWeights(M, W_out);

  for(n=WASHOUT+LEARNING;n<WASHOUT+LEARNING+EVAL;n++){
    updateInput(n,dt, T_input, input, in_num, in_p ,N, force_x, w_in,x);
     updateNarma(input, o_nrm2, o_nrm10, o_nrm20,o_nrm30);
     rk4(N, dt, x, y, u, v, fixed_num, fixed_p, G, p2l_mat, force_x, m, M, k, gamma1, natu_l);
//      rk4(N, dt, x_d, y_d, u_d, v_d, fixed_num, fixed_p, G, p2l_mat, force_x, m, M, k, gamma1, natu_l);
      getSpringLength(l,x,y,N,G);
  //    getSpringLength(l_d,x_d,y_d,N,G);
      updateOutputsMS(M, o_ms, l, W_out);
    // updateErr
    squared_err[0] += pow(o_nrm2[0]-o_ms[0], 2.0);
    squared_err[1] += pow(o_nrm10[0]-o_ms[1], 2.0);
    squared_err[2] += pow(o_nrm20[0]-o_ms[2], 2.0);
    squared_err[3] += pow(o_nrm30[0]-o_ms[3], 2.0);
    normalizer[0] += pow(o_nrm2[0],2.0);
    normalizer[1] += pow(o_nrm10[0],2.0);
    normalizer[2] += pow(o_nrm20[0],2.0);
    normalizer[3] += pow(o_nrm30[0],2.0);
    //  lyapunov = updateLyapunovExponent(n,M,l,l_d, initial_d, &sum_log);
      if(debug_flag==1){
         exportCoordinates(n,dirname,fp1,filename1, N, x, y);
        exportLength(n, dt, dirname, fp2, filename2, N, G, M, l, l_d);
        exportOutputs(n, dt, dirname, fp3, filename3, o_ms, o_nrm2, o_nrm10, o_nrm20, o_nrm30);
        exportLyapunovExponent(n, dt, lyapunov, dirname, fp5, filename5);
        exportInput(n, dt, dirname, fp7, filename7, input);
      }

  }
  //printf("sum_log=%.12lf\n",sum_log);
  // getErr
  /*printf("lyapunov exponent = %lf\n", lyapunov);
  printf("NMSE:\nNARMA2    10    20   30\n");
  for(i=0;i<4;i++){
    err[i] = squared_err[i] / normalizer[i];
    printf("%.12lf    ",err[i]);
  } */

// result.txtに結果を残す
if(debug_flag==1){
  //2つの系の2normを計算
  // array1(l_d[M])を、tmpにコピー
    cblas_dcopy(M, l, 1, tmp_array, 1);
  // tmp = -l_d+tmp
    cblas_daxpy(M, -1.0, l_d, 1, tmp_array, 1);
  // tmpの2ノルムを求める
    norm2 = cblas_dnrm2(M, tmp_array, 1);
    // make results.dat file (parameters and results)
    sprintf(filename4,"%s/results/results.txt",dirname);
    fp4 = fopen(filename4,"w");
    if( fp4 == NULL ){
      printf("cannot open file %s\n",filename4);
    }
    else{
    fprintf(fp4, "lyapunov exponent = %lf\n", lyapunov);
    fprintf(fp4,"NMSE:\nNARMA2    10    20    30\n");
    for(i=0;i<4;i++){
      fprintf(fp4,"%f    ",err[i]);
    }
    fprintf(fp4,"\n---parameters---\n");
    fprintf(fp4,"WASHOUT = %d, LEARNING = %d, EVAL = %d\n",WASHOUT,LEARNING,EVAL);
  //  fclose(fp4);
  //  fp4 = fopen(filename4,"a");
    fprintf(fp4,"N = %d, M = %d\n",N,M);
    fprintf(fp4,"differential_point_index: %d\n",d_idx);
    fprintf(fp4,"%.12lf=differential_initial_value\n",delta);
    fprintf(fp4,"%.12lf=initial_distance\n",initial_d);
    fprintf(fp4,"%.12lf=2-norm\n",norm2);
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

    fprintf(fp4,"\nk: mu = %f, sigma = %f",mu_k, sigma_k);
    fprintf(fp4,"\ngamma: mu = %f, sigma = %f",mu_g, sigma_g);

    fprintf(fp4,"\n---index_of_l[i]  k  gamma---");
    for(i=0;i<M;i++){
      fprintf(fp4,"\n%d   %f   %f",i,k[i],gamma1[i]);
    }
    fprintf(fp4,"\n");
      fclose(fp4);
    }
  }

  // 初期値を変えたときのlyapunov指数の平均と近似誤差の平均
  le_avr = le_avr/N;
//  for(i=0;i<4;i++){  err_avr[i] = err_avr[i]/N;  }
// loop回数(初期値入れた点)　lyapunov errorNARMA2 10 20 30

  printf("average lyapunov exponent = %.12lf\n", le_avr);
  printf("NMSE:\nNARMA2      10       20      30\n");
  for(i=0;i<4;i++){
    err[i] = squared_err[i] / normalizer[i];
    printf("%.12lf    ",err[i]);
  }
  printf("\n");

  sprintf(filename6,"%s/data.txt",dirname);
  fp6 = fopen(filename6,"a");
  fprintf(fp6,"%lf  %lf  %lf  %lf  %.12lf   %.12lf   %.12lf   %.12lf   %.12lf\n",mu_k, sigma_k, mu_g, sigma_g, le_avr, err[0],err[1],err[2],err[3] );
  fclose(fp6);

  free(l);
  free(l_d);
  free(T);
  free(L);
  free(W_out);
  return (0);
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
    printf("%2.0d ",j);
  }
  printf("\n");

  printf("   ");
  for(j=0;j<size;j++){
    printf("－");
  }
  printf("\n");

  for(i=0;i<size;i++){
    //行のindexを表示
    printf("%2.0d | ",i);
    //行毎に要素を出力
    for(j=0;j<size;j++){
      printf("%2.0d ",table[i][j]);
    }
    printf("\n");
  }
}

void exportGraph(char *dirname, int size, int graph[size][size]){
  int i,j;
  FILE *fp;
  char filename[40];
  sprintf(filename,"%s/results/graph.dat",dirname);
  fp = fopen(filename,"w");
  for(i=0;i<size;i++){
    //行毎に要素を出力
    for(j=0;j<size;j++){
      fprintf(fp,"%d\n",graph[i][j]);
    }
  //  fprintf(fp,"\n");
  }
  fclose(fp);
}

void initArray1dim(int size, double array[size], double initial_value){
  int i;
  for(i=0;i<size;i++){
    array[i] = initial_value;
  }
}
// x_d[0],y_d[0]のみ初期値をずらす
void initCoordinates(double natu_l, int size, double array_x[size], double array_y[size], double array_x_d[size], double array_y_d[size], double delta, int d_idx, int fixed_num, int fixed_p[fixed_num])
{
  // init coordinates of mass points
  int root_size = sqrt(size);
  int i,j;
  for(i=0;i<root_size;i++){
    for(j=0;j<root_size;j++){
  //    if(root_size*i+j == fixed_p[0] || root_size*i+j == fixed_p[1]){
        array_x[root_size*i+j] = j*natu_l;
        array_y[root_size*i+j] = i*natu_l;
        array_x_d[root_size*i+j] = j*natu_l;
        array_y_d[root_size*i+j] = i*natu_l;
    //  }
    //  else{
    //  array_x[root_size*i+j] = j;
    //  array_y[root_size*i+j] = i;
    //  array_x_d[root_size*i+j] = j + (double)rand()/RAND_MAX*delta-delta/2;
    //  array_y_d[root_size*i+j] = i + (double)rand()/RAND_MAX*delta-delta/2;
    //}
    }
  }
    array_x_d[d_idx] = array_x_d[d_idx] + delta;
    array_y_d[d_idx] = array_y_d[d_idx] + delta;
}

void initSpringParameters(int size, double mu_k, double sigma_k,double mu_g, double sigma_g, double array_k[size], double array_g[size]){
  int i;
  for(i=0;i<size;i++){
    do { array_k[i] = rand_normal(mu_k,sigma_k); } while( array_k[i]<0 );
    do { array_g[i] = rand_normal(mu_g,sigma_g); } while( array_g[i]<0 );
  }
}
// 質点idx1とidx2と、バネのインデックスarrayl_idxを対応づける行列p2l_matの初期化
void initP2lMat(int size, int graph[size][size], int table[size][size]){
  int i,j;
  int looked_idx=0;
  int arrayl_idx=0;

  for(i=0;i<size;i++){
    for(j=0;j<size;j++){
      table[i][j] = -1;
    }
  }

  for(i=0;i<size;i++){
    for(j=looked_idx;j<size;j++){
      if(graph[i][j] == 1){
        table[i][j] = arrayl_idx;
        table[j][i] = arrayl_idx;
        arrayl_idx++;
      }
    }
    looked_idx++;
  }
}

void updateInput(int time_steps, double dt, int T_input, double input[30], int in_num, int in_p[in_num],int N, double force_x[N],double w_in[in_num],double x[N]){
  int i=0;
  int idx=0;
  double t = dt*time_steps/T_input;
  for(i=0;i<29;i++){
    input[29-i] = input[28-i];
  }
  input[0] = 0.2*sin(2.0*M_PI*2.11*t)*sin(2.0*M_PI*3.73*t)*sin(2.0*M_PI*4.33*t);
  // -0.5-0.5のランダムな実数
  //input[0] = (double)rand()/RAND_MAX*-0.5;
  for(i=0;i<in_num;i++){
     idx = in_p[i];
     force_x[idx] = w_in[i]*input[0];
  }
}

void updateNarma(double input[30], double o_nrm2[3], double o_nrm10[11], double o_nrm20[21], double o_nrm30[31]){
  int i=0;
  double tmp=0.0;
  // update NARMA2 model
  // o_nrm2 = {o(t), o(t-1), o(t-2)}
  o_nrm2[2] = o_nrm2[1];
  o_nrm2[1] = o_nrm2[0];
  o_nrm2[0] = 0.4 * o_nrm2[1] + 0.4 * o_nrm2[1] * o_nrm2[2] + 0.6 * pow(input[0], 3.0) + 0.1;
  // update NARMA10 model
  for(i=0;i<10;i++){
    o_nrm10[10-i] = o_nrm10[9-i];
    tmp += o_nrm10[10-i];
  }
  o_nrm10[0] = 0.3*o_nrm10[1] + 0.05*o_nrm10[1]*tmp + 1.5*input[9]*input[0] + 0.1;
  //o_nrm10[0] = 0.2*o_nrm10[1] + 0.004*o_nrm10[1]*tmp + 1.5*input[9]*input[0] + 0.001;
  // update NARMA20 model
  tmp = 0.0;
  for(i=0;i<20;i++){
    o_nrm20[20-i] = o_nrm20[19-i];
    tmp += o_nrm20[20-i];
  }
  o_nrm20[0] = 0.3*o_nrm20[1] + 0.05*o_nrm20[1]*tmp + 1.5*input[19]*input[0] + 0.1;
  //o_nrm20[0] = 0.2*o_nrm20[1] + 0.004*o_nrm20[1]*tmp + 1.5*input[19]*input[0] + 0.001;
//  printf("%lf\n",o_nrm20[0]);
tmp = 0.0;
  for(i=0;i<30;i++){
    o_nrm30[30-i] = o_nrm30[29-i];
    tmp += o_nrm30[30-i];
  }
  // overflow
//  o_nrm30[0] = 0.3*o_nrm30[1] + 0.05*o_nrm30[1]*tmp + 1.5*input[29]*input[0] + 0.1;
  o_nrm30[0] = 0.2*o_nrm30[1] + 0.004*o_nrm30[1]*tmp + 1.5*input[29]*input[0] + 0.001;
  //printf("%.12lf\n",o_nrm30[0]);
}

void rk4(int N, double dt, double x[N], double y[N], double u[N], double v[N], int fixed_num, int fixed_p[fixed_num], int G[N][N], int p2l_mat[N][N],double force_x[N], double m[N], int M, double k[M], double gamma1[M], double natu_l){
  int fixed_flag = 0; // 0:false 1:true
  int i=0;
  int j=0;
  double tmp_u[N];
  double tmp_x[N];
  double tmp_v[N];
  double tmp_y[N];
  double ku1[N], ku2[N], ku3[N], ku4[N], kx1[N], kx2[N], kx3[N], kx4[N];
  double kv1[N], kv2[N], kv3[N], kv4[N], ky1[N], ky2[N], ky3[N], ky4[N];

  // skip calculation of fixed points
  for(i=0;i<N;i++){
    tmp_u[i] = 0.0;
    tmp_x[i] = 0.0;
    tmp_v[i] = 0.0;
    tmp_y[i] = 0.0;
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
  // update k1 vectors
  for(i=0;i<N;i++){
    for(j=0;j<fixed_num;j++){ if(i==fixed_p[j]){ fixed_flag = 1; break; }}
    if(fixed_flag == 0){
      ku1[i] = dt*Fx(x,y,u,i, N, G, p2l_mat, force_x, m, M, k, gamma1, natu_l);
      kx1[i] = dt*u[i];
      kv1[i] = dt*Fy(y,x,v,i, N, G, p2l_mat, force_x, m, M, k, gamma1, natu_l);
      ky1[i] = dt*v[i];
    }
    fixed_flag = 0;
  }
  // update k2 vectors
 for(i=0;i<N;i++){
    tmp_u[i] = u[i] + ku1[i]/2;
    tmp_x[i] = x[i] + kx1[i]/2;
    tmp_v[i] = v[i] + kv1[i]/2;
    tmp_y[i] = y[i] + ky1[i]/2;
  }
/*
   // uを、tmp_uにコピー
  cblas_dcopy(N, u, 1, tmp_u, 1);
  cblas_dcopy(N, x, 1, tmp_x, 1);
  cblas_dcopy(N, v, 1, tmp_v, 1);
  cblas_dcopy(N, y, 1, tmp_y, 1);
  // tmp_u <- 0.5*ku1 + tmp_u(=u)
  cblas_daxpy(N, 0.5, ku1, 1, tmp_u, 1);
  cblas_daxpy(N, 0.5, kx1, 1, tmp_x, 1);
  cblas_daxpy(N, 0.5, kv1, 1, tmp_v, 1);
  cblas_daxpy(N, 0.5, ky1, 1, tmp_y, 1);
*/
  for(i=0;i<N;i++){
   for(j=0;j<fixed_num;j++){
      if(i==fixed_p[j]){ fixed_flag = 1; break; }
    }
    if(fixed_flag == 0){
      ku2[i] = dt*Fx(tmp_x,tmp_y,tmp_u,i, N, G, p2l_mat, force_x, m, M, k, gamma1, natu_l);
      kx2[i] = dt*tmp_u[i];
      kv2[i] = dt*Fy(tmp_y,tmp_x,tmp_v,i, N, G, p2l_mat, force_x, m, M, k, gamma1, natu_l);
      ky2[i] = dt*tmp_v[i];
    }
    fixed_flag = 0;
  }
  // update k3 vectors
  for(i=0;i<N;i++){
    tmp_u[i] = u[i] + ku2[i]/2;
    tmp_x[i] = x[i] + kx2[i]/2;
    tmp_v[i] = v[i] + kv2[i]/2;
    tmp_y[i] = y[i] + ky2[i]/2;
  }
/*
  // uを、tmp_uにコピー
 cblas_dcopy(N, u, 1, tmp_u, 1);
 cblas_dcopy(N, x, 1, tmp_x, 1);
 cblas_dcopy(N, v, 1, tmp_v, 1);
 cblas_dcopy(N, y, 1, tmp_y, 1);
 // tmp_u <- 0.5*ku1 + tmp_u(=u)
 cblas_daxpy(N, 0.5, ku2, 1, tmp_u, 1);
 cblas_daxpy(N, 0.5, kx2, 1, tmp_x, 1);
 cblas_daxpy(N, 0.5, kv2, 1, tmp_v, 1);
 cblas_daxpy(N, 0.5, ky2, 1, tmp_y, 1);
*/
  for(i=0;i<N;i++){
    for(j=0;j<fixed_num;j++){
      if(i==fixed_p[j]){ fixed_flag = 1; break; }
    }
    if(fixed_flag == 0){
      ku3[i] = dt*Fx(tmp_x,tmp_y,tmp_u,i, N, G, p2l_mat, force_x, m, M, k, gamma1, natu_l);
      kx3[i] = dt*tmp_u[i];
      kv3[i] = dt*Fy(tmp_y,tmp_x,tmp_v,i, N, G, p2l_mat, force_x, m, M, k, gamma1, natu_l);
      ky3[i] = dt*tmp_v[i];
    }
    fixed_flag = 0;
  }
  // update k4 vectors
  for(i=0;i<N;i++){
    tmp_u[i] = u[i] + ku3[i];
    tmp_x[i] = x[i] + kx3[i];
    tmp_v[i] = v[i] + kv3[i];
    tmp_y[i] = y[i] + ky3[i];
  }
/*
  // uを、tmp_uにコピー
 cblas_dcopy(N, u, 1, tmp_u, 1);
 cblas_dcopy(N, x, 1, tmp_x, 1);
 cblas_dcopy(N, v, 1, tmp_v, 1);
 cblas_dcopy(N, y, 1, tmp_y, 1);
 // tmp_u <- 0.5*ku1 + tmp_u(=u)
 cblas_daxpy(N, 0.5, ku3, 1, tmp_u, 1);
 cblas_daxpy(N, 0.5, kx3, 1, tmp_x, 1);
 cblas_daxpy(N, 0.5, kv3, 1, tmp_v, 1);
 cblas_daxpy(N, 0.5, ky3, 1, tmp_y, 1);
*/
  for(i=0;i<N;i++){
    for(j=0;j<fixed_num;j++){
      if(i==fixed_p[j]){ fixed_flag = 1; break; }
    }
    if(fixed_flag == 0){
      ku4[i] = dt*Fx(tmp_x,tmp_y,tmp_u,i, N, G, p2l_mat, force_x, m, M, k, gamma1, natu_l);
      kx4[i] = dt*tmp_u[i];
      kv4[i] = dt*Fy(tmp_y,tmp_x,tmp_v,i, N, G, p2l_mat, force_x, m, M, k, gamma1, natu_l);
      ky4[i] = dt*tmp_v[i];
    }
    fixed_flag = 0;
  }
  // update u,x,v,y ( time step n -> (n+1) )
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
double Fx(double *array1, double *array2, double *array3, int idx1, int N, int G[N][N], int p2l_mat[N][N],double force_x[N], double m[N], int M, double k[M], double gamma1[M], double natu_l){
  int idx2=0;
  double ans = 0;
  // get sum of elastic forces to array1[idx1]
  double sum = 0;
  int arrayl_idx=-1;
  for(idx2=0;idx2<N;idx2++){
    if(G[idx1][idx2] == 1){
      arrayl_idx = p2l_mat[idx1][idx2];
      sum += f(array1,array2,idx1,idx2,N,p2l_mat,M,k, natu_l);
    }
  }
  if(arrayl_idx!=-1){
    ans = ( sum - force_x[idx1] - gamma1[arrayl_idx] * array3[idx1] ) / m[idx1];
  }
  else { ans = 0.0; }
  return ans;
}

// force_xをx方向にのみ入力するために分ける
double Fy(double *array1, double *array2, double *array3, int idx1, int N, int G[N][N], int p2l_mat[N][N],double force_x[N], double m[N], int M, double k[M], double gamma1[M], double natu_l){
  int idx2=0;
  double ans = 0;
  // get sum of elastic forces to array1[idx1]
  double sum = 0;
  int arrayl_idx = -1;
  for(idx2=0;idx2<N;idx2++){
    if(G[idx1][idx2] == 1){
      arrayl_idx = p2l_mat[idx1][idx2];
      sum += f(array1,array2,idx1,idx2,N,p2l_mat,M,k, natu_l);
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

double f(double *array1, double *array2, int idx1, int idx2, int N, int p2l_mat[N][N], int M, double k[M], double natu_l){
  double ans = 0;
  int arrayl_idx = p2l_mat[idx1][idx2];
  double distance = sqrt( pow((array1[idx1]-array1[idx2]), 2.0) + pow((array2[idx1]-array2[idx2]), 2.0) );
  ans = k[arrayl_idx] * ( natu_l* (array1[idx1]-array1[idx2]) / distance + array1[idx2] - array1[idx1]);
  return ans;
}

void getSpringLength(double *array_l, double *array_x, double *array_y, int N, int G[N][N]){
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

// 2つの出力の差　引数(i,l,l_d)
double updateLyapunovExponent(int time_steps,int M, double array1[M+1], double array2[M], double initial_d, double *sum_log){
  int i=0;
  double tmp[M];
//  double tmp=0.0;
  double norm2=0.0;
  double d_timesteps = (double)time_steps;
  int inc_array1=1;
  int inc_array2=1;
  int inc_tmp=1;
  /*for(i=0;i<M;i++){
      tmp = tmp + pow( (array1[i]-array2[i]), 2.0 );
      //printf("%14.12lf %14.12lf\n", array1[i],array2[i]);
  }
  norm2 = sqrt(tmp); */

  // array1(l_d[M])を、tmpにコピー
    cblas_dcopy(M, array1, inc_array1, tmp, inc_tmp);
  // tmp = -array2+tmp
    cblas_daxpy(M, -1.0, array2, inc_array2, tmp, inc_tmp);
  // tmpの2ノルムを求める
    norm2 = cblas_dnrm2(M, tmp, 1);

  // printf("norm2=%14.12lf\n", norm2);
  // overflowを防ぐ
  // for(i=0;i<M;i++){ array2[i] = array1[i] + (initial_d/norm2)*(array2[i]-array1[i]); }
  // printf("nomrm2=%14.12lf   log(norm2/initial_d)=%14.12lf\n",norm2, log(norm2/initial_d));
  *sum_log = *sum_log + log(norm2/initial_d);
//  printf("%d %.12lf %.12lf %.12lf %.12lf\n",time_steps, array1[0],array2[0], sum_log, norm2);
  // if(time_steps==1000) printf("%.12lf\n",norm2);
  // printf("%14.12lf %14.12lf %14.12lf %14.12lf\n",norm2,initial_d,sum_log,sum_log/d_timesteps);

  return(*sum_log/d_timesteps);
}

void exportLyapunovExponent(int time_steps, double dt, double lyapunov, char *dirname, FILE *fp5, char filename5[40]){
  double real_time = time_steps*dt;
  sprintf(filename5,"%s/data/le.dat",dirname);
  fp5 = fopen(filename5,"a");
  if( fp5 == NULL ){
    printf("loop count n=%d : cannot open file %s\n",time_steps,filename5);
  }
  fprintf(fp5,"%d  %lf %.12lf\n",time_steps, real_time, lyapunov);
  fclose(fp5);
}

void updateLearnigData(int time_steps, int M, int WASHOUT, int LEARNING, double l[M+1],double o_nrm2[3], double o_nrm10[11], double o_nrm20[21],double o_nrm30[31], double *T, double *L){
  int i=0;
  T[time_steps-WASHOUT] = o_nrm2[0];
  T[time_steps-WASHOUT+LEARNING] = o_nrm10[0];
  T[time_steps-WASHOUT+LEARNING*2] = o_nrm20[0];
  T[time_steps-WASHOUT+LEARNING*3] = o_nrm30[0];
   // printf("%lf %lf %lf\n",T[time_steps-WASHOUT], T[time_steps-WASHOUT+LEARNING], T[time_steps-WASHOUT+LEARNING*2]);
  for(i=0;i<(M+1);i++){
    L[time_steps-WASHOUT+LEARNING*i] = l[i];
  //    printf("%lf ",L[time_steps-WASHOUT+LEARNING*i]);
  }
  // printf("\n");
}

void getWeights(int M, int LEARNING, double *L, double *T, double *W_out){
  lapack_int nrhs, ldL, ldT, info;
  int i,j;
  nrhs = 4;
  ldL = LEARNING;
  ldT = LEARNING;
  //for(j=0;j<3;j++){
    //for(i=0;i<ldT;i++){
  //printf("%lf\n",T[i+ldT*j]);
//}
//}
  info = LAPACKE_dgels(LAPACK_COL_MAJOR,'N',LEARNING,M+1,nrhs,L,ldL,T,ldT);

  for(j=0;j<4;j++){
    for(i=0;i<(M+1);i++){
      W_out[i+(M+1)*j] = T[i+ldT*j];
    //  printf("%lf\n",W_out[i+M*j]);
    }
  }
}

void printWeights(int M, double *W_out){
  int i,j;
  printf("W_out:\n");
  printf("NARMA2  NARMA10 NARMA20 NARMA30\n");
  for(i=0;i<(M+1);i++){
    for(j=0;j<4;j++){
      printf("%15lf  ",W_out[i+(M+1)*j]);
    }
    printf("\n");
  }
}

void updateOutputsMS(int M, double o_ms[4], double l[M+1], double *W_out){
  int ldW_out = M+1; //W_outの行数
  int incl = 1;
  int inco_ms = 1;
//  double tmp[M+1];
//  cblas_dcopy(M+1, l, 1, tmp, 1);
//  tmp[M]=1.0;
  // 行列W_outの行数, 列数（転置前）
  cblas_dgemv(CblasColMajor,CblasTrans,M+1,4,
    1.0, W_out, ldW_out, l, incl,
      0.0, o_ms, inco_ms);
}

void exportCoordinates(int time_steps,char *dirname, FILE *fp1,  char *filename1, int N, double x[N], double y[N]){
  int i=0;
  for(i=0;i<N;i++){
    sprintf(filename1,"%s/data/points/point%d.dat",dirname,i);
    fp1 = fopen(filename1,"a");
    if( fp1 == NULL ){
      printf("loop count n=%d : cannot open file %s\n",time_steps,filename1);
    }
    fprintf(fp1,"%.12lf %.12lf\n",x[i],y[i]);
    fclose(fp1);
  }
}
void exportLength(int time_steps, double dt, char *dirname, FILE *fp2, char *filename2, int N, int G[N][N], int M, double l[M+1], double l_d[M]){
  int looked_idx = 0;
  int arrayl_idx = 0;
  int i=0,j=0;
  double real_time=time_steps*dt;
  for(i=0;i<N;i++){
    for(j=looked_idx;j<N;j++){
      if(G[i][j] == 1){
        sprintf(filename2,"%s/data/springs/length%d.dat",dirname,arrayl_idx);
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
  sprintf(filename2,"%s/data/springs/bias.dat",dirname);
  fp2 = fopen(filename2,"a");
  fprintf(fp2,"%d %lf %.12lf\n",time_steps,real_time,l[arrayl_idx]);
  fclose(fp2);
}

void exportOutputs(int time_steps, double dt, char *dirname, FILE *fp3, char *filename3, double o_ms[4], double o_nrm2[3], double o_nrm10[11], double o_nrm20[21],double o_nrm30[31] ){
  int i;
  double real_time = time_steps*dt;
  sprintf(filename3,"%s/data/outputs.dat",dirname);
  fp3 = fopen(filename3,"a");
  if( fp3 == NULL ){
    printf("loop count n=%d : cannot open file %s\n",time_steps,filename3);
  }
  // time_step, 質点バネ系による近似, NARMAモデルの出力値を1行に出力
  fprintf(fp3,"%d  %lf  ",time_steps,real_time);
  for(i=0;i<4;i++){
    fprintf(fp3,"%.12lf  ",o_ms[i]);
  }
  fprintf(fp3,"%.12lf  %.12lf  %.12lf  %.12lf",o_nrm2[0],o_nrm10[0],o_nrm20[0],o_nrm30[0]);
  fprintf(fp3,"\n");
  fclose(fp3);
}
void exportInput(int time_steps, double dt, char *dirname, FILE *fp7, char *filename7, double input[30]){
  int i;
  double real_time = time_steps*dt;
  sprintf(filename7,"%s/data/input.dat",dirname);
  fp7 = fopen(filename7,"a");
  if( fp7 == NULL ){
    printf("loop count n=%d : cannot open file %s\n",time_steps,filename7);
  }
  fprintf(fp7,"%d  %lf  %.12lf\n",time_steps,real_time,input[0]);
  fclose(fp7);
}
double rand_normal(double mu, double sigma){
      double z=sqrt( -2.0*log(Uniform()) ) * sin( 2.0*M_PI*Uniform() );
      return mu + sigma*z;
}
double Uniform( void ){
    return ((double)rand()+1.0)/((double)RAND_MAX+2.0);
}
