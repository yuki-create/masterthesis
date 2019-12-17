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
void initCoordinates(int size, double array_x[size], double array_y[size], double array_x_d[size], double array_y_d[size]);
void initSpringParameters(int size, double mu_k, double sigma_k,double mu_g, double sigma_g, double array_k[size], double array_g[size]);
void initP2lMat(int size, int graph[size][size], int table[size][size]);
void updateInput(int time_steps, double dt, int T_input, double input[20], int in_num, int in_p[in_num],int N, double force_x[N],double w_in[in_num], double x[N]);
void updateNarma(double input[20], double o_nrm2[3], double o_nrm10[11], double o_nrm20[21]);
void rk4(int N, double dt, double x[N], double y[N], double u[N], double v[N], int fixed_num, int fixed_p[fixed_num], int G[N][N], int p2l_mat[N][N],double force_x[N], double m[N], int M, double k[M], double gamma1[M], double natu_l);
double Fx(double *array1, double *array2, double *array3, int idx1, int N, int G[N][N], int p2l_mat[N][N],double force_x[N], double m[N], int M, double k[M], double gamma1[M], double natu_l );
double Fy(double *array1, double *array2, double *array3, int idx1, int N, int G[N][N], int p2l_mat[N][N],double force_x[N], double m[N], int M, double k[M], double gamma1[M], double natu_l );
double f(double *array1, double *array2, int idx1, int idx2, int N, int p2l_mat[N][N], int M, double k[M], double natu_l);
void getSpringLength(double *array_l, double *array_x, double *array_y, int N, int G[N][N]);
double updateLyapunovExponent(int time_steps, double *array1, double *array2, int M, double initial_d);
void exportLyapunovExponent(int time_steps, double dt, double lyapunov, char *dirname, FILE *fp5, char filename5[40]);
void updateLearnigData(int time_steps, int M, int WASHOUT, int LEARNING, double l[M],double o_nrm2[3], double o_nrm10[11], double o_nrm20[21], double *T, double *L);
void getWeights(int M, int LEARNING, double *L, double *T, double *W_out);
void printWeights(int M, double *W_out);
void updateOutputsMS(int M, double o_ms[3], double l[M], double *W_out);
void exportCoordinates(int time_steps,char *dirname, FILE *fp1,  char *filename1, int N, double x[N], double y[N]);
void exportLength(int time_steps, double dt, char *dirname, FILE *fp2, char *filename2, int N, int G[N][N], int M, double l[M]);
void exportOutputs(int time_steps, double dt, char *dirname, FILE *fp3, char *filename3, double o_ms[3], double o_nrm2[3], double o_nrm10[11], double o_nrm20[21]);

int main(int argc, char *argv[]){
  //////// parameters to be costomized ////////
  int seed_flag = 0;
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
  double w_in[] = {1.0};
  int fixed_p[] = {}; // index array of fixed points
  int in_p[] = {0}; // index array of input points

  //////// variables dosen't need to be costomized ////////
  double k[M];
  double gamma1[M];
  double mu_k = strtod(argv[3],NULL);
  double sigma_k = strtod(argv[4],NULL);
  double mu_g = strtod(argv[5],NULL);
  double sigma_g = strtod(argv[6],NULL);
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
  //double ku1[N], ku2[N], ku3[N], ku4[N], kx1[N], kx2[N], kx3[N], kx4[N];
  //double kv1[N], kv2[N], kv3[N], kv4[N], ky1[N], ky2[N], ky3[N], ky4[N];

  // variables for mass-spring system with minute initial state
  double x_d[N];
  double u_d[N]; // dx/dt
  double y_d[N];
  double v_d[N]; // dy/dt
  //double ku1_d[N], ku2_d[N], ku3_d[N], ku4_d[N], kx1_d[N], kx2_d[N], kx3_d[N], kx4_d[N];
  //double kv1_d[N], kv2_d[N], kv3_d[N], kv4_d[N], ky1_d[N], ky2_d[N], ky3_d[N], ky4_d[N];
  double l_d[M];

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
  double initial_d = 0.0; // initial distance of two systems
  // double norm2_pre = 0.0;

  // parameters for file I/O
  char *dirname = argv[7];
  FILE *fp1; // for export coodinates
  FILE *fp2; // for export length of springs
  FILE *fp3; // for outputs of MS and NARMA after learning
  FILE *fp4; // for parameters and results file
  FILE *fp5; // for export lyapunov exponent
  FILE *fp6; // for export many results
  char filename1[40],filename2[40],filename3[40],filename4[40],filename5[40],filename6[40];
  int i,j,n;
  double tmp = 0.0;
  // 初期化始まり　init start
  fixed_num = sizeof fixed_p / sizeof fixed_p[0];
  in_num = sizeof in_p / sizeof in_p[0];
  if(seed_flag == 1) srand((unsigned)time(NULL));
  initArray1dim(N,u,0.0);
  initArray1dim(N,v,0.0);
  initArray1dim(N,u_d,0.0);
  initArray1dim(N,v_d,0.0);
  initArray1dim(N,m,1.0);
  initArray1dim(N,force_x,0.0);
  //initArray1dim(M,l,natu_l);
  //initArray1dim(M,l_d,natu_l);
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      p2l_mat[i][j] = -1;
    }
  }
  genGraph(N,G);
  initP2lMat(N,G,p2l_mat);
  initCoordinates(N,x,y,x_d,y_d);
  initSpringParameters(M,mu_k,sigma_k,mu_g,sigma_g,k,gamma1);
  getSpringLength(l,x,y,N,G);
  getSpringLength(l_d,x_d,y_d,N,G);
  for(i=0;i<M;i++){
      initial_d += pow( (l[i]-l_d[i]), 2.0 );
  }
  initial_d = sqrt(initial_d);
//  printf("%14.12lf\n",initial_d);
  // lapackのライブラリを使う配列の初期化。行列の列ごとに一次元配列に格納。
  T = malloc(LEARNING*3*sizeof(double));
  L = malloc(LEARNING*M*sizeof(double));
  W_out = malloc(M*3*sizeof(double));
  // 考えている行列の、列ごとに要素を代入する
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
  // 出力波形を見たりする時、ファイルを初期化
  if(debug_flag == 1){
    int looked_idx = 0;
    int arrayl_idx = 0;
    int i,j;
    // make coodinates data files
    for(i=0;i<N;i++){
      sprintf(filename1,"%s/points/point%d.dat",dirname,i);
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
          sprintf(filename2,"%s/springs/length%d.dat",dirname,arrayl_idx);
          fp2 = fopen(filename2,"w");
          if( fp2 == NULL ){
            printf("cannot open file %s\n",filename2);
          }
          else{
            //    printf("open file %s\n",filename2);
          }
          fprintf(fp2,"# l[%d]( l_%d(t) ) connection: point%d-point%d\n",arrayl_idx,arrayl_idx,i,j);
          fprintf(fp2,"n t l_%d(t)\n",arrayl_idx);
          fprintf(fp2,"0 0 %.12lf\n",l[arrayl_idx]);
          fclose(fp2);
          arrayl_idx++;
        }
      }
      looked_idx++;
    }
    // make outputs.dat file
    sprintf(filename3,"%s/results/outputs.dat",dirname);
    fp3 = fopen(filename3,"w");
    if( fp3 == NULL ){
      printf("cannot open file %s\n",filename3);
    }
    else{
      //    printf("open file %s\n",filename3);
    }
    fprintf(fp3,"#n  real_time  ms_narma2  ms_narma10  ms_narma20  narma2  narma10  narma20\n");
    fclose(fp3);

    // make le.dat file
    sprintf(filename5,"%s/results/le.dat",dirname);
    fp5 = fopen(filename5,"w");
    if( fp5 == NULL ){
      printf("cannot open file %s\n",filename5);
    }
    else{
      //    printf("open file %s\n",filename3);
    }
    fprintf(fp5,"#n  real_time  lyapunov_exponential\n");
    fclose(fp5);
    //printGraph(N,G,p2l_mat);
  }
  // 初期化終わり　init end
  // washout
  for(n=0;n<WASHOUT;n++){
    updateInput(n,dt, T_input, input, in_num, in_p ,N, force_x, w_in,x);
     updateNarma(input, o_nrm2, o_nrm10, o_nrm20);
     rk4(N, dt, x, y, u, v, fixed_num, fixed_p, G, p2l_mat, force_x, m, M, k, gamma1, natu_l);
      rk4(N, dt, x_d, y_d, u_d, v_d, fixed_num, fixed_p, G, p2l_mat, force_x, m, M, k, gamma1, natu_l);
     getSpringLength(l,x,y,N,G);
     getSpringLength(l_d,x_d,y_d,N,G);
     lyapunov = updateLyapunovExponent(n,l,l_d, M, initial_d);
     if(debug_flag==1){
       exportCoordinates(n,dirname,fp1,filename1, N, x, y);
       exportLength(n, dt, dirname, fp2, filename2, N, G, M, l);
       exportLyapunovExponent(n, dt, lyapunov, dirname, fp5, filename5);
     }
  }
  // 学習
  for(n=WASHOUT;n<WASHOUT+LEARNING;n++){
    updateInput(n,dt, T_input, input, in_num, in_p ,N, force_x, w_in,x);
     updateNarma(input, o_nrm2, o_nrm10, o_nrm20);
     rk4(N, dt, x, y, u, v, fixed_num, fixed_p, G, p2l_mat, force_x, m, M, k, gamma1, natu_l);
      rk4(N, dt, x_d, y_d, u_d, v_d, fixed_num, fixed_p, G, p2l_mat, force_x, m, M, k, gamma1, natu_l);
      getSpringLength(l,x,y,N,G);
      getSpringLength(l_d,x_d,y_d,N,G);
      updateLearnigData(n, M, WASHOUT, LEARNING, l, o_nrm2, o_nrm10, o_nrm20, T, L);
      lyapunov = updateLyapunovExponent(n,l,l_d, M, initial_d);
      if(debug_flag==1){
        exportCoordinates(n,dirname,fp1,filename1, N, x, y);
        exportLength(n, dt, dirname, fp2, filename2, N, G, M, l);
        exportLyapunovExponent(n, dt, lyapunov, dirname, fp5, filename5);
      }
  }
  getWeights(M, LEARNING, L, T, W_out);
  //printWeights(M, W_out);

  for(n=WASHOUT+LEARNING;n<WASHOUT+LEARNING+EVAL;n++){
    updateInput(n,dt, T_input, input, in_num, in_p ,N, force_x, w_in,x);
     updateNarma(input, o_nrm2, o_nrm10, o_nrm20);
     rk4(N, dt, x, y, u, v, fixed_num, fixed_p, G, p2l_mat, force_x, m, M, k, gamma1, natu_l);
      rk4(N, dt, x_d, y_d, u_d, v_d, fixed_num, fixed_p, G, p2l_mat, force_x, m, M, k, gamma1, natu_l);
      getSpringLength(l,x,y,N,G);
      getSpringLength(l_d,x_d,y_d,N,G);
      updateOutputsMS(M, o_ms, l, W_out);
    // updateErr
    squared_err[0] += pow(o_nrm2[0]-o_ms[0], 2.0);
    squared_err[1] += pow(o_nrm10[0]-o_ms[1], 2.0);
    squared_err[2] += pow(o_nrm20[0]-o_ms[2], 2.0);
    normalizer[0] += pow(o_nrm2[0],2.0);
    normalizer[1] += pow(o_nrm10[0],2.0);
    normalizer[2] += pow(o_nrm20[0],2.0);

      lyapunov = updateLyapunovExponent(n,l,l_d, M, initial_d);
      if(debug_flag==1){
        exportCoordinates(n,dirname,fp1,filename1, N, x, y);
        exportLength(n, dt, dirname, fp2, filename2, N, G, M, l);
        exportOutputs(n, dt, dirname, fp3, filename3, o_ms, o_nrm2, o_nrm10, o_nrm20);
        exportLyapunovExponent(n, dt, lyapunov, dirname, fp5, filename5);
      }

  }
  // getErr
  printf("lyapunov exponent = %lf\n", lyapunov);
  printf("NMSE:\nNARMA2    10    20\n");
  for(i=0;i<3;i++){
    err[i] = squared_err[i] / normalizer[i];
    printf("%f    ",err[i]);
  }
  printf("\n");

    // make results.dat file (parameters and results)
    sprintf(filename4,"%s/results/results.txt",dirname);
    fp4 = fopen(filename4,"w");
    if( fp4 == NULL ){
      printf("cannot open file %s\n",filename4);
    }
    else{
    fprintf(fp4, "lyapunov exponent = %lf\n", lyapunov);
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

    fprintf(fp4,"\nk: mu = %f, sigma = %f",mu_k, sigma_k);
    fprintf(fp4,"\ngamma: mu = %f, sigma = %f",mu_g, sigma_g);

    fprintf(fp4,"\n---index_of_l[i]  k  gamma---");
    for(i=0;i<M;i++){
      fprintf(fp4,"\n%d   %f   %f",i,k[i],gamma1[i]);
    }
    fprintf(fp4,"\n");
      fclose(fp4);
    }

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

void initArray1dim(int size, double array[size], double initial_value){
  int i;
  for(i=0;i<size;i++){
    array[i] = initial_value;
  }
}
// x_d[0],y_d[0]のみ初期値をずらす
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
    array_x_d[0] = array_x_d[0] + 0.00000001;
    array_y_d[0] = array_y_d[0] + 0.00000001;
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

void updateInput(int time_steps, double dt, int T_input, double input[20], int in_num, int in_p[in_num],int N, double force_x[N],double w_in[in_num],double x[N]){
  int i=0;
  int idx=0;
  double t = time_steps*dt/T_input;
  for(i=0;i<19;i++){
    input[19-i] = input[18-i];
  }
  input[0] = 0.2*sin(2.0*M_PI*2.11*t)*sin(2.0*M_PI*3.73*t)*sin(2.0*M_PI*4.33*t);

  for(i=0;i<in_num;i++){

    idx = in_p[i];
    force_x[idx] = w_in[i]*input[0];
    //x[idx] = x[idx] + w_in[i]*input[0];
    // printf("%d: %f * %f = %f\n",idx, w_in[i],input[0],force_x[idx]);

  }
}

void updateNarma(double input[20], double o_nrm2[3], double o_nrm10[11], double o_nrm20[21]){
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
  // update NARMA20 model
  tmp = 0.0;
  for(i=0;i<20;i++){
    o_nrm20[20-i] = o_nrm20[19-i];
    tmp += o_nrm20[20-i];
  }
  o_nrm20[0] = 0.3*o_nrm20[1] + 0.05*o_nrm20[1]*tmp + 1.5*input[19]*input[0] + 0.1;
//  printf("%lf\n",o_nrm20[0]);
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
    for(j=0;j<fixed_num;j++){
      if(i==fixed_p[j]){ fixed_flag = 1; break; }
    }
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
  int arrayl_idx = p2l_mat[idx1][idx2];
  for(idx2=0;idx2<N;idx2++){
    if(G[idx1][idx2] == 1){
      sum += f(array1,array2,idx1,idx2,N,p2l_mat,M,k, natu_l);
    }
  }
  return ( sum - force_x[idx1] - gamma1[arrayl_idx] * array3[idx1] ) / m[idx1];
}
// force_xをx方向にのみ入力するために分ける
double Fy(double *array1, double *array2, double *array3, int idx1, int N, int G[N][N], int p2l_mat[N][N],double force_x[N], double m[N], int M, double k[M], double gamma1[M], double natu_l){
  int idx2=0;
  double ans = 0;
  // get sum of elastic forces to array1[idx1]
  double sum = 0;
  int arrayl_idx = p2l_mat[idx1][idx2];
  for(idx2=0;idx2<N;idx2++){
    if(G[idx1][idx2] == 1){
      sum += f(array1,array2,idx1,idx2,N,p2l_mat,M,k, natu_l);
    }
  }
  return ( sum  -  gamma1[arrayl_idx] * array3[idx1] ) / m[idx1];
}

double f(double *array1, double *array2, int idx1, int idx2, int N, int p2l_mat[N][N], int M, double k[M], double natu_l){
  double ans = 0;
  int arrayl_idx = p2l_mat[idx1][idx2];
  // debug
   if(arrayl_idx == -1){
      printf("arrayl_idx = -1: idx1=%d, idx2=%d\n",idx1,idx2);
  }
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
double updateLyapunovExponent(int time_steps, double *array1, double *array2, int M,  double initial_d){
  static double sum_log = 0.0;
  int i=0;
  double tmp=0;
  double norm2=0.0;
  double d_timesteps = (double)time_steps;
  for(i=0;i<M;i++){
      tmp = tmp + pow( (array1[i]-array2[i]), 2.0 );
    //  printf("%14.12lf %14.12lf\n", array1[i],array2[i]);
  }
  norm2 = sqrt(tmp);
  // norm2がおかしい
  // printf("norm2=%14.12lf\n", norm2);
  // overflowを防ぐ
  // for(i=0;i<M;i++){ array2[i] = array1[i] + (initial_d/norm2)*(array2[i]-array1[i]); }
  // printf("nomrm2=%14.12lf   log(norm2/initial_d)=%14.12lf\n",norm2, log(norm2/initial_d));
  sum_log = sum_log + log(norm2/initial_d);
  // printf("%14.12lf %14.12lf %14.12lf %14.12lf\n",norm2,initial_d,sum_log,sum_log/d_timesteps);
  return(sum_log/d_timesteps);
}

void exportLyapunovExponent(int time_steps, double dt, double lyapunov, char *dirname, FILE *fp5, char filename5[40]){
  sprintf(filename5,"%s/results/le.dat",dirname);
  fp5 = fopen(filename5,"a");
  if( fp5 == NULL ){
    printf("loop count n=%d : cannot open file %s\n",time_steps,filename5);
  }
  fprintf(fp5,"%d  %lf %.12lf\n",time_steps, dt, lyapunov);
  fclose(fp5);
}

void updateLearnigData(int time_steps, int M, int WASHOUT, int LEARNING, double l[M],double o_nrm2[3], double o_nrm10[11], double o_nrm20[21], double *T, double *L){
  int i=0;
  T[time_steps-WASHOUT] = o_nrm2[0];
  T[time_steps-WASHOUT+LEARNING] = o_nrm10[0];
  T[time_steps-WASHOUT+LEARNING*2] = o_nrm20[0];
   // printf("%lf %lf %lf\n",T[time_steps-WASHOUT], T[time_steps-WASHOUT+LEARNING], T[time_steps-WASHOUT+LEARNING*2]);
  for(i=0;i<M;i++){
    L[time_steps-WASHOUT+LEARNING*i] = l[i];
  //    printf("%lf ",L[time_steps-WASHOUT+LEARNING*i]);
  }
  // printf("\n");
}

void getWeights(int M, int LEARNING, double *L, double *T, double *W_out){
  lapack_int nrhs, ldL, ldT, info;
  int i,j;
  nrhs = 3;
  ldL = LEARNING;
  ldT = LEARNING;
  //for(j=0;j<3;j++){
    //for(i=0;i<ldT;i++){
  //printf("%lf\n",T[i+ldT*j]);
//}
//}
  info = LAPACKE_dgels(LAPACK_COL_MAJOR,'N',LEARNING,M,nrhs,L,ldL,T,ldT);

  for(j=0;j<3;j++){
    for(i=0;i<M;i++){
      W_out[i+M*j] = T[i+ldT*j];
    //  printf("%lf\n",W_out[i+M*j]);
    }
  }
}

void printWeights(int M, double *W_out){
  int i,j;
  printf("W_out:\n");
  printf("NARMA2  NARMA10 NARMA20\n");
  for(i=0;i<M;i++){
    for(j=0;j<3;j++){
      printf("%15lf  ",W_out[i+M*j]);
    }
    printf("\n");
  }
}

void updateOutputsMS(int M, double o_ms[3], double l[M], double *W_out){
  int ldW_out = M; //W_outの行数
  int incl = 1;
  int inco_ms = 1;
  // 行列W_outの行数, 列数（転置前）
  cblas_dgemv(CblasColMajor,CblasTrans,M,3,
    1.0, W_out, ldW_out, l, incl,
      0.0, o_ms, inco_ms);
}

void exportCoordinates(int time_steps,char *dirname, FILE *fp1,  char *filename1, int N, double x[N], double y[N]){
  int i=0;
  for(i=0;i<N;i++){
    sprintf(filename1,"%s/points/point%d.dat",dirname,i);
    fp1 = fopen(filename1,"a");
    if( fp1 == NULL ){
      printf("loop count n=%d : cannot open file %s\n",time_steps,filename1);
    }
    //  printf("%f %f\n",x[0],y[0]);
    fprintf(fp1,"%.12lf %.12lf\n",x[i],y[i]);
    fclose(fp1);
  }
}
void exportLength(int time_steps, double dt, char *dirname, FILE *fp2, char *filename2, int N, int G[N][N], int M, double l[M]){
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
        fprintf(fp2,"%d %lf %.12lf\n",time_steps,real_time,l[arrayl_idx]);
        fclose(fp2);
        arrayl_idx++;
      }
    }
    looked_idx++;
  }
}
void exportOutputs(int time_steps, double dt, char *dirname, FILE *fp3, char *filename3, double o_ms[3], double o_nrm2[3], double o_nrm10[11], double o_nrm20[21] ){
  int i;
  double real_time = time_steps*dt;
  sprintf(filename3,"%s/results/outputs.dat",dirname);
  fp3 = fopen(filename3,"a");
  if( fp3 == NULL ){
    printf("loop count n=%d : cannot open file %s\n",time_steps,filename3);
  }
  // time_step, 質点バネ系による近似, NARMAモデルの出力値を1行に出力
  fprintf(fp3,"%d  %lf  ",time_steps,real_time);
  for(i=0;i<3;i++){
    fprintf(fp3,"%.12lf  ",o_ms[i]);
  }
  fprintf(fp3,"%.12lf  %.12lf  %.12lf",o_nrm2[0],o_nrm10[0],o_nrm20[0]);
  fprintf(fp3,"\n");
  fclose(fp3);
}

double rand_normal(double mu, double sigma){
      double z=sqrt( -2.0*log(Uniform()) ) * sin( 2.0*M_PI*Uniform() );
      return mu + sigma*z;
}
double Uniform( void ){
    return ((double)rand()+1.0)/((double)RAND_MAX+2.0);
}
