// 実行は ./ms N NSTEP
#include <stdio.h>
#include <math.h>
#include "lapacke.h"
#include "cblas.h"
//////// parameters to be costomized ////////
#define N 36 // number of mass points
#define M 60 // number of springs (root_N-1)*root_N*2
const char *dirname = ".";
/* washout, learning, evaluating term (time steps) */
const int WASHOUT = 5000;
const int LEARNING = 5000;
const int EVAL = 1000;
/* simulating time steps */
const int NSTEP = WASHOUT+LEARNING+EVAL;

const double dt = 0.0025;
const int T_input = 1; // adjust frequency of input signal
const double gamma1 = 0.001;
const double k = 100.0;
const double natu_l = 1.0;
const double w_in[] = {1.0,1.5,0.5,2.0};
int fixed_p[] = {}; // index array of fixed points
int in_p[] = {0,10,20,30}; // index array of input points

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
int G[N][N];
double force_x[N];

/* variables and parameters for RK4 */
double ku1[N], ku2[N], ku3[N], ku4[N], kx1[N], kx2[N], kx3[N], kx4[N];
double kv1[N], kv2[N], kv3[N], kv4[N], ky1[N], ky2[N], ky3[N], ky4[N];

/* variables for NARMA models */
double o_nrm2[3]; // NARMA2 output
double o_nrm10[11]; // NARMA10 output
double o_nrm20[21]; // NARMA20 output

/* variables for learning */
double *T, *W_out, *L;

/* variables for EVAL phase */
double o_ms[3]; // [O_ms(t) for narma2, narma10, narma20]

/* parameters for file I/O */
FILE *fp1; // for export coodinates
FILE *fp2; // for export length of springs
char filename1[40]; // for export coodinates
char filename2[40];
/* decralation of functions */
void genGraph();
void printGraph();
void init();
void initForLapack();
void updateInput(int time_steps);
void updateNarma2();
void updateNarma10();
void updateNarma20();
void updateLearnigData(int time_steps);
void getWeights();
void printWeights();
void updateOutputsMS();
void test_updateOutputsMS(int time_steps);
void rk4();
void getSpringLength();
void exportCoordinates(int time_steps);
void exportLength(int time_steps);
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

  init();
  initForLapack();
  //  test_getSpringLength();
  // printGraph();

  /* washout phase */
  for(n=0;n<WASHOUT;n++){
    updateInput(n); //inputノードに外力を入力
    updateNarma2();
    updateNarma10();
    updateNarma20();
    rk4(); //全ての質点の座標の更新
  }

  /* learning phase */
  for(n=WASHOUT;n<WASHOUT+LEARNING;n++){
    //  test_updateNarma(n);
    updateInput(n); //inputノードに外力を入力
    updateNarma2();
    updateNarma10();
    updateNarma20();
    rk4(); //全ての質点の座標の更新
    getSpringLength(); //系の出力となるばねの長さを求め、配列l[]を更新
    //    exportCoordinates(n); //座標のデータをファイル出力
    //    exportLength(n); //ばねの長さをファイル出力
    updateLearnigData(n);
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
    rk4(); //全ての質点の座標の更新
    getSpringLength(); //系の出力となるばねの長さを求め、配列l[]を更新
    //    exportCoordinates(n); //座標のデータをファイル出力
    //    exportLength(n); //ばねの長さをファイル出力
    updateOutputsMS();
    test_updateOutputsMS(n);
  }

  //  test_updateLearningData();
  free(T);
  free(L);
  free(W_out);
  return (0);
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
  int looked_idx = 0;
  int arrayl_idx = 0;
  int i=0,j=0;

  genGraph();

  for(i=0;i<root_N;i++){
    for(j=0;j<root_N;j++){
      x[root_N*i+j] = j;
      y[root_N*i+j] = i;
    }
  }
  /* init spring length array l[M] with natural length*/
  for(i=0;i<M;i++){
    l[i] = natu_l;
  }
  for(i=0;i<N;i++){
    m[i] = 1.0;
    u[i] = 0.0;
    v[i] = 0.0;
    force_x[i] = 0.0;
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
  }

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

  /* make length of springs data files */
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

        arrayl_idx++;
      }
    }
    looked_idx++;
  }
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
        fprintf(fp2,"%d %f %f\n",time_steps,real_time,l[arrayl_idx]);
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

void getSpringLength(){
  int looked_idx = 0;
  int arrayl_idx = 0;
  int i,j=0;
  for(i=0;i<N;i++){
    for(j=looked_idx;j<N;j++){
      if(G[i][j] == 1){
        l[arrayl_idx] = sqrt( pow( (x[i]-x[j]), 2.0 ) + pow( (y[i]-y[j]), 2.0 ) );
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

void test_updateOutputsMS(int time_steps){
  int i;
/*  printf("outputs of MS for\n");
  printf("NARMA2   NARMA10   NARMA20:\n"); */
  printf("%d  ",time_steps);
  for(i=0;i<3;i++){
    printf("%f  ",o_ms[i]);
  }
  printf("%f  %f  %f",o_nrm2[0],o_nrm10[0],o_nrm20[0]);
  printf("\n");
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

void rk4(){
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

double Fx(double *array1, double *array2, double *array3, int idx1){
  int idx2=0;
  double ans = 0;
  /* get sum of elastic forces to array1[idx1] */
  double sum = 0;
  for(idx2=0;idx2<N;idx2++){
    if(G[idx1][idx2] == 1){
      sum += f(array1,array2,idx1,idx2);
    }
  }
  return ( sum - force_x[idx1] - gamma1 * array3[idx1] ) / m[idx1];
}
// force_xをx方向にのみ入力するために分ける
double Fy(double *array1, double *array2, double *array3, int idx1){
  int idx2=0;
  double ans = 0;
  /* get sum of elastic forces to array1[idx1] */
  double sum = 0;
  for(idx2=0;idx2<N;idx2++){
    if(G[idx1][idx2] == 1){
      sum += f(array1,array2,idx1,idx2);
    }
  }
  return ( sum  -  gamma1 * array3[idx1] ) / m[idx1];
}

double f(double *array1, double *array2, int idx1, int idx2){
  double ans = 0;
  double distance = sqrt( pow((array1[idx1]-array1[idx2]), 2.0) + pow((array2[idx1]-array2[idx2]), 2.0) );
  ans = k * ( natu_l* (array1[idx1]-array1[idx2]) / distance + array1[idx2] - array1[idx1]);
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
