// 実行は ./ms N NSTEP
#include <stdio.h>
#include <math.h>

/* parameters for mass spring */
#define N 9 // number of mass points
#define M 12 // number of springs (root_N-1)*root_N*2
const int NSTEP = 10000; /* simulating time steps */
const double gamma1 = 0.001;
const double k = 100.0;
const double natu_l = 1.0;
const char *dirname = "test_narma";
const double w_in[] = {1.0};
const double dt = 0.0025;
const int T_input = 1; // adjust frequency of input signal
int fixed_p[] = {}; // index array of fixed points
int in_p[] = {0}; // index array of input points
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
/* parameters for file I/O */
FILE *fp;
char filename[40];
/* decralation of functions */
void genGraph();
void dammy_genGraph();
void printGraph();
void init();
void updateInput(int time_steps);
void updateNarma2();
void updateNarma10();
void updateNarma20();
void rk4();
void eular();
void getSpringLength();
void exportCoordinates();
double Fx(double *array1, double *array2, double *array3, int idx1);
double Fy(double *array1, double *array2, double *array3, int idx1);
double f(double *array1, double *array2, int idx1, int idx2);
void test_getSpringLength();
void test_updateNarma(int time_steps);

//main(コマンドライン引数の個数，引数を文字列として保存する配列)
// argv[] = {"./ms", "N", "NSTEP"}
int main(int argc, char *argv[]){
  int i,n;
  fixed_num = sizeof fixed_p / sizeof fixed_p[0];
  in_num = sizeof in_p / sizeof in_p[0];

  init();
  genGraph();
  // printGraph();
//  printf("simulating mass-spring system...\n");
  for(n=0;n<NSTEP;n++){
    test_updateNarma(n);
    updateInput(n); //inputノードに外力を入力
    updateNarma2();
    updateNarma10();
    updateNarma20();
    rk4(); //全ての質点の座標の更新
    getSpringLength(); //系の出力となるばねの長さを求め、配列l[]を更新
    exportCoordinates(n); //座標のデータをファイル出力
  }
  return (0);
}

void init(){
  /* init coordinates of mass points */
  int root_N = sqrt(N);
  int i=0,j=0;
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

  /* make data files */
  for(i=0;i<N;i++){
    sprintf(filename,"./%s/data_point%d.dat",dirname,i);
    fp = fopen(filename,"w");
    /* error handling */
    if( fp == NULL ){
      printf("cannot open file %s\n",filename);
    }
    else{
  //    printf("open file %s\n",filename);
    }
    fprintf(fp,"x[%d] y[%d]\n",i,i);
    fprintf(fp,"%f %f\n",x[i],y[i]); //coordinates at n=0 (t=0)初期値を書き込み
    //  printf("%f %f\n",x[i],y[i]); //debug
    fclose(fp);
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
void exportCoordinates(int n){
  int i=0;
  for(i=0;i<N;i++){
    sprintf(filename,"./%s/data_point%d.dat",dirname,i);
    fp = fopen(filename,"a"); // is option "w" ok?
    if( fp == NULL ){
      printf("loop count n=%d : cannot open file %s\n",n,filename);
    }
    //  printf("%f %f\n",x[0],y[0]);
    fprintf(fp,"%f %f\n",x[i],y[i]);
    fclose(fp);
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
  // diverged
  o_nrm2[0] = 0.4 * o_nrm2[1] + 0.4 * o_nrm2[1] * o_nrm2[2] + 0.6 * pow(input[0], 3.0) + 0.1;
  //not diverged
  //o_nrm2[0] = 0.1 * o_nrm2[1] + 0.4 * o_nrm2[1] * o_nrm2[2] + 0.6 * pow(input[0], 3.0) + 0.1;
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
        arrayl_idx++;
        printf("l[%d] connection x[i]-x[j]: %d-%d\n",arrayl_idx,i,j);
      }
    }
    looked_idx++;
  }
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
