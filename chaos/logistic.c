#include <stdio.h>
#include <math.h>
double lyapunov=0.0;
double initial_d = 0.0;
double sum = 0.0;
double accurate_le = 0.0;
double sum_log = 0.0;
void updateLyapunovExponent(int time_steps, double x, double x_d);
void accurateLyapunovExponent(int time_steps, double a, double x);

int main(void){
  double a = 4.0; // chaos
  double x=0.01;
  double x_d=0.010000000001;
  int n=0;

//  initial_d = sqrt( pow( (x-x_d), 2.0 ) );
  initial_d = fabs(x-x_d);
  for(n=0;n<1000;n++){
    x=a*x*(1.0-x);
    x_d=a*x_d*(1.0-x_d);
    //printf("%d %f %f\n",n,x,x_d);
    updateLyapunovExponent(n,x,x_d);
    accurateLyapunovExponent(n,a,x);
    printf("%d %f %f\n",n, lyapunov, accurate_le);
  }
  return 0;
}

/* 2つの出力の差　引数(l,l_d) */
void updateLyapunovExponent(int time_steps, double x, double x_d){
  double norm2=0.0;
  norm2 = sqrt( pow( (x-x_d), 2.0 ) );
  sum += log(norm2/initial_d);
//  sum += log(norm2/initial_d);
  lyapunov = sum / time_steps;
}

// f'(x) = -2ax+a
void accurateLyapunovExponent(int time_steps, double a, double x){
  sum_log += log( fabs(-2.0*a*x+a) );
  accurate_le = sum_log/time_steps;
}

/*int i,j;
for(a=2.0;a<4;a+=0.002){
  for(i=0;i<100;i++){
    x=a*x*(1-x);
  }
  for(i=0;i<100;i++){
    x=a*x*(1-x);
    printf("%f %f\n",a,x);
  }
} */
