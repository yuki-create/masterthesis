#include <stdio.h>
#include <math.h>
double lyapunov=0.0;
double initial_d = 0.0;
double sum = 0.0;
double accurate_le = 0.0;
double sum_log = 0.0;
double a = 4.0; // chaos
double x=0.01;
double x_d=0.0100000001;
double sum_d = 0.0;
double norm2_pre = 0.0;

void updateLyapunovExponent(int time_steps);
void wolf(int time_steps);
void accurateLyapunovExponent(int time_steps);

int main(void){

  int n=0;

//  initial_d = sqrt( pow( (x-x_d), 2.0 ) );
  initial_d = fabs(x-x_d);
  norm2_pre = initial_d;
  for(n=0;n<1000;n++){
    x=a*x*(1.0-x);
    x_d=a*x_d*(1.0-x_d);
    //printf("%d %f %f\n",n,x,x_d);
    updateLyapunovExponent(n);
  //  wolf(n);
    accurateLyapunovExponent(n);
    printf("%d %f %f\n",n, lyapunov, accurate_le);
  }
  return 0;
}

/* 2つの出力の差　引数(l,l_d) */
 /*void updateLyapunovExponent(int time_steps, double x, double x_d){
  double norm2=0.0;
  norm2 =  fabs(x-x_d);
//  sum += log(norm2/initial_d);
  lyapunov = log(norm2/initial_d) / time_steps;
} */
void updateLyapunovExponent(int time_steps){
  double norm2=0.0;
  norm2 =  fabs(x-x_d); // gamma_k
//  sum += log(norm2/initial_d);
  x_d = x + (initial_d/norm2)*(x_d-x);
  sum_d = sum_d + log(norm2/initial_d);
  lyapunov = sum_d / time_steps;
}

void wolf(int time_steps){
  double norm2=0.0;
  norm2 =  fabs(x-x_d); // gamma_k
//  sum += log(norm2/initial_d);
  sum_d = sum_d + log(norm2/norm2_pre)/log(2.0);
  lyapunov = sum_d / time_steps;
  norm2_pre = norm2;
}


// f'(x) = -2ax+a
void accurateLyapunovExponent(int time_steps){
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
