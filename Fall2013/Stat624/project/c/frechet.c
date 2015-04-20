#include <stdio.h>             // standard input/output
#include <stdlib.h>            // malloc
#include <math.h>              // fabs, sqrt, etc.
#include <time.h>              // time
#include <unistd.h>            // getpid
#include <gsl/gsl_rng.h>       // GNU Scientific Library
#include <gsl/gsl_cdf.h>       // GNU Scientific Library
#include <gsl/gsl_randist.h>   // GNU Scientific Library

#define pi 3.14159265358979323846

double runif(){ // use a different library for this!
  return (double)rand()/(double)RAND_MAX;
}

double Min(double* x, int n){
  double m; m = 1000000;
  for (int i=0; i<n; i++){
    if (x[i] < m) {m = x[i];}
  }
  return m;
}

void dfrechet(double* x, double* p){
  // check: a = shape > 0; m = min; s = scale
  // x > m
  double a = p[0]; double m = p[1]; double s = p[2];
  *x = a/s * pow(((*x-m)/s),-1-a) * exp(-pow((*x-m)/s,-a));
}

void rfrechet(int* n, double* p, double* x){
  // check: a > 0
  double a = p[0]; double m = p[1]; double s = p[2];
  for (int i=0; i<*n; i++){
    double u = runif();
    x[i] = pow(-log(u),-1/a) * s + m;
  }  
}

double pfrechet(double x, double* p){
  // check: x > m
  double a = p[0]; double m = p[1]; double s = p[2];
  return exp(-pow((x-m)/s,-a));
}

double mean(double* x, int n){
  double sum = 0;
  for(int i=0; i<n; i++){
    sum = sum + x[i];
  }
  return sum/n;
}

double sd(double* x, int n){
  double sumSq = 0;
  double m = mean(x,n);
  for(int i=0; i<n; i++){
    sumSq += pow(x[i] - m, 2);
  }
  return sqrt(sumSq/n);

}

double lL(double* x, int n, double* p){
  double a = p[0]; double m = p[1]; double s = p[2];
  double z;
  double sum = 0;
  //printf("%s\t %f\t %f\t %f\t %f\t \n","LIKELIHOOD:", a, m, s);
  for (int i=0; i<n; i++){
    z = (x[i] - m) / s;
    //printf("%f\n",x[i]-m);
    sum += (1+a) * log(z) + pow(z,-a);
    //printf("%s\t %f\t %f\t %f\t %f\t \n","sum x m z:",sum,x[i],m,z);
  }
  return n*(log(a) - log(s)) - sum;
}
double lfa(double* x, int n, double* p){
  double a = p[0]; double m = p[1]; double s = p[2];
  double sum = 0;
  double z;
  for (int i=0; i<n; i++){
    z = (x[i] - m) / s;
    sum += log(z) * (pow(z,-a) - 1);
  }
  return n/a + sum;
}
double lfm(double* x, int n, double* p){
  double a = p[0]; double m = p[1]; double s = p[2];
  double sum = 0;
  double z;
  for (int i=0; i<n; i++){
    z = (x[i] - m) / s;
    sum += (a + 1 - a*pow(z,-a)) / (x[i]-m);
  }
  return sum;
}
double lfs(double* x, int n, double* p){
  double a = p[0]; double m = p[1]; double s = p[2];
  double sum = 0;
  double z;
  for (int i=0; i<n; i++){
    z = (x[i] - m) / s;
    sum += a+1 - a*pow(z,-a);
  }
  return -n/s + sum/s;
}

/* double derivatives: 6 unique ones */
double lfaa(double* x, int n, double* p){
  double a = p[0]; double m = p[1]; double s = p[2];
  double sum = 0;
  double z;
  for (int i=0; i<n; i++){
    z = (x[i] - m) / s;
    sum += log(z)*log(z) * pow(z,-a);
  }
  return -n/(a*a) - sum;
}

double lfmm(double* x, int n, double* p){
  double a = p[0]; double m = p[1]; double s = p[2];
  double sum = 0;
  double z;
  for (int i=0; i<n; i++){
    z = (x[i] - m) / s;
    //sum += (a+1) * (pow(z,a)-a) * pow(z*s,-2) * pow(z,-a);
    sum += (a+1)/pow(z*s,2) * (1-a*pow(z,-a));
  }
  return sum;
}

double lfss(double* x, int n, double* p){
  double a = p[0]; double m = p[1]; double s = p[2];
  double sum = 0;
  double z;
  for (int i=0; i<n; i++){
    z = (x[i] - m) / s;
    sum += (-(a+1)) - a*(a-1)*pow(z,-a);
  }
  //return n * (1+a + 1/(s*s)) - a*(a-1)*pow(s,a-2)*sum;
  return (n + sum)/(s*s);
}

double lfam(double* x, int n, double* p){
  double a = p[0]; double m = p[1]; double s = p[2];
  double sum = 0;
  double z;
  for (int i=0; i<n; i++){
    z = (x[i] - m) / s;
    sum += ((a * log(z) - 1) * pow(z,-a) + 1) / (z*s);
  }
  return sum;
}

double lfms(double* x, int n, double* p){
  double a = p[0]; double m = p[1]; double s = p[2];
  double sum = 0;
  double z;
  for (int i=0; i<n; i++){
    z = (x[i] - m) / s;
    sum += pow(z*s,-a-1);
  }
  return -a*a*pow(s,a-1) * sum; 
}

double lfas(double* x, int n, double* p){
  double a = p[0]; double m = p[1]; double s = p[2];
  double sum = 0;
  double z;
  for (int i=0; i<n; i++){
    z = (x[i] - m) / s;
    sum += (a*log(z) + pow(z,a) - 1) / (s*pow(z,a));
  }
  return sum;
}

void solve(double H[3][3], double g[3], double* Hig){
  double D;
  double Hi[3][3];

  //D = -H[0][0]*H[0][1]*H[0][1] + H[0][0]*H[0][0]*H[1][1] +
  //     H[0][2]*H[0][2]*H[1][1] + 2 * H[0][1]*H[0][2]*H[1][2] -
  //     H[0][0]*H[1][2]*H[1][2];

  D = H[0][0] * (H[2][2]*H[1][1] - H[2][1]*H[1][2]) -
      H[1][0] * (H[2][2]*H[0][1] - H[2][1]*H[0][2]) +
      H[2][0] * (H[1][2]*H[0][1] - H[1][1]*H[0][2]) ;
  
  
  Hi[0][0] = H[2][2] * H[1][1] - H[2][1] * H[1][2];
  Hi[0][1] = H[2][1] * H[0][2] - H[2][2] * H[0][1];
  Hi[0][2] = H[1][2] * H[0][1] - H[1][1] * H[0][2];
  Hi[1][0] = Hi[0][1];
  Hi[1][1] = H[2][2] * H[0][0] - H[2][0] * H[0][2];
  Hi[1][2] = H[1][0] * H[0][2] - H[1][2] * H[0][0];
  Hi[2][0] = Hi[0][2];
  Hi[2][1] = Hi[1][2];
  Hi[2][2] = H[1][1] * H[0][0] - H[1][0] * H[0][1];

  Hig[0] = (Hi[0][0]*g[0] + Hi[0][1]*g[1] + Hi[0][2]*g[2])/D;
  Hig[1] = (Hi[1][0]*g[0] + Hi[1][1]*g[1] + Hi[1][2]*g[2])/D;
  Hig[2] = (Hi[2][0]*g[0] + Hi[2][1]*g[1] + Hi[2][2]*g[2])/D; 
  
}

void newton(double* x, int* n1, double* p){
  int maxIt = 200;
  double n = *n1;
  double H[3][3] = { {lfaa(x,n,p),lfam(x,n,p),lfas(x,n,p)},
                     {lfam(x,n,p),lfmm(x,n,p),lfms(x,n,p)},
                     {lfas(x,n,p),lfms(x,n,p),lfss(x,n,p)} };

  double g[3] =      {lfa(x,n,p), lfm(x,n,p), lfs(x,n,p)}; 
  
  int N = 1;
  double tol = .1;
  double Hi[3][3];
  double* Hig = (double*) malloc(3*sizeof(double*));
  double D;

  while ( (N < maxIt) & 
          (sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]) > tol) ){

    solve(H,g,Hig);

    p[0] = p[0] - Hig[0];
    p[1] = p[1] - Hig[1];
    p[2] = p[2] - Hig[2];
    
    H[0][0]=lfaa(x,n,p); H[0][1]=lfam(x,n,p); H[0][2]=lfas(x,n,p);
    H[1][0]=H[0][1];     H[1][1]=lfmm(x,n,p); H[1][2]=lfms(x,n,p);
    H[2][0]=H[0][2];     H[2][1]=H[1][2];     H[2][2]=lfss(x,n,p);
    

    g[0] = lfa(x,n,p);
    g[1] = lfm(x,n,p);
    g[2] = lfs(x,n,p);

    N += 1;
    //printf("%s%d\n","Iteration: ",N);
  } //End of while loop
  
  free(Hig);
}

double ldgamma(double x, double a, double b){
  return (a-1) * log(x) -x/b;
}

double ldnorm(double x, double mu, double sig){
  return -pow((x-mu)/sig, 2)/2;
}

double logPo(double t, double* x, int n, double min, double* p, int par){
  double logPrior = 0;
  double logLik,z;
  p[par] = t;
  double a = p[0]; double m = p[1]; double s = p[2];
  double sum = 0;

  if (par==0){
    for (int i=0; i<n; i++){
      z = (x[i]-m)/s;
      sum += a*log(z) + pow(z,-a);
    }
    logLik = n*log(a) -sum;
    //logPrior = ldgamma(a,mean(x,n),1);
  }
  else if (par==1){
    for (int i=0; i<n; i++){
      z = (x[i]-m)/s;
      sum += (1+a)*log(z) + pow(z,-a);
    }
    logLik = -sum;
    //logPrior = ldnorm(m,*min * .8,1);
  }
  else {
    for (int i=0; i<n; i++){
      z = (x[i]-m)/s;
      sum += (1+a)*log(z) + pow(z,-a);
    }
    logLik = -n*log(s) -sum;
    //logPrior = ldgamma(s,sd(x,n),1);
  }

  return logPrior + logLik;
}

double rnorm(double mu, double sd){
  double ss = 2; 
  double u,v;
  while (ss > 1){
    u = runif()*2-1;
    v = runif()*2-1;
    ss = u*u+v*v;
  }
  double w = sqrt(-2*log(ss)/ss);

  return u*w * sd + mu;
}

void  mig(double* x,int* n1,int* N1,
          double* outA,double* outM,double* outS,
          double* init,double* cs){
  int n = *n1; int N = *N1;

  double min = Min(x,n);

  outA[0] = init[0];
  outM[0] = init[1];
  outS[0] = init[2];

  int cnta = 0; int cntm = 0; int cnts = 0;
  double csa = cs[0]; double csm = cs[1]; double css = cs[2];
  double cand,r,u;
  double* p = (double*) malloc(3*sizeof(double*)); 
  p[0] = outA[0];
  p[1] = outM[0];
  p[2] = outS[0];
  //printf("\n %s %f \n","MIN: ",*min);
  for (int i=1; i<N; i++){
    
    // NEED to fix my likelihoods!
    // Don't need to because my likelihoods behave poorly

    //UNIVERSAL UPDATES:
    outA[i] = outA[i-1]; outM[i] = outM[i-1]; outS[i] = outS[i-1];
    //p[0] = outA[i]; p[1] = outM[i]; p[2] = outS[i];


    cand = rnorm(outA[i], csa);  
    if (cand > 0){
      r = logPo(cand,x,n,min,p,0) - logPo(outA[i],x,n,min,p,0); 
      u = runif();
      if (r > log(u)){
        outA[i] = cand;
        p[0] = cand;
        cnta += 1;
      } 
    }

    cand = rnorm(outM[i], csm);
    if (cand < min){ // my new m's cannot be greater than min(x)
      r = logPo(cand,x,n,min,p,1) - logPo(outM[i],x,n,min,p,1); 
      u = runif();
      if (r > log(u)){
        outM[i] = cand;
        p[1] = cand;
        cntm += 1;
      } 
    } 

    cand = rnorm(outS[i], css);  
    if (cand > 0){
      r = logPo(cand,x,n,min,p,2) - logPo(outS[i],x,n,min,p,2); 
      u = runif();
      if (r > log(u)){
        outS[i] = cand;
        p[2] = cand;
        cnts += 1;
      } 
    }
    
  } // End of Metropolis Loop

  printf("%s %f\t %f\t %f\n","Acceptance%: ",
          cnta*1.0/N,cntm*1.0/N,cnts*1.0/N);
}

void myEst(double* x, int* n, double* param){
  //need to find quantiles in c 
  double m = Min(x,*n); 
  double s = sd(x,*n); 
  double a = m+mean(x,*n); 
  param[0] = a; param[1] = m; param[2] = s;
}

void simME(int* N, int* n, double* param, double* M){
  double* x = (double*) malloc(*n*sizeof(double*));
  double* init = (double*) malloc(3*sizeof(double*));
  //printf("%s\t %s\t\t %s\t\t %s\n","MLE","a","m","s");
  for (int i=0; i<*N; i++){
    rfrechet(n,param,x);
    for (int j=0; j<3; j++) {init[j] = param[j];}
    myEst(x,n,init); 
    //printf("%d\t %f\t %f\t %f\n",i,init[0],init[1],init[2]);
    M[i*3+0] = init[0];
    M[i*3+1] = init[1];
    M[i*3+2] = init[2];
  }
  free(x);
  free(init);
}

void simMLE(int* N, int* n, double* param, double* M){
  double* x = (double*) malloc(*n*sizeof(double*));
  double* init = (double*) malloc(3*sizeof(double*));
  //printf("%s\t %s\t\t %s\t\t %s\n","MLE","a","m","s");
  for (int i=0; i<*N; i++){
    rfrechet(n,param,x);
    for (int j=0; j<3; j++) {init[j] = param[j];}
    newton(x,n,init); 
    //printf("%d\t %f\t %f\t %f\n",i,init[0],init[1],init[2]);
    M[i*3+0] = init[0];
    M[i*3+1] = init[1];
    M[i*3+2] = init[2];
  }
  free(x);
  free(init);
}


void simBE(int* N, int* n, double* param, double* cs, double* M){
  double* x = (double*) malloc(*n*sizeof(double*));
  double* init = (double*) malloc(3*sizeof(double*));
  int* B = (int*) malloc(sizeof(int*));
  *B = 10000;

  double* outA = (double*) malloc(*B*sizeof(double*));
  double* outM = (double*) malloc(*B*sizeof(double*));
  double* outS = (double*) malloc(*B*sizeof(double*));

  //printf("%s\t %s\t\t %s\t\t %s\n","MLE","a","m","s");
  for (int i=0; i<*N; i++){
    rfrechet(n,param,x);
    for (int j=0; j<3; j++) {init[j] = param[j];}
    mig(x,n,B,outA,outM,outS,init,cs); 
    //printf("%d\t %f\t %f\t %f\n",i,init[0],init[1],init[2]);
    M[i*3+0] = mean(outA,*B);
    M[i*3+1] = mean(outM,*B);
    M[i*3+2] = mean(outS,*B);
  }
  free(x);
  free(B);
  free(init);
  free(outA);
  free(outM);
  free(outS);
}

int main(){ 
  time_t t; srand((unsigned) time(&t));
  int* N        = (int*)    malloc(sizeof(int*)); 
  int* n        = (int*)    malloc(sizeof(int*)); 
  double* param = (double*) malloc(3*sizeof(double*)); 
  
  param[0] = 8;
  param[1] = 4;
  param[2] = 1;

  //simMLE(N,n,param);
  free(N);
  free(n);
  free(param);

  return 0;
}

// http://www.emeraldinsight.com/journals.htm?articleid=871493&show=html#0830160102006.png
