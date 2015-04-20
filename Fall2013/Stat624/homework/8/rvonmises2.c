// rvonmises.c

#include <stdio.h>             // standard input/output
#include <stdlib.h>            // malloc
#include <math.h>              // fabs, sqrt, etc.
#include <time.h>              // time
#include <unistd.h>            // getpid
#include <gsl/gsl_rng.h>       // GNU Scientific Library
#include <gsl/gsl_cdf.h>       // GNU Scientific Library
#include <gsl/gsl_randist.h>   // GNU Scientific Library

#define pi 3.14159265358979323846

double runif(){
  return (double)rand()/(double)RAND_MAX;
}

// Returns the log of dvonmises, without normalizing constant
double logDvonmises(double* x, double mu, double nu, double kappa1, 
                    double kappa2, double lambda){

       double phi = x[0];
       double psi = x[1];

       return
         kappa1*cos(phi-mu) + kappa2*cos(psi-nu) +
         lambda*sin(phi-mu)*sin(psi-nu); 
}

void rvonmises(int* n, double* mu, double* nu, 
               double* k1, double* k2, double* l,
               double* phi, double* psi, int* PRINT){
  
  double mode[] = {0.0,0.0};
  
  if (*k1 * *k2 > *l * *l){ // UNIMODAL
    mode[0] = *mu; mode[1] = *nu;
  }
  else { // BIMODAL
    double phi0 = acos(*k1 / fabs(*l) * sqrt((*l * *l + *k2 * *k2) / 
                  (*l * *l + *k1 * *k1))); 
    double psi0 = acos(*k2 / fabs(*l) * sqrt((*l * *l + *k1 * *k1) / 
                  (*l * *l + * k2 * *k2)));
    
    if (l>0){
      mode[0] = phi0 + *mu; mode[1] = psi0 + *nu;
    }
    else{
      mode[0] = -phi0 + *mu; mode[1] = psi0 + *nu;
    }
    double nmode[] = {-mode[0],-mode[1]};
    if (logDvonmises( mode,*mu,*nu,*k1,*k2,*l) < 
        logDvonmises(nmode,*mu,*nu,*k1,*k2,*l)){
      mode[0] = nmode[0]; mode[1] = nmode[1];
    }  
  }

  double logF = logDvonmises(mode,*mu,*nu,*k1,*k2,*l);
  double logG = log(1/(4*pi*pi));
  double logAlpha = logG-logF;

  int i = 0;
  while(i<*n){
    double x[] = {runif()*2*pi-pi,runif()*2*pi-pi};
    logF = logDvonmises(x,*mu,*nu,*k1,*k2,*l);
    if( logF-logG+logAlpha > log(runif())){
      phi[i] = x[0]; psi[i] = x[1];
      if (*PRINT==1){
        printf("%f\t%f\n",x[0],x[1]);
      }  
      i++;
    }
  }

}

int main(int argc, char* argv[]){

   int    n   = atoi(argv[1]);
   double mu  = atof(argv[2]);
   double nu  = atof(argv[3]);
   double k1  = atof(argv[4]);
   double k2  = atof(argv[5]);
   double l   = atof(argv[6]);
   int    PRINT = 1;
   time_t t; srand((unsigned) time(&t));

   //double** x = rvonmises(n,mu,nu,k1,k2,l);  
   double* phi = (double*) malloc(n*sizeof(double*));
   double* psi = (double*) malloc(n*sizeof(double*));

   rvonmises(&n,&mu,&nu,&k1,&k2,&l,phi,psi,&PRINT);  

   free(phi);
   free(psi);
   
   return 0;
}

