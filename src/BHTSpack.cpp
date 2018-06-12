// Copyright (C) 2017, 2018 by Ivo D. Shterev
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <iostream>

using namespace std;

extern "C"{
SEXP hat_pai(SEXP z, SEXP ph1, SEXP ph0, SEXP mu1, SEXP mu0, SEXP sigma1, SEXP sigma0, SEXP pai, SEXP N, SEXP HH)
{
  int M = length(N);
  int H = INTEGER(coerceVector(HH, INTSXP))[0];
  double xpai = REAL(coerceVector(pai, REALSXP))[0];

  double *xz = REAL(PROTECT(coerceVector(z, REALSXP)));

  double *xmu1 = REAL(PROTECT(coerceVector(mu1, REALSXP)));
  double *xmu0 = REAL(PROTECT(coerceVector(mu0, REALSXP)));

  double *xsigma1 = REAL(PROTECT(coerceVector(sigma1, REALSXP)));
  double *xsigma0 = REAL(PROTECT(coerceVector(sigma0, REALSXP)));

  double *xph1 = REAL(PROTECT(coerceVector(ph1, REALSXP)));
  double *xph0 = REAL(PROTECT(coerceVector(ph0, REALSXP)));

  int *xN = INTEGER(PROTECT(coerceVector(N, INTSXP)));

  SEXP hatpai = PROTECT(allocVector(REALSXP, length(z)));
  double *xhatpai = REAL(hatpai);

  int it = 0;
  for (int m = 0; m < M; m++){
    for (int n = 0; n < xN[m]; n++){
      double valz = xz[it+n];
      double val0 = 0.0; 
      double val1 = 0.0;

      for (int h = 0; h < H; h++){
       val0 += xph0[m*H+h] * exp(-0.5*(valz-xmu0[m*H+h])*(valz-xmu0[m*H+h])/xsigma0[m*H+h]) / sqrt(xsigma0[m*H+h]);
       val1 += xph1[m*H+h] * exp(-0.5*(valz-xmu1[m*H+h])*(valz-xmu1[m*H+h])/xsigma1[m*H+h]) / sqrt(xsigma1[m*H+h]); 
      }

      // compute hatpai
      double r = val0/val1;

      if (R_FINITE(r))
        xhatpai[it+n] = 1.0 / (1.0+(1.0-xpai)/xpai*r);
      else
        xhatpai[it+n] = 0.5;
    }
    it += xN[m];
  }

  UNPROTECT(9);
  return hatpai;
}

SEXP stick_multnorm_z(SEXP z, SEXP ph, SEXP sigma, SEXP mu, SEXP N, SEXP HH)
{
  int M = length(N);
  int H = INTEGER(coerceVector(HH, INTSXP))[0];

  double *xz = REAL(PROTECT(coerceVector(z, REALSXP)));
  double *xmu = REAL(PROTECT(coerceVector(mu, REALSXP)));
  double *xsigma = REAL(PROTECT(coerceVector(sigma, REALSXP)));
  double *xph = REAL(PROTECT(coerceVector(ph, REALSXP)));
  int *xN = INTEGER(PROTECT(coerceVector(N, INTSXP)));

  SEXP pr = PROTECT(allocVector(REALSXP, length(z)*H));
  double *xpr = REAL(pr);

  for (int h = 0; h < H; h++){
    int it = 0;
    for (int m = 0; m < M; m++){
      for (int n = 0; n < xN[m]; n++){ 
        double val = (xz[it+n]-xmu[m*H+h]) * (xz[it+n]-xmu[m*H+h]) / xsigma[m*H+h];
    
        // compute pr=pdf*ph
        xpr[(it+n)*H+h] = xph[m*H+h]/sqrt(xsigma[m*H+h])/exp(0.5*val);
      }
      it += xN[m];
    }
  }

  // normalize pr
  int it = 0;
  for (int m = 0; m < M; m++){
    for (int n = 0; n < xN[m]; n++){
      double val = 0.0;
      for (int h = 0; h < H; h++)
        val += xpr[(it+n)*H+h];
      if (val != 0.0){
        for (int h = 0; h < H; h++)
          xpr[(it+n)*H+h] /= val;
      }
      else{
        for (int h = 0; h < H; h++)
          xpr[(it+n)*H+h] = 1.0 / H;
      }
    }
    it += xN[m];
  }

  UNPROTECT(6);
  return pr;
} 

SEXP stick_multnorm_h(SEXP z, SEXP ih, SEXP pk, SEXP sigma, SEXP mu, SEXP N, SEXP HH)
{
  double *xz = REAL(PROTECT(coerceVector(z, REALSXP)));
  double *xmu = REAL(PROTECT(coerceVector(mu, REALSXP)));
  double *xsigma = REAL(PROTECT(coerceVector(sigma, REALSXP)));
  double *xpk = REAL(PROTECT(coerceVector(pk, REALSXP)));
  int *xN = INTEGER(PROTECT(coerceVector(N, INTSXP)));
  int *xih = INTEGER(PROTECT(coerceVector(ih, INTSXP)));

  int M = length(N);
  int K = length(pk);
  int H = INTEGER(coerceVector(HH, INTSXP))[0];
  int MKH = M * K * H;

  SEXP pr = PROTECT(allocVector(REALSXP, MKH));
  double *xpr = REAL(pr);

  for (int mkh = 0; mkh < MKH; mkh++)
    xpr[mkh] = 0.0;

  int it = 0;
  for (int m = 0; m < M; m++){
    for (int h = 0; h < H; h++){
      for (int k = 0; k < K; k++){
        double val = 0.0;
        int c = 0;
        for (int n = 0; n < xN[m]; n++){ 
          if (xih[it+n]==h){
            val -= 0.5*(xz[it+n]-xmu[k])*(xz[it+n]-xmu[k])/xsigma[k] + 0.5*log(xsigma[k]);
            c++;
          }
        }
        if (c != 0){
          xpr[m*H*K+h*K+k] = xpk[k] * exp(val);
        }
      }
    }
    it += xN[m];
  }

  // normalize pr
  for (int m = 0; m < M; m++){
    for (int h = 0; h < H; h++){
      // normalize
      double val = 0.0;
      for (int k = 0; k < K; k++)
        val += xpr[m*H*K+h*K+k];
      if (val != 0.0){
        for (int k = 0; k < K; k++)
          xpr[m*H*K+h*K+k] /= val;
      }
      else{
        for (int k = 0; k < K; k++)
          xpr[m*H*K+h*K+k] = 1.0 / K;
      }
    }
  }

  UNPROTECT(7);
  return pr;
} 

void multinomind(const double *probs, const int *NN, const int *KK, int *res){
  int N = *NN;
  int K = *KK;
  double *prob;      
  prob = (double*)malloc(K*sizeof(double));
  int *rn;      
  rn = (int*)malloc(K*sizeof(int));

  for (int n=0; n<N; n++){
    for (int k=0; k<K; k++)
      prob[k] = probs[n*K+k];

    rmultinom(1, prob, K, rn);

    int val = 0;
    for (int k=0; k<K; k++){
      if (rn[k] == 1)
        val = k+1;
    }

    res[n] = val;
  }

  // cleanup
  free(prob);
  free(rn);
}

void abfun (const int *ind, const int *NN, const double *tau, const int *HH, double *nu){
  int N = *NN;
  int H = *HH;

  for (int h=0; h<H; h++){
    double vala = 1.0;
    double valb = *tau;
    
    for (int n=0; n<N; n++){
      if (ind[n] == h)
        vala += 1.0;
      if (ind[n] > h)
        valb += 1.0;
    }

    nu[h] = rbeta(vala, valb);
  }
} 

void lambda (const double *nu, const int *HH, double *ph){  
  int H = *HH;
  
  for (int h=0; h<H-1; h++){
    double val = nu[h];
    for (int g=0; g<h; g++)
      val *= 1.0 - nu[g];
    
    ph[h] = val;
  }

  double val = 1.0;
  for (int h=0; h<H-1; h++)
    val -= ph[h];
  ph[H-1] = max(val, 0.0);
}
}
