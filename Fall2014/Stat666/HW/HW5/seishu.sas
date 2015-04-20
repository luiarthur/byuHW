DATA words;
  INFILE 'words.dat';
  INPUT y1 y2 y3 y4 y5;
RUN;
*proc print data=words;

/*Q5:*/

data seishu;
  infile 'seishu.dat';
  input taste odor ph acid1 acid2 sake drsugar totsugar alcohol fn;
run;

title '#4 Seishu: PROC CALIS--Evaluating a specific hypothesis about the structure';
proc calis method=ml cov maxiter=5000 maxfunc=5000;
  lineqs
    taste    =        f1                                     + e01,
    ph       =                    f2                         + e03,
    totsugar =                                f3             + e08,
    alcohol  =                                            f4 + e09,
    odor     = lam021 f1 + lam022 f2 + lam023 f3 + lam024 f4 + e02,
    acid1    = lam041 f1 + lam042 f2 + lam043 f3 + lam044 f4 + e04,
    acid2    = lam051 f1 + lam052 f2 + lam053 f3 + lam054 f4 + e05,
    sake     = lam061 f1 + lam062 f2 + lam063 f3 + lam064 f4 + e06,
    drsugar  = lam071 f1 + lam072 f2 + lam073 f3 + lam074 f4 + e07,
    fn       = lam101 f1 + lam102 f2 + lam103 f3 + lam104 f4 + e10;
  std
    e01-e10 = psi01-psi10,
    f1-f4   = phi1-phi4;
  cov
    f1-f4= phi11-phi16;
  bounds
    0 <= phi1-phi4,
    0 <= psi01-psi10;
run;

/*
Absolute Index:
    Chi^2 = 14.8718
       df = 2
  p-value = .2485

5a) Is the model a good fit?
      CFI = .9824: Bentler Comparative Fit Index (want: > .9)            * Good fit
    RMSEA = .0908 (want: < .06)                                            Poor fit
     SRMR = 0.0308: Standardized Root Mean Square Residual (want: < .08) * Good fit
*/


title '#4 Seishu: Set lam022, 023, 041, 051, 053, 101, and 103 
      (all that had |t|<1) to 0 due to non-significant t-value';
proc calis method=ml cov maxiter=5000 maxfunc=5000;
  *Removed Lambda's with |t| < 1;
  lineqs
    taste    =        f1                                     + e01,
    ph       =                    f2                         + e03,
    totsugar =                                f3             + e08,
    alcohol  =                                            f4 + e09,
    odor     = lam021 f1 +                         lam024 f4 + e02,
    acid1    =             lam042 f2 + lam043 f3 + lam044 f4 + e04,
    acid2    =             lam052 f2 +             lam054 f4 + e05,
    sake     = lam061 f1 + lam062 f2 + lam063 f3 + lam064 f4 + e06,
    drsugar  = lam071 f1 + lam072 f2 + lam073 f3 + lam074 f4 + e07,
    fn       =             lam102 f2 +             lam104 f4 + e10;
  std
    e01-e10 = psi01-psi10,
    f1-f4   = phi1-phi4;
  cov
    f1-f4= phi11-phi16;
  bounds
    0 <= phi1-phi4,
    0 <= psi01-psi10;
run;

title '#4 Seishu: Additionally, set lam024, 042, 043, 061, 064, 071, 072, and 074 
      (all that had |t|<2) to 0 due to non-significant t-value';
proc calis method=ml cov maxiter=5000 maxfunc=5000;
  lineqs
     taste    =        f1                                     + e01,
     odor     = lam021 f1                                     + e02,
     ph       =                    f2                         + e03,
     acid1    =                                   + lam044 f4 + e04,
     acid2    =             lam052 f2             + lam054 f4 + e05,
     sake     =             lam062 f2 + lam063 f3             + e06,
     drsugar  =                         lam073 f3             + e07,
     totsugar =                                f3             + e08,
     alcohol  =                                            f4 + e09,
     fn       =             lam102 f2             + lam104 f4 + e10;
  std
    e01-e10= psi01-psi10,
    f1-f4= phi1-phi4;
  cov
    f1-f4= phi11-phi16;
  bounds
    0 <= phi1-phi4,
    0 <= psi01-psi10;
run;

