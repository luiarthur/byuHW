/*Q6:*/

data football;
  infile 'football.dat';
  input gp wdim circum fbeye eyehd earhd jaw;
run;
*proc print data=football;

title 'PROC CALIS--Football 14.6a';
proc calis method=ml cov maxiter=5000 maxfunc=5000;
  lineqs
    circum =       f1            + e1,
    eyehd  =                  f2 + e2,
    fbeye  = lam31 f1 + lam32 f2 + e3,
    wdim   = lam41 f1 + lam42 f2 + e4,
    earhd  = lam51 f1 + lam52 f2 + e5,
    jaw    = lam61 f1 + lam62 f2 + e6;
  std
    e1-e6 = psi1-psi6,
    f1-f2 = phi1-phi2;
  cov
    f1-f2 = phi12;
  bounds
    0 <= psi1-psi6,
    0 <= phi1-phi2;
run;

/*
6a) Is the model a good fit?
      CFI = .9401 Bentler Comparative Fit Index (want: > .9)             * Good fit
    RMSEA = .1690 (want: < .06)                                            Poor fit
     SRMR = 0.0630: Standardized Root Mean Square Residual (want: < .08) * Good fit
*/

title 'PROC CALIS--Football 14.6b';
proc calis method=ml cov maxiter=5000 maxfunc=5000;
  lineqs
    circum =       f1            + e1,
    eyehd  =                  f2 + e2,
    fbeye  = lam31 f1            + e3,
    wdim   = lam41 f1            + e4,
    earhd  =            lam52 f2 + e5,
    jaw    = lam61 f1            + e6;
  std
    e1-e6 = psi1-psi6,
    f1-f2 = phi1-phi2;
  cov
    f1-f2 = phi12;
  bounds
    0 <= psi1-psi6,
    0 <= phi1-phi2;
run;

/*
6b) Is the model a good fit?
    Chi^2 = 34.4404
       df = 9
      CFI = .8806 Bentler Comparative Fit Index (want: > .9)             * Good fit
    RMSEA = .1778 (want: < .06)                                            Poor fit
     SRMR = 0.0874: Standardized Root Mean Square Residual (want: < .08) * Good fit
*/
*All the statistics are worse for the simple model, so model a is the better model;
