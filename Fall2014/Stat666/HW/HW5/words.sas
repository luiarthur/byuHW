/*Q6:*/

data words;
  infile 'words.dat';
  input y1 y2 y3 y4 y5;
run;
*proc print data=words;

title 'PROC CALIS--Words 14.4a';
proc calis method=ml;* maxiter=5000 maxfunc=5000;
  lineqs
    y1 =       f1 + e1,
    y2 = lam21 f1 + e2,
    y3 = lam31 f1 + e3,
    y4 = lam41 f1 + e4,
    y5 = lam51 f1 + e5;
  std
    e1-e5 = psi1-psi5,
    f1 = phi1;
  bounds
    0 <= psi1-psi5,
    0 <= phi1;
run;

/*
6a) Is the model a good fit?
      CFI = .9280 Bentler Comparative Fit Index (want: > .9)             * Good fit
    RMSEA = .1778 (want: < .06)                                            Poor fit
     SRMR = 0.0711: Standardized Root Mean Square Residual (want: < .08) * Good fit
*/
/*
6b) Significant lambda's?
    |z_ij| > 1.96, for all lambda_ij. So, model should not be simplified.
*/

title 'PROC CALIS--Words 14.4c';
proc calis method=ml;* maxiter=5000 maxfunc=5000;
  lineqs
    y1 = lam11 f1 + e1,
    y2 =       f1 + e2,
    y3 = lam31 f1 + e3,
    y4 = lam41 f1 + e4,
    y5 = lam51 f1 + e5;
  std
    e1-e5 = psi1-psi5,
    f1 = phi1;
  bounds
    0 <= psi1-psi5,
    0 <= phi1;
run;
/*
6c) Is the model a good fit?
      CFI = .9280 Bentler Comparative Fit Index (want: > .9)             * Good fit
    RMSEA = .1778 (want: < .06)                                            Poor fit
     SRMR = 0.0711: Standardized Root Mean Square Residual (want: < .08) * Good fit
*/
