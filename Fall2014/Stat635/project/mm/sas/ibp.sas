data sim;
  infile 'mmDat.txt';
  input y x1 Z;
run;

proc MIXED data=sim;
  class z;
  model y = x1 / s;
  random intercept / subject = z G V SOLUTION;
run;
