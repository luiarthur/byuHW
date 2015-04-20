data cleft;
  infile 'cleft.dat';
  input tmt rep micenum cleftnum;
  lit = _N_;
run;
*proc print;

proc glimmix;
   class lit tmt;
   model cleftnum / micenum = tmt / dist=binomial;
   random intercept / subject=lit; 
run;
