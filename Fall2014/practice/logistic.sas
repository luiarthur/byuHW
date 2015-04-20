data sim;
  infile 'sim.dat';
  input y x1-x4;
run;

proc print;

proc logistic data=sim;
  model y(event='1') = x1-x4;
  *default is modelling P(y=0|X). y(event='1') models P(y=1|X);
run;  
