options nocenter nonumber nodate;
*Note: the herco variable was created: hercoc=1 or 3 then herco=1 else herco=hercoc;

*Getting the data into SAS.;
*data nurse; 
*  infile 'nursing.csv' dlm=',' firstobs=2;
*  input 
*run;

filename nursing '/data/arthurll/School/Winter2014/Stat538/Homework/4/nursing.csv'; *termstr=cr;
proc import datafile=nursing dbms=csv out=nursing1 replace;
run;

*Looking at the first 10 obs of the data.;
proc print data=nursing1 (obs=10);
run;

*Cox model;
proc phreg data=nursing1;
  class race(ref="1");
  model duration*completion(0) = race poverty smoked education;
run;

quit;
