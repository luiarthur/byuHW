options nocenter nonumber nodate;
*Note: the herco variable was created: hercoc=1 or 3 then herco=1 else herco=hercoc;

*Getting the data into SAS.;
data uis;
  set 'c:\survival\uis_small';
run;

*Looking at the first 10 obs of the data.;
proc print data=uis (obs=10);
run;

*K-M curve;
proc lifetest data=uis plots=(s);
  time time*censor(0);
run;

quit;
