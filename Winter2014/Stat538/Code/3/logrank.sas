options nocenter nonumber nodate;
*Note: the herco variable was created: hercoc=1 or 3 then herco=1 else herco=hercoc;

*Getting the data into SAS.;
data melanoma34;
  set 'c:\Documents and Settings\Administrator\Desktop\538\melanoma34';
run;

*Looking at the first 10 obs of the data.;
proc print data=melanoma34 (obs=10);
run;

*K-M curve;
proc lifetest data=melanoma34 plots=(s);
  time SurTime*SurStatus(0);
run;

* logrank test;
proc lifetest data=melanoma34 plots=(s);
  time SurTime*SurStatus(0);
  strata Stage;
run;

quit;
