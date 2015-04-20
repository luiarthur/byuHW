options nocenter nonumber nodate;
*Note: the herco variable was created: hercoc=1 or 3 then herco=1 else herco=hercoc;

data mel;
  infile '/data/arthurll/School/Winter2014/Stat538/Homework/2/melanoma.csv' dlm=',' firstobs=2;
  input Sur Age Gender InitStge $ Trtmnt Rmssn rcens SurvTime survcens;
run;


proc print data=mel (obs=30);
run;


*K-M curve;
ods graphics on / imagename="remission";
  proc lifetest data=mel plots=(s);
    time Rmssn*rcens(0);
    strata Trtmnt;
  run;
ods graphics off;

ods graphics on / imagename="survival";
  proc lifetest data=mel plots=(s);
    time SurvTime*survcens(0);
    strata Trtmnt;
  run;
ods graphics off;

quit;
