data carbon;
  set 'carbon.sas7bdat';
run;  
proc print;run;

*1;
proc MIXED data=carbon method=ml;
  where day <> 18;
  class T;
  model logwt = day T day*T/ ddfm=kr s;
  random intercept day / subject = ID type=vc;
run;

*2;
proc MIXED data=carbon method=ml;
  where day <> 18;
  class T;
  model logwt = day T day*T/ ddfm=kr s;
  random intercept / subject = ID type=vc;
run;

*3;
proc MIXED data=carbon method=ml;
  where day <> 18;
  class T;
  model logwt = day T day*T/ ddfm=kr s;
  random day / subject = ID type=vc;
run;

*4;
proc MIXED data=carbon method=ml;
  where day <> 18;
  class T;
  model logwt = day T day*T/ ddfm=kr s;
  random intercept day day*T/ subject = ID type=vc;
run;

