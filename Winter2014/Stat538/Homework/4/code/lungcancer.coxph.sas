options nocenter nonumber nodate;
*Note: the herco variable was created: hercoc=1 or 3 then herco=1 else herco=hercoc;

*Getting the data into SAS.;
data lungcancer;
  set 'c:\Documents and Settings\Administrator\Desktop\538\lungcancer';
run;

*Looking at the first 10 obs of the data.;
proc print data=lungcancer (obs=10);
run;

*Cox model;
proc phreg data=lungcancer;
  model time*status(0) = sex ph_ecog ph_karno pat_karno wt_loss;
run;

quit;
