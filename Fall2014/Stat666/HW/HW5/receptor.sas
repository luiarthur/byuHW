data receptor;
  infile 'receptor.dat';
  input Ethylene Propane nButane x2Mp x3Mp Benzene Cyclohexane x2Mh x224Tmp Acetylene Ethane;
run;


title '#2a (Receptor Data): Principal Factor Method with Varimax Rotation';
proc factor method=principal 
            priors=smc   /* COMMENT: smc tells Proc Factor to use the squared multiple 
                            correlation with all the other vars as the prior 
                      communality estimate. */
            mineigen=1
            rotate=promax   
            reorder 
            plot;
  var Ethylene Propane nButane x2Mp x3Mp Benzene Cyclohexane x2Mh x224Tmp Acetylene Ethane;
run;
/*
      Factor Pattern

                    Factor1

x224Tmp             0.98705
x2Mh                0.98486
Cyclohexane         0.98187
Benzene             0.98073
x2Mp                0.95607
Ethylene            0.94238
x3Mp                0.91900
Acetylene           0.91491
nButane             0.90806
Propane             0.70691
Ethane              0.58698
*/


title '#2b (Receptor Data): Maximum Likelihood Method (k=1 factor)';
proc factor method=ml
            heywood   /* sets to 1 any communality greater than var{x_i}=1  */
      n=1;   /* use 1 factor */
  var Ethylene Propane nButane x2Mp x3Mp Benzene Cyclohexane x2Mh x224Tmp Acetylene Ethane;
run;


title '#2b (Receptor Data): Maximum Likelihood Method (k=2 factor)';
proc factor method=ml
            heywood   /* sets to 1 any communality greater than var{x_i}=1  */
      n=2;   /* use 2 factors */
  var Ethylene Propane nButane x2Mp x3Mp Benzene Cyclohexane x2Mh x224Tmp Acetylene Ethane;
run;


title '#2b (Receptor Data): Maximum Likelihood Method (k=3 factor)';
proc factor method=ml
            heywood   /* sets to 1 any communality greater than var{x_i}=1  */
      n=3;   /* use 3 factors */
  var Ethylene Propane nButane x2Mp x3Mp Benzene Cyclohexane x2Mh x224Tmp Acetylene Ethane;
run;


title '#2b (Receptor Data): Maximum Likelihood Method (k=2 factor) with Promax Rotation';
proc factor method=ml
            heywood 
      n=2   /* use 2 factors */
            rotate = promax
            plot;
  var Ethylene Propane nButane x2Mp x3Mp Benzene Cyclohexane x2Mh x224Tmp Acetylene Ethane;
run;

