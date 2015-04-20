data child;
  infile '../child.dat';
  input kid sex $ age len;
run;
*proc print;

proc mixed data=child method=ml;
  class kid sex;
  model len = sex sex*age/noint s;
  random intercept/subject=kid G; *type=un;
run;


proc mixed data=child method=ml;
  class kid sex;
  model len = sex sex*age/noint s;
  repeated /subject=kid; *type=ar(1);
run;
