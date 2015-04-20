/*
data prim;
  infile 'primate.dat' delimiter=',' firstobs=2;
  input monkey $ treatment $ wk2 wk4 wk8 wk12 wk16;
run;

data prim2;
  set prim;
  drop w2--w16;
  w=w2; week=2; output;
  w=w4; week=4; output;
  w=w8; week=8; output;
  w=w12; week=12; output;
  w=w16; week=16; output;
run;
*/

data prim;
  infile 'prim2.txt' firstobs=2;
  input monkey $ treatment $ week score;
run;

data corn;
  set 'sasDat/corn.sas7bdat';
run;

data carbon;
  set 'sasDat/carbon.sas7bdat';
run;

proc sgpanel data=prim;
  title 'Primates';
  panelby treatment;
  series y=score x=week/group=monkey;
run;

proc sgplot data=corn;
  title 'Corn';
  series y=yld x=square/group=var;
run;

proc sgpanel data=carbon;
  title 'Carbon';
  panelby T;
  series y=wt x=day/group=ID;
run;

proc print data=carbon; run;
