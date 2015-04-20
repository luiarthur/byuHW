*Arthur Lui;

data mud;
  infile 'mud.dat';
  input mDate $ tank tmt total baet lest other;
        lchiro=log(baet);
        week = 0;
  if mDate = '6Jul10' then week=1;
  if mDate = '13Jul10' then week=2; 
  if mDate = '20Jul10' then week=3; 
  if mDate = '27Jul10' then week=4; 
  if mDate = '3Aug10' then week=5; 
  if mDate = '10Aug10' then week=6; 
  if mDate = '17Aug10' then week=7; 
  if mDate = '24Aug10' then week=8;
run;
*proc print;

*proc glimmix;
*   class week tmt;
*   model lchiro = week tmt total;
*   random intercept / subject=tmt; 
*run;

proc glimmix data=mud;
  week2 = week;
  class tank tmt;
  model lchiro =  week tmt*week / s;
  random week2 / subject=tank(tmt) type=rsmooth;
  random intercept / subject=tank(tmt);
  output out=mudresult pred(blup)=pred2;
run;

proc sgpanel data=mudresult;
  panelby tmt;
  series x=week y=pred2 / group=tank;
run;
