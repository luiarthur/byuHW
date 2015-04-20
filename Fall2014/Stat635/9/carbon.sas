data carbon;
  set 'carbon.sas7bdat';
run;  
proc print;run;

/*   variance - correlation table */
proc glm data=carbon;
  where day=4;
  class T;
  by T;
  model l4--l16 = /nouni; *nothing on RHS => anova;
  manova h=intercept/printe;
run;


