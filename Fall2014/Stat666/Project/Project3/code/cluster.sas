DATA collin;
  INFILE 'cleanData.txt';
  INPUT text $
        frstPrsn inrThnk thnkPstv thnkNgtv thnkAhd thnkBck reason socTies
        drctAct intract notify lnrGuide wrdPict spceInt motion pastEvnt timeInt
        shftEvnt txtCvrg genre counter corpus corpGen;
RUN;  
*proc print;

*proc factor corr mineigin=1 scree preplot rotate=varimax reorder plot;
*  var frstPrsn inrThnk thnkPstv thnkNgtv thnkAhd thnkBck reason socTies
*  drctAct intract notify lnrGuide wrdPict spceInt motion pastEvnt timeInt
*  shftEvnt;
*run;  

PROC DISCRIM method=NORMAL pool=yes LIST crosslist manova;
   CLASS group;
   VAR wdim circum fbeye eyehd earhd jaw;
   TITLE 'LINEAR Classification Analysis';
RUN;
