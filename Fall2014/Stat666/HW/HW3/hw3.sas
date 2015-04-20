DATA snap;
  INFILE '/home/luiarthur/666HW3/snap.dat';
  INPUT S V rep x1 x2 x3 x4;
RUN;
PROC PRINT data=snap;

*3A;
PROC GLM;
  CLASS S V;
  MODEL x1 x2 x3 x4 = S | V / NOUNI;
  MANOVA H = S | V/ PRINTE PRINTH;
RUN;

*3B;
PROC GLM;
  CLASS S V;
  MODEL x1 x2 x3 x4 = S | V/ nouni;
  CONTRAST 'Comparison of Variety 2 to Varieties  and 3 means' V 1 -2 1;
  MANOVA H = V / PRINTE PRINTH;
RUN;

*3c;
PROC GLM;
  CLASS S V;
  MODEL x1 x2 x3 x4 = S | V/ nouni;
  CONTRAST 'Linear Sowdate Contrast' S -3 -1 1 3;
  CONTRAST 'Quadratic Sowdate Contrast' S 1 -1 -1 1;
  CONTRAST 'Cubic Sowdate Contrast' S 1 -3 3 -1;
  MANOVA H = S / PRINTE PRINTH;
RUN;

*3d;
PROC GLM; /*Full Model*/
  CLASS S V;
  MODEL x1 x2 x3 x4 = S | V /NOUNI;
  MANOVA H = S | V / PRINTE PRINTH;
RUN;

PROC GLM; /*Reduced Model*/
  CLASS S V;
  MODEL x1 x2 = S | V /NOUNI;
  MANOVA H = S | V / PRINTE PRINTH;
RUN;

*3f;
PROC GLM;
  CLASS S V;
  MODEL x1 x2 x3 x4 = S | V /nouni;
  MANOVA H = V / PRINTE PRINTH;
  MEANS V / BON ALPHA = .0125; * = .05 / 4;
RUN;




/*------------------------------------------------------------------*/
*4;
DATA mand;
  INFILE '/home/luiarthur/666HW3/mandible.dat';
  INPUT S G a1t1 a1t2 a1t3 a2t1 a2t2 a2t3 a3t1 a3t2 a3t3;
RUN;
PROC PRINT DATA = mand;

*4a);
PROC GLM ;
  CLASS G;
  MODEL a1t1 a1t2 a1t3 a2t1 a2t2 a2t3 a3t1 a3t2 a3t3 = G / NOUNI;
  repeated ac 3, ti 3 / PRINTE SUMMARY;
RUN;


*4b);
PROC GLM ;
  CLASS G;
  MODEL a1t1 a1t2 a1t3 a2t1 a2t2 a2t3 a3t1 a3t2 a3t3 = G / NOUNI;
  REPEATED ac 3, ti 3 / PRINTE SUMMARY;
  MANOVA H = G / PRINTE PRINTH;
RUN;


*4c);
PROC GLM ;
  CLASS G;
  MODEL a1t1 a1t2 a1t3 a2t1 a2t2 a2t3 a3t1 a3t2 a3t3 = G / NOUNI;
  REPEATED ac 3, ti 3 POLYNOMIAL / PRINTE SUMMARY;
  MANOVA H = G / PRINTE PRINTH;
RUN;




