\documentclass{article}
\usepackage{fullpage}
\usepackage{amssymb}
\usepackage{Sweave}
\usepackage{bm}
\usepackage{mathtools}
\def\wl{\par \vspace{\baselineskip}}

\begin{document}
\SweaveOpts{concordance=TRUE}
\title{Stat535 HW 8}
\author{Arthur Lui}
\maketitle

<<results=hide>>=
rm(list=ls())
options('scipen'=8)
options(continue=" ")
@

\subsection*{Data:}
<<design,echo=T>>=
#TWO FACTOR ANOVA
#DATA:
###################### Factor B ###
##############    0     25     50 #
# Factor A   0  mu11  mu12   mu13 #
#          100  mu21  mu22   mu23 #
#          150  mu31  mu32   mu33 #
###################################
y11 <- c(-.04,-.02,-.09,.05)
y12 <- c(.06,-.18,0.0,-.04)
y13 <- c(-.36,-.26,-.27)
y21 <- c(0,-.05,-.26,.1)
y22 <- c(-.21,-.06,-.16,-.1)
y23 <- c(-.13,-.03,-.19,-.3)
y31 <- c(-.06,-.17,-.07,.15)
y32 <- c(-.26,.05,-.11,-.09)
#y33 <- c()

# Create Y
y <- list(y11,y12,y13,y21,y22,y23,y31,y32)
Y <- matrix(unlist(y))
n <- nrow(Y)

# Create X
k <- length(y)
W <- matrix(0,nrow(Y),k)
b <- 1; e <- length(y[[1]])
j <- 1
while(j <= k){
  W[b:e, j] <- 1
  j <- j + 1
  b <- e + 1
  e <- ifelse(j<=k,e + length(y[[j]]),0)
}
@
\newpage
\noindent
\textbf{Q1:} Estimate $\bm{\hat{\mu}}$ and provide a description of each element
             of $\hat{\mu}$ for the pharmaceutical collaborator.
<<design,echo=T>>=
ans1 <- mu <- solve(t(W) %*% W) %*% t(W) %*% Y #same as: lapply(y,mean)
@
\[ \hat{\mu}= \begin{pmatrix}
                 \hat{\mu}_{11}\\
                 \hat{\mu}_{12}\\
                 \hat{\mu}_{13}\\
                 \hat{\mu}_{21}\\
                 \hat{\mu}_{22}\\
                 \hat{\mu}_{23}\\
                 \hat{\mu}_{31}\\
                 \hat{\mu}_{32}\\
               \end{pmatrix} =
<<label=tab1,echo=F,results=tex>>=
  library(xtable)
  x=xtable(mu,align=rep("",ncol(mu)+1),digits=4)
  print(x,floating=FALSE,tabular.environment="pmatrix",hline.after=NULL,
        include.rownames=FALSE,include.colnames=FALSE)
@
\]
\wl\noindent
$\hat{\mu}_{11}=$The estimated difference in VAS (post-pre) for the patients 
                 given the placebo treatment.\wl\noindent
$\hat{\mu}_{12}=$The estimated difference in VAS (post-pre) for the patients 
                 given ONLY 25mg of treatment B.\wl\noindent
$\hat{\mu}_{13}=$The estimated difference in VAS (post-pre) for the patients
                 given ONLY 50mg of treatment B.\wl\noindent
$\hat{\mu}_{21}=$The estimated difference in VAS (post-pre) for the patients
                 given ONLY 100mg of treatment A.\wl\noindent
$\hat{\mu}_{22}=$The estimated difference in VAS (post-pre) for the patients
                 given the 100mg of A and 25mg of B treatment.\wl\noindent
$\hat{\mu}_{23}=$The estimated difference in VAS (post-pre) for the patients
                 given the 100mg of A and 50mg of B treatment.\wl\noindent
$\hat{\mu}_{31}=$The estimated difference in VAS (post-pre) for the patients 
                 given the 150mg of A and 25mg of B treatment.\wl\noindent
$\hat{\mu}_{32}=$The estimated difference in VAS (post-pre) for the patients
                 given the 150mg of A and 50mg of B treatment.
\newpage
\textbf{Q2:} Compute the test statistic and p-value for $H_0: \mu_{ij}$ 
             are all equal. Provide an explaination of what this test 
             corresponds to for the phamaceutical collaborator.
<<design,echo=T>>=
F.test <- function(C){
  q <- nrow(C)
  SSH <- t(C%*%mu) %*% solve(C %*% solve(t(W)%*%W) %*% t(C)) %*% (C%*%mu)
  SSE <- t(Y-W%*%mu) %*% (Y-W%*%mu)
  F.stat <- (SSH/q) / (SSE/(n-k))
  p.val  <- pf(F.stat,q,n-k,lower.tail=FALSE)
  out <- data.frame(F.stat=F.stat,p.val=p.val)
}
C2 <- cbind(diag(k-1),0) + cbind(0,-1*diag(k-1))
ans2 <- F.test(C2)
@

\[ C = 
<<label=tab1,echo=F,results=tex>>=
  library(xtable)
  x=xtable(C2,align=rep("",ncol(C2)+1),digits=0)
  print(x,floating=FALSE,tabular.environment="pmatrix",hline.after=NULL,
        include.rownames=FALSE,include.colnames=FALSE)
@
\]
\wl\noindent
\textbf{The F-statistic is: \Sexpr{ans2$F.stat} (p-value = \Sexpr{ans2$p.val})}
\wl\noindent
\textbf{We have tested whether the effects of each tested drug 
        combinations are the same. That is we have tested whether any one of
        the drugs have a different effect on VAS than that of the others. Since 
        the p-value obtained from the test is greater than 0.05, we do not have
        significant evidance to conclude that any one drug combination of drugs 
        and dosages have a different effect on VAS than that of others, and we 
        conclude that the drugs all have the same effect on VAS.}
\wl\wl\noindent 
\textbf{Q3:} The collaborator asks for a test of whether or not Drug A 
             is effective in combatting pain. Define a hypothesis of the 
             form $H_0: \bm{C}\mu=0$ and compute the test statistic and p-value.
             Explain to another statistician why the C matrix in your 
             $H_0$ is a valid test of the effectiveness of Drug A in 
             combatting pain.
<<design,echo=T>>=
#C3 <- matrix(c(-1,0,0,1,0,0,0,0,
#               -1,0,0,0,0,0,1,0),nrow=2,byrow=TRUE) What I had before
C3 <- matrix(c(1/3,1/3,1/3,-1/3,-1/3,-1/3,0, 0,
               1/3,1/3,1/3,0,0,0,-1/2,-1/2),nrow=2,byrow=TRUE)
ans3 <- F.test(C3)
@
$H_0: C\mu=0, where$
\[ C = 
<<label=tab1,echo=F,results=tex>>=
  library(xtable)
  x=xtable(C3,align=rep("",ncol(C3)+1),digits=2)
  print(x,floating=FALSE,tabular.environment="pmatrix",hline.after=NULL,
        include.rownames=FALSE,include.colnames=FALSE)
@
\]
\wl\noindent
\textbf{The F-statistic is: \Sexpr{ans3$F.stat} (p-value = \Sexpr{ans3$p.val})}
\wl\noindent
\textbf{The C matrix I have chosen is constructed such that I can
        compare the average effect of drug A 
        at dose 0mg (not administered) to the average effect of drug A 
        at dose 100mg, and the average effect of drug A at dose 0mg
        (not administered) to the average effect of drug A at dose 150mg.
        If increasing the dosage of A does not significantly alter the 
        change in VAS, we can conclude that Drug A is not effective 
        in combatting pain.}
\wl\wl\noindent
\textbf{Q4:} The collaborator asks for a test of whether or not Drug B 
             is effective in combatting pain. Define a hypothesis of the
             form $H_0: \bm{C}\mu=0$ and compute the test statistic and p-value. 
             Explain to another statistician why the C matrix in your 
             $H_0$ is a valid test of the effectiveness of Drug B in 
             combatting pain.
<<design,echo=T>>=
#C4 <- matrix(c(-1,1/2,1/2,0,0,0,0,0),nrow=1) I decided NOT to do this
C4 <- matrix(c(-1/3,1/3,0,-1/3,1/3,0,-1/3,1/3,
               -1/3,0,1/2,-1/3,0,1/2,-1/3,0),nrow=2,byrow=TRUE)
ans4 <- F.test(C4)
@
$H_0: C\mu=0, where$
\[ C = 
<<label=tab1,echo=F,results=tex>>=
  library(xtable)
  x=xtable(C4,align=rep("",ncol(C4)+1),digits=2)
  print(x,floating=FALSE,tabular.environment="pmatrix",hline.after=NULL,
        include.rownames=FALSE,include.colnames=FALSE)
@
\]
\wl\noindent
\textbf{The F-statistic is: \Sexpr{ans4$F.stat} (p-value = \Sexpr{ans4$p.val})}
\wl\noindent
\textbf{The C matrix I have chosen is constructed such that I can
        compare the average effect of drug B 
        at dose 0mg (not administered) to the average effect of drug B 
        at dose 25mg, and the average effect of drug B at dose 0mg
        (not administered) to the average effect of drug B at dose 50mg.
        If increasing the dosage of B does not significantly alter the change 
        in VAS, we can conclude that Drug B is not effective in combatting pain.}
\wl\wl\noindent
\textbf{Q5:} The collaborator asks for a test of whether or not there is 
             an interaction between Drug A and Drug B. Define a hypothesis
             of the form $H_0: \bm{C}\mu=0$ and compute the test statistic and 
             p-value. Explain to another statistician why the C matrix in 
             your $H_0$ is a valid test of the Drug A - Drug B interaction.
<<design,echo=T>>=
#C5 <- matrix(c(-1/5,-1/5,-1/5,-1/5,1/3,1/3,-1/5,1/3),nrow=1,byrow=T)
C5 <- matrix(c(1,-1,0,-1,1,0,0,0,
               0,1,-1,0,-1,1,0,0,
               0,0,0,1,-1,0,-1,1),nrow=3,byrow=TRUE)
ans5 <- F.test(C5)
@
\[ C = 
<<label=tab1,echo=F,results=tex>>=
  library(xtable)
  x=xtable(C5,align=rep("",ncol(C5)+1),digits=3)
  print(x,floating=FALSE,tabular.environment="pmatrix",hline.after=NULL,
        include.rownames=FALSE,include.colnames=FALSE)
@
\]
\wl\noindent
\textbf{The F-statistic is: \Sexpr{ans5$F.stat} (p-value = \Sexpr{ans5$p.val})}
\wl\noindent
\textbf{The C matrix is constructed such that I can compare the difference
        between the change in VAS for each dosage of A across different 
        dosages of B. If changing the dosage of B while holding the dosage of A 
        constant does not alter the change in VAS significantly, we can conclude 
        that there is not an interaction between Drug A and Drug B.}
\end{document}
