#3d

lam.to.F <- function(lam,vh,ve) {
  (ve-1)/vh * (1-sqrt(lam))/sqrt(lam)
}

F.no.sv <- lam.to.F(.6402317,6,46)
F.no.v  <- lam.to.F(.9336179,2,46)
F.no.s  <- lam.to.F(.122354,3,46)

F.no.sv
F.no.v
F.no.s

1-pf(F.no.sv,2*6,2*(46-1))
1-pf(F.no.v,2*2,2*(46-1))
1-pf(F.no.s,2*3,2*(46-1))

#3e)

Ed <- sqrt(diag(c(11.896,14.404,13.656,7245.6))) # diagonal elements of E
a.s <- c(.20895364,-.05317724,.17493015,-.00429399)
a.v <- c(.27,-.0267,.0556,.00170)
a.sv <- c(.2827,-.04609,.02655,.00025815)
ve <- 48

a.s.star  <- 1/sqrt(ve) * Ed %*% a.s
a.v.star  <- 1/sqrt(ve) * Ed %*% a.v
a.sv.star <- 1/sqrt(ve) * Ed %*% a.sv

a.s.star
a.v.star
a.sv.star


which(abs(a.s.star) > max(a.s.star)/2)
which(abs(a.v.star) > max(a.v.star)/2)
which(abs(a.sv.star) > max(a.sv.star)/2)


