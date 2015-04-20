rapply <- function(L,f) { # L is a list, f is a function to apply to L[[x]]
                          # L apply takes a list, applies a function, 
                          # and rbinds it. Assumes output is vector.
  n <- length(L)
  temp <- f(L[[1]])
  out <- matrix(0,n,length(temp))
  out[1,] <- temp

  for (i in 2:n) {
    out[i,] <- f(L[[i]]) 
  }
  out
}
