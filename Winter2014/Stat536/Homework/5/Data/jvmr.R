install.packages("jvmr")
library(jvmr)

scala <- scalaInterpreter()
interpret(scala,"x") <- matrix(1:30,5,6)
matrix(interpret(scala,"x"),5,6)
