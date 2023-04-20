library(tidyverse)
library(Rcpp)

temp_func <- function(x){
  return(x^2+1)
}

Rcpp::sourceCpp("src/Testing/Func_test.cpp")

cppF_using_RF(seq(1,1000,1))