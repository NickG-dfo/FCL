#The MP function can be defined using any number of variables you want
#However, all values and parameters must be called as a single vector argument

MP_function <- function(inputs){
  
  if(inputs[1] == 1){
    return( (inputs[2]-30)*.1 * (inputs[2] > 30) )
  }else if(inputs[1] == 2){
    return( exp(inputs[2]) / (1+exp(inputs[2])) )
  }
  
}