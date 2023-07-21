#Example custom recruitment function
Rec_Function <- function(SSB, parms){
  
  exp(parm[1] * SSB) / (1 + exp(parm[2] * SSB))
  
}