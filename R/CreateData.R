
Create_OM <- function(){

  OM_data <- list(
    #Integers
    'SR_lag' = NULL,
    'ymp' = NULL,
    #Numeric Atomics
    'terminalYield' = NULL,
    'Rstd' = NULL,
    'Max_Rec' = NULL,
    #Strings
    'Model_Name' = NULL,
    'Mortality_Name' = NULL,
    'Recruitment_Name' = NULL,
    'Rparm_type' = NULL,
    #Numeric Vectors
    'Catch_wt' = NULL,
    'Stock_wt' = NULL,
    'Maturity' = NULL,
    'N1' = NULL,
    'Rparameters' = NULL,
    'F_Mu' = NULL,
    'Selectivity' = NULL,
    'terminalF' = NULL,
    'terminalM' = NULL,
    'PE_Std' = NULL,
    #Integer Vectors
    'Fbar_ind' = NULL,
    'Mbar_ind' = NULL,
    #Numeric Matrices
    'RCoV' = NULL,
    'RKernel' = NULL,
    'Nrec' = NULL,
    #Numeric Vector or Matrix
    'F_sigma' = NULL,
    #Function
    'SR_function' = NULL
  )
  
  OM_data$Minfo <- list(
    #Integers
    'A' = NULL,
    #Numeric Vectors
    'Ages' = NULL,
    'base' = NULL,
    'std' = NULL,
    #Numeric Atomics
    'condition' = NULL, #size Y
    'age_effect' = NULL, #size A
    'year_effect' = NULL, #size Y
    #Numeric Matrix
    'correlates' = NULL #3x3 AR matrix
  )
  
  OM_data$Minfo$error_type <- c(
    #Logicals
    'base' = T,
    'condition' = F,
    'age_effect' = F,
    'year_effect' = F,
    'correlates' = F
  )
  
  OM_data$Rec_settings <- c(
    'Rparm_error' = T, #logical
    'R_error' = T, #logical
    'bias_correction' = F, #logical
    'Kernel_seed' = F #logical
  )
  
  OM_data$F_settings <- c(
    'F_error' = T, #logical
    'F_experror' = F #logical
  )
  
  OM_data$Mort_settings <- c(
    'condition' = T, #logical
    'age' = F, #logical
    'year' = F, #logical
    'AR' = F #logical
  )
  
  MP_data <- list(
    
    'Name' = NULL, #string
    'F_type' = NULL, #bool -- if True MP outputs F, otherwise TAC
    #'y0' = NULL, #int
    'delay' = NULL, #int
    'Cmin' = NULL, #double
    'TAC_Error' = NULL, #double
    'TAC0' = NULL, #double
    'Function' = NULL, #string -- name of function for MP in R
    'MP_inputs' = NULL #this can be any type, but must be the only argument to Function
                #e.g. list, vector, numeric, string, etc.
    
  )
  
  OM <<- OM_data
  MP <<- MP_data
  
}


  