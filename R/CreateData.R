OM_data <- list(
  #Integers
  'SR_lag',
  'y0',
  'ymp',
  #Numeric Atomics
  'PE_std',
  'terminalYield',
  'Min_Catch',
  'TAC_Error',
  'Max_Rec',
  #Strings
  'Model_Name',
  'Mortality_Name',
  'Recruitment_Name',
  'Rparm_type',
  #Numeric Vectors
  'Catch_wt',
  'Stock_wt',
  'Maturity',
  'N1',
  'Rparameters',
  'F_Mu',
  'Selectivity',
  'terminalF',
  'terminalM',
  'Mpar_est',
  'Mpar_sd',
  'PE_Std',
  'TAC0',
  #Integer Vectors
  'Fbar_ind',
  'Mbar_ind',
  #Numeric Matrices
  'RCoV',
  'RKernel',
  'Nrec',
  #Numeric Vector or Matrix
  'F_sigma'
)

OM_data$Minfo <- list(
  #Integers
  'A',
  #Numeric Vectors
  'Ages',
  'base',
  'std',
  #Numeric Atomics
  'condition', #size Y
  'age_effect', #size A
  'year_effect', #size Y
  #Numeric Matrix
  'correlates', #AR matrix
)

OM_data$Minfor$error_type <- c(
  #Logicals
  'base' = T,
  'condition' = F,
  'age_effect' = F,
  'year_effect' = F,
  'correlates' = F
)

OM_data$Rec_settings <- c(
  'Rparm_error', #logical
  'R_error', #logical
  'bias_correction', #logical
  'RKseed', #logical
)

OM_data$F_settings <- c(
  'F_error', #logical
  'F_experror' #logical
)

OM_data$Mort_settings <- c(
  'condition', #logical
  'age', #logical
  'year', #logical
  'AR' #logical
)