library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(ggplot2)

Rcpp::sourceCpp("src/v1.1/simulation.cpp")

source('MP_function.R')
source('R_function.R')

source('R/CreateData.R')
Create_OM()

#Modify for 3Ps 

load('data/m2021_assessment2021C.RData')
A <- m$tmb.data$A
Y <- m$tmb.data$Y
n_avg <- 3
avg_ind <- seq(Y-(n_avg-1),Y)
#Rec
num <- length(m$tmb.data$years)
sr.data <- list(
  years=m$tmb.data$years[1:(num-2)],
  recruits=m$rep$N_matrix[1,3:num],
  ssb=m$rep$ssb[1:(num-2)]
) %>% 
  as.data.frame() %>% 
  mutate(yearsR=years+2,colg=ifelse(yearsR>=1993,2,1))
sr.data.small <- sr.data %>% filter(yearsR>=1993)

bh <- log(recruits)~log((exp(log_a)*ssb)/(1+ exp(log_b)*ssb)) + 0
bhfit <- nls(bh,
             data=sr.data.small,
             start=list(log_a = log(max(sr.data.small$recruits)/max(sr.data.small$ssb)),
                        log_b = log(2/max(sr.data.small$ssb))),
             lower = list(log_a = -Inf, log_b = -Inf),
             upper = list(log_a = Inf, log_b = Inf),
             algorithm = 'port'
)
BH_parms <- as.vector( coef(bhfit) )
BH_cov <- as.matrix( vcov(bhfit) )
bherr <- sigma(bhfit)

#OM List
OM$Model_Name = 'Hybrid'
OM$Mortality_Name = 'Condition'
OM$Recruitment_Name = 'Test'
OM$N1 = m$rep$N_matrix[,Y]
OM$Stock_wt = stock_wt = rowSums(m$tmb.data$weight[,avg_ind])/n_avg
OM$Catch_wt = rowSums(m$tmb.data$Cw[,avg_ind])/n_avg
OM$Maturity = rowSums(m$tmb.data$mat[,avg_ind])/n_avg
OM$F_Mu = rep(0, A)
OM$Selectivity = rowSums(m$rep$F[,avg_ind])/n_avg
OM$terminalF = m$rep$F[,Y]
OM$terminalM = m$rep$M[,Y]
OM$terminalYield = m$rep$pred_landings[Y]
OM$PE_Std = exp(m$opt$par[names(m$opt$par) == 'log_std_s'])
OM$Fbar_ind = c(3:7) #ages 5-9
OM$Mbar_ind = c(3:7) #ages 5-9
OM$SR_lag = 0
OM$Rparameters = BH_parms
OM$Rstd = exp(m$opt$par['log_std_log_R'])
OM$Rparm_type = 'use_CoV' #or use_Kernel, or find_CoV
OM$RCoV = BH_cov
OM$RKernel = NULL
OM$Nrec = NULL
OM$F_sigma = m$rep$sigmaF2
OM$SR_function = 'Rec_function'

OM$Minfo$A = m$tmb.data$A
OM$Minfo$Ages = c(2:12)
OM$Minfo$base = rep(.3, A)
OM$Minfo$std = summary(m$sd.rep)%>%as_tibble(rownames="parm")%>%
  rename(estimate="Estimate",SE="Std. Error") %>% as.data.frame() %>% 
  filter(parm=="mpar") %>% dplyr::select(SE) %>% unlist() %>% 
  `[`((m$tmb.data$imest$keyM[1:A]+1)[1:A])
OM$Minfo$condition = 

#MP List
OM$ymp = 1

## Test Examples ##

load("data/OMassessment_MPno_fishing_Mtermyr_RSigBH_MSElite.Rdata")

obj <- Simulate(OM, MP, 10000, 40)

obj$All_summary %>% 
  filter(year > 2022 & year < 2048) %>%
  filter(var %in% c('Fbar', 'SSB', 'TAC', 'Rec')) %>% 
  dplyr::select(year, var, mp, p50, p05, p10, p90, p95) %>%
  rbind(MSElite$all_summary %>% 
          filter(year > 2022 & year < 2048) %>%
          filter(var  %in% c('Fbar', 'SSB', 'TAC', 'Recruits')) %>% 
          mutate(var = ifelse(var == 'Recruits', 'Rec', var)) %>% 
          dplyr::select(year, var, mp, p50, p5, p10, p90, p95) %>% 
          rename('p05' = p5)) %>% 
  mutate(mp = factor(mp, levels = c('no_fishing', 'No fishing'),
                     labels = c('R', 'C++'))) %>% 
  ggplot(aes(x = year, y = p50))+
  geom_line(aes(col = mp))+
  geom_ribbon(aes(ymax = p90, ymin = p10, fill = mp), alpha = .2)+
  guides(col = 'none', fill = guide_legend(title = 'Model'))+
  facet_wrap(~factor(var, levels = c('Fbar', 'SSB', 'TAC', 'Rec')), 
             ncol = 2, nrow = 2, scales = 'free_y')+
  labs(x = "Year", y = 'Median (80% CIs)')+
  theme_bw()

#Critical NDF

OM <- Create_OM('assessment', 'termyr', 'SigBH')
MP <- Create_MP('Critical 100')
load("data/OMassessment_MPprecautionary_Mtermyr_RSigBH_MSElite.Rdata")

obj <- Simulate(OM, MP, 10000, 40)

obj$All_summary %>% 
  filter(year > 2022 & year < 2048) %>%
  filter(var %in% c('Fbar', 'SSB', 'TAC', 'Rec')) %>% 
  dplyr::select(year, var, mp, p50, p05, p10, p90, p95) %>%
  rbind(MSElite$all_summary %>% 
          filter(year > 2022 & year < 2048) %>%
          filter(var  %in% c('Fbar', 'SSB', 'TAC', 'Recruits')) %>% 
          mutate(var = ifelse(var == 'Recruits', 'Rec', var)) %>% 
          dplyr::select(year, var, mp, p50, p5, p10, p90, p95) %>% 
          rename('p05' = p5)) %>% 
  mutate(mp = factor(mp, levels = c('precautionary', 'Critical 100'),
                     labels = c('R', 'C++'))) %>% 
  ggplot(aes(x = year, y = p50))+
  geom_line(aes(col = mp))+
  geom_ribbon(aes(ymax = p90, ymin = p10, fill = mp), alpha = .2)+
  guides(col = 'none', fill = guide_legend(title = 'Model'))+
  facet_wrap(~factor(var, levels = c('Fbar', 'SSB', 'TAC', 'Rec')), 
             ncol = 2, nrow = 2, scales = 'free_y')+
  labs(x = "Year", y = 'Median (80% CIs)')+
  theme_bw()