library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(ggplot2)

source('R/CreateData.R')

Rcpp::sourceCpp("src/v1.1/simulation.cpp")

load('data/m2021_assessment2021C.RData')

Create_OM <- function(om, mproj, rproj){
  
  A <- m$tmb.data$A
  Y <- m$tmb.data$Y
  n_avg <- 3
  avg_ind <- seq(Y-(n_avg-1),Y)
  
  mlowYr <- seq(1996,2005)
  mhighYr <- seq(2015,2021)
  
  mest <- summary(m$sd.rep)%>%as_tibble(rownames="parm")%>%
          rename(estimate="Estimate",SE="Std. Error") %>% as.data.frame() %>% 
          filter(parm=="mpar") %>% dplyr::select(estimate) %>% unlist()
  msd <- summary(m$sd.rep)%>%as_tibble(rownames="parm")%>%
          rename(estimate="Estimate",SE="Std. Error") %>% as.data.frame() %>% 
          filter(parm=="mpar") %>% dplyr::select(SE) %>% unlist()
  
  OM <- list(
        name = om,
        mproj = mproj,
        rproj = rproj,
        A = A,
        years = m$tmb.data$years,
        N1 = m$rep$N_matrix[,Y],
        N2020 =m$rep$N_matrix[,Y-1],
        std_pe = exp(m$opt$par[names(m$opt$par) == 'log_std_s']),
        logR = log(mean(m$rep$N_matrix[1,avg_ind])),
        log_std_log_r = m$opt$par['log_std_log_R'],
        Favg = rowSums(m$rep$F[,avg_ind])/n_avg,
        stock_wt = rowSums(m$tmb.data$weight[,avg_ind])/n_avg,
        catch_wt = rowSums(m$tmb.data$Cw[,avg_ind])/n_avg,
        maturity = rowSums(m$tmb.data$mat[,avg_ind])/n_avg,
        yieldterminal = m$rep$pred_landings[Y],
        Fterminal = m$rep$F[,Y],
        Mterminal = m$rep$M[,Y],
        F_err = m$rep$sigmaF2,
        keyF = m$tmb.data$keyF[1:A],
        mpar_est = mest,
        m_ind = m$tmb.data$imest$keyM[1:A]+1,
        mscale = m$tmb.data$imest$mscale[,Y]
  )
  
  OM$mscale_low <- apply(m$tmb.data$imest$mscale[, which(m$tmb.data$years%in%mlowYr)], 1, mean)
  OM$mscale_high <- apply(m$tmb.data$imest$mscale[, which(m$tmb.data$years%in%mhighYr)], 1, mean)
  OM$mpar_sd <- msd[OM$m_ind[1:A]]
  
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
  OM$BH_parms <- as.vector( coef(bhfit) )
  OM$BH_cov <- as.matrix( vcov(bhfit) )
  OM$bherr <- sigma(bhfit)
  
  OM$max_S = max(sr.data$ssb)
  OM$sigmoid_kernel <- as.matrix( data.table::fread('data/mcmc_sample_assessment.csv') )
  OM$bhsigerr <- median(OM$sigmoid_kernel[,4])
  
  return(OM)
  
}

Create_MP <- function(mp){
 
  MP <- list(
        # f = 0.035, 
        # fixedcatch = .8*1.346,
        name = mp,
        TAC = c(1.346,1.346),
        Fbound = 0.1,
        threshold=26.4,
        harvestrate=0.05,
        ssbrbracket = c(.4,.5,.6,.8,1,1.2,1.4,1.6,1.8,2),
        fbracket = c(.025,.03,0.03,.03,.05,.1,.125,.15,.18,.22,.26)
  )
  MP$y_hcr <- length(MP$TAC)
  
  return(MP)
  
}


## Test Examples ##

OM <- Create_OM('assessment', 'termyr', 'SigBH')
MP <- Create_MP('No fishing')
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