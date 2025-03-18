#Landsacpe portraits
library(dplyr)
library(doParallel)
library(infotheo)
library(ggplot2)
source("final_funs.R")

####################################
#Coexistence treatment space
####################################
p_length <- 5 #length of parameter values
par_space <- expand.grid(disp_rate = seq(from = 0, to = 0.4, length.out = p_length),
                         delta = seq(from = 0, to = 0.4, length.out = p_length),
                         aij = exp(1),
                         p = 0.001,
                         cost_c = 0,
                         shuffle_p = 0)%>%
  tibble::rowid_to_column(var = 'ID')%>%
  mutate(ID = as.character(ID))


####################################
#II. Initial condition and constants ----
#####################################
#Constant parameters
t <- 1500
t_steps <- 100 #time steps to record
L <- 64 #Landscape size (64 x 64)
K <- 50 #carrying capacity
no.spp <- 2
r <- rep(1.5, no.spp)
N0 <- array(10, dim = c(L, L, no.spp)) #initial condition


####################################
#Set up clusters -----------
####################################
n_cores <- 25
cl <- makeCluster(n_cores, type = "FORK")
registerDoParallel(cl)

####################################
#Simulate in parallel
####################################
#Simulate
sim_out <- foreach(i = 1:nrow(par_space)) %dopar%{
  set.seed(i);
  metaautomaton(N = N0, r = r, aij = par_space$aij[i], K = K,
                d_disp = par_space$disp_rate[i], delta = par_space$delta[i],
                cost_c = par_space$cost_c[i], ext_p = par_space$p[i],
                shuffle_p = par_space$shuffle_p[i], time = t,
                weight_fitness = FALSE, weight_overlap = TRUE,
                metrics = FALSE, record_every = t_steps)%>%
    reshape2::melt()%>%
    setNames(nm = c('row','col', 'Spp', 'time','N'))%>%
    mutate(time = time*t_steps)%>%
    mutate(Spp = plyr::mapvalues(Spp, from = c(1,2), to = c('Sp1','Sp2')))%>%
    tidyr::pivot_wider(names_from = 'Spp', values_from = "N")%>%
    mutate(State = case_when(Sp1 > 0 & Sp2 == 0 ~ "S1",
                             Sp1 > 0 & Sp2 >0 ~ "S1 + S2",
                             Sp1 == 0 & Sp2 > 0 ~ "S2",
                             TRUE ~ 'Extinct'))%>%
    mutate(State = factor(State, levels = c('S1','S1 + S2', 'S2', "Extinct")))%>%
    mutate(disp_rate = par_space$disp_rate[i],
           delta = par_space$delta[i])

}


stopCluster(cl)
####################################
####################################


####################################
#Export results
####################################
saveRDS(sim_out, file = 'landscape-portraits-10MAR25.rds')





