library(dplyr)
library(doParallel)
library(infotheo)
library(ggplot2)
source("final_funs.R")

####################################
#I. Dispersal-delta space ----------
####################################
reps <- 300
par_space <- expand.grid(disp_rate = seq(from = 0, to = 0.4, length.out = 5),
                         delta = seq(from = 0, to = 0.4, length.out = 5),
                         aij = exp(1),
                         p = 0.001,
                         cost_c = 0,
                         shuffle_p = 0)%>%
  slice(rep(1:n(), each = reps))%>%
  tibble::rowid_to_column(var = 'ID')%>%
  mutate(ID = as.character(ID))

####################################
#II. Initial condition and constants ----
#####################################
#Constant parameters
t <- 2000 #2000
t_steps <- 50 # record every 20 steps
L <- 64 #64 #Landscape size (64 x 64)
K <- 50 #carrying capacity
no.spp <- 2
r <- rep(1.5, no.spp)
N0 <- array(10, dim = c(L, L, no.spp)) #initial condition

####################################
#Set up clusters -----------
####################################
nodes <- 2
cores <- 35
n_cores <- nodes*cores
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
                metrics = TRUE, record_every = t_steps)

}


stopCluster(cl)
####################################
####################################



####################################
#Export results
####################################
#Combine
result <- lapply(sim_out, as.data.frame)%>%
  bind_rows(.id = 'ID')%>%
  left_join(par_space, by = 'ID')

#Export
saveRDS(result, file = 'stationary-states-28FEB25.rds')
