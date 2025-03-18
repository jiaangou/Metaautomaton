#Coxistence time phase space -----------
library(dplyr)
library(doParallel)
library(infotheo)
library(ggplot2)
#source("complete_functions.R")
source("final_funs.R")

####################################
#I. Simulation parameters  -------------
####################################
#Parameter space -------
par_range <- 20
replicates <- 15

par_space <- expand.grid(disp_rate = seq(from = 0, to = 0.4, length.out = par_range),
                         delta = seq(from = 0, to = 0.4, length.out = par_range),
                         aij  = exp(seq(from = -1, to = 1, by = 0.5)), #expoenntiatl sequence that extends to aij < 1
                         #aij = seq(from = 1, to = 2, length.out = 3),
                         p = c(0.001, 0.1),
                         cost_c = c(0, 3),
                         shuffle_p = c(0, 1))%>%
  slice(rep(1:n(), each = replicates))


#split up aij levels into separte simulations
par_space <- par_space%>%
  filter(!aij < 1)

# ggplot(aes(x = disp_rate, y = delta))+
#   geom_point()+
#   facet_grid(p~aij)



####################################
#II. Initial condition and constants ----
#####################################
#Constant parameters
t <- 10000 #timesteps (10K)
L <- 10 #Landscape size (10)
K <- 50 #carrying capacity
no.spp <- 2
r <- rep(1.5, no.spp)
#initial condition
N0 <- array(10, dim = c(L, L, no.spp))


####################################
#Set up clusters -----------
####################################
nodes <- 4
cores <- 30
n_cores <- nodes*cores
cl <- makeCluster(n_cores, type = "FORK")
registerDoParallel(cl)


####################################
#Simulate in parallel
####################################

#Simulate
#sim_out <- foreach(i = 1:nrow(par_space)) %dopar%{
sim_out <- foreach(i = 1:nrow(par_space)) %dopar%{
  set.seed(i);
  metaautomaton(N = N0, r = r, aij = par_space$aij[i], K = K,
                d_disp = par_space$disp_rate[i], delta = par_space$delta[i],
                cost_c = par_space$cost_c[i], ext_p = par_space$p[i],
                shuffle_p = par_space$shuffle_p[i], time = t,
                weight_fitness = FALSE, weight_overlap = FALSE,
                metrics = TRUE, record_every = NULL)
}

stopCluster(cl)


#Export
coexist_dat <- par_space%>%
  bind_cols(., bind_rows(sim_out))%>%
  mutate(r_freq = r_N1 / (r_N1  + r_N2))%>%
  mutate(total_N = r_N1  + r_N2)


saveRDS(coexist_dat, 'ASSsim_28FEB25.rds')


#############################################
#============  Make figure   ===============
#############################################

#############################################

#Average over replicates
coexist_summ <- coexist_dat%>%
  group_by(disp_rate, delta, aij, p, cost_c, shuffle_p)%>%
  summarise_all(list(mean), na.rm = TRUE)%>%
  mutate(aij = round(aij, 3))%>% #round aij so its less clunky
  mutate_at(c('aij', 'p','cost_c', 'shuffle_p'), as.factor)

#Export summary
saveRDS(coexist_summ, 'ASSsim_summary_28FEB25.rds')


#############################################
#Plot labels -----------
# Disturbance
levels(coexist_summ$p) <- c(expression("Disturbance:~'0.1%'"),
                            expression("Disturbance:~'10%'"))
# Dispersal cost
levels(coexist_summ$cost_c) <- c("Cost: 0", "Cost: 3")
#Aij
levels(coexist_summ$aij) <- paste0('alpha[ij] ~ ":" ~ ', round(unique(par_space$aij), 3))

#Shuffle
levels(coexist_summ$shuffle_p) <- c(expression("Shuffle:~'0'"),
                                    expression("Shuffle:~'1'"))

#############################################
#Plot 1: Coexistnece time  -------------
# #coexist_summ%>%
# coexist_summ%>%
#   ggplot(aes(x = disp_rate, y = delta))+
#   geom_tile(col = 'black',fill = 'white', size = 0.5)+
#   scale_fill_gradient(low='#bdbdbd', high="#252525", name = 'Coexistence time')+
#   facet_grid(p+cost_c ~ aij+shuffle_p,labeller = label_parsed)+
#   labs(x = "Dispersal rate (d)", y = expression("Dispersal assortativity ("~delta~")"))+
#   theme_bw()


co_plot <- coexist_summ%>%
  ggplot(aes(x = disp_rate, y = delta))+
  geom_tile(aes(fill = time), col = 'black', size = 0.5)+
  scale_fill_gradient(trans = "log",low = 'grey', high = 'black', breaks = c(500, 1000, 2000, 4000, 8000), name = "Coexistence time")+
  labs(x = "Dispersal rate (d)", y = expression("Dispersal assortativity ("~delta~")"))+
  facet_grid(p+cost_c ~ aij+shuffle_p,labeller = label_parsed)+
  theme_bw()


#Save plot
ggsave('ass_28FEB25.png', plot = co_plot, width = 16, height = 8)



#############################################
#Plot 2: Spatial structure effect plot  -------------
# #spat_struc_plot <-
#   coexist_summ%>%
#   select(time)
#   mutate(shuffle_p = ifelse(shuffle_p == 0, 'no_shuffle', 'shuffle')%>%as.factor)%>%
#   tidyr::pivot_wider(names_from = shuffle_p, values_from = time)
#   mutate(shuffle_effect = no_shuffle - shuffle)
#   ggplot(aes(x = disp_rate, y = delta))+
#   geom_tile(aes(fill = shuffle_effect))+
#   scale_fill_gradient2(low = '#f03b20', high = '#1c9099', mid = 0, name = "Effect of spatial structure")+
#   facet_grid(aij ~ p + cost_c,
#              labeller = labeller(cost_c = cost_labs, p = disturb_labs, aij = aij_labs))+
#   labs(x = "Dispersal rate (d)", y = expression("Dispersal assortativity ("~delta~")"))+
#   theme_bw()
#
# #Save plot
# ggsave('spat_str_12MAR24.png', plot = spat_struc_plot, width = 15, height = 12)
#
#

