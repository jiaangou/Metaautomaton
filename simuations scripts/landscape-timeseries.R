library(dplyr)
library(doParallel)
library(infotheo)
library(ggplot2)
source("final_funs.R")

####################################
#I. Dispersal-delta space ----------
####################################
par_space <- expand.grid(disp_rate = seq(from = 0, to = 0.4, length.out = 5),
                         delta = seq(from = 0, to = 0.4, length.out = 5),
                         aij = exp(1),
                         p = c(0.001, 0.1),
                         cost_c = c(0,3))%>%
  tibble::rowid_to_column(var = 'ID')%>%
  mutate(ID = as.character(ID))

####################################
#II. Initial condition and constants ----
#####################################
#Constant parameters
t <- 7000
t_steps <- 20 #time steps to record
L <- 64 #Landscape size (64 x 64)
K <- 50 #carrying capacity
no.spp <- 2
r <- rep(1.5, no.spp)
N0 <- array(10, dim = c(L, L, no.spp)) #initial condition

####################################
#Set up clusters -----------
####################################
n_nodes <- 5
n_cores <- 20 * n_nodes #100 cores
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
                shuffle_p = 0, time = t,
                weight_fitness = FALSE, weight_overlap = FALSE,
                metrics = TRUE, record_every = t_steps)
}

stopCluster(cl)
####################################
####################################



####################################
#Export results
####################################
#Combine
#sim_out
result <- lapply(sim_out, as.data.frame)%>%
  bind_rows(.id = 'ID')%>%
  left_join(par_space, by = 'ID')


#Check visual
result%>%
  ggplot(aes(x = time, y = cooccurrence))+
  geom_path(aes(col = delta))+
  scale_fill_gradient(high = '#3182bd', low = '#c6dbef', name = expression(delta))+
  facet_grid(disp_rate~p+cost_c)


#Export
saveRDS(result, file = 'macroscopic-timeseries-28FEB25.rds')




####################################
#Plot results --------
####################################
#
# landscape_ts <- here::here('simulated data/macroscopic-timeseries-28FEB25.rds')%>%
#   readRDS()
#
# #Plot labels -----------
# # Disturbance
# disturb_labs <- c("Disturbance: 0.1%", "Disturbance: 10%")
# names(disturb_labs) <- c(0.001, 0.1)
# # Dispersal cost
# cost_labs <- c("Cost: 0", "Cost: 1")
# names(cost_labs) <- c(0, 1)
#
# #NMI
# nmi_gg <- landscape_ts%>%
#   mutate(nmi = I/Hx)%>%
#   #mutate(t_scaled = time / 10^3)%>%
#   mutate(d = plyr::mapvalues(disp_rate, from = seq(from = 0, to = 0.4, length.out = 5),
#                              to = c('d = 0','d = 0.1','d = 0.2','d = 0.3','d = 0.4')))%>%
#   ggplot(aes(x = time, y = nmi, group = ID))+
#   #geom_point()+
#   geom_path(aes(col = delta))+
#   scale_color_gradient(high = '#3182bd', low = '#c6dbef', name = expression(delta))+
#   scale_x_continuous(breaks = c(0, 4000, 8000))+
#   labs(x = 'Time', y = 'Normalized mutual information (%)')+
#   facet_grid(p+cost_c
#              ~d, labeller = labeller(cost_c = cost_labs, p = disturb_labs))+
#   theme_bw()+
#   theme(strip.text = element_text(size = 6), axis.text = element_text(size = 7),
#         axis.title = element_text(size = 8))
#
#
# #Overlap
# overlap_gg <- landscape_ts%>%
#   mutate(d = plyr::mapvalues(disp_rate, from = seq(from = 0, to = 0.4, length.out = 5),
#                              to = c('d = 0','d = 0.1','d = 0.2','d = 0.3','d = 0.4')))%>%
#   ggplot(aes(x = time, y = overlap, group = ID))+
#   #geom_point()+
#   geom_path(aes(col = delta))+
#   scale_color_gradient(high = '#3182bd', low = '#c6dbef', name = expression(delta))+
#   scale_x_continuous(breaks = c(0, 4000, 8000))+
#   labs(x = 'Time', y = 'Overlap')+
#   facet_grid(p+cost_c~d, labeller = labeller(cost_c = cost_labs, p = disturb_labs))+
#   theme_bw()+
#   theme(strip.text = element_text(size = 6), axis.text = element_text(size = 7),
#         axis.title = element_text(size = 8))
#
#
# #co-occurrence
# co_gg <- landscape_ts%>%
#   mutate(d = plyr::mapvalues(disp_rate, from = seq(from = 0, to = 0.4, length.out = 5),
#                              to = c('d = 0','d = 0.1','d = 0.2','d = 0.3','d = 0.4')))%>%
#   ggplot(aes(x = time, y = cooccurrence, group = ID))+
#   #geom_point()+
#   geom_path(aes(col = delta))+
#   labs(x = 'Time', y = 'Co-occurrence (%)')+
#   scale_color_gradient(high = '#3182bd', low = '#c6dbef', name = expression(delta))+
#   scale_x_continuous(breaks = c(0, 4000, 8000))+
#   facet_grid(p+cost_c~d, labeller = labeller(cost_c = cost_labs, p = disturb_labs))+
#   theme_bw()+
#   theme(strip.text = element_text(size = 6), axis.text = element_text(size = 7),
#         axis.title = element_text(size = 8))
#
#
# #Patch
# library(patchwork)
# nmi_gg + overlap_gg + co_gg + plot_layout(guides = 'collect') +
#   plot_annotation(tag_levels = 'A') &  theme(plot.tag = element_text(size = 20))
#


