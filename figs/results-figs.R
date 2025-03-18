library(here)
library(patchwork)


state_cols <- c('#1979a9', '#fed976', '#e07b39', 'black')

################################################################
#Figure 5: Landscape portraits -----------------------------------
################################################################

#Load data and set parameters
p_length <- 5
L <- 64
portraits_data <- here::here('R/data/landscape-portraits-10MAR25.rds')%>%
  readRDS()%>%
  bind_rows(.id = 'ID')


#Axis labels
gap_size <- 5 #this is the amount space between replicate landscapes
label_pos <- L/2 + ((L+gap_size)*(0:4)) #position of labels (halfway points)

x_label <- portraits_data$disp_rate%>%
  unique
y_label <- portraits_data$delta%>%
  unique
axis_x <- data.frame(x = label_pos, label = x_label)
axis_y <- data.frame(y = label_pos, label = y_label)

#Portrait data
portraits_df <- portraits_data%>%
  mutate(s_disp = plyr::mapvalues(disp_rate, from = x_label, to = 0:4))%>%
  mutate(s_delta = plyr::mapvalues(delta, from = y_label, to = 0:4))%>%
  mutate(scale_x = col + ((L+gap_size)*s_disp))%>% #rescale x values
  mutate(scale_y = row + ((L+gap_size)*s_delta)) #rescale y values


#Plot time steps separately
time_steps <- seq(from = 500, to = 1500, by = 500)


portrait_time <- vector(mode = 'list', length(time_steps))
for(t in 1:length(time_steps)){

  portrait_time[[t]] <- portraits_df%>%
    filter(time == time_steps[t])%>%
    ggplot(aes(x = scale_x, y = scale_y))+
    geom_tile(aes(fill = State))+
    scale_fill_manual(values = state_cols)+
    geom_text(data = axis_x, aes(x = x, y = -15, label = label), size = 3.5)+
    geom_text(data = axis_y, aes(x = -18, y = y, label = label), size = 3.5)+
    annotate("segment", x = -5, xend = (L*p_length) + 30, y = -5, yend = -5,
             arrow = arrow(length = unit(0.1, "cm")))+
    annotate("segment", x = -5, xend = -5,  y = -5, yend = (L*p_length) + 30,
             arrow = arrow(length = unit(0.1, "cm")))+
    annotate("text", x = -40, y = ((L+gap_size)*p_length)/2,
             label = 'paste("Dispersal assortativity ("~delta~")")', parse = TRUE,
             size = 3.5, angle = 90)+
    annotate("text", y = -40, x = ((L+gap_size)*p_length)/2,
             label = 'Dispersal rate (d)', size = 3.5)+
    ggtitle(label = paste0("Time: ", time_steps[t]))+
    theme_void()+
    theme(plot.title = element_text(hjust = 0.1))  # Adjust hjust as needed
}


#Combine and Export
portraits_p <- portrait_time[[1]] / portrait_time[[2]] / portrait_time[[3]]
fig5 <- portraits_p + plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'A') &  theme(plot.tag = element_text(size = 20))

ggsave(fig5, file = 'fig5.png',
       path = here::here("figs"),
       width = 5, height = 11, units = 'in')


################################################################
# Figure S3: Time-series -----------------------------------
################################################################
grad_pal <- c('#3182bd', '#bdd7e7')

#Load data
landscape_ts <- here::here('R/data/macroscopic-timeseries-28FEB25.rds')%>%
  readRDS()


# #Plot labels -----------
# Disturbance
disturb_labs <- c("Disturbance: 0.1%", "Disturbance: 10%")
names(disturb_labs) <- c(0.001, 0.1)
# Dispersal cost
cost_labs <- c("Cost: 0", "Cost: 3")
names(cost_labs) <- c(0, 3)

#NMI
nmi_gg <- landscape_ts%>%
  mutate(d = plyr::mapvalues(disp_rate, from = seq(from = 0, to = 0.4, length.out = 5),
                             to = c('d = 0','d = 0.1','d = 0.2','d = 0.3','d = 0.4')))%>%
  ggplot(aes(x = time, y = NMI, group = ID))+
  #geom_point()+
  geom_path(aes(col = delta))+
  scale_color_gradient(high = '#3182bd', low = '#c6dbef', name = expression(delta))+
  scale_x_continuous(breaks = c(0, 4000, 8000))+
  labs(x = 'Time', y = 'Normalized mutual information (%)')+
  facet_grid(p+cost_c
             ~d, labeller = labeller(cost_c = cost_labs, p = disturb_labs))+
  theme_bw()+
  theme(strip.text = element_text(size = 9), axis.text = element_text(size = 7),
        axis.title = element_text(size = 10))

#Overlap
overlap_gg <- landscape_ts%>%
  mutate(d = plyr::mapvalues(disp_rate, from = seq(from = 0, to = 0.4, length.out = 5),
                             to = c('d = 0','d = 0.1','d = 0.2','d = 0.3','d = 0.4')))%>%
  ggplot(aes(x = time, y = rho, group = ID))+
  geom_path(aes(col = delta))+
  scale_color_gradient(high = '#3182bd', low = '#c6dbef', name = expression(delta))+
  scale_x_continuous(breaks = c(0, 4000, 8000))+
  labs(x = 'Time', y = 'Overlap')+
  facet_grid(p+cost_c~d, labeller = labeller(cost_c = cost_labs, p = disturb_labs))+
  theme_bw()+
  theme(strip.text = element_text(size = 9), axis.text = element_text(size = 7),
        axis.title = element_text(size = 10))


#co-occurrence
co_gg <- landscape_ts%>%
  mutate(d = plyr::mapvalues(disp_rate, from = seq(from = 0, to = 0.4, length.out = 5),
                             to = c('d = 0','d = 0.1','d = 0.2','d = 0.3','d = 0.4')))%>%
  ggplot(aes(x = time, y = cooccurrence, group = ID))+
  geom_path(aes(col = delta))+
  labs(x = 'Time', y = 'Co-occurrence (%)')+
  scale_color_gradient(high = '#3182bd', low = '#c6dbef', name = expression(delta))+
  scale_x_continuous(breaks = c(0, 4000, 8000))+
  facet_grid(p+cost_c~d, labeller = labeller(cost_c = cost_labs, p = disturb_labs))+
  theme_bw()+
  theme(strip.text = element_text(size = 9), axis.text = element_text(size = 7),
        axis.title = element_text(size = 10))


#Combine
figS3 <- nmi_gg + overlap_gg + co_gg + plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'A') &  theme(plot.tag = element_text(size = 20))


ggsave(figS3, file = 'figS3.png',
       path = here::here("figs"),
       width = 15, height = 6, units = 'in')







################################################################
# Stationary distributions -----------------------------------
################################################################
stationary_data <-  here::here('R/data/stationary-states-28FEB25.rds')%>%
  readRDS()

t_intervals <- c(seq(from = 50, to = 150, by = 50), 1000)

#Function to calculate the statistical mode
calc_mode <- function(x) {
  # Remove NA values if any
  x <- na.omit(x)
  # Compute the kernel density estimate
  dens <- density(x)
  # Return the x-value corresponding to the maximum density estimate
  dens$x[which.max(dens$y)]
}

#Function for facet labels
custom_labeller <- function(variable, value) {
  if (variable == "Metric") {
    # Parse the expression so that, e.g., "Niche~overlap~rho" becomes a mathematical label.
    return(parse(text = as.character(value)))
  } else if (variable == "time") {
    # Parse the time labels to show as "T = <value>"
    # Here we use '==' in the expression, which renders as an equals sign.
    return(parse(text = paste0("T==", value)))
  }
}

#Modify data for plotting
stationary_plot_data <- stationary_data%>%
  filter(time %in% t_intervals)%>%
  mutate(d = paste0('d = ', disp_rate))%>%
  mutate(t = paste0('T = ', time))%>%
  mutate(t = factor(t, levels = paste0('T = ', t_intervals)))

#Facet themes
facet_tm <- theme(strip.text = element_text(size = 12),
      strip.text.y = element_text(size = 8),
      axis.title = element_text(size = 10))


################################################################
#Figure S2: Statistical modes
#################################################################
figs2 <- stationary_data%>%
  filter(time %in% t_intervals)%>%
  tidyr::pivot_longer(cols = c('NMI','rho','cooccurrence'), names_to = 'Metric')%>%
  select(disp_rate, delta, time, Metric, value)%>%
  group_by(disp_rate, delta, time, Metric)%>%
  summarise(mode = calc_mode(value))%>%
  mutate(Metric = case_when(
    Metric == 'NMI' ~ "NMI~'%'",
    Metric == 'rho' ~ "Niche~overlap~rho",
    TRUE ~ "Co-occurrence~'%'"
  )) %>%
  mutate(t = paste0('T = ', time))%>%
  mutate(t = factor(t, levels = paste0('T = ', t_intervals)))%>%
  ggplot(aes(x = disp_rate, y = mode, group = delta, col = delta))+
  geom_point()+
  geom_line()+
  scale_color_gradient(low = grad_pal[2], high = grad_pal[1], name = expression(delta))+
  facet_grid(Metric ~ t, scales = "free_y",
             labeller = labeller(Metric = label_parsed))+
  labs(x = "Dispersal rate (d)", y = "Statistical mode")+
  theme_bw()+
  theme(strip.text = element_text(size = 12),
        strip.text.y = element_text(size = 12), axis.title = element_text(size = 14))


ggsave(figs2, file = 'figS2.png',
       path = here::here("figs"),
       width = 9, height = 6, units = 'in')




################################################################
#Figure 6: Stationary distributions
#################################################################
#Overlap
overlap_gg <- stationary_plot_data%>%
  ggplot(aes(x = rho, group = delta))+
  geom_density(aes(fill = delta),col = 'black', alpha = 0.8)+
  scale_fill_gradient(high = '#3182bd', low = '#c6dbef', name = expression(delta))+
  scale_color_gradient(high = '#3182bd', low = '#c6dbef', name = expression(delta))+
  #labs(y = 'Count', x = "Mutual information (bits)")+
  #xlim(c(0,0.11))+
  labs(y = 'Count', x = expression("Niche overlap ("~rho~")"))+
  facet_grid(d~t)+
  theme_bw()+
  facet_tm


#Scale of local exclusion (normalized mutual information)
nmi_gg <- stationary_plot_data%>%
  ggplot(aes(x = NMI, group = delta))+
  geom_density(aes(fill = delta),col = 'black', alpha = 0.8)+
  scale_fill_gradient(high = '#3182bd', low = '#c6dbef', name = expression(delta))+
  scale_color_gradient(high = '#3182bd', low = '#c6dbef', name = expression(delta))+
  labs(y = 'Count', x = "Normalized mutual information (%)")+
  facet_grid(d~t)+
  theme_bw()+
  facet_tm

#Co-occurence
co_gg <- stationary_plot_data%>%
  ggplot(aes(x = cooccurrence, group = delta))+
  geom_density(aes(fill = delta),col = 'black', alpha = 0.8)+
  scale_fill_gradient(high = '#3182bd', low = '#c6dbef', name = expression(delta))+
  scale_color_gradient(high = '#3182bd', low = '#c6dbef', name = expression(delta))+
  labs(y = 'Count', x = "Co-occurrence (%)")+
  facet_grid(d~t)+
  theme_bw()+
  facet_tm


fig6 <- overlap_gg / nmi_gg / co_gg + plot_annotation(tag_levels = 'A') + plot_layout(guides = "collect") &  theme(plot.tag = element_text(size = 20))


ggsave(fig6, file = 'fig6.png',
       path = here::here("figs"),
       width = 7, height = 10, units = 'in')


################################################################
# Coexistence space  -----------------------------------
################################################################
#Load data
coexist_summ <- here::here('R/data/COEXsim_summary_28FEB25.rds')%>%
  readRDS()

# Disturbance
levels(coexist_summ$p) <- c(expression("Disturbance:~'0.1%'"),
                            expression("Disturbance:~'10%'"))
# Dispersal cost
levels(coexist_summ$cost_c) <- c("Cost: 1", "Cost: 3")

#Aij
levels(coexist_summ$aij) <- paste0('alpha[ij] ~ ":" ~ ', unique(coexist_summ$aij))

#Shuffle
levels(coexist_summ$shuffle_p) <- c(expression("Shuffle:~'0'"),
                                    expression("Shuffle:~'1'"))

############################################################
# Supplement S1: Full parameter space -----------------
############################################################
#Load data
coexist_summ <- here::here('R/data/COEXsim_summary_28FEB25.rds')%>%
  readRDS()%>%
  bind_rows(here::here('R/data/ASSsim_summary_28FEB25.rds')%>%
              readRDS())


#Combine and add labels to variables
summ_data <- coexist_summ%>%
  #bind_rows(ass_summ)%>%
  mutate(aij_lab = aij)%>%
  mutate(p_lab = p)%>%
  mutate(cost_lab = cost_c)%>%
  mutate(shuffle_lab = shuffle_p) #create a separate aij variable for plot labels


# Disturbance
levels(summ_data$p_lab) <- c(expression("Disturbance:~'0.1%'"),
                            expression("Disturbance:~'10%'"))
# Dispersal cost
levels(summ_data$cost_lab) <- c("Cost: 0", "Cost: 3")

#aij
levels(summ_data$aij_lab) <- paste0('alpha[ij] ~ ":" ~ ', unique(summ_data$aij_lab))

#Shuffle
levels(summ_data$shuffle_lab) <- c(expression("Shuffle:~'0'"),
                                    expression("Shuffle:~'1'"))



# Figure S1: Full parmeter space -------
figs1 <- summ_data%>%
  ggplot(aes(x = disp_rate, y = delta))+
  geom_tile(aes(fill = time), col = 'black', size = 0.5)+
  scale_fill_gradient(trans = "log",low = 'grey', high = 'black', breaks = c(500, 1000, 2000, 4000, 8000), name = "Coexistence time")+
  labs(x = "Dispersal rate (d)", y = expression("Dispersal assortativity ("~delta~")"))+
  facet_grid(p_lab+cost_lab ~ aij_lab+shuffle_lab,labeller = label_parsed)+
  theme_bw()


ggsave(figs1, file = 'figS1.png',
       path = here::here("figs"),
       width = 14, height = 7, units = 'in')

################################################################
# Alternative stable states  -----------------------------------
################################################################

#Generate reference line by fitting a spline ------------
# 1. subset the coexistence times between interval C1 and C2
# 2. Fit a gam to predict delta with disp_rate using the subsetted data
# NOTE: This prediction corresponds to the region of delta~d where coexistence is equal mean(C1, C2)

#Subset only ass summaries
ass_summ <- summ_data%>%
  mutate(aij = as.character(aij)%>%as.numeric)%>%
  filter(aij > 1)

#Define the range of coexistnce time to reference and subset data
time_range <- c(4000, 4500) # t = 4000  - 4500
#Subset data within that range
sub_time <- ass_summ%>%
  filter(shuffle_p == 0, cost_c == 0, p == 0.001, aij == 2.718)%>%
  ungroup%>%
  select(delta, disp_rate, time)%>%
  filter(time > min(time_range) & time < max(time_range))%>%
  filter(disp_rate > 0)

#fit GAM
cotime_gam <- mgcv::gam(delta ~ s(disp_rate, k = 1, fx = TRUE), data = sub_time)

#generate prediction
co4000_prediction <- data.frame(disp_rate = seq(from = 0, to = 0.4, by =0.001))%>%
  mutate(delta = predict(cotime_gam, newdata = data.frame(disp_rate = disp_rate)))%>%
  filter(delta > 0 & delta < 0.4)

################################################################
#Figure 1: Coexistnece heatmap ------------------
################################################################
main_gg <- ass_summ%>%
  filter(shuffle_p == 0, cost_c == 0, p == 0.001, aij == '2.718')%>%
  #filter(shuffle_p == 0, p == 0.001, aij == '2.718')%>%
  ggplot(aes(x = disp_rate, y = delta))+
  geom_tile(aes(fill = time), col = 'black', size = 0.5)+
  geom_line(data = co4000_prediction, aes(x = disp_rate, y = delta), size = 2, col = '#cb181d')+
  scale_fill_gradient(trans = "log",low = 'grey', high = 'black', breaks = c(500, 1000, 2000, 4000, 8000), name = "Time")+
  labs(x = "Dispersal rate (d)", y = expression("Dispersal assortativity ("~delta~")"))+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 15),
        axis.text.x = element_text(margin = margin(-15,0,0,0)),
        axis.text.y = element_text(margin = margin(0,-15,0,0)),
        axis.text = element_text(size = 15), axis.title = element_text(size = 17),
        legend.position= "left",
        legend.text.position = "left")

#Dispersal marginal plot
disp_gg <- ass_summ%>%
  filter(shuffle_p == 0, cost_c == 0, p == 0.001, aij == '2.718')%>%
  #filter(shuffle_p == 0, cost_c == 0, p == 0.001, aij == 2)%>%
  group_by(disp_rate)%>%
  summarise(time = mean(time))%>%
  ggplot(aes(x = disp_rate, y = time))+
  geom_point(size = 2)+
  geom_line(size = 1)+
  labs(y = "Time")+theme_classic()+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 13), axis.title.y = element_text(size = 15),
        axis.line =element_line(size = 1))

#Delta marginal plot
delta_gg <- ass_summ%>%
  filter(shuffle_p == 0, cost_c == 0, p == 0.001, aij == '2.718')%>%
  group_by(delta)%>%
  summarise(time = mean(time))%>%
  ggplot(aes(x = delta, y = time))+
  geom_point(size = 2)+
  geom_line(size = 1)+
  theme_bw()+
  labs(y = "Time")+
  coord_flip()+
  theme_classic()+
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12, angle = -70), axis.title.x = element_text(size = 15), axis.line = element_line(size = 1))


fig1 <- disp_gg + plot_spacer() + main_gg + delta_gg + plot_layout(ncol = 2, heights = c(0.8, 3), widths = c(3, 0.8))
ggsave(fig1, file = 'fig1.png',
       path = here::here("figs"),
       width = 9, height = 8, units = 'in')


################################################################
#Figure 2: Marginal effects of competition, disturbance, and dispersal costs ------------
#################################################################
aij_gg <- ass_summ%>%
  mutate(Shuffle = ifelse(shuffle_p == 1, 'yes', 'no'))%>%
  group_by(aij, Shuffle)%>%
  summarise(time = mean(time))%>%
  ggplot(aes(x = aij, y = time, group = Shuffle))+
  geom_path(position = position_dodge(0.3))+
  geom_point(size = 8, aes(fill = Shuffle), shape = 21, position = position_dodge(0.3))+
  geom_text(aes(y = time + 250, label = round(time)), position = position_dodge(0.3), col = 'red', fontface = 'bold')+
  scale_fill_manual(values = c('white','black'))+
  labs(x = expression("Interspecific competition"~"("~alpha['ij']~")"), y = "Coexistence time")+
  ylim(c(0, 3000))+
  theme_bw()+
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 13))


dist_gg <- ass_summ%>%
  mutate(Shuffle = ifelse(shuffle_p == 1, 'yes', 'no'))%>%
  group_by(p, Shuffle)%>%
  summarise(time = mean(time))%>%
  ggplot(aes(x = p, y = time, group = Shuffle))+
  geom_path(position = position_dodge(0.3))+
  geom_point(size = 8, aes(fill = Shuffle), shape = 21, position = position_dodge(0.3))+
  geom_text(aes(y = time + 250, label = round(time)), position = position_dodge(0.3), col = 'red', fontface = 'bold')+
  scale_fill_manual(values = c('white','black'))+
  #scale_shape_manual(values = c(16, 21))+
  labs(x = 'Disturbance (p)', y = "")+
  ylim(c(0, 3000))+
  #scale_y_log10()+
  theme_bw()+
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 13))


cost_gg <- ass_summ%>%
  mutate(Shuffle = ifelse(shuffle_p == 1, 'yes', 'no'))%>%
  group_by(cost_c, Shuffle)%>%
  summarise(time = mean(time))%>%
  ggplot(aes(x = cost_c, y = time, group = Shuffle))+
  geom_path(position = position_dodge(0.3))+
  geom_point(size = 8, aes(fill = Shuffle), shape = 21, position = position_dodge(0.3))+
  geom_text(aes(y = time + 250, label = round(time)), position = position_dodge(0.3), col = 'red', fontface = 'bold')+
  scale_fill_manual(values = c('white','black'))+
  ylim(c(0, 3000))+
  #scale_y_log10()+
  labs(x = 'Cost of Assortative dispersal (c)', y = "")+
  theme_bw()+
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 13))


fig2 <- aij_gg + dist_gg + cost_gg + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 23))
ggsave(fig2, file = 'fig2.png',
       path = here::here("figs"),
       width = 16, height = 5, units = 'in')


################################################################
#Figure 3: Effects of spatial structure ------------
################################################################

fig3 <-ass_summ%>%
  filter(aij > 1)%>%
  #ungroup()%>%
  #select(disp_rate, delta, aij, p, cost_c, shuffle_p, time)%>%
  #tidyr::pivot_wider(names_from = shuffle_p, values_from = time)%>%
  select(disp_rate, shuffle_lab, aij_lab, p_lab, cost_lab, time)%>%
  #filter(p != 0.01, aij != 1.5, cost_c != 1)%>%
  mutate(shuffle_lab = ifelse(shuffle_lab == "Shuffle:~'0'", 'no_shuffle', 'shuffle')%>%
           as.factor)%>%
  tidyr::pivot_wider(names_from = shuffle_lab, values_from = time)%>%
  #mutate(shuffle_effect = no_shuffle - shuffle)%>%
  mutate(shuffle_log_ratio = log(no_shuffle/shuffle))%>%
  ggplot(aes(x = disp_rate, y = delta))+
  geom_tile(aes(fill = shuffle_log_ratio), col = 'black')+
  scale_fill_gradient2(low = '#f03b20', high = '#1c9099', mid = 0, name = expression(ln(frac(T[c], T[s]))))+
  facet_grid(aij_lab ~ p_lab + cost_lab ,labeller = label_parsed)+
  labs(x = "Dispersal rate (d)", y = expression("Dispersal assortativity ("~delta~")"))+
  theme_bw()+
  theme(panel.background = element_blank(),
        #panel.background = element_rect(fill = "black",colour = "black", linetype = "solid"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        axis.text = element_text(size = 10), axis.title = element_text(size = 13),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))


ggsave(fig3, file = 'fig3.png',
       path = here::here("figs"),
       width = 10, height = 6, units = 'in')

################################################################
#Figure 4: auxillary variable contrasts --------------------
################################################################
#Remove shuffle from averages
ass_raw <- here::here('R/data/ASSsim_28FEB25.rds')%>%
    readRDS()


ass_aux <- ass_raw%>%
  filter(shuffle_p == 0)%>%
  #filter(aij %in% c(1, exp(1)))%>%
  filter(aij >= 1)%>%
  mutate(aij = round(aij, 3)%>%as.factor)%>%
  group_by(disp_rate, delta, cost_c, aij, p)%>%
  summarise(time = mean(time))%>%
  mutate(p = factor(p, levels = c(0.001, 0.1),
               labels = c(expression("Disturbance:~'0.1%'"),
                          expression("Disturbance:~'10%'"))),
         cost_c = factor(cost_c, levels = c(0, 3),
                    labels = c("Cost: 0", "Cost: 3")) )%>%
  mutate(aij = factor(aij, levels = c(1, 1.649, 2.718),
                      labels = paste0('alpha[ij] ~ ":" ~ ', c(1, 1.649, 2.718))))

#1. Disturbance
dist_p <- ass_aux%>%
  filter(aij %in% c('alpha[ij] ~ ":" ~ 1.649', 'alpha[ij] ~ ":" ~ 2.718'))%>%
  mutate(p = ifelse(p == "Disturbance:~'0.1%'", 'p0.1', 'p10'))%>%
  tidyr::pivot_wider(names_from = p, values_from = time)%>%
  mutate(p_effect = log(p10/p0.1))%>%
  ggplot(aes(x = disp_rate, y = delta))+
  geom_tile(aes(fill = p_effect), col = 'black')+
  scale_fill_gradient2(low = '#f03b20', high = '#1c9099', mid = 0, name = expression(ln(frac(T['p = 0.1'], T['p = 0.001']))))+
  facet_grid(aij ~ cost_c, labeller = label_parsed)+
  labs(x = "Dispersal rate (d)", y = expression("Dispersal assortativity ("~delta~")"))+
  theme_bw()+
  theme(panel.background = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10), axis.title = element_text(size = 12),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10))

#2. Costs
cost_p <- ass_aux%>%
  filter(aij %in% c('alpha[ij] ~ ":" ~ 1.649', 'alpha[ij] ~ ":" ~ 2.718'))%>%
  mutate(cost_c = ifelse(cost_c == "Cost: 0", 'cost_0', 'cost_3'))%>%
  tidyr::pivot_wider(names_from = cost_c, values_from = time)%>%
  mutate(cost3_effect = log(cost_3/cost_0))%>%
  ggplot(aes(x = disp_rate, y = delta))+
  geom_tile(aes(fill = cost3_effect), col = 'black')+
  scale_fill_gradient2(low = '#f03b20', high = '#1c9099', mid = 0, name = expression(ln(frac(T['c = 3'], T['c = 0']))))+
  facet_grid(aij ~ p,labeller = label_parsed)+
  labs(x = "Dispersal rate (d)", y = expression("Dispersal assortativity ("~delta~")"))+
  theme_bw()+
  theme(panel.background = element_blank(),
        #panel.background = element_rect(fill = "black",colour = "black", linetype = "solid"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10), axis.title = element_text(size = 12),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10))

#3. Interspecific interaction
int_p <- ass_aux%>%
  #filter(aij > 1)%>%
  #mutate(aij = ifelse(aij == 'alpha[ij] ~ ":" ~ 1', 'aij_1', 'aij_2.7'))%>%
  mutate(aij = case_when(aij == 'alpha[ij] ~ ":" ~ 1.649' ~ "aij_1.6",
                         aij == 'alpha[ij] ~ ":" ~ 1' ~ "aij_1",
                         TRUE ~ "aij_2.7"))%>%
  tidyr::pivot_wider(names_from = aij, values_from = time)%>%
  #mutate(aij_effect = log(aij_2.7/aij_1))%>%
  mutate(aij_effect = log(aij_2.7/aij_1.6))%>%
  ggplot(aes(x = disp_rate, y = delta))+
  geom_tile(aes(fill = aij_effect), col = 'black')+
  scale_fill_gradient2(low = '#f03b20', high = '#1c9099', mid = 0,
                       name = expression(ln(frac(T[alpha[ij] == 2.71], T[alpha[ij] == 1.65]))))+
  facet_grid(cost_c ~ p, labeller = label_parsed)+
  labs(x = "Dispersal rate (d)", y = expression("Dispersal assortativity ("~delta~")"))+
  theme_bw()+
  theme(panel.background = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10), axis.title = element_text(size = 12),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10))


fig4 <- int_p / dist_p / cost_p + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 20))


ggsave(fig4, file = 'fig4.png',
       path = here::here("figs"),
       width = 7, height = 12, units = 'in')








