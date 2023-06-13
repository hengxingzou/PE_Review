# Simulating Priority Effects in a Consumer-Resource Model
# Hengxing Zou

# packages
library(tidyverse)
library(magrittr)
library(patchwork)
library(doParallel)
registerDoParallel(cores = 4)

# color scheme from https://jfly.uni-koeln.de/color/
color_scheme = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##### Functions #####

# The Model

model_kl = function(time, state, pars, deltap, type) {
  
  all_pop = tibble()
  
  # type = 1: trait-mediated priority effects (k_i changes with \Delta p)
  # type = 0: numeric priority effects (k_i constant)
  
  k1 = type*get_k(2*k1, scal, xmid, deltap) + (1-type)*k1
  k2 = type*get_k(2*k2, -scal, xmid, deltap) + (1-type)*k2
  
  for (t in 1:time) {
    
    new_state = with(as.list(c(state, pars)), {
      
      new_R = r*R*(1-R/K)+R-Q1*lambda1*N1*R^2/(k1+R^2)-Q2*lambda2*N2*R/(k2+R+R^2/ki2)
      new_N1 = lambda1*N1*R^2/(k1+R^2)+(1-mu1)*N1
      new_N2 = lambda2*N2*R/(k2+R+R^2/ki2)+(1-mu2)*N2
      
      c(N1 = new_N1, N2 = new_N2, R = new_R)
    })
    
    all_pop = rbind(all_pop, tibble(Time = t, 
                                    N1 = new_state["N1"], N2 = new_state["N2"], R = new_state["R"]))
    state = new_state
    
  }
  
  return(all_pop)
}

# Trait-mediated Priority Effects

get_k = function(base_k, scal, xmid, deltap) {
  return(base_k/(1+exp((xmid-deltap)/scal)))
}

# Simulating Early/Late Arrival over Generations

arrival_time = function(time_total, time_early, init_pop, deltap, type) {
  
  time_both = time_total-abs(time_early)
  
  if (time_early == 0) {
    out = model_kl(time_both, init_pop, pars = NULL, deltap, type)
    return(out)
  } 
  
  else {
    
    late_spp = if_else(time_early < 0, 2, 1)
    late_spp_init = init_pop[late_spp]
    init_pop[late_spp] = 0
    
    out1 = model_kl(abs(time_early), init_pop, pars = NULL, deltap, type)
    
    init_pop["N1"] = as.numeric(tail(out1, 1)[, 2])
    init_pop["N2"] = as.numeric(tail(out1, 1)[, 3])
    init_pop["R"] = as.numeric(tail(out1, 1)[, 4])
    init_pop[late_spp] = late_spp_init
    
    out2 = model_kl(time_both, init_pop, pars = NULL, deltap, type)
    
    return(rbind(out1, mutate(out2, Time = Time+abs(time_early))))
    
  }
}

##### Parameters #####

lambda1 = 0.029
lambda2 = 0.2
k1 = 0.02
k2 = 3
Q1 = 0.01
Q2 = 0.01
mu1 = mu2 = 0.01
ki2 = 1

r = 0.5
K = 3

scal = 7.5
xmid = 0

init_pop = c(N1 = 5, N2 = 5, R = K)

##### Simulation I. Numeric and Trait-mediated Priority Effects #####

# Species 2 arrives late by 100 generations (Species 1 wins)
spp1wins = arrival_time(3000, -100, init_pop, 0, type = 0)

# Species 2 arrives late by 100 generations but early by 4 time steps within a generation (\Delta p=4; Species 2 wins)
spp2wins = arrival_time(3000, -100, init_pop, 4, type = 1)

##### Simulation II. Priority Effects Both within and across Seasons #####

# List of arrival times in generations

early_gens = seq(-100, 20, 5)

# List of relative arrival times (\Delta p) within a generation

deltap_lst = seq(-4, 2, 1)

## Simulation II.i. Numeric Priority Effects ##

# Increase run time to 12000 to get rid of all transients; this takes a while, even with parallelization
# Run through both lists, collect final populations

num_results = data.frame()
for (dp in deltap_lst) {
  endpop = foreach(t = early_gens, .combine = rbind) %dopar% {
   arrival_time(12000, t, init_pop, dp, type = 0) %>% 
      tail(1) %>% 
      mutate(early = t, deltap = dp)
  }
  num_results = rbind(endpop, num_results)
}

# Categorize by outcomes of competition; 
# consumer is excluded if population of consumers < 1e-2
# 1, 2, and 3 codes for species 1 wins, 2 wins, and coexistence

num_results %<>% mutate(outcome = if_else(N1 > 1e-2 & N2 > 1e-2, 3, 
                                          if_else(N1 > 1e-2 & N2 < 1e-2, 1, 2)))

## Simulation II.ii. Trait-mediated Priority Effects ##

tra_results = data.frame()
for (dp in deltap_lst) {
  endpop = foreach(t = early_gens, .combine = rbind) %dopar% {
    arrival_time(12000, t, init_pop, dp, type = 1) %>% 
      tail(1) %>% 
      mutate(early = t, deltap = dp)
  }
  tra_results = rbind(endpop, tra_results)
}

# Categorize by outcomes of competition; 
# consumer is excluded if population of consumers < 1e-2
# 1, 2, and 3 codes for species 1 wins, 2 wins, and coexistence

tra_results %<>% mutate(outcome = if_else(N1 > 1e-2 & N2 > 1e-2, 3, 
                                          if_else(N1 > 1e-2 & N2 < 1e-2, 1, 2)))

##### Write Data #####

dta_dir = "Data"
if (!dir.exists(dta_dir)) {
  dir.create(dta_dir)
}

write_csv(spp1wins, "Data/1wins.csv")
write_csv(spp2wins, "Data/2wins.csv")

write_csv(num_results, "Data/num_results.csv")
write_csv(tra_results, "Data/tra_results.csv")

##### Figures #####

# Read Data

spp1wins = read_csv("Data/1wins.csv")
spp2wins = read_csv("Data/2wins.csv")

num_results = read_csv("Data/num_results.csv")
tra_results = read_csv("Data/tra_results.csv")

## Simulation I. Numeric and Trait-mediated Priority Effects

p_1wins =
  spp1wins %>% 
  pivot_longer(cols = 2:3, names_to = "Species", values_to = "Population") %>% 
  ggplot(aes(x = Time, y = Population, color = Species)) + 
  geom_line(size = 1.5) + 
  xlab("Years") +
  scale_y_continuous(name = "Density", limits = c(0, 1300), breaks = c(500, 1000)) + 
  # scale_y_log10() + 
  scale_color_manual(values = color_scheme, labels = c("1", "2")) + 
  ggtitle("A") + 
  theme_bw() + 
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 17), legend.title = element_text(size = 20), 
        plot.title = element_text(size = 20), panel.grid.minor = element_blank())

p_2wins =
  spp2wins %>% 
  pivot_longer(cols = 2:3, names_to = "Species", values_to = "Population") %>% 
  ggplot(aes(x = Time, y = Population, color = Species)) + 
  geom_line(size = 1.5) + 
  xlab("Years") +
  scale_y_continuous(name = "Density", limits = c(0, 1300), breaks = c(500, 1000)) + 
  # scale_y_log10() +
  scale_color_manual(values = color_scheme, labels = c("1", "2")) + 
  ggtitle("B") + 
  theme_bw() + 
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20), axis.title.y = element_blank(),
        legend.text = element_text(size = 17), legend.title = element_text(size = 20), 
        plot.title = element_text(size = 20), panel.grid.minor = element_blank())

# Simulation II: Priority Effects Both within and across Seasons

num_phase = 
  num_results %>% 
  ggplot(aes(x = deltap, y = early, fill = as.factor(outcome))) + 
  geom_tile() + 
  # geom_rect(color = "black", alpha = 0, 
  #           aes(xmin = -0.5, xmax = 0.5, ymin = min(early_gens)-2.5, ymax = max(early_gens)+2.5)) + 
  # geom_rect(color = "black", alpha = 0, linetype = 2, 
  #           aes(xmin = min(deltap_lst)-0.5, xmax = max(deltap_lst)+0.5, ymin = -2.5, ymax = 2.5)) +
  scale_fill_manual(values = color_scheme, labels = c("1 wins", "2 wins"), name = "Outcome") + 
  scale_x_continuous(name = expression(paste("Relative arrival times (", Delta, p, ")")), 
                     breaks = deltap_lst, label = deltap_lst) +
  scale_y_continuous(name = "Years early") +
  ggtitle("C") + 
  theme_bw() + 
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 17), legend.title = element_text(size = 20), legend.position = "left",
        plot.title = element_text(size = 20), panel.grid = element_blank())

tra_phase = 
  tra_results %>% 
  ggplot(aes(x = deltap, y = early, fill = as.factor(outcome))) + 
  geom_tile() + 
  # geom_rect(color = "black", alpha = 0, 
  #           aes(xmin = -0.5, xmax = 0.5, ymin = min(early_gens)-2.5, ymax = max(early_gens)+2.5)) + 
  # geom_rect(color = "black", alpha = 0, linetype = 2, 
  #           aes(xmin = min(deltap_lst)-0.5, xmax = max(deltap_lst)+0.5, ymin = -2.5, ymax = 2.5)) + 
  scale_fill_manual(values = color_scheme, labels = c("1 wins", "2 wins", "Coexistence"), name = "Outcome", guide = "none") + 
  scale_x_continuous(name = expression(paste("Relative arrival times (", Delta, p, ")")), 
                     breaks = deltap_lst, label = deltap_lst) +
  scale_y_continuous(name = "Years early") +
  theme_bw() + 
  ggtitle("D") +
  theme(axis.text = element_text(size = 15), axis.title.x = element_text(size = 20), axis.title.y = element_blank(),
        legend.text = element_text(size = 17), legend.title = element_text(size = 20), legend.position = "left",
        plot.title = element_text(size = 20), panel.grid = element_blank())

# Create directory

fig_dir = "Figures"
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir)
}

# Composite phase diagram

ggsave("Figures/Fig3I.pdf", 
       plot = (p_1wins | p_2wins) / (num_phase | tra_phase) + plot_layout(guides = "collect"), 
       width = 4000, height = 3000, units = "px")
  
