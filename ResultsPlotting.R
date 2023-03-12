# Loading tidyverse
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
# loading gridExtra
if (!require('gridExtra')) install.packages('gridExtra'); library('gridExtra')

### Loading the results of the analysis:
# Real networks
real_net_autoreg <- read_csv("real_net_autoreg.csv")
real_net_FFL <- read_csv("real_net_FFL.csv")
real_net_path <- scan("real_net_path.csv", sep = ",") # csv of only one line
names(real_net_path) <- c("all", "wsc", "sc")

# Simulated networks
# all evidences (null, Weak, Strong, Confirmed):
all_autoreg <- read_csv("all_random_autoreg.csv")
all_FFL <- read_csv("all_random_FFL.csv")
all_path <- scan("all_random_path.csv") # only one value

# wsc evidences (Weak, Strong, Confirmed):
wsc_autoreg <- read_csv("wsc_random_autoreg.csv")
wsc_FFL <- read_csv("wsc_random_FFL.csv")
wsc_path <- scan("wsc_random_path.csv", sep = ",") # csv of only one line

# sc evidences (Strong, Confirmed):
sc_autoreg <- read_csv("sc_random_autoreg.csv")
sc_FFL <- read_csv("sc_random_FFL.csv")
sc_path <- scan("sc_random_path.csv", sep = ",") # csv of only one line

##### Plots:

#### FFL

# All:

real_all_FFL_stacked <- real_net_FFL[1, ] %>% 
                            stack() %>% 
                                cbind(network = rep("real", 12))

rnd_all_FFL_stacked <- all_FFL %>% 
                            stack() %>% 
                                 cbind(network = rep("random", 12))

png("all_FFL_real_vs_rnd.png", height = 360, width = 720)
rbind(real_all_FFL_stacked, rnd_all_FFL_stacked) %>%
  ggplot(aes(fill = network, y=values, x=ind)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) +
  labs(x = "Type of FFL",
       y = "Amount",
       title = "FFL comparison - real vs random networks (all evidence)") +
  scale_y_log10()
dev.off()


# Wsc and Sc:

# only the random networks, bubble plots to represent discrete data
wsc_FFL_plot <- stack(wsc_FFL) %>%
  ggplot(aes(x=ind, y=values)) + 
  geom_count(aes(size = ..n..)) +
  scale_size(range=c(5, 15)) +
  labs(x = "Type of FFL",
       y = "Amount of FFLs",
       title = "Feed Forward Loops - Evidences Weak, Strong and Confirmed") +
  ylim(0, 5) # so it is comparable with the sc plot

# adding the number of simulations (n) to each circle
full_wsc_FFL_plot <- wsc_FFL_plot + 
  geom_text(data = ggplot_build(wsc_FFL_plot)$data[[1]], 
            aes(x, y, label = n), color = "#ffffff")

png("wsc_rnd_FFL.png", height = 360, width = 600)
full_wsc_FFL_plot
dev.off()

sc_FFL_plot <- stack(sc_FFL) %>%
  ggplot(aes(x=ind, y=values)) + 
  geom_count(aes(size = ..n..)) +
  scale_size(range=c(5, 15)) +
  labs(x = "Type of FFL",
       y = "Amount of FFLs",
       title = "Feed Forward Loops - Evidences Strong and Confirmed")

# adding the number of simulations (n) to each circle
full_sc_FFL_plot <- sc_FFL_plot + 
  geom_text(data = ggplot_build(sc_FFL_plot)$data[[1]], 
            aes(x, y, label = n), color = "#ffffff")

png("sc_rnd_FFL.png", height = 360, width = 600)
full_sc_FFL_plot
dev.off()

# comparison between real and random networks
# getting the real df in an appropriate format and adding the network = "real" column
real_wsc_FFL_stacked <- real_net_FFL[2, ] %>% 
                            stack() %>% 
                                cbind(network = rep("real", 12))
# computing the mean of the results of the random network according to the type,
# and getting it in an appropriate format and adding the network = "random" column
rnd_wsc_FFL_mean_stacked <- wsc_FFL %>%
                              stack() %>%
                                group_by(ind) %>%
                                  summarise(across(where(is.numeric), mean)) %>%
                                    cbind(network = rep("random", 12))

# grouped bar plots
png("wsc_FFL_real_vs_rnd.png", height = 360, width = 720)
rbind(real_wsc_FFL_stacked, rnd_wsc_FFL_mean_stacked) %>%
  ggplot(aes(fill = network, y=values, x=ind)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) +
  labs(x = "Type of FFL",
       y = "Amount",
       title = "FFL comparison - real vs random networks (Weak and strong evidence)") +
  ylim(0, 4)
dev.off()



# getting the real df in an appropriate format and adding the network = "real" column
real_sc_FFL_stacked <- real_net_FFL[3, ] %>% 
                          stack() %>% 
                              cbind(network = rep("real", 12))
# computing the mean of the results of the random network according to the type,
# and getting it in an appropriate format and adding the network = "random" column
rnd_sc_FFL_mean_stacked <- sc_FFL %>%
                             stack() %>%
                               group_by(ind) %>%
                                 summarise(across(where(is.numeric), mean)) %>%
                                   cbind(network = rep("random", 12))

# grouped bar plots
png("sc_FFL_real_vs_rnd.png", height = 360, width = 720)
rbind(real_sc_FFL_stacked, rnd_sc_FFL_mean_stacked) %>%
  ggplot(aes(fill = network, y=values, x=ind)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) +
  labs(x = "Type of FFL",
       y = "Amount",
       title = "FFL comparison - real vs random networks (Strong evidence)") +
  ylim(0, 4)
dev.off()



#### Autoregulations

# All:
real_all_auto_stacked <- real_net_autoreg[1, ] %>% 
  stack() %>% 
  cbind(network = rep("real", 12))

rnd_all_auto_stacked <- all_autoreg %>% 
  stack() %>% 
  cbind(network = rep("random", 12))

png("all_autoreg_real_vs_rnd.png", height = 360, width = 480)
rbind(real_all_auto_stacked, rnd_all_auto_stacked) %>%
  ggplot(aes(fill = network, y=values, x=ind)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) +
  labs(x = "Type of autoregulation",
       y = "Amount",
       title = "Autoregulation comparison - real vs random networks (all evidence)")
dev.off()

# Wsc and Sc:

# only the random networks, bubble plots to represent discrete data
# wsc
wsc_AR_plot <- stack(wsc_autoreg) %>%
  ggplot(aes(x=ind, y=values)) + 
  geom_count(aes(size = ..n..)) +
  scale_size(range=c(5, 15)) +
  labs(x = "Type of autoregulation",
       y = "Amount of motifs",
       title = "Autoregulation - Weak and strong evidence")

# adding the number of simulations (n) to each circle
full_wsc_AR_plot <- wsc_AR_plot + 
  geom_text(data = ggplot_build(wsc_AR_plot)$data[[1]], 
            aes(x, y, label = n), color = "#ffffff")

png("wsc_rnd_autoreg.png", height = 360, width = 480)
full_wsc_AR_plot
dev.off()

# sc
sc_AR_plot <- stack(sc_autoreg) %>%
  ggplot(aes(x=ind, y=values)) + 
  geom_count(aes(size = ..n..)) +
  scale_size(range=c(5, 15)) +
  labs(x = "Type of autoregulation",
       y = "Amount of motifs",
       title = "Autoregulation - Strong evidence")

# adding the number of simulations (n) to each circle
full_sc_AR_plot <- sc_AR_plot + 
  geom_text(data = ggplot_build(sc_AR_plot)$data[[1]], 
            aes(x, y, label = n), color = "#ffffff")

png("sc_rnd_autoreg.png", height = 360, width = 480)
full_sc_AR_plot
dev.off()

# comparison between real and random networks
# getting the real df in an appropriate format and adding the network = "real" column
real_wsc_AR_stacked <- real_net_autoreg[2, ] %>% 
  stack() %>% 
  cbind(network = rep("real", 12))
# computing the mean of the results of the random network according to the type,
# and getting it in an appropriate format and adding the network = "random" column
rnd_wsc_AR_mean_stacked <- wsc_autoreg %>%
  stack() %>%
  group_by(ind) %>%
  summarise(across(where(is.numeric), mean)) %>%
  cbind(network = rep("random", 12))

# grouped bar plots
png("wsc_autoreg_real_vs_rnd.png", height = 360, width = 480)
rbind(real_wsc_AR_stacked, rnd_wsc_AR_mean_stacked) %>%
  ggplot(aes(fill = network, y=values, x=ind)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) +
  labs(x = "Type of autoregulation",
       y = "Amount",
       title = "Autoregulation comparison (Weak and strong evidence)") +
  ylim(0, 10)
dev.off()



# getting the real df in an appropriate format and adding the network = "real" column
real_sc_AR_stacked <- real_net_autoreg[3, ] %>% 
  stack() %>% 
  cbind(network = rep("real", 12))
# computing the mean of the results of the random network according to the type,
# and getting it in an appropriate format and adding the network = "random" column
rnd_sc_AR_mean_stacked <- sc_autoreg %>%
  stack() %>%
  group_by(ind) %>%
  summarise(across(where(is.numeric), mean)) %>%
  cbind(network = rep("random", 12))

# grouped bar plots
png("sc_autoreg_real_vs_rnd.png", height = 360, width = 480)
rbind(real_sc_AR_stacked, rnd_sc_AR_mean_stacked) %>%
  ggplot(aes(fill = network, y=values, x=ind)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) +
  labs(x = "Type of autoregulation",
       y = "Amount",
       title = "Autoregulation comparison (Strong evidence)") +
  ylim(0, 10)
dev.off()



#### Longest cascade

# All:
all_path_table = as.data.frame(list(Real = real_net_path[1], Random = all_path))
png("all_path_table.png", height = 50*nrow(all_path_table), width = 200*ncol(all_path_table))
grid.table(all_path_table)
dev.off()

# Wsc and Sc:
wsc_path_plot <- as.data.frame(wsc_path) %>%
  ggplot(aes(wsc_path)) + 
  geom_bar() +
  labs(x = "Maximum path length",
       y = "Number of simulations",
       title = "Longest cascade - Weak and Strong evidence") +
  xlim(0, 40) +
  ylim(0, 18) +
  geom_segment(aes(x = real_net_path[2], y = 0, 
                   xend = real_net_path[2], yend = 18), 
               lwd = 1.5, lty = "dashed", color = "darkblue")

sc_path_plot <- as.data.frame(sc_path) %>%
  ggplot(aes(sc_path)) + 
  geom_bar() +
  labs(x = "Maximum path length",
       y = "Number of simulations",
       title = "Longest cascade - Strong evidence") +
  xlim(0, 40) +
  ylim(0, 18) +
  geom_segment(aes(x = real_net_path[3], y = 0, 
                   xend = real_net_path[3], yend = 18), 
               lwd = 1.5, lty = "dashed", color = "darkblue")

png("wsc_sc_path_barplot.png", height = 854, width = 640)
grid.arrange(wsc_path_plot, sc_path_plot, nrow = 2)
dev.off()
