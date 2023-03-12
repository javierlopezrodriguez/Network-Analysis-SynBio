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


# Z-score = (frequency_real - mean_frequency_random) / sqrt(variance_random)
# Statistical significance: Z > 2.0 (equivalent to p < 0.01)
# Abundance = (frequency_real - mean_frequency_random) / (frequency_real + mean_frequency_random + small_positive_value)
# Abundance = -1: underrepresented
# Abundance = 1: overrepresented
# doi: 10.1093/bib/bbr033

# Autoregulation, WSC:

real_freqs <- real_net_autoreg[2, ] %>% 
  stack() %>% 
  rename(type = ind, f_real = values) %>% 
  as_tibble()

mean_rnd_freqs <- wsc_autoreg %>% 
  stack() %>% 
  group_by(ind) %>%
  summarise(across(where(is.numeric), mean)) %>% 
  rename(type = ind, f_mean_rnd = values)

variance_rnd <- wsc_autoreg %>% 
  stack() %>% 
  group_by(ind) %>%
  summarise(across(where(is.numeric), var)) %>% 
  rename(type = ind, var_rnd = values)

all_data <- mean_rnd_freqs %>% 
  inner_join(real_freqs) %>% 
  inner_join(variance_rnd) %>%
  relocate(f_real, .after = type)

stats_wsc_autoreg <- all_data %>% 
  mutate(Z_score = (f_real - f_mean_rnd) / sqrt(var_rnd)) %>%
  mutate(Abundance = (f_real - f_mean_rnd) / (f_real + f_mean_rnd + 0.001)) %>%
  rename(mean.freq.random = f_mean_rnd,
         real.freq = f_real,
         var.random = var_rnd,
         Autoregulation.wsc = type) %>%
  mutate(Significance = ifelse(Z_score > 2.0, 'Yes', 'No'))

# Autoregulation, SC:

real_freqs <- real_net_autoreg[3, ] %>% 
  stack() %>% 
  rename(type = ind, f_real = values) %>% 
  as_tibble()

mean_rnd_freqs <- sc_autoreg %>% 
  stack() %>% 
  group_by(ind) %>%
  summarise(across(where(is.numeric), mean)) %>% 
  rename(type = ind, f_mean_rnd = values)

variance_rnd <- sc_autoreg %>% 
  stack() %>% 
  group_by(ind) %>%
  summarise(across(where(is.numeric), var)) %>% 
  rename(type = ind, var_rnd = values)

all_data <- mean_rnd_freqs %>% 
  inner_join(real_freqs) %>% 
  inner_join(variance_rnd) %>%
  relocate(f_real, .after = type)

stats_sc_autoreg <- all_data %>% 
  mutate(Z_score = (f_real - f_mean_rnd) / sqrt(var_rnd)) %>%
  mutate(Abundance = (f_real - f_mean_rnd) / (f_real + f_mean_rnd + 0.001)) %>%
  rename(mean.freq.random = f_mean_rnd,
         real.freq = f_real,
         var.random = var_rnd,
         Autoregulation.sc = type) %>%
  mutate(Significance = ifelse(Z_score > 2.0, 'Yes', 'No'))

# FFL, WSC:

real_freqs <- real_net_FFL[2, ] %>% 
  stack() %>% 
  rename(type = ind, f_real = values) %>% 
  as_tibble()

mean_rnd_freqs <- wsc_FFL %>% 
  stack() %>% 
  group_by(ind) %>%
  summarise(across(where(is.numeric), mean)) %>% 
  rename(type = ind, f_mean_rnd = values)

variance_rnd <- wsc_FFL %>% 
  stack() %>% 
  group_by(ind) %>%
  summarise(across(where(is.numeric), var)) %>% 
  rename(type = ind, var_rnd = values)

all_data <- mean_rnd_freqs %>% 
  inner_join(real_freqs) %>% 
  inner_join(variance_rnd) %>%
  relocate(f_real, .after = type)

stats_wsc_FFL <- all_data %>% 
  mutate(Z_score = (f_real - f_mean_rnd) / sqrt(var_rnd)) %>%
  mutate(Abundance = (f_real - f_mean_rnd) / (f_real + f_mean_rnd + 0.001)) %>%
  rename(mean.freq.random = f_mean_rnd,
         real.freq = f_real,
         var.random = var_rnd,
         FFL.wsc = type) %>%
  mutate(Significance = ifelse(Z_score > 2.0, 'Yes', 'No'))

# FFL, SC:

real_freqs <- real_net_FFL[3, ] %>% 
  stack() %>% 
  rename(type = ind, f_real = values) %>% 
  as_tibble()

mean_rnd_freqs <- sc_FFL %>% 
  stack() %>% 
  group_by(ind) %>%
  summarise(across(where(is.numeric), mean)) %>% 
  rename(type = ind, f_mean_rnd = values)

variance_rnd <- sc_FFL %>% 
  stack() %>% 
  group_by(ind) %>%
  summarise(across(where(is.numeric), var)) %>% 
  rename(type = ind, var_rnd = values)

all_data <- mean_rnd_freqs %>% 
  inner_join(real_freqs) %>% 
  inner_join(variance_rnd) %>%
  relocate(f_real, .after = type)

stats_sc_FFL <- all_data %>% 
  mutate(Z_score = (f_real - f_mean_rnd) / sqrt(var_rnd)) %>%
  mutate(Abundance = (f_real - f_mean_rnd) / (f_real + f_mean_rnd + 0.001)) %>%
  rename(mean.freq.random = f_mean_rnd,
         real.freq = f_real,
         var.random = var_rnd,
         FFL.sc = type) %>%
  mutate(Significance = ifelse(Z_score > 2.0, 'Yes', 'No'))

### Tables to image:
stats_table <- stats_wsc_autoreg
png("stats_wsc_AR.png", height = 50*nrow(stats_table), width = 200*ncol(stats_table))
grid.table(stats_table)
dev.off()

stats_table <- stats_sc_autoreg
png("stats_sc_AR.png", height = 50*nrow(stats_table), width = 200*ncol(stats_table))
grid.table(stats_table)
dev.off()

stats_table <- stats_wsc_FFL
png("stats_wsc_FFL.png", height = 50*nrow(stats_table), width = 200*ncol(stats_table))
grid.table(stats_table)
dev.off()

stats_table <- stats_sc_FFL
png("stats_sc_FFL.png", height = 50*nrow(stats_table), width = 200*ncol(stats_table))
grid.table(stats_table)
dev.off()
