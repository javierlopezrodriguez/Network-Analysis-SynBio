# installing the packages if they aren't installed and loading them
if (!require('igraph')) install.packages('igraph'); library('igraph')
source('NetworkFunctions.R')

# reading the file 'network_tf_gene.txt'
tf_gene_net_df <- read.csv('network_tf_gene.txt', 
                  header = FALSE, 
                  sep = "\t",  
                  comment.char = "#") # removes comments starting with #

# removing the last column (it contains no information)
tf_gene_net_df <- tf_gene_net_df[ , -8]

# naming the columns
colnames(tf_gene_net_df) <- c("TF.id", "TF.name", "Gene.id", "Gene.name", 
                              "Effect", "Evidence", "Evidence.type")

# adding a field with the lowercase TF and Gene names
tf_gene_net_df$Source.name <- tolower(tf_gene_net_df$TF.name)
tf_gene_net_df$Target.name <- tolower(tf_gene_net_df$Gene.name)

### Creating an empty adjacency matrix:
# getting a list of the IDs, removing duplicates
source_names <- unique(tf_gene_net_df$Source.name)
target_names <- unique(tf_gene_net_df$Target.name)
all_names <- unique(c(source_names, target_names))
# the number of nodes (dimensions of the matrix) is the number of unique ids.
num_nodes <- length(all_names)
# creating a square matrix full of zeros with dimensions num_nodes x num_nodes
# with all_names as row and column names
tf_gene_adj_mat <- matrix(0, nrow = num_nodes, ncol = num_nodes,
                          dimnames = list(all_names, all_names))

### Populating the matrix:
# To populate the matrix, I'm going to have different cases:
# - The whole network, including interactions for which there is null evidence
all.hash <- list(null = 1, Weak = 1, Strong = 1, Confirmed = 1)
# - One in which I consider all of the evidence (weak, strong, confirmed)
w.s.c.hash <- list(null = 0, Weak = 1, Strong = 1, Confirmed = 1)
# - Another in which I only consider the strong and confirmed evidence
s.c.hash <- list(null = 0, Weak = 0, Strong = 1, Confirmed = 1) # Weak -> not interaction
# - The last one in which I only consider confirmed evidence
conf.hash <- list(null = 0, Weak = 0, Strong = 0, Confirmed = 1) # Weak or Strong -> not interaction

# Adjacency matrix for all types of evidence, including null (all)
tf_gene_all_adj <- populate(adj_matrix = tf_gene_adj_mat,
                            hash = all.hash,
                            network_df = tf_gene_net_df)
all_adj <- simplify_matrix(tf_gene_all_adj) # removing 
# Adjacency matrix for Weak, Strong and Confirmed evidences (wsc)
tf_gene_wsc_adj <- populate(adj_matrix = tf_gene_adj_mat,
                            hash = w.s.c.hash,
                            network_df = tf_gene_net_df)
wsc_adj <- simplify_matrix(tf_gene_wsc_adj)
# Adjacency matrix for Strong and Confirmed evidences (sc)
tf_gene_sc_adj <- populate(adj_matrix = tf_gene_adj_mat,
                           hash = s.c.hash,
                           network_df = tf_gene_net_df)
sc_adj <- simplify_matrix(tf_gene_sc_adj)
# Adjacency matrix for Confirmed evidences (conf)
tf_gene_conf_adj <- populate(adj_matrix = tf_gene_adj_mat,
                             hash = conf.hash,
                             network_df = tf_gene_net_df)
conf_adj <- simplify_matrix(tf_gene_conf_adj)
# conf_adj is empty, there are no Confirmed evidences, only Weak and Strong


### Visualizing the networks:
wsc_graph <- graph.adjacency(wsc_adj, mode="directed")
sc_graph <- graph.adjacency(sc_adj, mode="directed")

# Check if the node is a TF (True) or a non-TF gene (False):
node_is_TF <- function(vertices, source_names) {
  result_vector <- logical(length(vertices)) # default FALSE
  names(result_vector) <- vertices
  # for each node
  for (node in vertices) {
    # if it is a TF (it is in source_names)
    if (node %in% source_names) {
      result_vector[node] <- TRUE
    }
  }
  return(result_vector)
}

# Coloring the TFs as purple and the non-TF genes as orange
V(wsc_graph)$color <- ifelse(node_is_TF(V(wsc_graph)$name, source_names), "purple", "orange")
V(sc_graph)$color <- ifelse(node_is_TF(V(sc_graph)$name, source_names), "purple", "orange")

# Plotting the graphs:
# I'm not plotting all_graph because it is too complex to be seen properly
png("wsc_graph.png", width = 1080, height = 1080, units = "px")
plot(wsc_graph, edge.arrow.size = 0.7, vertex.label = NA, vertex.size = 5, margin = 0.3)
dev.off()
png("sc_graph.png", width = 1080, height = 1080, units = "px")
plot(sc_graph, edge.arrow.size = 0.7, vertex.label = NA, vertex.size = 5, margin = 0.3)
dev.off()

########## - Applying the algorithms - ##########

### Evidence: null, Weak, Strong and Confirmed
# Characteristics of the natural network:
all_size <- nrow(all_adj) # number of nodes
all_n_edges <- sum(all_adj) # number of edges

# number of each effect: positive ('+'),  negative ('-') or other (different from '+' or '-')
all_effect_vec <- tf_gene_net_df[tf_gene_net_df$Evidence.type == "null" |
                               tf_gene_net_df$Evidence.type == "Weak" |
                               tf_gene_net_df$Evidence.type == "Strong" |
                               tf_gene_net_df$Evidence.type == "Confirmed", ]$Effect

all_pos_effect <- sum(all_effect_vec == "+")
all_neg_effect <- sum(all_effect_vec == "-")
all_total_effect <- length(all_effect_vec)
all_other_effect <- all_total_effect - (all_pos_effect + all_neg_effect)

# Finding the number of autoregulations
all_autoreg <- find_autoreg(all_adj)
all_num_autoreg <- count_autoreg(all_autoreg, tf_gene_net_df)
# Finding the number of feed-forward loops
all_FFL <- find_feed_forw(all_adj)
all_num_FFL <- count_feed_forw(all_FFL, tf_gene_net_df)
# Finding the length of the longest path
all_longestpath <- find_longest_path(all_adj)

# Random networks:
all_results = simulate_and_analyze_networks(num_sim = 1,
                                            adj_matrix = all_adj,
                                            pos_freq = all_pos_effect / all_total_effect,
                                            neg_freq = all_neg_effect / all_total_effect,
                                            other_freq = all_other_effect / all_total_effect)

### Evidence: Weak, Strong and Confirmed
# Characteristics of the natural network:
wsc_size <- nrow(wsc_adj) # number of nodes
wsc_n_edges <- sum(wsc_adj) # number of edges

# number of each effect: positive ('+'),  negative ('-') or other (different from '+' or '-')
wsc_effect_vec <- tf_gene_net_df[tf_gene_net_df$Evidence.type == "Weak" | 
                                tf_gene_net_df$Evidence.type == "Strong" | 
                                tf_gene_net_df$Evidence.type == "Confirmed", ]$Effect

wsc_pos_effect <- sum(wsc_effect_vec == "+")
wsc_neg_effect <- sum(wsc_effect_vec == "-")
wsc_total_effect <- length(wsc_effect_vec)
wsc_other_effect <- wsc_total_effect - (wsc_pos_effect + wsc_neg_effect)

# Finding the number of autoregulations
wsc_autoreg <- find_autoreg(wsc_adj)
wsc_num_autoreg <- count_autoreg(wsc_autoreg, tf_gene_net_df)
# Finding the number of feed-forward loops
wsc_FFL <- find_feed_forw(wsc_adj)
wsc_num_FFL <- count_feed_forw(wsc_FFL, tf_gene_net_df)
# Finding the length of the longest path
wsc_longestpath <- find_longest_path(wsc_adj)

# Random networks:
wsc_results = simulate_and_analyze_networks(num_sim = 100,
                                            adj_matrix = wsc_adj,
                                            pos_freq = wsc_pos_effect / wsc_total_effect,
                                            neg_freq = wsc_neg_effect / wsc_total_effect,
                                            other_freq = wsc_other_effect / wsc_total_effect,
                                            plot_name = "wsc_random.png",
                                            plot_graph = TRUE)

### Evidence: Strong and Confirmed
# Characteristics of the natural network:
sc_size <- nrow(sc_adj) # number of nodes
sc_n_edges <- sum(sc_adj) # number of edges

# number of each effect: positive ('+'),  negative ('-') or other (different from '+' or '-')
sc_effect_vec <- tf_gene_net_df[tf_gene_net_df$Evidence.type == "Strong" | 
                                   tf_gene_net_df$Evidence.type == "Confirmed", ]$Effect

sc_pos_effect <- sum(sc_effect_vec == "+")
sc_neg_effect <- sum(sc_effect_vec == "-")
sc_total_effect <- length(sc_effect_vec)
sc_other_effect <- sc_total_effect - (sc_pos_effect + sc_neg_effect)

# Finding the number of autoregulations
sc_autoreg <- find_autoreg(sc_adj)
sc_num_autoreg <- count_autoreg(sc_autoreg, tf_gene_net_df)
# Finding the number of feed-forward loops
sc_FFL <- find_feed_forw(sc_adj)
sc_num_FFL <- count_feed_forw(sc_FFL, tf_gene_net_df)
# Finding the length of the longest path
sc_longestpath <- find_longest_path(sc_adj)

# Random networks:
sc_results = simulate_and_analyze_networks(num_sim = 100,
                                            adj_matrix = sc_adj,
                                            pos_freq = sc_pos_effect / sc_total_effect,
                                            neg_freq = sc_neg_effect / sc_total_effect,
                                            other_freq = sc_other_effect / sc_total_effect,
                                            plot_name = "sc_random.png",
                                            plot_graph = TRUE)

### Save the results to a file
if (!require('data.table')) install.packages('data.table'); library('data.table')

# Real networks, I'm concatenating the information of each feature into a single file
# in order: all, wsc, sc
real_net_autoreg <- rbind(all_num_autoreg, wsc_num_autoreg, sc_num_autoreg)
real_net_FFL <- rbind(all_num_FFL, wsc_num_FFL, sc_num_FFL)
real_net_path <- rbind(all_longestpath, wsc_longestpath, sc_longestpath)

fwrite(real_net_autoreg, "real_net_autoreg.csv")
fwrite(real_net_FFL, "real_net_FFL.csv")
fwrite(as.list(real_net_path), "real_net_path.csv")

# Simulated networks, each result in a different file
# Autoregulation:
fwrite(all_results$results_autoreg, "all_random_autoreg.csv")
fwrite(wsc_results$results_autoreg, "wsc_random_autoreg.csv")
fwrite(sc_results$results_autoreg, "sc_random_autoreg.csv")
# FFL:
fwrite(all_results$results_FFL, "all_random_FFL.csv")
fwrite(wsc_results$results_FFL, "wsc_random_FFL.csv")
fwrite(sc_results$results_FFL, "sc_random_FFL.csv")
# Longest path:
fwrite(as.list(all_results$results_path), "all_random_path.csv")
fwrite(as.list(wsc_results$results_path), "wsc_random_path.csv")
fwrite(as.list(sc_results$results_path), "sc_random_path.csv")



