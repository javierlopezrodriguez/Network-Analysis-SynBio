########## Helper functions for the main program ##########

### Function to populate the matrix:
# We pass it the adjacency matrix, the hash which codes which evidence we'll keep,
# and the network dataframe
populate <- function(adj_matrix, hash, network_df) {
  # for each line in the network dataframe
  for (i in 1:nrow(network_df)) {
    line <- network_df[i, ]
    # we replace the corresponding position in the adj_matrix with the value
    # defined by the evidence type and our hash,
    # if that interaction is not already 1 (to avoid overwriting a 1 with a 0)
    if (adj_matrix[line$Source.name, line$Target.name] == 0) {
      adj_matrix[line$Source.name, line$Target.name] <- hash[[line$Evidence.type]]
    }
  }
  return (adj_matrix)
}

### Simplifying an adjacency matrix:
# Takes an adjacency matrix and removes the ids for which both column and row
# are completely empty
simplify_matrix <- function(adj_matrix) {
  # boolean vector that marks if a row is empty (if all values are 0)
  empty_rows = apply(X = adj_matrix, MARGIN = 1, FUN = function(x) all(x == 0))
  # boolean vector that marks if a column is empty (if all values are 0)
  empty_cols = apply(X = adj_matrix, MARGIN = 2, FUN = function(x) all(x == 0))
  # boolean vector that marks the ids whose rows and columns are empty
  empty_ids = empty_rows & empty_cols
  # indexes vector corresponding to the TRUE elements of empty_ids
  indexes = (1:length(empty_ids))[empty_ids]
  # if there are indexes to remove:
  if (length(indexes) > 0) {
    # removing the corresponding rows and columns from the adj_matrix
    result_matrix <- adj_matrix[-indexes, -indexes]
    # warning if result_matrix is empty
    if(nrow(result_matrix) == 0 | ncol(result_matrix) == 0) {
      print("WARNING: matrix is now empty")
    }
    # returning the matrix
    return (result_matrix)
  } else {
    return(adj_matrix) # returning the unchanged matrix
  }
}



########## - Creating the three algorithms used - ##########

### Checking for autoregulation motifs is simple, we have to see if the diagonal
### (where a node meets itself) has an interaction or not.
# Instead of simply counting, I want to first return all of the autoregulation
# instances, so that I can later check if it is positive or negative.

find_autoreg <- function(adj_matrix) {
  # creating an empty vector to store the results
  autoregulation <- character()
  # for each row of the adjacency matrix
  for (node in rownames(adj_matrix)) {
    # if the node interacts with itself
    if (adj_matrix[node, node] != 0) {
      autoregulation <- c(autoregulation, node) # adding it to the vector
    }
  }
  return(autoregulation)
}

# Function that, given the result of the previous function, counts the total, 
# positive, negative, (and other, if any) autoregulations.
# I'm including in 'other' the cases where, at the same time, 
# A induces itself and A represses itself ('+' and '-'),
# or the cases in which there is an effect different from '+' or '-'.
count_autoreg <- function(autoreg, net_df) {
  # initializing the counters
  total <- length(autoreg)
  positive <- 0
  negative <- 0
  other <- 0
  # for each node
  for (node in autoreg) {
    # getting the type of Effect from the network dataframe (and removing duplicates)
    reg_effect <- unique(net_df[net_df$Source.name == node & net_df$Target.name == node, "Effect"])
    # incrementing the counters
    if (length(reg_effect) == 1 && reg_effect == "+") { # only 1 effect, and that effect is '+'
      positive <- positive + 1
    } else if (length(reg_effect) == 1 && reg_effect == "-") { # only 1 effect, and that effect is '-'
      negative <- negative + 1
    } else { # more than 1 different effects, or effects that aren't '+' or '-'
      other <- other + 1
    }
  }
  # returning a list of the counts
  return (list(total = total, 
              positive = positive, 
              negative = negative, 
              other = other))
}

### Creating an algorithm for finding feed forward loops:

# feed forward loops (FFL): when, at the same time, A->B->C and A->C exist
# and none of the reverse interactions (C->B->A, C->A) exist (it would not be a FFL)
# A: start position
# B: mid position
# C: end position

# for each row of the adjacency matrix (start)
# store the column numbers/names that have interactions (possible mid or end)
# for each row equivalent to those columns (mid)
# store the column numbers/names that have interactions (possible end)
# check that start isn't among these interactions (avoid mid -> start)
# if the start and the mid connect both to the same column
# and this column isn't start or mid
# then we have a feedforward loop: start->mid->end and start->end
# (if end->mid and end->start aren't interactions)

find_feed_forw <- function(adj_matrix) {
  # creating an empty dataframe to store the results.
  feed_forward <- data.frame(Start = character(),
                             Mid = character(),
                             End = character())
  
  # for each row of the adjacency matrix (start)
  for (start_name in rownames(adj_matrix)) {
    
    # getting a boolean vector of the columns that have interactions
    start_int <- adj_matrix[start_name, ] != 0
    # getting the names of the TRUE elements of the vector
    mid_or_end_names <- names(adj_matrix[start_name, ][start_int])
    
    # if there is more than one element (can't have feed forward loop with 0 or 1)
    if (length(mid_or_end_names) > 1) {
      
      # removing start_name if it was in mid_or_end_names
      # (we don't care about autoregulation)
      mid_or_end_names <- mid_or_end_names[mid_or_end_names != start_name]
      
      # for each name in the possible mid or end names
      for(mid_name in mid_or_end_names) {
        
        # getting a boolean vector of the columns that have interactions
        mid_int <- adj_matrix[mid_name, ] != 0
        # getting the names of the TRUE elements of the vector
        end_names <- names(adj_matrix[mid_name, ][mid_int])
        
        # if end_names is not empty
        # and end_names doesn't contain start_name (don't want Mid -> Start interaction)
        if (!is.null(end_names) && !any(end_names == start_name)) {
          
          # removing mid_name if it was in end_names
          end_names <- end_names[end_names != mid_name]
          # if there are still elements
          if (length(end_names) != 0) {
            
            # for each end_name, there can be one feed forward loop
            for (end_name in end_names) {
              
              # if start_name interacts with end_name, we have a feed forward loop
              # if end_name doesn't interact with mid_name nor start_name
              # (don't want End -> Mid and End -> Start interactions)
              if (adj_matrix[start_name, end_name] != 0 & 
                  adj_matrix[end_name, mid_name] == 0 &
                  adj_matrix[end_name, start_name] == 0) {
                
                # adding start_name, mid_name, end_name to the feed_forward df
                position <- nrow(feed_forward) + 1
                feed_forward[position, ] <- list(Start = start_name, 
                                                 Mid = mid_name, 
                                                 End = end_name)
              }
            }
          }
        }
      }
    }
  }
  return (feed_forward) # at the end, returns the data frame with the results
}

# Function that, given the result of the previous function, counts the total, 
# coherent (each of the four types and the total), 
# incoherent (each of the four types and the total), 
# and other (if any) feed-forward loops.

# I'm including in 'other' the cases where, for two genes X and Y,
# the effect is different than '+' or '-', or there are more than one different
# effects at the same time.

count_feed_forw <- function(feed_forw, net_df) {
  
  # Start -> Mid -> End, Start -> End
  # Interactions: SM (Start -> Mid), ME (Mid -> End), SE (Start -> End)
  # Coherent (C), Incoherent (I)
  # C1: SM +, ME +, SE +
  # C2: SM -, ME +, SE -
  # C3: SM +, ME -, SE -
  # C4: SM -, ME -, SE +
  # I1: SM +, ME -, SE +
  # I2: SM -, ME -, SE -
  # I3: SM +, ME +, SE -
  # I4: SM -, ME +, SE +
  
  # initializing the counters
  total <- nrow(feed_forw)
  C1 <- 0
  C2 <- 0
  C3 <- 0
  C4 <- 0
  I1 <- 0
  I2 <- 0
  I3 <- 0
  I4 <- 0
  coherent <- 0
  incoherent <- 0
  other <- 0
  
  if (total > 0) {
    # for each feed forward loop
    for (i in 1:nrow(feed_forw)) {
      row <- feed_forw[i, ]
      
      # getting the effect of the interactions (removing duplicates if any)
      SM <- unique(net_df[net_df$Source.name == row$Start & net_df$Target.name == row$Mid, "Effect"])
      ME <- unique(net_df[net_df$Source.name == row$Mid & net_df$Target.name == row$End, "Effect"])
      SE <- unique(net_df[net_df$Source.name == row$Start & net_df$Target.name == row$End, "Effect"])
      
      # If any of them have more than 1 different effect, classify the FFL as other:
      if (length(SM) > 1 | length(ME) > 1 | length(SE) > 1) { other <- other + 1 }
      else {
        # C1: SM +, ME +, SE +
        if (SM == '+' & ME == '+' & SE == '+') {C1 <- C1 + 1}
        # C2: SM -, ME +, SE -
        else if (SM == '-' & ME == '+' & SE == '-') {C2 <- C2 + 1}
        # C3: SM +, ME -, SE -
        else if (SM == '+' & ME == '-' & SE == '-') {C3 <- C3 + 1}
        # C4: SM -, ME -, SE +
        else if (SM == '-' & ME == '-' & SE == '+') {C4 <- C4 + 1}
        # I1: SM +, ME -, SE +
        else if (SM == '+' & ME == '-' & SE == '+') {I1 <- I1 + 1}
        # I2: SM -, ME -, SE -
        else if (SM == '-' & ME == '-' & SE == '-') {I2 <- I2 + 1}
        # I3: SM +, ME +, SE -
        else if (SM == '+' & ME == '+' & SE == '-') {I3 <- I3 + 1}
        # I4: SM -, ME +, SE +
        else if (SM == '-' & ME == '+' & SE == '+') {I4 <- I4 + 1}
        # else (at least 1 effect is different from '+' or '-')
        else {other <- other + 1}
      }
    }
    # total coherent and incoherent loops
    coherent <- C1 + C2 + C3 + C4
    incoherent <- I1 + I2 + I3 + I4
  }
  
  # returning a list with the counters for every possibility
  return (list(coherent = coherent, C1 = C1, C2 = C2, C3 = C3, C4 = C4,
               incoherent = incoherent, I1 = I1, I2 = I2, I3 = I3, I4 = I4,
               other = other, total = total))
}

### Algorithm for finding the longest path in an adjacency matrix of a directed graph

# If the graph is not acyclic (i.e. there is A->B->C->A), then the problem
# becomes harder because technically the longest path has infinite length.
# To check the longest path with unique nodes, we have to check every path
# (which is computationally expensive).
# (https://cs.stackexchange.com/questions/14991/a-to-find-the-longest-path-in-a-directed-cyclic-graph)


# This function recursively looks for the longest simple path (crossing each node once at most)
# I'm basically copying this pseudocode, adapting it to an adjacency matrix and to R code
# (https://stackoverflow.com/questions/21880419/how-to-find-the-longest-simple-path-in-a-graph?rq=1)
path_fn <- function(node_i, input_list) {
  #input_list contains: adj_matrix, visited, current_path_length, best_path_length
  
  # Setting the starting node (node_i) to "visited" so that we don't cross it again
  input_list$visited[node_i] <- "visited"
  
  # to get the not visited child nodes immediately adjacent to node_i:
  # getting a boolean vector with the interactions of node_i
  row_i <- input_list$adj_matrix[node_i, ] != 0 # interactions
  # getting a boolean vector with the not_visited nodes
  not_visited <- input_list$visited == "not_visited" # not_visited
  # getting the names of the nodes that are both interactions of node_i and "not_visited"
  child_nodes <- names(row_i)[row_i & not_visited]
  
  # If there are possible child nodes
  if (length(child_nodes) != 0) {
    # for each child node node_j
    for (node_j in child_nodes) {
      
      # we're increasing the current path length by 1 (we've moved 1 node forwards)
      input_list$current_path_length <- input_list$current_path_length + 1
      
      # storing the maximum between the current length and the stored max length
      input_list$best_path_length <- max(input_list$current_path_length, 
                                         input_list$best_path_length)
      
      # recursively calling path_fn with the child node node_j, which will do all of the above again
      input_list <- path_fn(node_i = node_j, input_list = input_list)
      
      # after all of the recursive calls have been resolved, we reduce the current length by 1
      # because we have finished with this child node and we will pass to the next child node
      input_list$current_path_length <- input_list$current_path_length - 1
    }
    # at the end of the for loop, we set the starting node to "not_visited" so that
    # we can start again from the next starting node
    input_list$visited[node_i] <- "not_visited"
  }
  # returns all of the parameters as a list
  return (input_list)
}

# The function that calls the previous function
find_longest_path <- function(adj_matrix) {
  # setting the initial best path length to 0
  best_path_length <- 0
  # creating a vector that keeps track of the nodes that have been visited 
  visited <- rep_len("not_visited", nrow(adj_matrix))
  names(visited) <- rownames(adj_matrix)
  # for each starting node:
  for (start_node in rownames(adj_matrix)) {
    # creating input list:
    input_list <- list(adj_matrix = adj_matrix,
                       visited = visited,
                       current_path_length = 0,
                       best_path_length = best_path_length)
    # calling the recursive function defined above
    input_list <- path_fn(node_i = start_node, input_list = input_list)
    best_path_length <- input_list$best_path_length
  }
  return (best_path_length)
}

########## Simulate random networks ##########
simulate_and_analyze_networks <- function(num_sim, adj_matrix, pos_freq, neg_freq, other_freq, 
                                          plot_name = "default.png", plot_graph = FALSE) {
  # num_sim = number of simulations
  # adj_matrix = the adjacency matrix of the real network
  # pos_freq = the frequency of '+' effects in the real network
  # neg_freq = the frequency of '-' effects in the real network
  # other_freq = the frequency of other effects in the real network
  # plot_name = the name of the png file that will be created if the graph is plotted
  # plot_graph = plot a graph of the first random network if TRUE (default FALSE)
  
  # creating the empty result dataframes
  results_autoreg <- data.frame(total = numeric(num_sim), 
                                positive = numeric(num_sim), 
                                negative = numeric(num_sim), 
                                other = numeric(num_sim))
  
  results_FFL <- data.frame(coherent = numeric(num_sim), 
                            C1 = numeric(num_sim), 
                            C2 = numeric(num_sim), 
                            C3 = numeric(num_sim), 
                            C4 = numeric(num_sim),
                            incoherent = numeric(num_sim),
                            I1 = numeric(num_sim),
                            I2 = numeric(num_sim), 
                            I3 = numeric(num_sim),
                            I4 = numeric(num_sim),
                            other = numeric(num_sim),
                            total = numeric(num_sim))
  
  results_path <- numeric(num_sim)
  
  # getting the number of nodes and number of edges of the adjacency matrix
  size = nrow(adj_matrix)
  n_edges = sum(adj_matrix)
  
  # for each simulation
  for (i in 1:num_sim) {
    # random empty (full of 0s) adjacency matrix of the same size and same names as adj_matrix
    rnd_adj <- matrix(0, size, size, dimnames = list(rownames(adj_matrix), colnames(adj_matrix)))
    
    # getting the random row coordinates and random column coordinates that will define the interactions
    random_row_coords <- sample(1:size, n_edges, replace = TRUE)
    random_col_coords <- sample(1:size, n_edges, replace = TRUE)
    
    # joining the row and column coordinates, using that matrix to index rnd_adj,
    # and setting those entries to 1 (interaction)
    rnd_adj[cbind(random_row_coords, random_col_coords)] <- 1
    
    # number of total interactions of the random matrix
    rnd_total_int <- sum(rnd_adj)
    
    # simulating the effect of each interaction ('+', '-' or 'other')
    # with the same frequencies as the ones in the original network
    rnd_effect <- sample(c("+", "-", "other"), 
                         size = rnd_total_int, 
                         prob = c(pos_freq, neg_freq, other_freq), 
                         replace = TRUE)
    # creating an empty dataframe that will act as the random network dataframe
    rnd_effect_df <- data.frame(Source.name = character(),
                                Target.name = character(),
                                Effect = character())
    # input each interaction in the rnd_adj matrix to the random network dataframe,
    # including the randomly generated effect
    for (row_name in rownames(rnd_adj)) {
      for (col_name in colnames(rnd_adj)) {
        if (rnd_adj[row_name, col_name] == 1) {
          index <- nrow(rnd_effect_df) + 1 # 1-indexed
          rnd_effect_df[index, ] <- list(Source.name = row_name, 
                                         Target.name = col_name,
                                         Effect = rnd_effect[index])
        }
      }
    }
    # plotting the graph of the first random network if plot_graph == TRUE
    if (plot_graph == TRUE && i == 1) {
      if (!require('igraph')) install.packages('igraph'); library('igraph')
      rnd_graph <- graph.adjacency(rnd_adj, mode = "directed")
      png(plot_name, width = 1080, height = 1080, units = "px")
      plot(rnd_graph, edge.arrow.size = 0.7, vertex.label = NA, vertex.size = 5, margin = 0.3, vertex.color = "salmon")
      dev.off()
    }
    
    # finding the autoregulations, feed forward loops, and longest path
    #print("autoreg")
    autoreg <- find_autoreg(rnd_adj)
    #print("count_autoreg")
    num_autoreg <- count_autoreg(autoreg, rnd_effect_df)
    #print("FFL")
    FFL <- find_feed_forw(rnd_adj)
    #print("count_FFL")
    num_FFL <- count_feed_forw(FFL, rnd_effect_df)
    #print("path")
    longestpath <- find_longest_path(rnd_adj)
    
    # inputing the results into the results dataframes
    results_autoreg[i, ] <- num_autoreg
    results_FFL[i, ] <- num_FFL
    results_path[i] <- longestpath
  }
  
  # returns a list containing the results dataframes
  return (list(results_autoreg = results_autoreg,
               results_FFL = results_FFL,
               results_path = results_path))
}

