# install required packages
list.of.packages <- c("ggraph","igraph","readr","robustbase","statnet")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(readr)
library(ggraph)
library(igraph)

# load in data
body_df <-readr::read_tsv("soc-redditHyperlinks-body.tsv")
title_df <-readr::read_tsv("soc-redditHyperlinks-title.tsv")

# set working directory
setwd("C:/Users/17145/Documents/IEMS_341")

# -------------------------------------------------------------------------------------------------
######################## Finding Top 1000 Nodes ########################
# -------------------------------------------------------------------------------------------------

# Finding the top 1000
# let's find the 1000 most active communities and how they engage with one another 
by_source <- body_df %>% group_by(SOURCE_SUBREDDIT)
by_receiver <- body_df %>% group_by(TARGET_SUBREDDIT)

# Indegree 
# 1 post can impact multiple subreddits 
send_count <- by_source %>% summarise(
  sends = n()
)

# OutDegree
receive_count <- by_receiver %>%  summarise(
  receptions = n()
)

# merge the two right_join
engagement_subreddits <- merge(send_count, receive_count, by.x = "SOURCE_SUBREDDIT", by.y = "TARGET_SUBREDDIT", all = TRUE)

# replace all NA with 0
engagement_subreddits[is.na(engagement_subreddits)] <- 0

# rename subreddit
colnames(engagement_subreddits)[1] = "subreddit"

# add up sends and receptions to get engagement
engagement_subreddits$total_engagement <- engagement_subreddits$sends + engagement_subreddits$receptions

# reorder by most total_engagement
engagement_subreddits <- engagement_subreddits[order(-engagement_subreddits$total_engagement), ]

# create 3 different subreddit cutoffs 
top1000 <- engagement_subreddits[1:1000, ]
top100 <- engagement_subreddits[1:100, ]
top50 <- engagement_subreddits[1:50, ]

# only keep most engaged n-subreddits
subs_1000 <- body_df %>%  filter(SOURCE_SUBREDDIT %in% c(top1000[, 1]) & TARGET_SUBREDDIT %in% c(top1000[, 1]))
subs_100 <- body_df %>%  filter(SOURCE_SUBREDDIT %in% c(top100[, 1]) & TARGET_SUBREDDIT %in% c(top100[, 1]))
subs_50 <- body_df %>%  filter(SOURCE_SUBREDDIT %in% c(top50[, 1]) & TARGET_SUBREDDIT %in% c(top50[, 1]))

# count number of actually unique nodes 
# 998 for subs_1000
length(unique(c(unique(subs_1000$SOURCE_SUBREDDIT), unique(subs_1000$TARGET_SUBREDDIT))))
# 100 for subs_100
length(unique(c(unique(subs_100$SOURCE_SUBREDDIT), unique(subs_100$TARGET_SUBREDDIT))))
# 49 for subs_50
length(unique(c(unique(subs_50$SOURCE_SUBREDDIT), unique(subs_50$TARGET_SUBREDDIT))))

# lowest 50 of 1000
# Lowest Total Engagement is 96
top1000[950:1000, ]

# -------------------------------------------------------------------------------------------------
######################## Data Aggregating ########################
# -------------------------------------------------------------------------------------------------

# Data aggregating

# INPUT IN SUB CUTOFF YOU WANT
# group by unique combos of source/target subs
unique_combos_1000 <- subs_1000 %>% group_by(SOURCE_SUBREDDIT, TARGET_SUBREDDIT)
unique_combos_100 <- subs_100 %>% group_by(SOURCE_SUBREDDIT, TARGET_SUBREDDIT)
unique_combos_50 <- subs_50 %>% group_by(SOURCE_SUBREDDIT, TARGET_SUBREDDIT)

# relevant stats = average sentiment sent over
# 1. # of contacts 
# 2. average sentiment (mean of all sentiments combined)

relevant_stats_1000 <- unique_combos_1000 %>% summarise(
  contacts = n(),
  average_sent = mean(LINK_SENTIMENT)
)

relevant_stats_100 <- unique_combos_100 %>% summarise(
  contacts = n(),
  average_sent = mean(LINK_SENTIMENT)
)

relevant_stats_50 <- unique_combos_50 %>% summarise(
  contacts = n(),
  average_sent = mean(LINK_SENTIMENT)
)

# create dataframe of subreddit combos
# relevant_stats <- relevant_stats[relevant_stats$average_sent <= 0, ]
sub_combos_1000 <- relevant_stats_1000[,1:2]
sub_combos_100 <- relevant_stats_100[,1:2]
sub_combos_50 <- relevant_stats_50[,1:2]

print(relevant_stats_1000)


# -------------------------------------------------------------------------------------------------
######################## Plotting Networks ########################
# -------------------------------------------------------------------------------------------------

# With Top 1000
# create directed, weighted network graphs (lab2)
network_graph_1000 <- as.network.matrix(sub_combos_1000, matrix.type = "edgelist") 
set.edge.attribute(network_graph_1000, "contacts", relevant_stats_1000$contacts)
set.edge.attribute(network_graph_1000, "sentiment", relevant_stats_1000$average_sent)

# Set default plot options
igraph_options(vertex.size = 12, # vertex.size changes the size of nodes; 
               vertex.color = 'grey', #vertex.color changes the color of nodes
               edge.color='gray80', # edge.color changes the color of ties;
               edge.arrow.size=.1,  # edge.arrow.size changes the size of tie arrow heads
               vertex.label = NA)   # vertex.label = NA specifies not to display vertex labels in the plot

# Plot the Subreddit network
network_igraph_1000 <- graph.adjacency(as.matrix.network(network_graph_1000)) # make an igraph network object from statnet network object

network_igraph_1000 <- set_edge_attr(graph = network_igraph_1000, name = "contacts", value = relevant_stats_1000$contacts)
network_igraph_1000 <- set_edge_attr(graph = network_igraph_1000, name = "sentiment", value = relevant_stats_1000$average_sent)


net_layout_1000 <- layout_with_fr(network_igraph_1000) # Calculates and stores a spring-embedded layout

# Network with Top 1000 Nodes [Plot 1]
plot(network_igraph_1000,
     layout=net_layout_1000,
     edge.color='black', # edge.color changes the color of ties;
     vertex.size = 2, # vertex.size changes the size of nodes; 
     vertex.color = 'purple', #vertex.color changes the color of nodes
     edge.arrow.size=.1,  # edge.arrow.size changes the size of tie arrow heads
     vertex.label = NA  # vertex.label = NA specifies not to display vertex labels in the plot
     )

# Plot the Advice network with node coloring based on average sentiment 
E(network_igraph_1000)$color = ifelse (E(network_igraph_1000)$sentiment < 0, " red ", " darkseagreen ")

# Network with Top 1000 Nodes (Average Sentiment) [Plot 2]
plot(network_igraph_1000,
     layout=net_layout_1000,
     vertex.size = 2, # vertex.size changes the size of nodes; 
     vertex.color = 'purple', #vertex.color changes the color of nodes
     edge.arrow.size=.1,  # edge.arrow.size changes the size of tie arrow heads
     vertex.label = NA  # vertex.label = NA specifies not to display vertex labels in the plot
)

###############################################################################################

# With Top 100
# create directed, weighted network graphs (lab2)
network_graph_100 <- as.network.matrix(sub_combos_100, matrix.type = "edgelist") 
set.edge.attribute(network_graph_100, "contacts", relevant_stats_100$contacts)
set.edge.attribute(network_graph_100, "sentiment", relevant_stats_100$average_sent)

# Set default plot options
igraph_options(vertex.size = 12, # vertex.size changes the size of nodes; 
               vertex.color = 'grey', #vertex.color changes the color of nodes
               edge.color='gray80', # edge.color changes the color of ties;
               edge.arrow.size=.1,  # edge.arrow.size changes the size of tie arrow heads
               vertex.label = NA)   # vertex.label = NA specifies not to display vertex labels in the plot

# Plot the Subreddit network
network_igraph_100 <- graph.adjacency(as.matrix.network(network_graph_100)) # make an igraph network object from statnet network object

network_igraph_100 <- set_edge_attr(graph = network_igraph_100, name = "contacts", value = relevant_stats_100$contacts)
network_igraph_100 <- set_edge_attr(graph = network_igraph_100, name = "sentiment", value = relevant_stats_100$average_sent)


net_layout_100 <- layout_with_fr(network_igraph_100) # Calculates and stores a spring-embedded layout

# Unused
plot(network_igraph_100,
     layout=net_layout_100,
     edge.color='black', # edge.color changes the color of ties;
     vertex.size = 3, # vertex.size changes the size of nodes; 
     vertex.color = 'purple', #vertex.color changes the color of nodes
     edge.arrow.size=.1,  # edge.arrow.size changes the size of tie arrow heads
     vertex.label = NA  # vertex.label = NA specifies not to display vertex labels in the plot
)

# Plot the Advice network with node coloring based on average sentiment 
E(network_igraph_100)$color = ifelse (E(network_igraph_100)$sentiment < 0, " red ", " darkseagreen ")

# Network with Top 100 Nodes only (Average Sentiment) [Plot 3]
plot(network_igraph_100,
     layout=net_layout_100,
     vertex.size = 3, # vertex.size changes the size of nodes; 
     vertex.color = 'purple', #vertex.color changes the color of nodes
     edge.arrow.size=.1,  # edge.arrow.size changes the size of tie arrow heads
     vertex.label = NA  # vertex.label = NA specifies not to display vertex labels in the plot
)


####################################################################################

# With Top 50
# create directed, weighted network graphs (lab2)
network_graph_50 <- as.network.matrix(sub_combos_50, matrix.type = "edgelist") 
set.edge.attribute(network_graph_50, "contacts", relevant_stats_50$contacts)
set.edge.attribute(network_graph_50, "sentiment", relevant_stats_50$average_sent)

# Set default plot options
igraph_options(vertex.size = 12, # vertex.size changes the size of nodes; 
               vertex.color = 'grey', #vertex.color changes the color of nodes
               edge.color='gray80', # edge.color changes the color of ties;
               edge.arrow.size=.1,  # edge.arrow.size changes the size of tie arrow heads
               vertex.label = NA)   # vertex.label = NA specifies not to display vertex labels in the plot

# Plot the Subreddit network
network_igraph_50 <- graph.adjacency(as.matrix.network(network_graph_50)) # make an igraph network object from statnet network object

network_igraph_50 <- set_edge_attr(graph = network_igraph_50, name = "contacts", value = relevant_stats_50$contacts)
network_igraph_50 <- set_edge_attr(graph = network_igraph_50, name = "sentiment", value = relevant_stats_50$average_sent)

net_layout_50 <- layout_with_fr(network_igraph_50) # Calculates and stores a spring-embedded layout


# Unused
plot(network_igraph_50,
     layout=net_layout_50,
     edge.color='black', # edge.color changes the color of ties;
     vertex.size = 3, # vertex.size changes the size of nodes; 
     vertex.color = 'purple', #vertex.color changes the color of nodes
     edge.arrow.size=.1,  # edge.arrow.size changes the size of tie arrow heads
     vertex.label = NA  # vertex.label = NA specifies not to display vertex labels in the plot
)

# Plot the Advice network with node coloring based on average sentiment 
E(network_igraph_50)$color = ifelse (E(network_igraph_50)$sentiment < 0, " red ", " darkseagreen ")

# Network with Top 50 Nodes only (Average Sentiment) [Plot 4]
plot(network_igraph_50,
     layout=net_layout,
     vertex.size = 3, # vertex.size changes the size of nodes; 
     vertex.color = 'purple', #vertex.color changes the color of nodes
     edge.arrow.size=.1,  # edge.arrow.size changes the size of tie arrow heads
     vertex.label = NA  # vertex.label = NA specifies not to display vertex labels in the plot
)


# -------------------------------------------------------------------------------------------------
######################## Analyzing Centrality ########################
# -------------------------------------------------------------------------------------------------

# calculate centrality measures (lab1b)

# convert to "sna" graph object 
sub_graph_sna <- igraph::get.adjacency(network_igraph_1000, sparse=FALSE) %>% network::as.network.matrix()

# indegree measurement 
ideg_subs <- degree(sub_graph_sna, cmode = 'indegree')

# Store the information in centralities_sub dataframe
centralities_subs <- data.frame('node_name' = as.character(network.vertex.names(sub_graph_sna)),
                                'in_degree' = ideg_subs)

# Calculate out-degree centrality and store it in the data.frame called 'centralities'
centralities_subs$out_degree <- degree(sub_graph_sna, cmode = 'outdegree')

# Calculate betweenness centrality and store it in the data.frame called 'centralities'
centralities_subs$betweeness <- betweenness(sub_graph_sna)

#igraph calculations (centrality)

# Calculate closeness centrality and store it in the data.frame called 'centralities'
centralities_subs$incloseness <- igraph::closeness(network_igraph_1000, mode = 'in')
centralities_subs$outcloseness <- igraph::closeness(network_igraph_1000, mode = 'out')


# Calculate eigenvector centrality and store it in the data.frame called 'centralities'
# using 'igraph' because the code implemented in 'sna' is unreliable
centralities_subs$eigen <- igraph::eigen_centrality(network_igraph_1000)$vector

# Calculate Burt's network constraint and store it in the data.frame called 'centralities'
# using 'igraph' because 'sna' doesn't have the function
centralities_subs$netconstraint <- igraph::constraint(network_igraph_1000)

# Calculate authority and store it in the data.frame called 'centralities'
# using 'igraph' because 'sna' doesn't have the function
# 'igraph::' allows calling for any igraph function without loading the package
centralities_subs$authority <- igraph::authority_score(network_igraph_1000, scale = TRUE)$`vector`

# Calculate hub and store it in the data.frame called 'centralities'
# using 'igraph' because 'sna' doesn't have the function
centralities_subs$hub <- igraph::hub_score(network_igraph_1000, scale = TRUE)$`vector`

# read centrality table
View(centralities_subs)

#install.packages("writexl")
library("writexl")
write_xlsx(centralities_subs,"C:/Users/17145/Documents/IEMS_341\\full_dataset_centralities.xlsx")

# -------------------------------------------------------------------------------------------------
######################## ERGM models ########################
# -------------------------------------------------------------------------------------------------


