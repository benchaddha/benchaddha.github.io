# Exponential Random Graph Models (ERGMs)
#
# Install packages below if needed:
# -------------------------------------------------
if (!"statnet" %in% installed.packages()) install.packages("statnet") # For fitting ERGMs
if (!"igraph" %in% installed.packages()) install.packages("igraph") # For network plotting
if (!"texreg" %in% installed.packages()) install.packages("texreg") # For printing "nicer" model output

library(statnet)
?statnet # About the package
list.files()

# ----------------------------------------------------------------------------------------------------
######################## PART I: Building and Visualizing the Networks ########################
# ----------------------------------------------------------------------------------------------------
# Dependent variable:
# Responses to the question:
adviceEdgelist <- read.csv("adviceEdgelist.csv")
# View the first rows of the edgelist to make sure it imported correctly:
head(adviceEdgelist)
# Convert the edgelist to a network object in statnet format:
advice <- as.network.matrix(adviceEdgelist, matrix.type = "edgelist") 
advice # View a summary of the network object

# Independent variables:
# Load node attributes, and store them in the advice network object we have created
set.vertex.attribute(advice, "department",read.csv("departmentNode.csv",stringsAsFactors=FALSE)$department) # Categorical variable for department
set.vertex.attribute(advice, "leader",read.csv("leaderNode.csv")$leader) # Indicator variable for department leader
set.vertex.attribute(advice, "tenure",read.csv("tenureNode.csv")$tenure) # Years tenure
set.vertex.attribute(advice, "office",read.csv("officeNode.csv")$office) # Indicator variable for whether they are located in the main or secondary office
set.vertex.attribute(advice, "female",read.csv("femaleNode.csv")$female) # Indicator variable for female vs. male
advice # These five variables should now be listed as vertex attributes when viewing the summary of the network

# Double-check the values for all of the node variables
get.vertex.attribute(advice,"department")
get.vertex.attribute(advice,"leader")
get.vertex.attribute(advice,"tenure")
get.vertex.attribute(advice,"office")
get.vertex.attribute(advice,"female")


# Finally, we will import data on the counts of direct messages sent between two employees
# The file "messageEdgelist.csv" contains a messaging edgelist, that we will convert to a matrix
# Statnet allows adjacency matrices with valued ties to be used as predictors (edge covariates) in ERGMs
messageEdgelist <- read.csv("messageEdgelist.csv")
head(messageEdgelist) # Check the first five rows of the edgelist. Third column is the message count
messages <- matrix(nrow = 66, ncol = 66) # Number of Direct messages sent from i to j
for (i in 1:nrow(messageEdgelist)) { # Read from edgelist
  messages[messageEdgelist$SenderId[i], messageEdgelist$ReceiverId[i] ] <- as.numeric(messageEdgelist$MessagesSent[i])
}
for (i in 1:66) { # Remove self-ties (messages sent to self)
  messages[i,i] <- as.numeric(0)
}
hundreds_messages <- messages / 100 # Change weights to represent hundred of messages sent
                                   # This will make viewing/interpreting ergm coefficients easier

# ---------------------------------------------------------------------------------------
# Basic descriptive information
# ---------------------------------------------------------------------------------------
summary(advice, print.adj = FALSE)           # summarize the advice From You network

# ---------------------------------------------------------------------------------------
# Visualize networks
# ---------------------------------------------------------------------------------------
library('igraph') # Ignore messages on any objects that are masked

# Set default plot options
igraph_options(vertex.size = 9, vertex.color = 'grey', # vertex.size changes the size of nodes; vertex.color changes the color of nodes
               edge.color='gray80', edge.arrow.size=.1, # edge.color changes the color of ties; edge.arrow.size changes the size of tie arrow heads
               vertex.label = NA)                       # vertex.label = NA specifies not to display vertex labels in the plot

# Plot the Advice network
advice_igraph <- graph.adjacency(as.matrix.network(advice)) # make an igraph network object from statnet network object
advice_igraph <- set_vertex_attr(advice_igraph,"female",value = read.csv("femaleNode.csv")$female)
net_layout <- layout_with_fr(advice_igraph) # Calculates and stores a spring-embedded layout
                                           # We will re-use this same layout for each plot, so nodes are always in the same place
plot(advice_igraph, layout=net_layout, edge.color='black', vertex.label = V(advice_igraph))

# Plot the Advice network with node coloring based on sex
V(advice_igraph)$color = ifelse (V(advice_igraph)$female ==1, " orange ", " grey ")
plot(advice_igraph, layout=net_layout, edge.color='black', vertex.label = V(advice_igraph))

# Plot the network of who messages whom
# One unit of weight is 100 messages
messages_igraph <- graph.adjacency(messages, weighted = TRUE) ## weighted = TRUE creates edge weight in the igraph object
## below because it is added not only edge weight, but also changed the transparency of edges, it will take more time to plot
## also, if you use zoom in RStudio, it may take about 1 minute until you see the plot.
plot(messages_igraph, layout=net_layout, edge.color = adjustcolor('blue',alpha=.2), vertex.label = V(advice_igraph), edge.width=log(E(messages_igraph)$weight), edge.arrow.width =.2)

# -------------------------------------------------------------------------------------------------
######################## PART II: Build the ERGM models ########################
#
# -------------------------------------------------------------------------------------------------
detach(package:igraph) # Remove the 'igraph' package from your environment. 
library(statnet)
options(ergm.loglik.warn_dyads=FALSE) #Whether or not a warning should be issued when sample space constraints render the observed number of dyads ill-defined

# Ergm Terms are statistics: They are some deterministic function of the ties, node attributes, and edge covariates of a network.
help("ergm-terms",package = "ergm") # Documentation that contains definitions for all of the terms we are using
                                    # ex. what does "mutual" test and how is it calculated
# We will use the ergm-terms to perform hypothesis testing using ERGMs
# But we can note that any of the ERGM terms can also be examined directly for your observed network, by creating a formula in R

# Look at Endogenous statistics: terms based on only ties in the advice network
summary(advice ~ edges)                     # Number of edges (ties)
summary(advice ~ mutual)                    # Number of pairs of reciprocated ties
summary(advice ~ odegree(0:5))              # Outdegree distribution. (# of nodes with outdegree of 0, # nodes outdegree of 1, etc.)
                                           # Remember, respondents could nominate at most five employees in our survey
summary(advice ~ idegree(0:65))             # Indegree distribution.
summary(advice ~ gwodegree(log(2),fixed=T)) # One parameter summarizing outdegree distribution - tendency against outdegree hubs
summary(advice ~ gwidegree(log(2),fixed=T)) # One parameters summarizing indegree distribution - tendency against indegree hubs
summary(advice ~ desp(1:5))                 # Pairs of nodes with one shared partner, two shared partners, etc.
summary(advice ~ dgwesp(log(2),fixed = T))  # One parameter summarizing 

# Look at Exogenous statistics: terms based on advice ties AND other ties / node attributes
summary(advice ~ nodeicov("office"))             # Ties directed towards employees at the main office (as opposed to secondary office)
summary(advice ~ nodeocov("office"))             # Ties originating from employees at the main office (as opposed to secondary office)
summary(advice ~ nodematch("female"))            # Number of ties between people of the same sex
summary(advice ~ nodematch("department"))        # Number of ties between people working in the same department
summary(advice ~ nodemix("leader",levels2=NULL)) # Number of ties between different combinations of leaders(1) and non-leaders(0)
summary(advice ~ diff("tenure"))                 # Total difference in tenure: sum of (sender's tenure - receivers-tenure) for all ties
summary(advice ~ edgecov(hundreds_messages))     # Total messages sent: sum of (messages sent from sender to receiver)/100 for all Advice ties
                                                # e.g., a total of 5669 messages were sent from employees to those they go to for advice during the observed period

# The above are statistics - counts of these patterns for our networks
# What fitting the ERGM model will tell is whether these counts are relatively high/low
# in comparison to what we would expect based on random chance, controlling for the other effects in our model.
# This type of analysis can be helpful for understanding your network, as well as troubleshooting issues with ERGM regression


# The following commands do model estimation for ERGMs.
# This may take a second. Text will print in-console to update you on progress in model estimation.
model1 <- ergm(advice ~ edges                 # This is  a tendency towards a greater number of advice ties existing. Based on a statistic counting the number of ties.
               # Structural patterns
               + mutual                      # This is a tendency towards reciprocity for the advice ties. Based on a statistic counting the number of reciprocated ties.
               + edgecov(hundreds_messages)  # This is the effect of every 100 messages sent from i->j on likelihood of an advice tie. Based on a weighted sum of advice ties x 100s of messages sent
               + nodemix("leader",base = 3)
               # Model constraints
               , constraints =~ bd(maxout=5) # This constraint enforces the maximum outdegree is 5
) 
summary(model1) 

## Convert a log-odds (e.g., -0.68122) ratio to an odds ratio
exp(-0.68122)


model2 <- ergm(advice ~  # This model will be slower to estimate than model 1
                        # Expect roughly 2-7 minutes. If it gets stuck for longer than that, try hitting "stop" and re-running it
               # Structural patterns
               # edges
               mutual
               + gwidegree(log(2), fixed = T)                 # Inverted preferential attachment (indegree)
               + gwodegree(2, fixed = T, cutoff = 5)              # Inverted preferential attachment (outdegree)
               + dgwesp(log(2), type = "OTP", fixed = T, cutoff =5)    # A modified version of Outgoing Two Path(i->j + i->k->j) structures. Geometrically weighted version of transitivity
               # Node attribute effects
               + nodematch("female")                                   # Homophily on a categorical variable 
               + nodemix("leader", base = 3)                           # Mixing matrix of all different combinations of node attributes (ex. A -> A ties, A-> B ties, B -> A ties, B -> B ties). 
               + nodematch("department") 
               + nodeicov("office")                                    # Covariance between in-degree of nodes and attributes of nodes
               + nodeocov("office")                                    # Covariance between out-degree of nodes and attributes of nodes
               + diff("tenure")                                        # Difference is computed as (tenure_i - tenure_j) i: sending node, j: receiving node
               + edgecov(hundreds_messages)                            # Covariance between edges of two networks (predictor can be continous)
               # Constraints on network
               , constraints =~ bd(maxout=5)                           # This constraint enforces the maximum outdegree is 5
               # Control settings for MCMC-MLE algorithm
               , control = control.ergm(MCMC.effectiveSize = 50)
) 
summary(model2) 

# Easy side-by-side model comparison:
library(texreg)
screenreg(list("model1"=model1,"model2"=model2))
# Screenshot this output and include it in your report

# Ergm uses a random algorithm, so every time you estimate the model you may get slightly different results
# Thus, you may want to save your R environment, into a .Rdata file:
# save.image("Lab2_files.RData")
# And then you can re-load your results into R later by clicking to open the file, or by running the line of code below:
# load("Lab2_files.RData")

# -------------------------------------------------------------------------------------------------
######################## PART III: Model diagnostics ########################
# -------------------------------------------------------------------------------------------------
# MCMC diagnostics - will save to a pdf in the current working directory

pdf('model1diagnostics.pdf')              # Open a pdf file to save to
mcmc.diagnostics(model1) # Run the markov chain monte carlo diagnostics
dev.off()                # Closes and saves the pdf

pdf('model2diagnostics.pdf')              # Open a pdf file to save to
mcmc.diagnostics(model2) # Run the markov chain monte carlo diagnostics
dev.off()                # Closes and saves the pdf

# -------------------------------------------------------------------------------------------------
# Goodness of fit test - will display in RStudio
# Check how well the estimated model captures certain features of the observed network, for example triangles in the network.
# -------------------------------------------------------------------------------------------------
# Look at networks simulated according to model 1
# This first command simulates 100 networks.
# These networks, if we use sufficient burnin steps in the markov chain used to generate them,
# may be thought of as random samples from the joint probability distribution that is our fitted ERGM.
sim1 <- simulate(model1, burnin=100000, interval=100000, nsim=100, verbose=T)  # Uses the ergm model to simulate a null model
# Plot the first of the simulated networks
sim1_net1 <- igraph::graph.adjacency(as.matrix.network(sim1[[1]]))
igraph::plot.igraph(sim1_net1,layout=net_layout,edge.color="brown",  
                    vertex.color = 'grey',edge.arrow.size=.1)                                                               
# Plot the 10th simulated network
sim1_net10 <- igraph::graph.adjacency(as.matrix.network(sim1[[10]]))
igraph::plot.igraph(sim1_net10,layout=net_layout,edge.color="red",  
                    vertex.color = 'grey',edge.arrow.size=.1)

# Repeat, now looking at networks simulated according to model 2
sim2 <- simulate(model2, burnin=100000, interval=100000, nsim=100, verbose=T)  # Uses the ergm model to simulate a null model
# Plot the first of the simulated networks
sim2_net1 <- igraph::graph.adjacency(as.matrix.network(sim2[[1]]))
igraph::plot.igraph(sim2_net1,layout=net_layout,edge.color="grey",  
                    vertex.color = 'grey',edge.arrow.size=.1)                                                               
# Plot the 10th simulated network
sim2_net10 <- igraph::graph.adjacency(as.matrix.network(sim2[[10]]))
igraph::plot.igraph(sim2_net10,layout=net_layout,edge.color="purple",  
                    vertex.color = 'grey',edge.arrow.size=.1)

# -------------------------------------------------------------------------------------------------
# Extract the number of triangles from each of the 100 samples and
# compare the distribution of triangles in the sampled networks with the observed network
# -------------------------------------------------------------------------------------------------
# Model 1:
model1.tridist <- sapply(1:100, function(x) summary(sim1[[x]] ~triangle)) # Extracts the triangle data from the simulated networks
hist(model1.tridist,xlim=c(0,1000),breaks=10)                             # Plots that triangle distribution as a histogram, change xlim to change the x-axis range if necessary
advice.tri <- summary(advice ~ triangle)                                    # Stores the number of observed triangles
advice.tri
arrows(advice.tri,20, advice.tri, 0.5, col="red", lwd=3)                      # Adds an arrow to the plotted histogram
c(obs=advice.tri,mean=mean(model1.tridist),sd=sd(model1.tridist),
  tstat=abs(mean(model1.tridist)-advice.tri)/sd(model1.tridist))

# Model 2:
model2.tridist <- sapply(1:100, function(x) summary(sim2[[x]] ~triangle)) # Extracts the triangle data from the simulated networks
hist(model2.tridist,xlim=c(0,1000),breaks=10)                             # Plots that triangle distribution as a histogram, change xlim to change the x-axis range if necessary
arrows(advice.tri,20, advice.tri, 0.5, col="red", lwd=3)                    # Adds an arrow to the plotted histogram
c(obs=advice.tri,mean=mean(model2.tridist),sd=sd(model2.tridist),
  tstat=abs(mean(model2.tridist)-advice.tri)/sd(model2.tridist))

# -------------------------------------------------------------------------------------------------
# Test the goodness of fit of the model
# Compiles statistics for these simulations as well as the observed network, and calculates p-values 
# -------------------------------------------------------------------------------------------------

# Model 1:
# It may take a second for this command to run.
gof1 <- gof(model1, verbose=T, burnin=1e+5, interval=1e+5, control = control.gof.ergm(nsim = 200))
# If you run below and then wouldn't see the plot, trypar(mar=c(2,2,2,2))
dev.off()           # Clear any other plots from the plot window
plot(gof1)          # Plot the goodness of fit
                    # Note: This should produce five separate plots that you should look through.
                    #       In RStudio, scroll between the plots using the arrow buttons
gof1                # Display the goodness of fit info in the console

# Model 2:
# It may take a second for this command to run.
gof2 <- gof(model2, verbose=T, burnin=1e+5, interval=1e+5, control = control.gof.ergm(nsim = 200))
dev.off()
plot(gof2)
gof2
