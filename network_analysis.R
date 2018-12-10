library(ggplot2)
library(gridExtra)
library(dplyr)
library(igraph)
library(qgraph)
library(RColorBrewer)
library(relaimpo)

# Read data
imdata <- read.table("imitation_data.csv", sep=',', stringsAsFactors=T, header=T)
imdata <- imdata[imdata$nasality!='ref',] # Items used for tongue height reference will not be used

# Unique identities for each token
imdata$ident <- as.numeric(as.factor(paste0(imdata$imitator, imdata$order)))

# Select random seed for continuity throughout all functions
myseed <- 105

# Run the function for randomized nasal-oral matched pairing
diffdata <- oral.nasal.matching(imdata, myseed)

# Run the function for the relative importance analysis (RIA)
RIA.table <- RIA(diffdata, myseed)

# Run the function to perform vowel-wise cross validation of linear models for each speaker
cross.val   <- lm.cross.validation(diffdata, myseed) 
cross.val.R <- cross.val[1] # Get the correlation coefficients
cross.val.P <- cross.val[2] # Get the p-values

# Summarize articulatory variable rankings
RIA.table %>% group_by(vars, ranks) %>% summarise(n = n()) %>% mutate(freq = 100*n / sum(n))

# Run the function for determining the number of network clusters in 100 random seeds
clust.dat <- network.clusters(RIA.table)

# Print summaries for each vowel pair to determine the cluster frequency for each vowel pair
summary(clust.dat$clusters[clust.dat$vowel=='a'])
summary(clust.dat$clusters[clust.dat$vowel=='e'])
summary(clust.dat$clusters[clust.dat$vowel=='o'])

# Get the most frequent number of network clusters for the three vowel pairs
a.clust <- as.numeric(names(which.max(summary(clust.dat$clusters[clust.dat$vowel=='a']))[1]))
e.clust <- as.numeric(names(which.max(summary(clust.dat$clusters[clust.dat$vowel=='e']))[1]))
o.clust <- as.numeric(names(which.max(summary(clust.dat$clusters[clust.dat$vowel=='o']))[1]))

# Get the separate vowel pair seed results
a.clust.dat <- clust.dat[clust.dat$vowel=='a',]
e.clust.dat <- clust.dat[clust.dat$vowel=='e',]
o.clust.dat <- clust.dat[clust.dat$vowel=='o',]

# Determine which random seeds generate the most frequent network clusters for each of the vowel pairs, individually
a.seeds <- a.clust.dat$seed[a.clust.dat$clusters==a.clust]
e.seeds <- e.clust.dat$seed[e.clust.dat$clusters==e.clust]
o.seeds <- o.clust.dat$seed[o.clust.dat$clusters==o.clust]

# Determine which unique seeds generate the most frequent network clusters for ALL three vowel pairs
matches <- e.seeds[e.seeds %in% a.seeds]
seeds   <- o.seeds[o.seeds %in% matches]


# Pick a new seed for clustering!
myseed <- 3

# Run the function for performing RIA coefficient Euclidean distance measurements and getting average F1 values
ED.output <- ED.measurements(RIA.table, diffdata, myseed)

RIA.table.new <- ED.output[[1]] # New RIA table with clustering info and average F1 values
sim.mats      <- ED.output[[2]] # Distance-based similarity score matricies
sim.mats.sc   <- ED.output[[3]] # Scaled similarity matrices
groupings     <- ED.output[[4]] # Cluster groupings for plotting


# Plot results
mycols  <- brewer.pal(12, "Set3") # Create qualitative color palette for plotting
npars   <- length(unique(diffdata$participant)) # The number of participants is needed for various network plotting variables

# Plot a network graph
# Example given here for the /a/-/~a/ vowel pair
qgraph(sim.mats.sc$a, layout="groups", graph="cor", sampleSize=npars, groups=groupings$a,
       posCol="orange", color=mycols, legend=T, GLratio=10, legend.cex=1,
       layout.par = list(init=matrix(rnorm(npars*2), npars, 2)),
       vsize=7, cut=0, border.width=1.5)


# Plot the average RIA coefficients and average acoustic/articulatory variables by group/cluster
# Example given here for the /a/-/~a/ vowel pair
plotdat <- RIA.table.new$a

# Average RIA coefficients, by group
p1 <- ggplot(plotdat[plotdat$vars!='f1',], aes(x=vars, y=coeffs, fill=cluster2)) + 
  geom_hline(yintercept=0, lty=2) + geom_boxplot(notch=F) + labs(fill="Group") +
  geom_vline(xintercept=1.5) + geom_vline(xintercept=2.5) +
  xlab('Articulatory variable') + ylab('Relative importance coefficient') +
  scale_x_discrete(breaks=c("cq","height","nasalance"), labels=c("Contact quotient","Tongue height","Nasalance")) +
  scale_fill_manual(values = mycols) + theme_bw() + theme(legend.position="none",axis.text=element_text(size=14),
                                                            axis.title=element_text(size=15))

# Average acoustic/articulatory variables, by group
p2 <- ggplot(plotdat, aes(x=vars, y=means, fill=cluster2)) + 
  geom_hline(yintercept=0, lty=2) + geom_boxplot(notch=F) + labs(fill="Group") +
  geom_vline(xintercept=1.5) + geom_vline(xintercept=2.5) + geom_vline(xintercept=3.5) +
  xlab('Articulatory/acoustic variable') + ylab('Average nasal-oral z-score difference') +
  scale_x_discrete(breaks=c("cq","height","nasalance","f1"), labels=c("Contact quotient","Tongue height","Nasalance","F1")) +
  scale_fill_manual(values = mycols) + theme_bw() + theme(legend.position="none",axis.text=element_text(size=14),
                                                            axis.title=element_text(size=15))

# Combine both plots
myplot <- grid.arrange(p1, p2, ncol=1, nrow=2)