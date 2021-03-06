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


# Relative importance analysis on combined vowel data
RIA.table.allV <- c()

# Create lists of names for all participants, naive AE imitators, and native FR speakers
plist       <- levels(diffdata$participant) 
AElist      <- unique(sort(diffdata$participant[diffdata$group=='naive']))
FRlist      <- unique(sort(diffdata$participant[diffdata$group=='native']))

naive       <- diffdata[diffdata$group=='naive',] # Get naive AE imitator data 
naive$item  <- as.numeric(as.factor(paste0(naive$imitator,naive$speaker,naive$order))) # Unique identities for each token

native      <- diffdata[diffdata$group=='native',] # Get the native FR speaker data
native$item <- as.numeric(as.factor(paste0(native$imitator,native$speaker,native$order))) # Unique identities for each token


# Since there is an imbalance in the number of native FR speakers (4) and the number of naive AE imitators (9),
# there is a greater number of items for each FR speaker than for each AE imitator. In order to maintain a balanced
# item set for the Euclidean distances measurements, the item sets for the FR speakers need to be reduced to equal
# the size of the item sets for the AE imitators. This is performed in a randomized manner.

# Get the number of individual oral-nasal matches for each imitator
mylength    <- length(unique(naive$order[naive$participant==AElist[1]])) 

nativenew   <- c()
for (onset in c('b','p')){
  for (vowel in c('a','e','o')){
    for (participant in FRlist){
      subdata   <- native[native$participant==participant & native$vowel==vowel & native$onset==onset,] # Get the unique data set
      keep      <- sample(unique(subdata$item), mylength/6) # Get a randomized subset of the data (6: 2 onsets X 3 vowels)
      
      trimmed   <- native[native$item %in% keep,] # Select the randomized data subset
      nativenew <- rbind(nativenew,trimmed) # Add data to the new native data set
    }
  }
}
combdata <- rbind(naive,nativenew) # Combine new FR native speaker data set to naive AE imitator data set


for (participant in plist){
  subdata <- combdata[combdata$participant==participant,] # Get the unique data set
  
  group <- as.character(unique(subdata$group)) # Which language group?
  vars  <- c('nasalance','height','cq')
  
  lm1 <- lm(f1 ~ nasalance + height + cq, data=subdata) # Create a linear model to predict F1 from the articulatory variables
  rel.weights <- calc.relimp(lm1, type = c("lmg"), rela = F) # Use RIA to find the accurate coefficients for each variable
  
  coeffs <- as.vector(rel.weights$ave.coeffs[,1]) # Extract the RIA coefficients
  
  R2 <- rel.weights$R2 # Extract the average R^2 for each of the articulatory variables
  
  ranks <- as.vector(rel.weights$lmg.rank) # Extract the rankings for the articulatory variables
  
  # Get the average values for the three articulatory variables
  means <- as.vector(c(
    mean(subdata$nasalance),
    mean(subdata$height),
    mean(subdata$cq)
  ))
  
  binddat   <- cbind(group,participant,vowel,vars,means,coeffs,R2,ranks) # Combine data
  RIA.table.allV <- rbind(RIA.table.allV,binddat) # Add data to RIA table
}
RIA.table.allV         <- as.data.frame(RIA.table.allV)
RIA.table.allV$group   <- as.factor(RIA.table.allV$group)
RIA.table.allV$means   <- as.numeric(as.character(RIA.table.allV$means))
RIA.table.allV$coeffs  <- as.numeric(as.character(RIA.table.allV$coeffs))
RIA.table.allV$R2      <- as.numeric(as.character(RIA.table.allV$R2))
RIA.table.allV$ranks   <- as.numeric(as.character(RIA.table.allV$ranks))

# Reorder participant factor
RIA.table.allV$participant <- factor(RIA.table.allV$participant, level=levels(diffdata$participant))


# Summarize articulatory variable rankings
RIA.table.allV %>% group_by(vars, ranks) %>% summarise(n = n()) %>% mutate(freq = 100*n / sum(n))


# Determine the number of clusters generated by the network model from each of 100 different random seeds
cluster.table.allV <- c()
myrow <- 1

plist <- levels(RIA.table.allV$participant) # List of all participant names

for (myseed in 1:100){
  set.seed(myseed)
  
  dists <- matrix(nrow=length(plist), ncol=length(plist)) # Preallocate similarity matrix
  
  for (x in 1:length(plist)){
    # Get participant 1
    p1    <- plist[x]
    p1dat <- RIA.table.allV$coeffs[RIA.table.allV$participant==p1]
    
    for (y in 1:length(plist)){
      # Get participant 2
      p2    <- plist[y]
      p2dat <- RIA.table.allV$coeffs[RIA.table.allV$participant==p2]
      
      # Euclidean distance between participants 1 and 2
      dist  <- sqrt( (p1dat[1] - p2dat[1])^2 + (p1dat[2] - p2dat[2])^2 + (p1dat[3] - p2dat[3])^2 )
      
      # Distance-based similarity score
      dists[x,y] <- 1/(1+dist)
    }
  }
  
  # Scale the similarity matrix: 0-1
  sim.matx <- (dists - min(dists)) / (1 - min(dists))
  
  # Create a blank network graph for the spinglass algorithm
  graph1 <- qgraph(sim.matx, layout="groups", graph="cor", sampleSize=length(plist), DoNotPlot=T)
  
  g   <- as.igraph(graph1, attributes=T) # Convert to igraph
  sgc <- spinglass.community(g) # Perform spinglass algorithm to get groupings/clusters
  
  # Get groupings/clusters
  clusters  <- length(sgc$csize)
  cluster.table.allV$seed[myrow]      <- myseed
  cluster.table.allV$clusters[myrow]  <- clusters
  
  myrow <- myrow + 1 # Iterate for next seed
}

cluster.table.allV <- as.data.frame(cluster.table.allV)
cluster.table.allV$clusters <- as.factor(cluster.table.allV$clusters)


# Print summaries to determine the cluster frequencies
summary(cluster.table.allV$clusters)


# Pick a new seed for clustering! (doesn't matter in this case, because 2 groups were found for all 100 random seeds)
set.seed(2)

### Euclidean distance measures between speakers
plist <- levels(RIA.table.allV$participant) # List of all participant names

subdata <- RIA.table.allV

subdata$cluster <- c()
subdata$cluster2 <- c()

dists <- matrix(nrow=length(plist), ncol=length(plist)) # Preallocate similarity matrix

for (x in 1:length(plist)){
  # Get participant 1
  p1    <- plist[x]
  p1dat <- subdata$coeffs[subdata$participant==p1]
  
  for (y in 1:length(plist)){
    # Get participant 2
    p2    <- plist[y]
    p2dat <- subdata$coeffs[subdata$participant==p2]
    
    # Euclidean distance between participants 1 and 2
    dist  <- sqrt( (p1dat[1] - p2dat[1])^2 + (p1dat[2] - p2dat[2])^2 + (p1dat[3] - p2dat[3])^2 )
    
    # Distance-based similarity score
    dists[x,y] <- 1/(1+dist)
  }
}

# Add participant names to similarity matrix
colnames(dists)  <- plist
rownames(dists)  <- plist

# Scale the similarity matrix: 0-1
sim.matx <- (dists - min(dists))/(1 - min(dists))

# Create a blank network graph for the spinglass algorithm
graph1 <- qgraph(sim.matx, layout="groups", graph="cor", sampleSize=length(plist), DoNotPlot=T)

g   <- as.igraph(graph1, attributes=T) # Convert to igraph
sgc <- spinglass.community(g) # Perform spinglass algorithm to get groupings/clusters

# Preallocate arrays and variables for groupings/clusters
groupings <- list()
names     <- list()
clusters  <- length(sgc$csize)
letters   <- c('A','B','C','D','E','F','G','H','I','J','K','L','M')

# Get groupings/clusters and assign to arrays
for (clust in 1:clusters){
  groupings[clust]  <- sgc[clust]
  names[clust]      <- as.data.frame(rownames(sim.matx)[as.numeric(unlist(groupings[clust]))])
}

# Add groupings/clusters to subdata table
for (clust in 1:clusters){
  subdata$cluster[subdata$participant %in% unlist(names[clust])]  <- clust # Cluster number
  subdata$cluster2[subdata$participant %in% unlist(names[clust])] <- letters[clust] # Cluster letter
}
subdata$cluster  <- as.factor(subdata$cluster)
subdata$cluster2 <- as.factor(subdata$cluster2)

# Add groupings/clusters to the RIA table
RIA.table.allV$cluster   <- subdata$cluster
RIA.table.allV$cluster2  <- as.character(subdata$cluster2)

# Get average F1 values
f1.table <- c()
for (participant in unique(RIA.table.allV$participant)){
  subdata  <- RIA.table.allV[RIA.table.allV$participant==participant,c('group','participant','vowel','cluster2')]
  f1       <- diffdata[diffdata$participant==participant,] %>% summarize(f1 = mean(f1))
  f1.table <- rbind(f1.table,
                    cbind(subdata[,1:3],'f1',f1$f1,0,0,0,0,subdata[,4])
  )
}
colnames(f1.table) <- colnames(RIA.table.allV)
f1.table$cluster   <- as.factor(f1.table$cluster)

# Combine RIA table with F1 value table
comb.table <- rbind(RIA.table.allV,f1.table)

# Assign clusters to F1 data
for (participant in unique(RIA.table.allV$participant)){
  # Cluster numbers
  comb.table$cluster[comb.table$participant==participant & comb.table$vars=='f1'] <-
    comb.table$cluster[comb.table$participant==participant & comb.table$vars=='nasalance']
  # Cluster letters
  comb.table$cluster2[comb.table$participant==participant & comb.table$vars=='f1'] <-
    comb.table$cluster2[comb.table$participant==participant & comb.table$vars=='nasalance']
  
  comb.table$cluster  <- as.factor(comb.table$cluster)
  comb.table$cluster2 <- as.factor(comb.table$cluster2)
}


# Plot results
mycols  <- brewer.pal(12, "Set3") # Create qualitative color palette for plotting
npars   <- length(unique(diffdata$participant)) # The number of participants is needed for various network plotting variables

# Plot the network graph
qgraph(sim.matx, layout="groups", graph="cor", sampleSize=npars, groups=groupings,
       posCol="orange", color=mycols, legend=T, GLratio=10, legend.cex=1,
       layout.par = list(init=matrix(rnorm(npars*2), npars, 2)),
       vsize=7, cut=0, border.width=1.5)

# Plot the average RIA coefficients and average acoustic/articulatory variables by group/cluster
plotdat <- comb.table

# Average RIA coefficients, by group
p1 <- ggplot(plotdat[plotdat$vars!='f1',], aes(x=vars, y=coeffs, fill=cluster2)) + 
  geom_hline(yintercept=0, lty=2) + geom_boxplot(notch=F) + labs(fill="Group") +
  geom_vline(xintercept=1.5) + geom_vline(xintercept=2.5) +
  xlab('Articulatory variable') + ylab('Relative importance coefficient') +
  scale_x_discrete(breaks=c("cq","height","nasalance"), labels=c("Contact quotient","Tongue height","Nasalance")) +
  scale_fill_manual(values = mycols) +  theme_bw()  + theme(legend.position="none",axis.text=element_text(size=14),
                                                            axis.title=element_text(size=15))

# Average acoustic/articulatory variables, by group
p2 <- ggplot(plotdat, aes(x=vars, y=means, fill=cluster2)) + 
  geom_hline(yintercept=0, lty=2) + geom_boxplot(notch=F) + labs(fill="Group") +
  geom_vline(xintercept=1.5) + geom_vline(xintercept=2.5) + geom_vline(xintercept=3.5) +
  xlab('Articulatory/acoustic variable') + ylab('Average nasal-oral z-score difference') +
  scale_x_discrete(breaks=c("cq","height","nasalance","f1"), labels=c("Contact quotient","Tongue height","Nasalance","F1")) +
  scale_fill_manual(values = mycols) +  theme_bw()  + theme(legend.position="none",axis.text=element_text(size=14),
                                                            axis.title=element_text(size=15))

# Combine both plots
myplot <- grid.arrange(p1, p2, ncol=1, nrow=2)