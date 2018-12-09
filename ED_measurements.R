### Euclidean distance measures between speakers
ED.measurements <- function(RIA.table, diffdata, myseed) {
  # Preallocate arrays
  RIA.table$cluster   <- c()
  RIA.table$cluster2  <- c()
  comb.table  <- c()
  dists       <- c()
  sim.mats    <- c()
  groupings   <- c()
  names       <- c()
  
  plist <- levels(RIA.table$participant) # List of all participant names
  
  for (vowel in c('a','e','o')){
    set.seed(myseed)
    subdata <- RIA.table[RIA.table$vowel==vowel,] # Get unique data set
    
    subdata$cluster <- c()
    subdata$cluster2 <- c()
    
    dists[[vowel]] <- matrix(nrow=length(plist), ncol=length(plist)) # Preallocate similarity matrix
    
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
        dists[[vowel]][x,y] <- 1/(1+dist)
      }
    }
    
    # Add participant names to similarity matrix
    colnames(dists[[vowel]])  <- plist
    rownames(dists[[vowel]])  <- plist
    
    # Scale the similarity matrix: 0-1
    sim.mats[[vowel]] <- (dists[[vowel]] - min(dists[[vowel]]))/(1 - min(dists[[vowel]]))
    
    # Create a blank network graph for the spinglass algorithm
    graph1 <- qgraph(sim.mats[[vowel]], layout="groups", graph="cor", sampleSize=length(plist), DoNotPlot=T)
    
    g   <- as.igraph(graph1, attributes=T) # Convert to igraph
    sgc <- spinglass.community(g) # Perform spinglass algorithm to get groupings/clusters
    
    # Preallocate arrays and variables for groupings/clusters
    groupings[[vowel]] <- list()
    names[[vowel]]     <- list()
    clusters  <- length(sgc$csize)
    letters   <- c('A','B','C','D','E','F','G','H','I','J','K','L','M')
    
    # Get groupings/clusters and assign to arrays
    for (clust in 1:clusters){
      groupings[[vowel]][clust]  <- sgc[clust]
      names[[vowel]][clust]      <- as.data.frame(rownames(sim.mats[[vowel]])[as.numeric(unlist(groupings[[vowel]][clust]))])
    }
    
    # Add groupings/clusters to subdata table
    for (clust in 1:clusters){
      subdata$cluster[subdata$participant %in% unlist(names[[vowel]][clust])]  <- clust # Cluster number
      subdata$cluster2[subdata$participant %in% unlist(names[[vowel]][clust])] <- letters[clust] # Cluster letter
    }
    subdata$cluster  <- as.factor(subdata$cluster)
    subdata$cluster2 <- as.factor(subdata$cluster2)
    
    # Add groupings/clusters to the RIA table
    RIA.table$cluster[RIA.table$vowel==vowel]   <- subdata$cluster
    RIA.table$cluster2[RIA.table$vowel==vowel]  <- as.character(subdata$cluster2)
    
    # Get average F1 values
    f1.table <- c()
    for (participant in unique(RIA.table$participant)){
      subdata  <- RIA.table[RIA.table$participant==participant,c('group','participant','vowel','cluster2')]
      subdata  <- subdata[subdata$vowel==vowel,]
      f1       <- diffdata[diffdata$participant==participant & diffdata$vowel==vowel,] %>% summarize(f1 = mean(f1))
      f1.table <- rbind(f1.table,
                        cbind(subdata[,1:3],'f1',f1$f1,0,0,0,0,subdata[,4])
      )
    }
    colnames(f1.table) <- colnames(RIA.table)
    f1.table$cluster   <- as.factor(f1.table$cluster)
    
    # Combine RIA table with F1 value table
    comb.table[[vowel]] <- rbind(RIA.table[RIA.table$vowel==vowel,],f1.table)
    
    # Assign clusters to F1 data
    for (participant in unique(RIA.table$participant)){
      # Cluster numbers
      comb.table[[vowel]]$cluster[
        comb.table[[vowel]]$participant==participant & comb.table[[vowel]]$vars=='f1'
        ] <-
        comb.table[[vowel]]$cluster[
          comb.table[[vowel]]$participant==participant & comb.table[[vowel]]$vars=='nasalance'
          ]
      # Cluster letters
      comb.table[[vowel]]$cluster2[
        comb.table[[vowel]]$participant==participant & comb.table[[vowel]]$vars=='f1'
        ] <-
        comb.table[[vowel]]$cluster2[
          comb.table[[vowel]]$participant==participant & comb.table[[vowel]]$vars=='nasalance'
          ]
    }
    
    comb.table[[vowel]]$cluster  <- as.factor(comb.table[[vowel]]$cluster)
    comb.table[[vowel]]$cluster2 <- as.factor(comb.table[[vowel]]$cluster2)
  }
  
  return(list(comb.table, dists, sim.mats, groupings))
}