# Perform within-participant, between-vowel cross validation on linear models
lm_cross_validation <- function(diffdata, myseed) {
  set.seed(myseed)
  
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
  
  # Preallocate cross validation data frames
  cross.val.R <- data.frame(participant=plist, 
                          'a.pred.e'=rep(NA,length(plist)),
                          'e.pred.a'=rep(NA,length(plist)),
                          'a.pred.o'=rep(NA,length(plist)),
                          'o.pred.a'=rep(NA,length(plist)),
                          'e.pred.o'=rep(NA,length(plist)),
                          'o.pred.e'=rep(NA,length(plist)))
  
  cross.val.P <- data.frame(participant=plist, group=rep(NA,length(plist)),
                            'a.pred.e'=rep(NA,length(plist)),
                            'e.pred.a'=rep(NA,length(plist)),
                            'a.pred.o'=rep(NA,length(plist)),
                            'o.pred.a'=rep(NA,length(plist)),
                            'e.pred.o'=rep(NA,length(plist)),
                            'o.pred.e'=rep(NA,length(plist)))
  
  myrow <- 1
  
  for (participant in plist){
    # Get the unique data sets for each of the three vowel pairs
    a.data <- combdata[combdata$participant==participant & combdata$vowel=='a',] 
    e.data <- combdata[combdata$participant==participant & combdata$vowel=='e',] 
    o.data <- combdata[combdata$participant==participant & combdata$vowel=='o',] 
    
    group <- as.character(unique(a.data$group)) # Which language group?
    
    # Add participant and language group data to cross validation tables
    cross.val.R$participant[myrow]  <- participant
    cross.val.P$participant[myrow]  <- participant
    cross.val.R$group[myrow]        <- group
    cross.val.P$group[myrow]        <- group
    
    vars  <- c('nasalance','height','cq')
    
    # Create linear models for each vowel subset to predict F1 from the articulatory variables
    lm.a <- lm(f1 ~ nasalance + height + cq, data=a.data) 
    lm.e <- lm(f1 ~ nasalance + height + cq, data=e.data) 
    lm.o <- lm(f1 ~ nasalance + height + cq, data=o.data) 
    
    
    # Cross validation between a-pair and e-pair
    # a-pair model predicts e-pair F1 values, using e-pair articulatory variables
    a.pred.e <- predict(lm.a, newdata=e.data[,vars])
    cv <- cor.test(a.pred.e, e.data$f1)
    cross.val.R$a.pred.e[myrow] <- cv$estimate
    cross.val.P$a.pred.e[myrow] <- cv$p.value
    # e-pair model predicts a-pair F1 values, using a-pair articulatory variables
    e.pred.a <- predict(lm.e, newdata=a.data[,vars])
    cv <- cor.test(e.pred.a, a.data$f1)
    cross.val.R$e.pred.a[myrow] <- cv$estimate
    cross.val.P$e.pred.a[myrow] <- cv$p.value
    
    
    # Cross validation between a-pair and o-pair
    # a-pair model predicts o-pair F1 values, using o-pair articulatory variables
    a.pred.o <- predict(lm.a, newdata=o.data[,vars])
    cv <- cor.test(a.pred.o, o.data$f1)
    cross.val.R$a.pred.o[myrow] <- cv$estimate
    cross.val.P$a.pred.o[myrow] <- cv$p.value
    # o-pair model predicts a-pair F1 values, using a-pair articulatory variables
    o.pred.a <- predict(lm.o, newdata=a.data[,vars])
    cv <- cor.test(o.pred.a, a.data$f1)
    cross.val.R$o.pred.a[myrow] <- cv$estimate
    cross.val.P$o.pred.a[myrow] <- cv$p.value
    
    
    # Cross validation between e-pair and o-pair
    # e-pair model predicts o-pair F1 values, using o-pair articulatory variables
    e.pred.o <- predict(lm.e, newdata=o.data[,vars])
    cv <- cor.test(e.pred.o, o.data$f1)
    cross.val.R$e.pred.o[myrow] <- cv$estimate
    cross.val.P$e.pred.o[myrow] <- cv$p.value
    # o-pair model predicts e-pair F1 values, using e-pair articulatory variables
    o.pred.e <- predict(lm.o, newdata=e.data[,vars])
    cv <- cor.test(o.pred.e, e.data$f1)
    cross.val.R$o.pred.e[myrow] <- cv$estimate
    cross.val.P$o.pred.e[myrow] <- cv$p.value
    
    myrow <- myrow + 1
  }
  
  cross.val.R             <- as.data.frame(cross.val.R)
  cross.val.R$participant <- as.factor(cross.val.R$participant)
  cross.val.R$group       <- as.factor(cross.val.R$group)
  
  cross.val.P             <- as.data.frame(cross.val.P)
  cross.val.P$participant <- as.factor(cross.val.P$participant)
  cross.val.P$group       <- as.factor(cross.val.P$group)

  # Reorder participant factor
  cross.val.R$participant <- factor(cross.val.R$participant, level=levels(diffdata$participant))
  cross.val.P$participant <- factor(cross.val.P$participant, level=levels(diffdata$participant))
  
  return(list(cross.val.R,cross.val.P))
}