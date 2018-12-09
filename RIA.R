# Relative importance analysis (RIA)
RIA <- function(diffdata, myseed) {
  set.seed(myseed)
  
  RIA.table   <- c()
  
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
  
  
  for (vowel in c('a','e','o')){
    for (participant in plist){
      subdata <- combdata[combdata$participant==participant & combdata$vowel==vowel,] # Get the unique data set
      
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
      RIA.table <- rbind(RIA.table,binddat) # Add data to RIA table
    }
  }
  RIA.table         <- as.data.frame(RIA.table)
  RIA.table$group   <- as.factor(RIA.table$group)
  RIA.table$means   <- as.numeric(as.character(RIA.table$means))
  RIA.table$coeffs  <- as.numeric(as.character(RIA.table$coeffs))
  RIA.table$R2      <- as.numeric(as.character(RIA.table$R2))
  RIA.table$ranks   <- as.numeric(as.character(RIA.table$ranks))
  
  # Reorder participant factor
  RIA.table$participant <- factor(RIA.table$participant, level=levels(diffdata$participant))
  
  return(RIA.table)
}