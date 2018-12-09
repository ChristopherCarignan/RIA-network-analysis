# Randomized oral-nasal matching
oral.nasal.matching <- function(mydata, myseed) {
  set.seed(myseed)
  diffdata <- c()
  
  for (vowel in c('a','e','o')){
    
    voweldata <- mydata[mydata$vowel==vowel,] # Get vowel pair data
    
    for (imitator in unique(voweldata$imitator)){
      subdata <- voweldata[voweldata$imitator==imitator,] # Get imitator data for vowel pair
      
      orals   <- unique(subdata$ident[subdata$nasality=='oral']) # Oral vowel tokens
      nasals  <- unique(subdata$ident[subdata$nasality=='nasal']) # Nasal vowel tokens
      
      # Keep iterating until all pairs are matched
      while (nrow(subdata) > 0){
        this.oral   <- orals[sample(length(orals),1)] # Select a random oral token
        onset       <- unique(subdata[subdata$ident==this.oral,'onset']) # Determine the onset consonant
        speaker     <- unique(subdata[subdata$ident==this.oral,'speaker']) # Determine the native FR who produced the token
        
        # Find all nasal vowel tokens that were produced by the same FR speaker and have the same onset consonant
        matches     <- unique(subdata$ident[subdata$nasality=='nasal' & subdata$onset==onset & subdata$speaker==speaker])
        this.nasal  <- matches[sample(length(matches),1)] # Select a random nasal token from the matches
        
        rows.oral   <- which(subdata$ident==this.oral) # Row numbers for the oral token
        rows.nasal  <- which(subdata$ident==this.nasal) # Row numbers for the nasal token
        
        dat.oral    <- subdata[rows.oral,] # Get the selected oral token data
        dat.nasal   <- subdata[rows.nasal,] # Get the selected nasal token data
        
        # Remove these two tokens from the various data frames, so that they won't be selected again
        subdata     <- subdata[-c(rows.oral,rows.nasal),] 
        orals       <- orals[orals!=this.oral]
        nasals      <- nasals[nasals!=this.nasal]
        
        # Get the variable measurements from the native FR speaker data
        diffs.fr <- c()
        for (var in c('f1z_s','nasalancez_s','heightz_s','cqz_s')){
          diffs.fr[[var]]  <- dat.nasal[[var]] - dat.oral[[var]] 
        }
        diffs.fr <- as.data.frame(diffs.fr)
        colnames(diffs.fr) <- c('f1','nasalance','height','cq')
        
        # Get the variable measurements from the naive AE imitator data
        diffs.ae <- c()
        for (var in c('f1z','nasalancez','heightz','cqz')){
          diffs.ae[[var]]  <- dat.nasal[[var]] - dat.oral[[var]] 
        }
        diffs.ae <- as.data.frame(diffs.ae)
        colnames(diffs.ae) <- c('f1','nasalance','height','cq')
        
        vars <- dat.oral[,c('imitator','speaker','order','onset','vowel','time_seq')] # Data frame of variables
        
        # Combine the AE imitator data and add to the main data frame
        binddat <- cbind('naive',imitator,vars,diffs.ae)
        colnames(binddat)[1] <- 'group'
        colnames(binddat)[2] <- 'participant'
        diffdata <- rbind(diffdata,binddat)
        
        # Combine the FR imitator data and add to the main data frame
        binddat <- cbind('native',speaker,vars,diffs.fr)
        colnames(binddat)[1] <- 'group'
        colnames(binddat)[2] <- 'participant'
        diffdata <- rbind(diffdata,binddat)
      }
    }
  }
  
  diffdata$time_seq <- as.factor(diffdata$time_seq) # Time sequence as factor (1:25%, 2:50%, 3:75%)
  diffdata$group    <- factor(diffdata$group, levels=c("native","naive")) # Language group as factor
  
  # Reorder factor level for participant
  diffdata$participant <- factor(diffdata$participant, 
                                 levels=c('FR01','FR02','FR03','FR04',
                                          'AE01','AE02','AE03','AE04','AE05','AE06','AE07','AE08','AE09'))
  
  return(diffdata)
}