# checkPatterns seeks for two patterns in data: (a) a row of same values and 
# (b) "zigzag" or "tree-like" pattern (i.e., the value of a next item differs
# by 1 from the previous item (the direction can change).
# The patterns are assessed in a frame the width of which can be set manually (i.e., 
# the number of adjacent items that are assessed at a time). A narrower frame is more
# sensitive but will tend catch even non-problematic cases. It is set to 5 by default.
# The frame "moves" one by one item until the end of the scale.
# The function returns a "suspect score" for each case with values ranging from 0 to 1
# (0 = no problem, 1 = a pattern of A or B or their combination is present across all
# variables).
# Alternatively, you can use showPatterns to display problematic cases (see below).

# Examples of use:
# data$suspect.score <- checkPatterns(data)
# data$suspect.score <- checkPatterns(data, frame.length=4)

checkPatterns <- function(data, frame.length=5) {
  if( !is.data.frame(data) & !is.matrix(data) ) {
    warning( "Data must be a data frame or matrix." )
    return()
  }
  if( frame.length > ncol(data) ) {
    frame.length <- ncol(data)
    warning( paste0("Frame.length larger than the number of variables; value set to ", frame.length ) )
  }
  if( frame.length < 2 ) {
    frame.length <- 2
    warning( paste0("Minimum frame.length to detect patterns is 2; value set to ", frame.length ) )
  }
  suspect <- c()
  for( observation in 1:nrow(data) ) {
    suspect.score <- 0
    for( frames in 1:(ncol(data)-frame.length+1) ) {
      values <- data[ observation, c(frames:(frames+frame.length-1)) ]
      #Pattern A: all values are the same
      if ( all( sapply( values, function(x) x == values[1] ) ) )
        suspect.score <- suspect.score + 1
      #Pattern B: a tree-like patern (values in a subsequent variable differ by 1)
      difference <- c()
      for( i in 1:(length(values)-1) ) {
        difference[i] <- abs( values[i] - values[i+1] )
      }
      if( all( sapply( difference, function(x) x == 1 ) ) )
        suspect.score <- suspect.score + 1
    }
    suspect[observation] <- round( suspect.score / (ncol(data)-frame.length+1), 3 )
  }
  return(suspect)
}

# --
# showPatterns is a wrapper for checkPatterns. In addition, you can set a cutoff
# for the identification of potentiallz problematic cases. It is set to 0.7 by default.
# The function returns a list of problematic cases with variable values.

# Examples of use:
# showPatterns(data)
# showPatterns(data, frame.length=4, cut=.6)

showPatterns <- function(data, frame.length=5, cut=.7) {
  data$suspect.score <- checkPatterns(data, frame.length=frame.length)
  return( data[ data$suspect >= cut, ] )
}

# --
# Based on auto-correlations, checkPatterns2 evaluates the similarity of data 
# to several prototypical patterns: (a) 1111111111, (b) 0101010101, and (c) 1234321234. 
# First, generates prototypical sequences of values based the data (the same number of 
# variables and range of values) and computes their auto-correlations (lags 1 to 3).
# For pattern A, auto-correlations cannot be computed (it results in NaN because the 
# data has no variance). For the purpose of this function, auto-correlation is considered
# to be r=1 in this case. Next, the distance of each case from each of the prototypes is
# computed (using the "manhattan" method) and the least distance is reported as a measure
# of "suspiciousness" of the case. The less the distance, the higher similarity to some
# of the prototypes, ie. the higher the suspicion. Therefore, the scores are interpretated
# in the opposite way compared to CheckPatterns.

# Exmples of use:
# data$distances <- checkPatterns2(data)
# data[order(data$distances), c("id","distances")]
# data[data$distances < .5, c("id","distances")]

checkPatterns2 <- function(data) {
  if( !is.data.frame(data) & !is.matrix(data) ) {
    warning( "Data must be a data frame or matrix." )
    return()
  }
  distances <- c()
  min <- min(data)
  max <- max(data)
  # Find prototype patterns auto-correlations
  acors.prototypes <- matrix( nrow = 3, ncol = 3, dimnames = list( c(1:3), c("lag1","lag2","lag3") ) )
  # Prototype 1: 0101010101
  acors.prototypes[1,] <- acf( rep_len( c(0,1), ncol(data) ), lag.max = 3, plot = F )[[1]][2:4]
  # Prototype 2: 1234321234
  acors.prototypes[2,] <- acf( rep_len( c( seq(min,max), seq(max-1,min+1) ), ncol(data) ), lag.max = 3, plot = F )[[1]][2:4]
  # Prototype 3: 1111111111
  acors.prototypes[3,] <- c(1,1,1)
  # Compute distances between the actual measurement and each prototype
  # Then, select the least distance as a suspect score
  distances <- apply( data, 1, function(row){
    acors.row <- acf( row, lag.max = 3, plot = F )[[1]][2:4]
    if( is.nan(acors.row[1]) )
      acors.row <- c(1,1,1)
    d <- c()
    for( i in 1:nrow(acors.prototypes) ) {
      d <- append( d, dist( rbind( acors.prototypes[i,], acors.row ), method = "manhattan" ) )
    }
    min(d)
  })
  distances <- round( distances, 3 )
  return(distances)
}
