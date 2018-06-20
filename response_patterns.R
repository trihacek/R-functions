# Funkce checkPatterns v datech vyhledava dva typy vzorcu: (a) sloupec stejnych hodnot a 
# (b) "cik-cak" nebo "stromecky", tj. pripady, kdy se hodnota v nasledujici polozce lisi 
# o 1 od predchazejici polozky (smer se pritom muze menit).
# Vzorce se posuzuji v "okenku", jehoz sirku lze nastavit (tj. pocet sousedicich 
# promennych, ktere se berou v danou chvili v potaz). Mensi okenko je citlivejsi, ale 
# samozrejme bude casteji zachycovat i neproblematicke pripady. Defaultne je nastavena 
# sirka 5. Okenko se "posouva" po jedne polozce az do konce dotazniku.
# Vystupem funkce je "suspect score" s hodnotami od 0 do 1, ktery lze priradit 
# do promenne a pak filtrovat (0 = zadny problem, 1 = vzorec A nebo B nebo jejich
# kombinace jsou pritomny napric vsemi promennymi).
# Alternativne si lze nechat vypsat problematicka mereni (viz nize funkce showPatterns).
# (Funkce zatim neni zcela zabezpecena proti zadavani nesmyslnych hodnot.)

# Priklad pouziti:
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

# Fuknce showPatterns je wrapper pro funkci checkPatterns. Lze v ni navic zadat cutoff
# pro identifikaci potencialne problematickych mereni. Defaultne nastaveno na 0.7.
# Vystupem je vypis problematickych mereni i s hodnotami promennych.

# Priklady pouziti:
# showPatterns(data)
# showPatterns(data, frame.length=4, cut=.6)

showPatterns <- function(data, method="exact", frame.length=5, cut=.7) {
  data$suspect.score <- checkPatterns(data, method, frame.length=frame.length)
  return( data[ data$suspect >= cut, ] )
}


# Funkce checkPatterns2 pomocina zaklade autokorelaci vyhodnocuje podobnost dat 
# nekolika prototypickym vzorcum: (a) 1111111111, (b) 0101010101, (c) 1234321234, 
# (d) nahodne odpovedi. Nejprve na zaklade zadanych dat vygeneruje prototypy 
# (stejny pocet promennych, stejne rozpeti hodnot) a vypocita jejich autokorelace 
# (lag 1 az 3). Pak vypocita vzdalenost kazdeho respondenta od kazdeho z techto
# prototypu (metodou manhattan) a tu nejmensi z nich reportuje jako meritko 
# "podezrelosti" odpovedi. Cim mensi vzdalenost, tim vetsi shoda s nekterym 
# prototypem, tj. tim vetsi podezreni. Interpretace skoru je tedy opacna nez u 
# funkce checkPatterns.

# Priklad pouziti:
# data$distances <- checkPatterns2(data)
# data[order(data$distances), c("id","distances")]
# data[data$distances < .5, c("id","distances")]

checkPatterns2 <- function(data) {
  if( !is.data.frame(data) & !is.matrix(data) ) {
    warning( "Data must be a data frame or matrix." )
    return()
  }
  distances <- c()
  min <- min( data )
  max <- max( data )
  # Set prototype patterns
  prototypes <- data.frame( rbind(
    # 0101010101
    rep_len( c(0,1), ncol(data)),
    # 1234321234
    rep_len( c( seq(min,max), seq(max-1,min+1) ), ncol(data) ),
    # Random responses
    sample( c(min:max), ncol(data), replace=T)
  ), row.names = c(1:3) )
  # Prototype auto-correlations (row = prototype, col = lag)
  acors.prototypes <- cbind( apply( prototypes, 1, function(x){
    acf( x, lag.max = 3, plot = F )[[1]][2:4]
  }) )
  # Adding pattern 1111111111, for which autocorrelation cannot be computed
  acors.prototypes <- rbind( acors.prototypes, c(1,1,1))
  # Compute distances between the actual measurement and each prototype
  # Then, select the least distance as a suspect score
  distances <- apply( data, 1, function(row){
    acors.row <- acf( row, max.lag = 3, plot = F )[[1]][2:4]
    if( is.nan(acors.row[1]) )
      acors.row <- c(1,1,1)
    d <- c()
    for( i in 1:nrow(acors.prototypes) ) {
      d <- append( d, dist( rbind( acors.prototypes[i,], acors.row ), method = "manhattan" ) )
    }
    min(d)
  })
  distances <- round( distances, 3)
  return(distances)
}