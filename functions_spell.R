########################################
########################################
# Helper functions                     #
# 1) mean dry/wet days                 #
# 2) mean wet/dry spell                #
# 3) n-index applied to dry spells     #      
# Autors: Robert Monjo, Dominic Royé   #
# Email: robert@ficlima.org            #
########################################
########################################


## dry/wet days ##

#N: number of years
#thres: threshold for dry o wet days

mean_dw_days <- function(z, N = 38, thres = 0.1){
  
    z <- ifelse(z >= thres, 1, 0)
    
    wet <- sum(z==1)/N
    dry <- sum(z!=1)/N
    
  return(list(wet = wet, dry = dry))
  
}


## dry/wet spells ##
#thres: threshold for dry o wet days

mean_spells <- function(z, thres=0.1){
  
  z <- ifelse(z >= thres, 1, 0)
  
  temp <- data.frame(l = rle(z)$lengths, v = rle(z)$values)
    
 wet <- mean(temp[temp$v == 1, 1], na.rm = TRUE)
 dry <- mean(temp[temp$v != 1, 1], na.rm = TRUE)
    
    return(list(wet = wet, dry = dry))
  
}


## spell n-index ##

nindex_spell <- function(x){
  
    sp <- spell(x)
  
    #dry spell
    rdry <- sp$tdry
    
    rdry.sp <- spell(rdry, thres = 1)
    
    I0 <- mean(rdry.sp$I0)
    n <- mean(rdry.sp$nabs)
    
    return(list(I0 = I0, n = n))
}


