


#dry/wet days; N = 38

mean_dw_days <- function(z, N = 38){
  
    z <- ifelse(z >= 0.1, 1, 0)
    
    wet = sum(z==1)/N
    dry = sum(z!=1)/N
    
  return(list(wet = wet, dry=dry))
  
}

#####################

#dry/wet spells

mean_spells <- function(z){
  
  z <- ifelse(z >= 0.1, 1, 0)
  temp <- data.frame(l = rle(z)$lengths, v = rle(z)$values)
    
    wet <- mean(temp[temp$v == 1, 1], na.rm = TRUE)
    dry <- mean(temp[temp$v != 1, 1], na.rm = TRUE)
    
    return(list(wet = wet, dry=dry))
  
}


#####################
#spell n-index
nindex_spell <- function(x, I=FALSE){
  
  sp = spell(x)
  
    #racha seca
    rdry = sp$tdry
    
    #Calculamos las rachas de las rachas de sequia para las dos opciones
    rdry.sp = spell(rdry, thres=1)
    
    if(I==TRUE){
      #Por lo tanto el clima de sequias de p se resume con una "intensidad" y una "n":  I = I0*(1/i)^n
      return(mean(rdry.sp$I0))
      
    }else{
      #n-index
      return(mean(rdry.sp$nabs))
    }
  
}


