


#dry/wet days; N = 37

mean_spells <- function(z, N){
  
  z <- ifelse(z>=0.1,1,0)
    
    wet = sum(z==1)/N
    dry = sum(z!=1)/N
    
  return(list(wet = wet, dry=dry))
  
}

#####################

#dry/wet spells

mean_spells2 <- function(z){
  
  z <- ifelse(z>=0.1,1,0)
  temp <- data.frame(l=rle(z)$lengths,v=rle(z)$values)
    
    wet = mean(temp[temp$v==1,1],na.rm=TRUE)
    dry = mean(temp[temp$v!=1,1],na.rm=TRUE)
    
    return(list(wet = wet, dry=dry))
  
}


#spell n-index
fun_aux <- function(x,opcion=1,I=TRUE){
  
  sp = spell(x)
  
  if(opcion==1){
    
    
    #Opcion 1: racha seca
    rdry = sp$tdry
    
    #Calculamos las rachas de las rachas de sequia para las dos opciones
    rdry.sp = spell(rdry, thres=1)
    
    if(I==TRUE){
      #Por lo tanto el clima de sequias de p se resume con una "intensidad I0" y una "n":  I = I0*(1/i)^n
      return(mean(rdry.sp$I0))
    }else{
      return(mean(rdry.sp$nabs))
    }
  }else{
    
    #Opcion 2: balance seco - humedo
    n = min(c(length(sp$tdry),length(sp$twet)))
    bdry = sp$tdry[1:n]-sp$twet[1:n]
    #Calculamos las rachas de las rachas de sequia para las dos opciones
    bdry.sp = spell(bdry, thres=0)
    #Por lo tanto el clima de sequias de p se resume con una "intensidad I0" y una "n":  I = I0*(1/i)^n
    if(I==TRUE){
      return(mean(bdry.sp$I0))
    }else{
      return(mean(bdry.sp$nabs))
    }
  }
}



fun_auxh <- function(x,I=TRUE){
  
  sp = spell(x)
  
  #Opcion 1: racha seca
  rhum = sp$twet
  
  #Calculamos las rachas de las rachas de pr para las dos opciones
  rhum.sp = spell(rhum, thres=1)
  
  if(I==TRUE){
    #Por lo tanto el clima de sequias de p se resume con una "intensidad I0" y una "n":  I = I0*(1/i)^n
    return(mean(rhum.sp$I0))
  }else{
    return(mean(rhum.sp$nabs))
  }
}

