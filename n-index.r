#######################################################################################
#                                    ABSOLUTE N-INDEX   (CLASSIC)
#######################################################################################

Nindex = function(pr, profile=F, NL=NULL, linear=F, x=NULL)
{
  pr   = pr[is.na(pr)==F];        #Extracting the NA values
  if(is.null(NL)) NL = length(pr)
  NL = min(c(NL,length(pr)))
  if(NL==1)stop("The length of data must be greater than 1\n")

  I0   = NA;                      #Initial value of I0
  n    = NA;                      #Initial value of n
  m    = NULL;                    #Variational exponent of n depending on time
  prof = NULL                     #Initial value of prof
  xk   = NULL                     #Sequence to estimate profile
  Pmax = rep(NA,NL)               #Initial values of Pmax

  if(profile)
  {
    xk   = round( seq(1,  length(pr),  length.out=NL) )
    prof = rep(NA,NL)
  }

  if(length(pr)>1)
  {
    prck = 0;
    for(i in 1:NL)
    {
      prck = prck + c(rep(NA,i-1),pr[1:(length(pr)-(i-1))])
      Pmax[i] = max(prck,na.rm=T)[1]

      if(profile & i %in% xk & linear)  prof[xk %in% i]  = 1-lm(log(Pmax[1:i])~log(1:i))$coefficients[2];
      if(profile & i %in% xk & !linear) prof[xk %in% i]  = nlp(Pmax[1:i])$n;
    }

    if(linear)
    {
     if(is.null(x)) x = 1:length(Pmax)
     lm0 = try(lm(log(Pmax)~log(x)),silent=T);  if(class(lm0)=="try-error") lm0 = list(coefficients=c(NA,NA))
     n   = 1-lm0$coefficients[2];    names(n)=NULL
     I0  = exp(lm0$coefficients[1]); names(I0)=NULL

     prop_sigma = summary(lm0)$coefficients[,2]
     upper = lm0$coefficients+ 1.96*prop_sigma
     lower = lm0$coefficients- 1.96*prop_sigma
     Dn = 1-c(upper[2],lower[2])
     DI =  exp(c(upper[1],lower[1]))
    }
    if(!linear)
    {
     lm0 = try(nlp(Pmax, x=x),silent=T);  if(class(lm0)=="try-error") lm0 = list(n=NA,I0=NA)
     n   = lm0$n;  names(n)=NULL
     I0  = lm0$I0; names(I0)=NULL
     Dn  = lm0$Dn
     DI  = lm0$DI
    }

    if(profile)
    {
      prof= round(prof,4)
      m = lm(log(prof[prof>0])~log(xk[prof>0]))$coefficients[2]
    }
  }

   return(list(n=n,I0=I0,Dn=Dn,DI=DI,prof=cbind(xk,prof),m=m))
}

#Modified n-index (using a nugget 'tm' at time=0)
MNindex = function(pr, profile=F, NL=NULL, linear=F, tm=1)
{
  pr   = pr[is.na(pr)==F];        #Extracting the NA values
  if(is.null(NL)) NL = length(pr)
  NL = min(c(NL,length(pr)))
  if(NL==1)stop("The length of data must be greater than 1\n")
  
  I0   = NA;                      #Initial value of I0
  n    = NA;                      #Initial value of n
  m    = NULL;                    #Variational exponent of n depending on time
  prof = NULL                     #Initial value of prof
  xk   = NULL                     #Sequence to estimate profile
  Pmax = rep(NA,NL)               #Initial values of Pmax
  
  if(profile)
  {
    xk   = round( seq(1,  length(pr),  length.out=NL) )
    prof = rep(NA,NL)
  }
  
  if(length(pr)>1)
  {
    prck = 0;
    for(i in 1:NL)
    {
      prck = prck + c(rep(NA,i-1),pr[1:(length(pr)-(i-1))])
      Pmax[i] = max(prck,na.rm=T)[1]
      
      if(profile & i %in% xk & linear)  prof[xk %in% i]  = 1-lm(log(Pmax[1:i])~log(1:i))$coefficients[2];
      if(profile & i %in% xk & !linear) prof[xk %in% i]  = nlp(Pmax[1:i])$n;
    }
    
    if(linear)
    {
      lm0 = try(lm(log(Pmax/(1:length(Pmax)))~log((1:length(Pmax)+tm)/(1+tm))),silent=T);  if(class(lm0)=="try-error") lm0 = list(coefficients=c(NA,NA))
      n   = -lm0$coefficients[2];     names(n)=NULL
      I0  = exp(lm0$coefficients[1]); names(I0)=NULL
      
      prop_sigma = lm0$coefficients
      upper = lm0$coefficients+ 1.96*prop_sigma
      lower = lm0$coefficients- 1.96*prop_sigma
      Dn = -c(upper[2],lower[2])
      DI =  exp(c(upper[1],lower[1]))
    }
    if(!linear)
    {
      lm0 = try(nlp0(Pmax,tm),silent=T);  if(class(lm0)=="try-error") lm0 = list(n=NA,I0=NA)
      n   = lm0$n;  names(n)=NULL
      I0  = lm0$I0; names(I0)=NULL
      Dn  = lm0$Dn
      DI  = lm0$DI
    }
    
    if(profile)
    {
      prof= round(prof,4)
      m = lm(log(prof[prof>0])~log(xk[prof>0]))$coefficients[2]
    }
  }
  
  return(list(n=n,I0=I0,Dn=Dn,DI=DI,prof=cbind(xk,prof),m=m))
}

#######################################################################################
#                                    NESTED N-INDEX   (INTERPRETED)
#######################################################################################

Nested = function(pr, profile=F, NL=NULL, linear=F)
{
  pr   = pr[is.na(pr)==F];        #Extracting the NA values
  if(is.null(NL)) NL = length(pr)
  NL = min(c(NL,length(pr)))
  if(NL==1)stop("The length of data must be greater than 1\n")

  I0   = NA;                      #Initial value of I0
  n    = NA;                      #Initial value of n
  m    = NULL;                    #Variational exponent of n depending on time
  prof = NULL                     #Initial value of prof
  xk   = NULL                     #Sequence to estimate profile
  Pmax = rep(NA,NL)               #Initial values of Pmax

  if(profile)
  {
    xk   = round( seq(1,  length(pr),  length.out=NL) )
    prof = rep(NA,NL)
  }

  if(length(pr)>1)
  {
    Imax = max(pr,na.rm=T)[1];     #Maximum value
    i0 = which(pr==Imax)[1]        #position of the maximum
    i1 = i0-1                      #The left  position of the maximum
    i2 = i0+1                      #The right position of the maximum

    while(i1 >= 1 | i2 <= NL)
    {
      if(i1<=0)         { Imax = c(Imax,pr[i2]); i2=i2+1}
      if(i2> length(pr)){ Imax = c(Imax,pr[i1]); i1=i1-1}
      if(i1>=1 & i2<=length(pr))
      {
         if(pr[i1]> pr[i2])
         {
           Imax = c(Imax,pr[i1]); i1=i1-1
         }else{
           Imax = c(Imax,pr[i2]); i2=i2+1
         }
      }
      if(profile  &  length(Imax) %in% xk)
      {
        Pmax = Imax*NA; for(i in 1:length(Imax)) Pmax[i] = sum(Imax[1:i],na.rm=T)
        if(linear)  prof[xk %in% length(Imax)]  = 1-lm(log(Pmax)~log(1:length(Pmax)))$coefficients[2];
        if(!linear) prof[xk %in% length(Imax)]  = nlp(Pmax)$n;
      }
    }
    Pmax = Imax*NA
    for(i in 1:length(Imax)) Pmax[i] = sum(Imax[1:i],na.rm=T)

    if(linear)
    {
      lm0 = try(lm(log(Pmax)~log(1:length(Pmax))),silent=T);  if(class(lm0)=="try-error") lm0 = list(coefficients=c(NA,NA))
      n   = 1-lm0$coefficients[2];    names(n)=NULL
      I0  = exp(lm0$coefficients[1]); names(I0)=NULL
    }
    if(!linear)
    {
     lm0 = try(nlp(Pmax),silent=T);  if(class(lm0)=="try-error") lm0 = list(n=NA,I0=NA)
     n   = lm0$n;  names(n)=NULL
     I0  = lm0$I0; names(I0)=NULL
    }

    if(profile)
    {
      prof= round(prof,4)
      m = lm(log(prof[prof>0])~log(xk[prof>0]))$coefficients[2]
    }
  }
  return(list(n=n,I0=I0,prof=cbind(xk,prof),m=m))
}


#######################################################################################
#                                    ORDERED N-INDEX   (CLIMATIC ANALYSIS)
#######################################################################################

Norder = function(pr, profile=F, NL=NULL, linear=F)
{
  pr   = pr[is.na(pr)==F];        #Extracting the NA values
  pr   = pr[pr>0];                #Extracting the 0 values
  if(is.null(NL)) NL = length(pr)
  NL = min(c(NL,length(pr)))
  if(NL==1)stop("The length of data must be greater than 1\n")

  I0   = NA;                      #Initial value of I0
  n    = NA;                      #Initial value of n
  m    = NULL;                    #Variational exponent of n depending on time
  prof = NULL                     #Initial value of prof
  xk   = NULL                     #Sequence to estimate profile
  Pmax = rep(NA,NL)               #Initial values of Pmax

  if(profile)
  {
    xk   = round( seq(1,  length(pr),  length.out=NL) )
    prof = rep(NA,NL)
  }


  Imax = sort(pr,decreasing=T)
  for(i in 1:NL)
  {
   Pmax[i] = sum(Imax[1:i])
   if(profile & i %in% xk & linear)  prof[xk %in% i]  = 1-lm(log(Pmax[1:i])~log(1:i))$coefficients[2];
   if(profile & i %in% xk & !linear) prof[xk %in% i]  = nlp(Pmax[1:i])$n;
  }
  if(linear)
  {
   lm0 = try(lm(log(Pmax)~log(1:length(Pmax))),silent=T);  if(class(lm0)=="try-error") lm0 = list(coefficients=c(NA,NA))
   n   = 1-lm0$coefficients[2];    names(n)=NULL
   I0  = exp(lm0$coefficients[1]); names(I0)=NULL
  }
  if(!linear)
  {
   lm0 = try(nlp(Pmax),silent=T);  if(class(lm0)=="try-error") lm0 = list(n=NA,I0=NA)
   n   = lm0$n;  names(n)=NULL
   I0  = lm0$I0; names(I0)=NULL
  }

  if(profile)
  {
    prof= round(prof,4)
    m = lm(log(prof[prof>0])~log(xk[prof>0]))$coefficients[2]
  }


  return(list(n=n,I0=I0,prof=cbind(xk,prof),m=m))
}

#######################################################################################
#                                    FRACTAL N-INDEX
#######################################################################################

Nfract = function(pr, profile=F, NL=NULL, linear=F)
{
  pr   = pr[is.na(pr)==F];        #Extracting the NA values
  if(is.null(NL)) NL = length(pr)
  NL = min(c(NL,length(pr)))
  if(NL==1)stop("The length of data must be greater than 1\n")

  I0   = NA;                      #Initial value of I0
  n    = NA;                      #Initial value of n
  m    = NULL;                    #Variational exponent of n depending on time
  prof = NULL                     #Initial value of prof
  xk   = NULL                     #Sequence to estimate profile

  NL2  = round(log(NL)/log(2)-0.499)
  SE2  = 2^seq(0,NL2,1)
  Pmax = rep(NA,length(SE2))      #Initial values of Pmax

  if(profile)
  {
    xk   = round( seq(1,  length(pr),  length.out=NL) )
    prof = rep(NA,NL)
  }

  if(length(pr)>1)
  {
    prck = 0;

    for(i in 1:max(SE2))
    {
      prck  = prck + c(rep(NA,i-1),pr[1:(length(pr)-(i-1))])
      prck2 = prck[seq(i,length(pr),i)]
      if(i %in% SE2)
      Pmax[SE2 %in% i] = max(prck2,na.rm=T)[1];

      if(profile & i %in% xk & linear)  prof[xk %in% i]  = 1-lm(log(Pmax[1:i])~log(SE2[1:i]))$coefficients[2];
      if(profile & i %in% xk & !linear) prof[xk %in% i]  = nlp(Pmax[1:i],SE2[1:i])$n;
    }

    if(linear)
    {
     lm0 = try(lm(log(Pmax)~log(SE2)),silent=T);  if(class(lm0)=="try-error") lm0 = list(coefficients=c(NA,NA))
     n   = 1-lm0$coefficients[2];    names(n)=NULL
     I0  = exp(lm0$coefficients[1]); names(I0)=NULL
    }
    if(!linear)
    {
     lm0 = try(nlp(Pmax,SE2),silent=T);  if(class(lm0)=="try-error") lm0 = list(n=NA,I0=NA)
     n   = lm0$n;  names(n)=NULL
     I0  = lm0$I0; names(I0)=NULL
    }

    if(profile)
    {
      prof= round(prof,4)
      m = lm(log(prof[prof>0])~log(xk[prof>0]))$coefficients[2]
    }
  }

   return(list(n=n,I0=I0,prof=cbind(xk,prof),m=m))
}


#######################################################################################
#                        SIMPLIFIED GINI INDEX    (see Monjo & Martin-Vide, 2015)
#######################################################################################

SGI. = function(prc)
{
 if(class(prc)=="list" | class(prc)=="data.frame")  prc = prc$dat
 pr = prc[is.na(prc)==F]
 pr = sort(pr[pr>0])
 nt = length(pr);
 PR  = pr*0; for(i in 1:nt)PR[i] = sum(pr[1:i])/sum(pr)

 GI  = 1-2/nt*sum(PR)
 return(GI)
}

Ecdf = function(x)
{
  x = x[!is.na(x)]
  z = sort(unique(x))
  m = cbind(z, ecdf(x)(z))
  return(m)
}

#######################################################################################
#                                    SPELL CALCULATION
#######################################################################################
#prc: precipitacion
#thres: threshold to consider dry/wet spell
#allow: number of dry pauses allowed in a wet spell defined.

spell = function(prc,thres=NULL,steps=3,allow=0,calc=T,linear=F, tm=0, part=1)
{
  #Control over data class
  if(class(prc)!="list" & class(prc)!="data.frame")
    prc = data.frame(dat=prc)
  if("dat" %in% names(prc) == F)
    prc = data.frame(dat=c(cbind(prc)[,dim(cbind(prc))[2]]))
  
  #Control over threshold
  if(is.null(thres))
    thres = min(prc$dat,na.rm=T)
  
  rl = rle(ifelse(prc$dat > thres,1,0))
  tdry=rl$lengths[rl$values==0 & !is.na(rl$values)]+allow  #Length of each dry spell. Correcting the length by the "allow number"
  twet=rl$lengths[rl$values==1 & !is.na(rl$values)]-allow  #Length of each wet spell. Correcting the length by the "allow number"
  iwet = which(twet>=steps);    #Position where we can find the spells of at least 3 consecutive days (if steps=3)
  
  #VERY IMPORTANT ##############################
  #This modefication is very important because "part" reduces the calculus with this factor, omiting the rest of spells
  iwets = iwet[seq(1,length(iwet)+1,part)-1]
  ###tdry  = tdry[seq(1,length(tdry),part)]
  
  twets = twet[iwets]           #Length of each wet spell of at least 3 consecutive days (if steps=3)
  sgi  = rep(NA,length(iwets))  #Simplified Gini Index
  Ptot = rep(NA,length(iwets))  #Total precipitation of each wet spell
  nord = rep(NA,length(iwets))  #Ordered Index n of each wet spell
  nfra = rep(NA,length(iwets))  #Fractal Index n of each wet spell
  nabs = rep(NA,length(iwets))  #Index n of each wet spell
  Dnab = rep(NA,length(iwets))  #Interval of the index n of each wet spell
  I0   = rep(NA,length(iwets))  #Fitted MAI of the minimum time step (time resolution of data)
  nest = rep(NA,length(iwets))  #Estimated n using the nested method
  Iref = rep(NA,length(iwets))  #Fitted MAI using the nested method
  inis = rep(NA,length(iwets))
  ends = rep(NA,length(iwets))
  
  
  if(calc)
    if(max(twet)>=steps)
      for(k in 1:length(iwets))
      {
        i = which(rl$values==1 & rl$length >= steps)[k]
        i1 = sum(rl$length[1:i]) - rl$length[i] + 1
        i2 = sum(rl$length[1:i])
        pr = prc$dat[i1:i2]
        
        inis[k] = i1
        ends[k] = i2
        Ptot[k] = (mean(pr,na.rm=T)*length(pr) + sum(pr,na.rm=T))/2
        
        n1 = suppressWarnings(try(MNindex(pr,linear=linear,tm=tm),silent=T));
        n2 = suppressWarnings(try(Norder(pr,linear=linear),silent=T));
        n3 = suppressWarnings(try(Nfract(pr,linear=linear),silent=T));
        g1 = try(SGI.(pr),silent=T);
        if(class(n2)!="try-error") nord[k] = n2$n
        if(class(n3)!="try-error") nfra[k] = n3$n
        if(class(g1)!="try-error") sgi[k] = g1
        if(class(n1)!="try-error")
          if(!is.na(n1$n))
          {
            nabs[k] = n1$n
            Dnab[k] = (n1$Dn[2]- n1$Dn[1])/(2*n1$n)
            I0[k]   = n1$I0
          }
        n4  = try(Nested(pr),silent=T);
        if(class(n4)!="try-error")
        {
          nest[k] = n4$n
          Iref[k] = n4$I0
        }
      }
  mn    = mean(prc$dat,na.rm=T) 
  sd    = sd(prc$dat,na.rm=T) 
  p0    = min(prc$dat,na.rm=T)
  pecdf = Ecdf(prc$dat)
  
  return(list(tdry=tdry, twet=twet, twets=twets, sgi=sgi, Ptot=Ptot, nabs=nabs, nest=nest, nord=nord, nfra=nfra, Dnab=Dnab, I0=I0, Iref=Iref, inis=inis, ends=ends, e=pecdf, p0=p0, mn=mn, sd=sd))
}
#Results
#tdry:  Sequence of consecutive dry days: Length of each relaxed 'dry spell'
#twet:  Sequence of consecutive wet days: Length of each relaxed 'wet spell'
#twets: Sequence of wet spells (with at least 3 consecutive days, if steps=3)
#Ptot: Total precipitation of each wet spell
#nabs: Absolute index n of each wet spell
#I0:   Fitted MAI of each wet spell
#nest: Nested index n of each wet spell
#Iref: Fitted MAI linked to nest.


#######################################################################################
#                                    AUXILIAR FUNCTIONS
#######################################################################################

nlp = function(Pmax, x=NULL, ref=T)    #if ref=T, reference intensity is fitted, else it is equal to 1
{
  if(is.null(x)) x = c(1:length(Pmax))
  x = x/length(x)
  y = Pmax/max(Pmax)

  ajust_f <- function (p)
  {
    if(ref)  z = p[2]*x^p[1]
    if(!ref) z =      x^p[1]
    squdiff  = (z-y)^2; squdiff = squdiff[y*z>0]
    sigmasqu = mean(squdiff,na.rm=T);
    minloglik=(length(y)/2)*log(sigmasqu);# negativo de la verosimilitud profile     -n/2 ln Suma(dif)^2 -> (-)(-)
  }

  lik = try(nlm(ajust_f,p=c(1,1),hessian=TRUE), silent=T);
  if(class(lik)!="try-error")
  {
    n  = 1-lik$estimate[1]
    if(ref)  I0 = lik$estimate[2]*max(Pmax)*(1/length(Pmax))^(1-n)
    if(!ref) I0 = max(Pmax)*(1/length(Pmax))^(1-n)
    prop_sigma<-sqrt(diag(solve(lik$hessian)))
  }
  if(class(lik)=="try-error")
  {
    lm0 = lm(log(Pmax)~log(1:length(Pmax)));
    n   = 1-lm0$coefficients[2];    names(n)=NULL
    I0  = exp(lm0$coefficients[1]); names(I0)=NULL

    lik = list(estimate=lm0$coefficients)
    prop_sigma = summary(lm0)$coefficients[,2]
  }

  upper = lik$estimate + 1.96*prop_sigma
  lower = lik$estimate - 1.96*prop_sigma
  Dn = 1-c(upper[1],lower[1])
  DI = c(upper[2],lower[2])*max(Pmax)*(1/length(Pmax))^(1-n)

  return(list(n=n,I0=I0,Dn=Dn,DI=DI))
}

nlp0 = function(Pmax, tm=0, x=NULL, ref=T)    #if ref=T, reference intensity is fitted, else it is equal to 1
{
  if(is.null(x)) x = c(1:length(Pmax))
  #x = x/length(x)
  y = Pmax/max(Pmax)/x
  
  ajust_f <- function (p)
  {
    if(ref)  z = p[2]*((x+tm)/(1+tm))^p[1]
    if(!ref) z =      ((x+tm)/(1+tm))^p[1]
    squdiff  = (z-y)^2; squdiff = squdiff[y*z>0]
    sigmasqu = mean(squdiff,na.rm=T);
    minloglik=(length(y)/2)*log(sigmasqu);# negativo de la verosimilitud profile     -n/2 ln Suma(dif)^2 -> (-)(-)
  }
  
  lik = try(nlm(ajust_f,p=c(1,1),hessian=TRUE), silent=T);
  if(class(lik)!="try-error")
  {
    n  = -lik$estimate[1]
    if(ref)  I0 = lik$estimate[2]*max(Pmax)
    if(!ref) I0 = max(Pmax)
    prop_sigma<-sqrt(diag(solve(lik$hessian)))
  }
  if(class(lik)=="try-error")
  {
    lm0 = lm(log(Pmax/(1:length(Pmax)))~log((1:length(Pmax)+tm)/(1+tm)));
    n   = -lm0$coefficients[2];    names(n)=NULL
    I0  = exp(lm0$coefficients[1]); names(I0)=NULL
    
    lik = list(estimate=lm0$coefficients)
    prop_sigma = summary(lm0)$coefficients[,2]
  }
  
  upper = lik$estimate + 1.96*prop_sigma
  lower = lik$estimate - 1.96*prop_sigma
  Dn = -c(upper[1],lower[1])
  DI = c(upper[2],lower[2])*max(Pmax)*(1/length(Pmax))^(-n)
  
  return(list(n=n,I0=I0,Dn=Dn,DI=DI))
}

Pmax_ = function(pr,NL=NULL)
{
  if(is.null(NL)) NL = length(pr)
  Pmax = rep(NA,NL)               #Initial values of Pmax

  if(length(pr)>1)
  {
    prck = 0;
    for(i in 1:NL)
    {
      prck = prck + c(rep(NA,i-1),pr[1:(length(pr)-(i-1))])
      Pmax[i] = max(prck,na.rm=T)[1]
    }
   }
  return(Pmax)
}
