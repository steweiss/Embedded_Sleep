load_krawtchouk<-function(pow){
  krw_fun=data.frame(read.csv(paste0(getwd(),"/Data/Krawtchouk/K",2^pow,"D.txt"),header=F))
  cn=krw_fun[1:(nrow(krw_fun)/(2^pow+1))*(2^pow+1)-2^pow,1]
  krw_fun= data.frame(as.numeric(krw_fun[-(1:(nrow(krw_fun)/(2^pow+1))*(2^pow+1)-2^pow),1]))
  Krawfun<-list()
  n=0
  j=1
  while(j<5){
    l=1
    Krawfun[j][[1]]<-data.frame(matrix(ncol=9,nrow=2^pow))
    while(l<10){
      n=n+1
      Krawfun[j][[1]][,l]<-krw_fun[(((n-1)*(2^pow))+1):(n*(2)^pow),1]
      l=l+1
    }
    j=j+1
  }
  return(Krawfun)
}

krawtchouk<-function(ts,pow){
  ## load fourier functions according to power - this should be according to length of timeseries, but has not been implemented yet
  ## however power can be picked via knob in graphics implementation
  kraw_fou=load_kraw_fou(pow)
  krw_coef=list()
  krw_enrgy=data.frame(matrix(nrow=9,ncol=4))
  
  ## main calculations
  
  ## SLIDING WINDOW
  t=1
  check=1
  
  while(t<length(ts)-2^pow){
    j=0
    while(j<length(kraw_fou$X)){
      j=j+1
      l=0
      krw_coef[j][[1]]<-data.frame(matrix(ncol=9,nrow=2^pow))
      while(l<length(kraw_fou$X[[1]])){
        l=l+1
        y_fou=fft(
          as.numeric(
            window(ts,
                 start=index(ts)[t],
                 end=index(ts)[(t-1+2^pow)]
            )
          )
        )
        y_re=Re(y_fou)
        y_im=Im(y_fou)
        df_res_x=(y_re*kraw_fou$X[[j]][,l]+y_im*kraw_fou$Y[[j]][,l])*(2^pow)#Why times t1?#
        df_res_y=(y_re*kraw_fou$Y[[j]][,l]-y_im*kraw_fou$X[[j]][,l])*(2^pow)
        krw_coef[j][[1]][,l]=Re(fft(complex(real=df_res_x,imag=df_res_y),inverse=T))
        
        ## Power Spectrum Density
        krw_enrgy[l,j]=sum(krw_coef[j][[1]][,l]^2)
        
      }
    }
    mf=max(apply(krw_enrgy,1,max))
    lev=c(60000,50000,40000,20000,10000,8000,7000,5000,3000)
    if(all(mf<lev)){
        check=0
        print(t)
    }
    
    
    t=t+2^pow 
  }
  ## Treshold

  if(check!=1){
    return(0)
  }else{
    
    return(krw_coef)
  }
  
}

load_kraw_fou<-function(pow){
  krw_fun=data.frame(read.csv(paste0(getwd(),"/Data/Krawtchouk/K",2^pow,"F.txt"),header=F))
  cn=krw_fun[1:(nrow(krw_fun)/(2^pow+1))*(2^pow+1)-2^pow,1]
  krw_fun= data.frame(as.numeric(krw_fun[-(1:(nrow(krw_fun)/(2^pow+1))*(2^pow+1)-2^pow),1]))
  
  Kraw_x<-list()
  Kraw_y<-list()
  n=0
  j=1
  while(j<5){
    l=1
    Kraw_x[j][[1]]<-data.frame(matrix(ncol=9,nrow=2^pow))
    Kraw_y[j][[1]]<-data.frame(matrix(ncol=9,nrow=2^pow))
    while(l<10){
      n=n+1
      Kraw_x[j][[1]][,l]<-krw_fun[(((n-1)*(2^pow))+1):(n*(2)^pow),1]
      n=n+1
      Kraw_y[j][[1]][,l]<-krw_fun[(((n-1)*(2^pow))+1):(n*(2)^pow),1]
      l=l+1
    }
    j=j+1
  }
  return(list(X=Kraw_x,Y=Kraw_y))
}
