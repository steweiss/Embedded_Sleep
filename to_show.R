loadpackage<-function(x){
  for(i in x){
    if(!(i %in% rownames(installed.packages()))){
      if(i!="eegUtils"){
        eval(parse(text=paste0("install.packages(\"",i,"\")")))
        eval(parse(text=paste0("library(\"",i,"\")")))
      }else{
        eval(parse(text=paste0("remotes::install_github(","\'craddm/eegUtils\'",")")  )  )
        eval(parse(text=paste0("library(\"",i,"\")")))
      }
    }else{eval(parse(text=paste0("library(\"",i,"\")")))}
  }
}

loadpackage(c("xts","ggplotify","data.table","tidyr","ggplot2"))

mainDir<<-getwd()
reactiveValues<-function(...)list(...)

init.ExG<<-data.frame(fread(paste0(mainDir,"/Data/t6"), header=F))
init.srate=250
ExG<<-reactiveValues(
  ################## BASIC OBJECT 
  ts=as.xts(init.ExG,seq(as.POSIXct("2020-03-03 00:40:00"), by=1/init.srate, length=nrow(init.ExG)))
)
input=list()
input$pow=5
ExG$pow=input$pow

source(paste0(mainDir,"/krawtchouk.R"))

ExG$krw_fun=load_krawtchouk(ExG$pow)## incorrect should be chosen based on length of timeseries
ExG$krw_coef = krawtchouk(ExG$ts,ExG$pow)

for(i in 1:length(ExG$krw_fun)){
  df=data.frame("ind" = seq(1,length(ExG$krw_fun[[1]][,1]) ))
  df[,paste0("j=", seq(1,9) ) ] = ExG$krw_fun[[i]]
  df=gather(df,key="j",value = "Value",-ind)
  #df=df[df$Comp %in% paste0("Comp",unlist(ExG$gssa[[i]])),]
  #df$Comp=factor(df$Comp,levels = paste0("Comp",unlist(ExG$gssa[[i]])))
  df$Row=i
  #df$Column=input$stages_rows_selected
  ExG$krw_fun_res=rbind(ExG$krw_fun_res,df)
  # ### Fourier Coefficients
  df=data.frame("ind" = seq(1,length(ExG$krw_coef[[1]][,1]) ))
  df[,paste0("j=", seq(1,9) ) ] = ExG$krw_coef[[i]]
  df=gather(df,key="j",value = "Value",-ind)
  #df=df[df$Comp %in% paste0("Comp",unlist(ExG$gssa[[i]])),]
  #df$Comp=factor(df$Comp,levels = paste0("Comp",unlist(ExG$gssa[[i]])))
  df$Row=i
  ExG$krw_coef_res=rbind(ExG$krw_coef_res,df)
  
  # df[,paste0("Real",ExG$gssa[[i]])]=ExG$ssa$V[,ExG$gssa[[i]]]
  # df=gather(df,key="Real",value = "Value",-ind)
  # df=df[df$Real %in% paste0("Real",unlist(ExG$gssa[[i]])),]
  # df$Real=factor(df$Real,levels = paste0("Real",unlist(ExG$gssa[[i]])))
  # df$Row=i
  # df$Column=input$stages_rows_selected
  # ExG$df_res2=rbind(ExG$df_res2,df)
}
p=ggplot()+
  geom_line(
    data=ExG$krw_fun_res,
    aes(x=ind,y=Value,group=j,color=j))+
  theme_classic()+facet_wrap(c( vars(Row),vars(j)),nrow=4,ncol=9)+
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),legend.position = "none")+
  scale_y_continuous(name="",breaks=c(-0.15,0,0.15))+scale_x_continuous("Index")


p=ggplot()+
  geom_line(
    data=ExG$krw_coef_res,
    aes(x=ind,y=Value,group=j,color=j))+
  theme_classic()+facet_wrap(c( vars(Row),vars(j)),nrow=4,ncol=9)+
  theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),legend.position = "none")+
  scale_y_continuous(name="",breaks=c(-0.15,0,0.15))+scale_x_continuous("Index")
