
#Island biogeography model producing stochastic time series of species number
#Used in section 2 on duration in Gonzalez et al. submitted to Ecology
#This code estimates effect sizes for the trend in diversity. Adjusting the extinction and colonization probabilities
#allows us to create strong and weak trends in diversity. The aim is to see if we can get precise estimates of the trend
#using short time series. Then create a fake meta-analytic dataset to see how power it has to detecting a mean effect size different
#from zero. Note, we are assuming perfect sampling of the community, so no error in estimate of true local diversity.

rm(list=ls())
require(ggplot2)
require(nlme)
require(lattice)

IBG_fun<-function(SR.init=22,SR.pool=100,MAXGEN=50,Area=100,MAX.ext=0.5,MAX.col=0.3,area.exp=0.3,distance=100,distance.exp=0.1){
  
  e.prob<-rnorm(MAXGEN, mean=MAX.ext, sd=0.07)
  c.prob<-rnorm(MAXGEN, mean=MAX.col, sd=0.07)
  
  Time<-seq(from=1, to=MAXGEN)
  #vector for richness values
  SRtime<-matrix(0,MAXGEN,1)
  SRtime[1]<-SR.init

  
  #iterate equation over time 
  for(i in 2:MAXGEN){SRtime[i]<-SRtime[i-1]+(c.prob[i-1]*(SR.pool-SRtime[i-1]))/(distance.exp*distance)-(e.prob[i-1]*SRtime[i-1]/(Area^area.exp))}
  
  return(list(SRtime=SRtime,Time=Time))
}

MAXGEN<-50
IBG_output<-IBG_fun(MAXGEN = MAXGEN)
attach(IBG_output)

#########Sample time series shorter than MAXGEN#######
#Generate and store 100 time series from model. 
#Next take many short sample time series from these 100 series. These should start at random location in the series and extract continuous series 
#while respecting the time order of the species richness values. We create time series with durations of 10, 20, 30, 40, 50 years.
#We calculate the effect size of the trend (see below) for each of sample time series and store these values. 
#With the fake dataset we can calculate mean effect size using a hierarchical mixed model, with duration as a predictor.

reps<-100
sample_lengths<-seq(5,50,by=5)

sample_output<-data.frame(log_ratio=NA,slope=NA,rep=1:reps,sample_length=rep(sample_lengths,each=reps))
full_TS<-matrix(NA,50,reps)
sampled_TS<-data.frame(SR=NA, duration=NA,rep=NA, Time=NA)
for(r in 1:reps){
  output<-IBG_fun()
  Time<-output$Time
  SRtime<-output$SRtime
  full_TS[,r]<-SRtime
  for(i in 1:length(sample_lengths)){
    sample_length<-sample_lengths[i]
    sample_start<-sample(Time[1:(length(Time)-sample_length)],1)
    sample_TS<-SRtime[sample_start:(sample_start+sample_length-1)]
    sample_output$log_ratio[sample_output$sample_length==sample_length & sample_output$rep==r]<-log(sample_TS[length(sample_TS)]/sample_TS[1])
    sample_output$slope[sample_output$sample_length==sample_length & sample_output$rep==r]<-coef(lm(sample_TS~c(1:sample_length)))[2]
    sampled.df<-data.frame(SR=sample_TS,duration=rep(sample_length,sample_length),rep=r, Time=c(1:sample_length))
    sampled_TS<-rbind(sampled_TS,sampled.df)
    #sample_output$slope[sample_output$sample_length==sample_length & sample_output$rep==r]<-coef(gls(sample_TS~c(1:sample_length),cor=corAR1(0.7)))[2]
  }}
sampled_TS<-sampled_TS[-1,]

#matplot full series dataframe
colr <- rgb(0.1,0.1,0.1,0.2, names = NULL, maxColorValue = 1)
matplot(full_TS, type="l", col=colr, lty=1, lwd = 2, xlab="time", ylab="species richness")
lines(rowMeans(full_TS), type='l', lwd=3)

#plot(SRtime ~ Time, data = x.df, xlab='time (years)', ylab='species richness', type ='l' )
#abline(m2)

#plot summaries for full time series
#log ratio boxplot
ggplot(sample_output,aes(x=as.factor(sample_length),y=log_ratio))+geom_hline(yintercept = 0.0) + 
  geom_boxplot(notch = TRUE)+ labs(x = "duration (years)")+
  theme_grey()+theme(axis.text=element_text(size=16),axis.title=element_text(size=16))
#slope estimates boxplot
ggplot(sample_output,aes(x=as.factor(sample_length),y=slope))+geom_hline(yintercept = 0.0) + 
  geom_boxplot(notch = TRUE)+  labs(x = "duration (years)")+
  theme_grey()+theme(axis.text=element_text(size=16),axis.title=element_text(size=16))

#histogram effect size
hist(sample_output$log_ratio)
hist(sample_output$slope)
hist(SR)

#END  
