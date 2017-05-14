#####################################################################
# Title:Power and Effect Size
# Author: Steph de Silva.
# Email: steph@rex-analytics.com
# Date created: 12/05/17
# Date last altered: 13/05/17
# Attributions and acknowledgment of derivation:
# This script is due to information, tips and advice given in:
# Purpose: This script examines the power and size of a simple z statistic
# over different sample sizes and generating processes
#########################################################################
# Data Used: Generated
# Source: the computer
# Specifically: NA
# Translation by:NA
# Date Accessed: NA
# Gutenberg Number: NA
#########################################################################
# Script Outline:
# 1. Load Libraries/parameter settings
# 2. Open generation loops. Generate DGP
# 3. Calculate test statistic, implement decision rule
# 4. Close loops/calculate size and power
# 5. Plot differences
#########################################################################
# 1. Load libraries, define control parameters
#########################################################################

rm(list=ls(all=TRUE)) # clear the workspace of anything else
set.seed(170513) # we may want to do this again some day. 
library(ggplot2)
library(directlabels)
setwd("~/Documents/Rex Analytics/Blog/asymptotics")
sample_sizes = c(10, 30, 50, 100, 1000,10000,  100000)
mu <- c(0, 0.01, 0.1, 1, 10) # mu =0 is H0 (size), mu >0 is alternate/effect size
options = "cauchy" # generating error terms, see below
joption = 2 # degrees of freedom parameter for error generation, see below
iter = 5000 # iterations used

#########################################################################
# 2. Generate DGP Function, Open Loops, call function
#########################################################################

# Here, our true DGP is y(i)= mu +epsilon(i)
# epsilon(i) is iid with a mean of zero and variance constant
# there are several options for epsilon(i) we will explore
# (1) standard normal (Plain vanilla)
# (2) t(4) - all required moments for generalised estimation
# exist, but the tails are VERY fat. The distribution is symmmetric
# like normal.
# (3) centred and standardised chi squared (2). It's skewed, fat
# tailed and very, very difficult to estimate with.
# (4) Cauchy with central location zero and scale 1. This thing is a 
# statisticians nightmare. No moments >= first are defined, only fractionals.

# Define a function that will allow us to generate the data
# Generate y(i) = mu +epsilon(i)

DGP.process = function (n, mu, options, joption){
    # Generate an error term which is zero mean
    if (options=="normal"){
    epsilon = rnorm(n, mean=0, sd=1)
  } else if (options == "t"){
    epsilon = rt(n, joption, ncp=0)
    }else if (options == "chisq") {
    epsilon1 = rchisq(n, joption, ncp=0) # OK some backgroud. This does not have mean zero
    epsilon = (epsilon1 - joption) /sqrt(2*joption)  # the non central chi squared dist has moments
                                          # mu1 = joptions + ncp (non centrality parameter)
                                          # mu 2= (joptions +ncp)^2 + 2(joptions +2*ncp)
                                          # As ncp is zero here we have mean = joptions
                                          # and sd = sqrt(mu2-mu1^2)
                                          #        = sqrt(joptions ^2 + 2*joptions - joptions^2)
                                          #        = sqrt(2*joptions)
                                          # to standardise: epsilon = (epsilon1- mu1)/sd
                                          # More detail? C.f. https://en.wikipedia.org/wiki/Noncentral_chi-squared_distribution#Properties
    } else if (options =="cauchy") {
      epsilon = rcauchy(n, location=0, scale=1)
    }
  
  
  # Generate y
  
  y = mu + epsilon
  
  return (y)
}

rejections = matrix(data=0, nrow=length(sample_sizes), ncol=length(mu))
#########################################################################
# 3 Calculate test statistics
#########################################################################

for (i in 1:length(sample_sizes)) {
  n <- sample_sizes[i]
  for (j in 1: length(mu)) {
    mu_rep <- mu[j]
    for (r in 1:iter){
      y = DGP.process (n,mu_rep, options, joption)
      x_bar=mean(y)
      std_dev = sd(y)
      z = sqrt(n)*(x_bar)/std_dev
     
      if ((z)>=1.645){
        rejections[i,j] <- rejections [i,j] +1
      }
    }
#########################################################################
    # 4. Loops closing, calculate size/power
#########################################################################
    
  }
}
rejections<- rejections/iter
#########################################################################
# 5. Start Charting
#########################################################################
effect_size_list <-{}
for (k in 1:(length(mu))){
  effect_size_list <-c(effect_size_list, (paste("mu:", mu[k], sep="")))
}

size_output <- as.data.frame(cbind(sample_sizes, rejections[,1]))
size_output$alternative <- rep(effect_size_list[1], length(sample_sizes))
size_output <- as.data.frame(size_output)
size_output$Panel <- "1. Size" ### 
colnames (size_output) <- c("n", "rejections","alternative", "Panel")

power_output <- {}
for (k in 2:(length(mu))){
  power_output_tmp<-cbind(sample_sizes, rejections[,k], rep(effect_size_list[k], length(sample_sizes)))
  power_output<-rbind(power_output, power_output_tmp)
}  
power_output <- as.data.frame(power_output)
power_output$sample_sizes<-as.numeric(as.character(power_output$sample_sizes))
power_output$Panel <- "2. Power" ## 
colnames (power_output) <- c("n", "rejections", "alternative","Panel")


caption.data <- paste("Generating process: ", options, "distribution" )
if (options!="normal"){
  if(options!="cauchy"){
    caption.data <- paste (caption.data, " DF=", joption)
  }
}
chart.data <- as.data.frame(rbind(size_output, power_output))
# chart.data$n <-as.factor(chart.data$n)
chart.data$rejections <- as.numeric(chart.data$rejections)
col_vec<-c("slategray3","plum4","violetred", "steelblue4", "turquoise", "slateblue4")
p_ind<- ggplot(chart.data)+
  labs(x="Sample size", y="Rejection rates", caption=caption.data)+
  facet_grid(alternative~., scale="free_y")+
  geom_line(aes(n, rejections), group=1)+
  theme(plot.margin = unit(c(1,1,1,1), "lines"))+
  theme_light()+
 # scale_colour_manual(name="Index",values=col_vec)+
#  scale_linetype_manual(name="", values=c("solid","longdash", "dashed", "twodash", "dotted", "dotdash"))+
 # theme(legend.position="bottom")+
  scale_y_continuous(limits=c(0, 1))+
  scale_x_log10(breaks=sample_sizes, labels=sample_sizes)
  ggtitle("Size and power of the Z-test")
print(p_ind)

p_ind2 <-ggplot(chart.data)+
  labs(x="Sample size", y="Rejection rate", caption=caption.data)+
  facet_grid(Panel~., scale="free")+
  geom_line(aes(n, rejections, colour=alternative, linetype=alternative, group=alternative))+
  theme(plot.margin = unit(c(1,1,1,1), "lines"))+
  theme_light()+
  scale_colour_manual(name="Index",values=col_vec)+
  scale_linetype_manual(name="", values=c("solid","longdash", "dashed", "twodash", "dotted", "dotdash"))+
  theme(legend.position="left")+
  directlabels::geom_dl(aes(label=alternative), method="smart.grid")
  ggtitle("Size and power of the Z-test")
print(p_ind2)


p_ind3 <-ggplot(chart.data)+
  labs(x="Sample size", y="Rejection rate", caption=caption.data)+
  facet_grid(Panel~., scale="free")+
  geom_line(aes(n, rejections, colour=alternative, linetype=alternative, group=alternative))+
  theme(plot.margin = unit(c(1,1,1,1), "lines"))+
  theme_light()+
  scale_colour_manual(name="Index",values=col_vec)+
  scale_linetype_manual(name="", values=c("solid","longdash", "dashed", "twodash", "dotted", "dotdash"))+
  theme(legend.position="left")+
  ggtitle("Size and power of the Z-test")
print(p_ind3)

p_ind4 <-ggplot(chart.data)+
  labs(x="Sample size", y="Rejection rates", caption=caption.data)+
  facet_grid(Panel~., scale="free")+
  geom_line(aes(n, rejections, colour=alternative, group=alternative), show.legend=FALSE)+
  theme(plot.margin = unit(c(1,1,1,1), "lines"))+
  theme_light()+
  scale_colour_manual(name="Index",values=col_vec)+
  scale_x_log10(breaks=sample_sizes, labels=sample_sizes)
  ggtitle("Size and power of the Z-test")
direct.label(p_ind4, "last.points")
print(p_ind4)

#complicated <- list(dl.trans(x=x-2, y=y+0), rot=c(0), gapply.fun(d[-1:(-length(sample_sizes)+1),])) 
#direct.label(p_ind4,complicated)

direct.label(p_ind4, 
             list(rot=0, cex=0.75, dl.trans(x=x-1, y=y+0),
                  alpha=0.5,
                  gapply.fun(d[-1:(-length(sample_sizes)+1),])))


