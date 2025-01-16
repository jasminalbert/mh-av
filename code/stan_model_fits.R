#setwd("/Users/jasminalbert/Documents/HALLETTLAB/mh-av/code")
library(rstan)
source("./stan_model_fxn.R")
#option(mc.cores = parallel:detectCores()) #idk
#rstan_options(auto_write=TRUE) #idk

datLoc <- "../data/"
datname <- "modeldat.RDS"
data <- readRDS(paste0(datLoc,datname))

vars <-c("N","Fecundity", "intra","av","mh","P","Plot","mg","ag")

### TACA ### (mh)
## all treats ##
getdata(data,"TACA")

#dirty run
initials <- list(lambda=5, alpha_av=0.2, alpha_mh=0.0, epsilon=rep(1,P), sigma=10)
initials1 <- list(initials, initials, initials)

Sys.time()
mh_all_all <- stan(file="mh_bhmodel.stan", data=vars, iter=5000, chains=3, thin=3, control=list(adapt_delta = 0.9, max_treedepth=10), init=initials1)
Sys.time();beep(2)

#model fit
initials <- list(lambda=27, alpha_av=0.07, alpha_mh=0.03, epsilon=rep(1,P), sigma=10) #change? 
initials1 <- list(initials, initials, initials)

Sys.time()
mh_all_all <- stan(file="mh_bhmodel.stan", data=c("N","Fecundity","intra","av","mh","P","Plot","mg","ag"), iter=12000, chains=3, thin=3, control=list(adapt_delta = 0.9, max_treedepth=10), init=initials1)
Sys.time();beep(2)

traceplot(mh_all_all, pars=c("lambda","alpha_av","alpha_mh"))
pairs(mh_all_all, pars=c("lambda","alpha_av","alpha_mh"))

#save posterior distribution
save(mh_all_all, file="mh_all_all.rdata")

#look at resulting estimared parameter distributions
stan_dens(mh_all_all, pars=c("lambda","alpha_av","alpha_mh"))

#extract all parameter estimates
mh_all_all_ests <- rstan::extract(mh_all_all)
acf(mh_all_all_ests$lambda)

rm("N","Fecundity","intra","av","mh","P","Plot","mg","ag")

## xc xc ##
getdata(data,"TACA",nut="XC",ppt="XC")

#dirty run
initials <- list(lambda=5, alpha_av=0.2, alpha_mh=0.0, epsilon=rep(1,P), sigma=10)
initials1 <- list(initials, initials, initials)

Sys.time()
mh_xc_xc <- stan(file="mh_bhmodel.stan", data=vars, iter=5000, chains=3, thin=3, control=list(adapt_delta = 0.9, max_treedepth=10), init=initials1)
Sys.time();beep(2)

#model fit
initials <- list(lambda=27, alpha_av=0.07, alpha_mh=0.03, epsilon=rep(1,P), sigma=10) #change? 
initials1 <- list(initials, initials, initials)

Sys.time()
mh_xc_xc <- stan(file="mh_bhmodel.stan", data=vars, iter=12000, chains=3, thin=3, control=list(adapt_delta = 0.9, max_treedepth=10), init=initials1)
Sys.time();beep(2)

traceplot(mh_xc_xc, pars=c("lambda","alpha_av","alpha_mh"))
pairs(mh_xc_xc, pars=c("lambda","alpha_av","alpha_mh"))

#save posterior distribution
save(mh_xc_xc, file="mh_xc_xc.rdata")

#look at resulting estimared parameter distributions
stan_dens(mh_xc_xc, pars=c("lambda","alpha_av","alpha_mh"))

#extract all parameter estimates
mh_xc_xc_ests <- rstan::extract(mh_xc_xc)
acf(mh_xc_xc_ests$lambda)

rm("N","Fecundity","intra","av","mh","P","Plot","mg","ag")


### AVBA ### (av)
## all treats ##
getdata(data,"AVBA")

#dirty run
initials <- list(lambda=3, alpha_av=0, alpha_mh=0, epsilon=rep(1,P), sigma=10)
initials1 <- list(initials, initials, initials)

Sys.time()
av_all_all <- stan(file="av_bhmodel.stan", data=vars, iter=5000, chains=3, thin=3, control=list(adapt_delta = 0.9, max_treedepth=10), init=initials1)
Sys.time();beep(2)

#model fit
initials <- list(lambda=45, alpha_av=0.16, alpha_mh=0.01, epsilon=rep(1,P), sigma=10) #change? 
initials1 <- list(initials, initials, initials)

Sys.time()
av_all_all <- stan(file="av_bhmodel.stan", data=vars, iter=20000, chains=3, thin=3, control=list(adapt_delta = 0.9, max_treedepth=10), init=initials1)
Sys.time();beep(2)

traceplot(av_all_all, pars=c("lambda","alpha_av","alpha_mh"))
pairs(av_all_all, pars=c("lambda","alpha_av","alpha_mh"))

#save posterior distribution
save(av_all_all, file="av_all_all.rdata")

#look at resulting estimared parameter distributions
stan_dens(av_all_all, pars=c("lambda","alpha_av","alpha_mh"))

#extract all parameter estimates
av_all_all_ests <- rstan::extract(av_all_all)
acf(av_all_all_ests$lambda)
acf(av_all_all_ests$alpha_mh)
acf(av_all_all_ests$alpha_av)

rm("N","Fecundity","intra","av","mh","P","Plot","mg","ag")

## xc xc ##
getdata(data,"AVBA",nut="XC",ppt="XC")

#dirty run
initials <- list(lambda=30, alpha_av=0.0, alpha_mh=0.0, epsilon=rep(1,P), sigma=10)
initials1 <- list(initials, initials, initials)

Sys.time()
av_xc_xc <- stan(file="av_bhmodel.stan", data=vars, iter=5000, chains=3, thin=3, control=list(adapt_delta = 0.9, max_treedepth=10), init=initials1)
Sys.time();beep(2)

#model fit
initials <- list(lambda=30, alpha_av=0.1, alpha_mh=0.0, epsilon=rep(1,P), sigma=10) #change? 
initials1 <- list(initials, initials, initials)

Sys.time()
av_xc_xc <- stan(file="av_bhmodel.stan", data=vars, iter=20000, chains=3, thin=3, control=list(adapt_delta = 0.9, max_treedepth=10), init=initials1)
Sys.time();beep(2)

traceplot(av_xc_xc, pars=c("lambda","alpha_av","alpha_mh"))
pairs(av_xc_xc, pars=c("lambda","alpha_av","alpha_mh"))

#save posterior distribution
save(av_xc_xc, file="av_xc_xc.rdata")

#look at resulting estimared parameter distributions
stan_dens(av_xc_xc, pars=c("lambda","alpha_av","alpha_mh"))

#extract all parameter estimates
av_xc_xc_ests <- rstan::extract(av_xc_xc)
acf(av_xc_xc_ests$lambda)

rm("N","Fecundity","intra","av","mh","P","Plot","mg","ag")

load("av_xc_xc.rdata")
load("mh_xc_xc.rdata")
load("av_all_all.rdata")
load("mh_all_all.rdata")

str(av_xc_xc@sim$samples)
names(av_xc_xc@sim$samples[[1]])
av1 <- av_xc_xc@sim$samples[[1]]$lambda
av2 <- av_xc_xc@sim$samples[[2]]$lambda
av3 <- av_xc_xc@sim$samples[[3]]$lambda
hist(av1);abline(v=mean(av1),col='red');text(mean(av1),-100,paste(round(mean(av1),2)))

means <- function(out,parm,plot=T){
	chains <- length(out@sim$samples)
	parm1 <- out@sim$samples[[1]][[parm]]
	res <- array(NA,dim=c(length(parm1),chains))
	res[,1] <- parm1
	for (i in 2:chains){
		res[,i] <- out@sim$samples[[i]][[parm]]
	}
	mns <- colMeans(res)
	mn <- mean(mns)
	if (plot==T){
		for (i in 1:chains){
			histm(res[,i],main=parm)
		}
		histm(res,main=parm)
	}
	return(list(res,mns,mn))
}

histm <- function(av1,...){
	hist(av1,...);
	abline(v=mean(av1),col='red');
	text(mean(av1),-500,paste(round(mean(av1),2)),xpd=NA, col='darkred')
	abline(v=median(av1),col='green');
	text(median(av1),-700,paste(round(median(av1),2)),xpd=NA, col='darkgreen')
}

avxmn <- means(av_xc_xc, 'lambda')
resa <- list();resm <- resa
for (p in c('lambda','alpha_mh','alpha_av')){
	par(mfrow=c(2,2))
	resa[[p]] <- means(av_xc_xc,p)
	resm[[p]] <- means(mh_xc_xc,p)
}
res2 <- list(av=resa,mh=resm) #mg=ag=0.9
saveRDS(res1,"par_ests_1.rds")
saveRDS(res2,"par_ests_2.rds")

res1<-readRDS("par_ests_1.rds")
str(res1)
est1=sapply(res1,function(X){lapply(X,function(Y){Y[[3]]})})
est2=sapply(res2,function(X){lapply(X,function(Y){Y[[3]]})})
#effect of changing germination rate
matrix(est1)-matrix(est2)
diff(as.matrix(est1[,1]),as.matrix(est2[,1]))
str(as.matrix(est2))
est2 = matrix(unlist(est2),ncol=2)
est1 = matrix(unlist(est1),ncol=2)
est1-est2



