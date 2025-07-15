#popparms
# r1 <- 1:5
# r2 <- 5#1:5
# alpha11<- seq(0.01,0.05,0.01)
# alpha22<- alpha11
# alpha12<- 0.03
# alpha21<- seq(0.02,0.06,0.01)
# g1<- 0.5
# g2<- seq(0.1,0.5,0.1)
# popparmdf <- expand.grid(r1=r1,r2=r2,alpha11=alpha11,alpha22=alpha22, alpha12=alpha12, alpha21=alpha21, g1=g1,g2=g2)


r1 <- 1:5
r2 <- 5#1:5
alpha11<- 0.01
alpha22<- 0.01
alpha12<- 0.03
alpha21<- seq(0.02,0.04,0.005)
g1<- 0.5
g2<- 0.1
popparmdf <- expand.grid(r1=r1,r2=r2,alpha11=alpha11,alpha22=alpha22, alpha12=alpha12, alpha21=alpha21, g1=g1,g2=g2)



#litparms
l1 <- seq(0.01,0.1,length.out=5)
l2 <- l1 #AV's N to litter ratio
hali1 <- seq(0.5,4,length.out=4)#4 #half life
hali2 <- (1:3)/2
theta1 <- 0
theta2 <- seq(0.3,0.5,0.1) #AV's germination response to litter
rho1 <- seq(0.01,0.03,length.out=3)
rho2 <- 0 #AV's fecundity response to litter
#dispparms
d1 <- seq(0.01,0.5,length.out=3)
d2 <- d1
gd1 <- 0.001
gd2 <- 0.001

litparmdf <- expand.grid(l1=l1,l2=l2, hali1=hali1, hali2=hali2, theta1=theta1, theta2=theta2, rho1=rho1, rho2=rho2)
select <- litparmdf$hali1<litparmdf$hali2
litparmdf <- litparmdf[!select,]
litparmdf <- litparmdf[!select2,]
select2 <- litparmdf$l1>litparmdf$l2
dispparmdf <- expand.grid(d1=d1, d2=d2, gd1=gd1, gd2=gd2)
str(popparmdf)
str(litparmdf)
str(dispparmdf)
dispparmdf <- data.frame(t(dispparms))
litparmdf <- data.frame(t(litparms))
gridsizes=2
gridnames <- as.character(gridsizes)
poplist <- vector(mode="list",length=length(gridsizes))
names(poplist) <- gridnames 
grid <- gridnames
grid_size <- as.numeric(grid)


N1 <- matrix(c(100,5,10,80),nrow=grid_size)
N2 <- matrix(c(5,100,80,10),nrow=grid_size)
seepatches=T
timesteps=20
pdf("paramsim2.pdf")
for (p in 5:nrow(popparmdf)){
	pdf(paste0(p,"paramsim_.pdf"))
	for (l in 1:nrow(litparmdf)){
		for (d in 1:nrow(dispparmdf)){
			p.par.i <- unlist(popparmdf[p,])
			d.par.i <- unlist(dispparmdf[d,])
			l.par.i <- unlist(litparmdf[l,])
			cat("\n")
			poplist[[grid]] <- runsim(timesteps, N1, N2,p.par.i,d.par.i,l.par.i)
			cat("\n sim",p,l,d)			
			#if (seepatches){
			
			pop1end <- apply(poplist[[grid]]$sp1[,,(timesteps-5):timesteps],c(1,2),max)
			pop2end <- apply(poplist[[grid]]$sp2[,,(timesteps-5):timesteps],c(1,2),max)
			popsum <- sum(pop1end*pop2end)	
			#popsum <- sum(poplist[grid][[1]]$sp1[,,timesteps]*poplist[grid][[1]]$sp2[,,timesteps])
			if(popsum>0.5){	
				plotpatches(poplist[grid], ylim=200)
				cat("plotted",popsum)}
		}
	}
	dev.off()
}
dev.off()

#immediatelty crashed or no growth at all
	#2526-3125 g1=g2=0.5 & alpha22>0.01; otherwise sp2 grows when alpha22 is smaller. interesting sp2 can grow is alpha22 small
2525-3125
nrow(popparmdf[popparmdf$g1==0.5 &popparmdf$g2==0.5,])
	#1901-2500 g1=0.5,g2=0.4,alpha22>0.01
	#1276-1875 g1=0.5,g2=0.3,alpha22>0.01
	#651-1250 g1=0.5,g2=0.2,alpha22>0.01
	#606-625 g1=0.5,g2=0.1,alpha21=0.06,alpha12=0.03,alpha22=0.05,alpha11>0.1
	#581-600 g1=0.5,g2=0.1,alpha21=0.06,alpha12=0.03,alpha22=0.04,alpha11>0.1
	#556:575 g1=0.5,g2=0.1,alpha21=0.06,alpha12=0.03,alpha22=0.03,alpha11>0.1
	#531:550 g1=0.5,g2=0.1,alpha21=0.06,alpha12=0.03,alpha22=0.02,alpha11>0.1
	#506:525 g1=0.5,g2=0.1,alpha21=0.06,alpha12=0.03,alpha22=0.01,alpha11>0.1
	#481:500 g1=0.5,g2=0.1,alpha21=0.05,alpha12=0.03,alpha22=0.05,alpha11>0.1
	#456:475 g1=0.5,g2=0.1,alpha21=0.05,alpha12=0.03,alpha22=0.04,alpha11>0.1
	#431:450 g1=0.5,g2=0.1,alpha21=0.05,alpha12=0.03,alpha22=0.03,alpha11>0.1
	#406:425 g1=0.5,g2=0.1,alpha21=0.05,alpha12=0.03,alpha22=0.02,alpha11>0.1
	#381:400 g1=0.5,g2=0.1,alpha21=0.05,alpha12=0.03,alpha22=0.01,alpha11>0.1
	#356:375 g1=0.5,g2=0.1,alpha21=0.04,alpha12=0.03,alpha22=0.05,alpha11>0.1
	#331:350 g1=0.5,g2=0.1,alpha21=0.04,alpha12=0.03,alpha22=0.04,alpha11>0.1
	#306:325 g1=0.5,g2=0.1,alpha21=0.04,alpha12=0.03,alpha22=0.03,alpha11>0.1
	#281:300 g1=0.5,g2=0.1,alpha21=0.04,alpha12=0.03,alpha22=0.02,alpha11>0.1
	#231:250 g1=0.5,g2=0.1,alpha21=0.03,alpha12=0.03,alpha22=0.05,alpha11>0.1
	#206:225 g1=0.5,g2=0.1,alpha21=0.03,alpha12=0.03,alpha22=0.04,alpha11>0.1
	#181:200 g1=0.5,g2=0.1,alpha21=0.03,alpha12=0.03,alpha22=0.03,alpha11>0.1
	#156:175 g1=0.5,g2=0.1,alpha21=0.03,alpha12=0.03,alpha22=0.02,alpha11>0.1
	#106:125 g1=0.5,g2=0.1,alpha21=0.02,alpha12=0.03,alpha22=0.05,alpha11>0.1
	#81:100 g1=0.5,g2=0.1,alpha21=0.02,alpha12=0.03,alpha22=0.04,alpha11>0.1
	#56:75 g1=0.5,g2=0.1,alpha21=0.02,alpha12=0.03,alpha22=0.03,alpha11>0.1
	
no growth happens:
if g2 is more than 0.1 and alpha22 is more than 0.01 when alpha21<0.3
when alpha21>0.2, alpha22 can be less than 0.2 when g2>0.1
if g2=0.1, alpha11>0.01
if g2=0.1 and alpha21<0.05, alpha22>0.01  
if g2=0.1 and alpha21<0.03, alpha22>0.02

lowg2 <- popparmdf[popparmdf$g2==0.1 &lowg2$alpha11>0.01,]
nrow(lowg2)
check4 <- lowg2[lowg2$alpha21<0.03 & lowg2$alpha22>0.02,]
nrow(check4)
check3 <- lowg2[lowg2$alpha21<0.05 & lowg2$alpha21>0.02 &lowg2$alpha22>0.01,]
nrow(check3)
check2 <- lowg2[lowg2$alpha21>0.04,]
nrow(check2)
check1 <- popparmdf[popparmdf$g2>0.1,]
cut <- check1$alpha22<0.02 & check1$alpha21<0.03
check1 <- check1[!cut,]
nrow(check1)
(nogrowth[1:nrow(check4),-9]==check4)	#TRUE
(nogrowth[(nrow(check4)+1):(nrow(check4)+nrow(check3)),-9]==check3)	#TRUE
e <- nrow(check4)+nrow(check3)
(nogrowth[(e+1):(e+nrow(check2)),-9]==check2)	#TRUE
(nogrowth[(e+nrow(check2)+1):(e+nrow(check2)+nrow(check1)),-9]==check1)	#TRUE

nogrowth[520:530,]
check1[98:105,]

check <- rbind(popparmdf[popparmdf$alpha11>0.01 & popparmdf$g2==0.1,], popparmdf[popparmdf$alpha22>0.01 & popparmdf$g2>0.1,])
(nogrowth[1:500,-9]==check[1:500,])
head(check)
head(nogrowth,200)
select <-sort(c(2526:3125,1901:2500,1276:1875,651:1250,606:625,581:600,556:575,531:550,506:525,481:500,456:475,431:450,406:425,381:400,356:375,331:350,306:325,281:300,231:250,206:225,181:200,156:175,106:125,81:100,56:75))
nogrowth <- popparmdf[select,]
nogrowth <- cbind(nogrowth,select)

	
#sp1 grows
	#601-605
	#576-580	
	#551:555
	#526:530
	#501:505
	#476:480
	#451:455
	#426:430
	#401:405
	#376:380
	#351:355
	#326:330
	#301:305
	#276:280
	#226:230
	#201:205
	#176:180
	#151:155
	#101:105
	#76:80
	#51:55
		
#sp2 grows or wins all
	#2501-2525
	#1876-1900
	#1251-1275
	#626-650
	#256:275
	#131:150
	#31:50
	#6:25
	
#both grow
select <- sort(c(251:255,126:130,26:30,1:5))
bothgrow<- cbind(popparmdf[select,],select)
#always: g2=0.1, g1=0.5, alpha21>0.01, alpha12=0.03, alpha11=0.01, r1:[1:5]
#if alpha21=0.02, alpha22=0.02
#otherwise, alpha22=0.01
popparmdfx=popparmdf
popparmdf=bothgrow[,-9]	
	
popparmdf=popparmdfx
which(popparmdf==)
popparmdf[popparmdf$alpha21==0.03 & popparmdf$g2==0.1 & popparmdf$alpha11==0.01,]

#spatial coexistence 
#129 - the only one out of 3125 simulations that show coexistence
#and its exactly the parameters I've been using

#look at sims around those parameters and maybe simulate finer parms



p;l;d