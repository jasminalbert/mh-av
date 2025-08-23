
litdifa <- array(0,dim=dim(lit))
litdifr <- litdifa
for(t in 1:timesteps){
	#sp2[,,"seeds",t]
	litdifa[,,t] <- litDif(lit[,,t])
	litdifr[,,t] <- rank(-litdifa[,,t])
	#sp2[,,"dispered",t] <- biased_diffuse(sp2[,,"seeds",t],litdif,D=D0,beta=5)$N_new
}
#sp2[,,"dispered",1:5]
timesteps=18
dd <- rep(D1,timesteps)
dd[1:(tD-1)] <- 0
y <- seq(0,1,length.out=gridsize^2)
percent<-F;t1<-10
source("plotmovementfunction.R")
pdf("../../figures/movementtheta1_long.pdf",width=13)
par(mfrow=c(1,1),xpd=NA)
plotmovement(sparray=sp2,theta=1,pall="Dark 3",trim=0,w=1)
dev.off()

percent=F
pdf("../../figures/movsp1theta1grid3.pdf",width=13)
par(mfrow=c(1,1),xpd=NA)
plotmovement(sparray=sp1,theta=0,pall="Dark 3",trim=0,w=1)
dev.off()

#t=6

v <- hcl.colors(100,"Viridis",rev=T)
h<-1:100
plot(y=rep(1,100),x=log(h)/50,type="p",pch=19,col=v,cex=5)#,ylim=c(1,5))
plot(y=rep(2,100),x=log(h),type="p",pch=19,col=v,cex=5)#,ylim=range(P))
plot(y=rep(3,100),x=h,type="p",pch=19,col=v,cex=5)#,ylim=range(P))
#/max(log(h))
#/max(log(h/5))
plot(x=log(h/5),type="p",pch=19,col=v,cex=5,ylim=c(-2,10))
P<-ceiling(out$Pout*100)
rank(P,T,"min")
order(P,T,F,"min")
pv[pv[,2]%in%log(P),1]
pv <- data.frame(v,round(log(1:100),2))
pv<-pv[!duplicated(pv[,2]),]

plot(pv[,2],col=pv[,1],pch=19,cex=5)
Pvec <- c(1,55,88,89,91,100,100);i=1
P <- seq(1,100,20)
P=Pvec#[i];i=i+1;P;log(P)
plot(P,log(P),pch=19,col=pv[pv[,2]%in%round(log(P),1),1][as.factor(P)],cex=5,xlim=c(1,100))
plot(P,pch=19,col=cols[P],cex=5,xlim=c(1,100))

i=i+1;out <- trackOr(dif);out=out[i,];out
#if bigger difference, should be closer to x
points(x=x[1:nrow(out)]+.3,y=out$yhlf+(out$y0-out$y1)/3, pch=15,cex=6,col="white")

# Example: emphasize high-end values

vals <- 1:100
ncol <- 100
pal  <- colorRampPalette(c("navy","cyan","yellow","red"))(ncol)
pal <-v#[-96:-100]
# Color map with high-end emphasis (power = 2)
scaled <- stretch_high(vals, power = 3)
cols   <- pal[1 + floor(scaled*(ncol-1))]

barplot(rep(1, 98), col = cols, border = NA, space = 0, axes = FALSE)


yd <- seq(-1,1,1/3)[-4]
exp(1/abs(yd)*.08)/4.25 #for x align
yd/abs(yd)*exp(abs(yd)*.75)*35
yd/abs(yd)*exp(abs(yd)*.5)*45
text(6.3,.1,"hello",srt=45)
yd/abs(yd)*exp(abs(yd))*35
ya <- abs(yd)
xa <- ya/abs(ya)*exp(abs(ya)*.5)*45;xa
pdf('../../figures/angles.pdf');plot(ya,xa,type="b");dev.off()
weight <- seq(0,0.25,0.02)
exp(weight*5)*3-2.5 

50/exp(abs(yd))[5]
y[2]
y1-y0
sum(orgins[orgins$y==y[2],]$Nout,na.rm=T)
originsall = orgins
orgins=originsall
orgins=orgins[orgins$y==y[2],]
xy <- expand.grid(i=1:nr,j=1:nc)
xy$y<-y
((y1-y0)*80)[f]
atan(1/3)/(pi)*180
90-atan((y1-y0)[f])/pi*180
xy[xy$i==orgins[1,]$out2_i & xy$j==orgins[1,]$out2_j,]$y
xy[orgins$out2_i%in%xy$i & xy$j%in%orgins$out2_j,]$y
i=orgins$out2_i;j=orgins$out2_j

orgins
# orgins <- expand.grid(i=1:nr, j=1:nc,nb=1:8)
# orgins$Nout <- NA
# orgins$out2_i <- NA 
# orgins$out2_j <- NA 
# orgins$dx<-NA;orgins$dy<-NA
# orgins[orgins$i%in%which(N_origin>0,T)[,1] & orgins$j%in%which(N_origin>0,T)[,2],]

# rows_src <- (1:nr) - dx
# cols_src <- (1:nc) - dy
# nb=1
# i=1;d=nbrs[[i]]
# for (d in nbrs){
	# dx=-d[1];dy=-d[2]
	# N_origin <- shift(N, dx, dy)
    # sumw_origin <- shift(sum_w_origin, dx, -dy)
    # sumw_origin[sumw_origin == 0] <- 1
    # inflowd <- D * N_origin * (w / sumw_origin)
    # inflow <- inflow + inflowd
    
    # select_o <- (orgins$i%in%(which(N_origin>0,T)[,1] -dx) & orgins$j%in%(which(N_origin>0,T)[,2] -dy)) 
    # select_x <- (orgins$i%in%which(N_origin>0,T)[,1] & orgins$j%in%which(N_origin>0,T)[,2]) 
    # orgins[select_o & orgins$nb==nb,]$Nout<- inflowd[inflowd>0]
    # orgins[select_o & orgins$nb==nb,c("out2_i","out2_j")] <- orgins[select_x & orgins$nb==nb,c("i","j")] 
    # orgins$dy[orgins$nb==nb]<-dy; orgins$dx[orgins$nb==nb]<-dx
    # nb=nb+1
# }
# i=i+1;d=nbrs[[i]]

# sum(orgins[orgins$y==0,]$Nout,na.rm=T)
# inflow
# N*D
# t=1






