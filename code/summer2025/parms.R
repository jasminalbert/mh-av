r1 <- 5
r2 <- 5# 0.05
alpha11<- 0.01#0.02
alpha22<- 0.01#2
alpha12<- 0.03
alpha21<- 0.03#6
g1<- 0.5
g2<- 0.3#1
lambda1<-r1
lambda2<-r2
popparms <- c(r1=r1,r2=r2,alpha11=alpha11,alpha22=alpha22, alpha12=alpha12,alpha21=alpha21,g1=g1,g2=g2)


l1 <- 0.1
l2 <- 0.1#2 #AV's N to litter ratio
hali1 <- 4#4 #half life
hali2 <- 0.5
theta1 <- 0
theta2 <- 0.1 #AV's germination response to litter
rho1 <- 0.03
rho2 <- 0 #AV's fecundity response to litter
litparms <-c(l1=l1,l2=l2,hali1=hali1,hali2=hali2,theta1=theta1,theta2=theta2,rho1=rho1,rho2=rho2)
#litparms <-c(l1=l1,l2=l2,hali1=hali1,hali2=hali2,theta1=0,theta2=0,rho1=0,rho2=0)

d1 <- 0.1
d2 <- 0.1
gd1 <- 0.001
gd2 <- 0.001
#dispparms <- c(d1=d1,d2=d2,gd1=gd1,gd2=gd2)
#dispparms <- c(d1=0,d2=0,gd1=0,gd2=0)
dispparms <- c(d1=d1,d2=d2,gd1=0.001,gd2=0.001)
grid_size <- 5
timesteps <- 20
