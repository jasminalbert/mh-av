#popparms
r1 <- 1:5
r2 <- 1:5
alpha11<- seq(0.01,0.05,0.01)
alpha22<- alpha11
alpha12<- seq(0.02,0.06,0.01)
alpha21<- alpha12
g1<- seq(0.1,0.5,0.1)
g2<- g2
#litparms
l1 <- seq(0.01,0.1,length.out=5)
l2 <- l1 #AV's N to litter ratio
hali1 <- seq(0.5,4,length.out=4)#4 #half life
hali2 <- (1:3)/2
theta1 <- 0
theta2 <- seq(0.1,0.5,0.1) #AV's germination response to litter
rho1 <- seq(0.01,0.1,length.out=5)
rho2 <- 0 #AV's fecundity response to litter
#dispparms
d1 <- seq(0.1,0.5,0.1)
d2 <- d1
gd1 <- 0.001
gd2 <- 0.001

popparms <- expand.grid(r1=r1,r2=r2,alpha11=alpha11,alpha22=alpha22, alpha12=alpha12, alpha21=alpha21, g1=g1,g2=g2)
litparms <- expand.grid(l1=l1,l2=l2, hali1=hali1, hali2=hali2, theta1=theta1, theta2=theta2, rho1=rho1, rho2=rho2)
dispparms <- expand.grid(d1=d1, d2=d2, gd1=gd1, gd2=gd2)
str(popparms)
str(litparms)




