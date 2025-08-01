# Parameters
lambda1 <- 1.5
lambda2 <- 1.4
beta12  <- 0.5
beta21  <- 0.6

# Discrete-time update functions
f1 <- function(x, y) x * lambda1 * (1 - x - beta12 * y)
f2 <- function(x, y) y * lambda2 * (1 - y - beta21 * x)

# Isoclines
isocline_x <- function(y) 1 - (1/lambda1) - beta12 * y
isocline_y <- function(x) 1 - (1/lambda2)- beta21 * x


equils <- function(){
  #beta12 <- alpha
  #beta21 <- beta
  # Solve for y*
  ystarx0 <- isocline_y(x=0)
  numerator <- 1 - beta21 * (1 - 1/lambda1) - 1/lambda2
  denominator <- 1 - beta21 * beta12
  y_star <- numerator / denominator
  # Solve for x* using x = 1 - 1/lambda1 - beta12 * y
  x_star <- isocline_x(y=y_star)
    #Alternative
  numerator <- (lambda1-1)/lambda1 - beta12*(lambda2-1)/lambda2
  xstar <- numerator/denominator
  ystar <- isocline_y(xstar)
  xstary0 <- isocline_x(y=0)
  if((x_star-xstar)>1e-10){print("Unequal calculations in x")}
  if((y_star-ystar)>1e-10){print("Unequal calculations in y")}
  return(c(xstary0=xstary0, ystarx0=ystarx0, xstar=x_star, ystar=y_star))
}
eq <- equils()
#eq <- equils(lambda1, lambda2, alpha=beta12, beta=beta21)
abline(v=xstar, col=2)
abline(v=x_star)
abline(h=y_star)
abline(h=ystar, col=2)


coexist_eq <- c(x_star, y_star)
cat("Coexistence equilibrium: (x*, y*) =", round(coexist_eq, 3), "\n")
library(ggplot2)

vecField <- function(x_vals,y_vals){
  # Define grid
  grid <- expand.grid(x = x_vals, y = y_vals)
  xran <- range(x_vals)
  yran <- range(y_vals)
  # Compute vector field
  grid$dx <- with(grid, f1(x, y) - x)
  grid$dy <- with(grid, f2(x, y) - y)
  # Normalize arrows for direction field
  mag <- sqrt(grid$dx^2 + grid$dy^2)
  grid$dxn <- grid$dx / mag
  grid$dyn <- grid$dy / mag
  
  eq <- equils()
  par(mgp=c(2,0.1,0), tcl=-0.3)
  plot(isocline_x(y_vals),y_vals,type='n', 
       xlim=range(x_vals),ylim=range(x_vals),
       xlab="N1",ylab="N2")
  arrows(x0=grid$x, y0=grid$y, x1=grid$dxn*0.05+grid$x, 
         y1=grid$dyn*0.05+grid$y, length=0.03, lwd=0.5)
  lines(x_vals,isocline_y(x_vals), col=2,lty=2)
  lines(isocline_x(y_vals),y_vals, col=3,lty=2)
  #xstar label
  points(eq['xstary0'],0)
  segments(x0=eq['xstary0'],y0=yran[1],y1=yran[1]-yran[2]*.06,xpd=NA)
  labelx <- bquote(.(round(eq['xstary0'],2)) == 1-frac(1,lambda[1]))
  text(x=eq['xstary0'],y=yran[1]-yran[2]*0.15,labels=labelx, xpd=NA, cex=0.8, adj=0.9)
  #ystar label
  points(0,eq['ystarx0'])
  segments(y0=eq['ystarx0'],x0=xran[1],x1=xran[1]-xran[2]*.06,xpd=NA)
  labely <- bquote(.(round(eq['ystarx0'],2)) == 1-frac(1,lambda[2]))
  text(y=eq['ystarx0'],x=xran[1]-xran[2]*0.15,labels=labely, xpd=NA, cex=0.8, adj=0.7)
  #title
  parnames <- c("lambda1","lambda2","alpha","beta")
  pars <- c(lambda1,lambda2,beta12,beta21)
  title(main=paste(parnames,pars,sep="=",collapse=', '),
        font.main=1, cex.main=0.8, line=0.1)
  #coexistence star label
  points(equils()["xstar"],equils()["ystar"],pch=8,col=4)
  points(xran[1],yran[1]-yran[2]*0.2,pch=8, xpd=NA,col=4)
  text(xran[2]*0.12,yran[1]-yran[2]*0.2, 
       labels=bquote( "("* .(round(eq['xstar'],2))*","*
        .(round(eq['ystar'],2))*")" ),xpd=NA, cex=0.8)
}
x_vals <- seq(0, 0.8, length.out = 20)
y_vals <- seq(0, 0.8, length.out = 20)
x_vals <- seq(0, 0.4, length.out = 20)
y_vals <- seq(0, 0.4, length.out = 20)
vecField(x_vals,y_vals)

