# Biased dispersal on 2D grid (4-neighbor). 

#N = matrix(c(0.4,0.2,0.05,0.3),ncol=2)
#E = (1-N)/3
biased_diffuse <- function(N, E, D = 0.2, beta = 1, wrap = 0, n_nbs=8, origins=F) {
  stopifnot(is.matrix(N), is.matrix(E), all(dim(N) == dim(E)))
  nr <- nrow(N); nc <- ncol(N)
  
  # helper: shift matrix by (dx,dy). dx positive -> shift down (origin moves up -> dest)
  if (wrap=="tord"){
  shift <- function(mat, dx, dy) {
    rows <- ((1:nr) - dx - 1) %% nr + 1
    cols <- ((1:nc) - dy - 1) %% nc + 1
    mat[rows, cols]
  }}
  if (wrap==0){
  shift <- function(mat, dx, dy) {
    out <- matrix(0, nr, nc)
    rows_src <- (1:nr) - dx
    cols_src <- (1:nc) - dy
    valid_r <- which(rows_src >= 1 & rows_src <= nr)
    valid_c <- which(cols_src >= 1 & cols_src <= nc)
    out[valid_r, valid_c] <- mat[rows_src[valid_r], cols_src[valid_c]]
    out
  }}
  
  # weights based on environment at destination
  w <- exp(beta * E)
  # List of neighbor offsets (8 directions)
  nbrs <- list(
    up=c(-1,  0), # up
    down=c( 1,  0), # down
    left=c( 0, -1), # left
    right=c( 0,  1), # right
    upLeft=c(-1, -1), # up-left
    upRight=c(-1,  1), # up-right
    downLeft=c( 1, -1), # down-left
    downRight=c( 1,  1))  # down-right
  #if (n_nbs==4){nbrs<-nbrs[1:4]}
  nbrs <- nbrs[1:n_nbs]

    # for each origin, compute sum of neighbor weights (denominator for probabilities)  
  # Sum of destination weights for each origin
  #sum_w_origin <- matrix(0, nr, nc)
  #for (d in nbrs) {
  #  sum_w_origin <- sum_w_origin + shift(w, -d[1], -d[2])
  #}
  
  # Sum of neighbor weights for each origin
  sum_w_origin <- matrix(0, nr, nc)
  neighbor_count <- matrix(0, nr, nc) # to track actual number of neighbors
  
  for (d in nbrs) {
    w_dest <- shift(w, -d[1], -d[2])
    sum_w_origin <- sum_w_origin + w_dest
    neighbor_count <- neighbor_count + (w_dest > 0)  # counts non-zero-weight neighbors
  }
  
  
  # If no neighbor has weight > 0, fall back to equal split among actual neighbors
  sum_w_origin[sum_w_origin == 0 & neighbor_count > 0] <- neighbor_count[sum_w_origin == 0 & neighbor_count > 0]
  
  # if some origin has zero total weight (rare), replace with number of neighbors (4) -> uniform
  # Avoid division by zero (isolated cell → uniform)
  #sum_w_origin[sum_w_origin == 0] <- length(nbrs)
  
  
  # contributions to a destination cell from each neighbor origin:
  # for each direction, get N at the origin (i.e., shift N), and denominator is sum_w_origin at that origin
  orgins <-NA
  if (origins==T){
  	orgins <- expand.grid(i=1:nr, j=1:nc)
	orgins$Nout <- NA
	orgins$out2_i <- NA 
	orgins$out2_j <- NA 
	orgins$dx<-NA;orgins$dy<-NA
	dimnames <- list(NULL,names(orgins),names(nbrs))
	orgins <- array(rep(as.matrix(orgins),8),c(dim(orgins),8),dimnames)
  }

  
  # Inflow from all neighbors
  inflow <- matrix(0, nr, nc)
  dirname <- names(nbrs)
  inflowd <- array(0, c(nr,nc,n_nbs),dimnames=list(i=NULL,j=NULL,dir=dirname))
  nb=1
  for (n in 1:n_nbs) {
  	d <- nbrs[[n]]#;names(d)<-names(nbrs[n])
  	# Abundance at origin based on directional shift
  	#up, [1,] is 0 because none is above - so from d (1 is up)
  	#choose origins for all x based on directional shift
  	#indexed aat x 
  	dirname <- list(i=NULL,j=NULL,dir=names(nbrs[n])[1])
    N_origin <- array(shift(N, -d[1], -d[2]),c(nr,nc,1),dirname)
    # Sum of weights at origin
    sumw_origin <- shift(sum_w_origin, -d[1], -d[2])#;sumw_origin
    # 🔹 FIX: avoid NaNs by replacing zeros in per-origin sums
    sumw_origin[sumw_origin == 0] <- 1
    # Add contribution to current cell
    inflowd[,,n] <- D * N_origin[,,1] * (w / sumw_origin)
    inflow <- inflow + inflowd[,,n]
    if (origins==T & D>0){
    	dx=-d[1];dy=-d[2]
    	xy <- matrix(1:nc^2,nc=nc)
    	fromO<-diag(xy[(which(N_origin>0,T)[,1] -dx),(which(N_origin>0,T)[,2] -dy)])
    	toX <- diag(xy[which(N_origin>0,T)[,1],which(N_origin>0,T)[,2]])
    	orgins[fromO,"Nout",n] <-  inflowd[,,n][inflowd[,,n]>0]
    	orgins[fromO,c("out2_i","out2_j"),n] <- orgins[toX,c("i","j"),n]
    	orgins[,"dy",n]<-dy; orgins[,"dx",n]<-dx
    }
  }
  # remaining residents that didn't disperse
  stay <- (1 - D) * N
    #check 
    if(origins==T){
  if (sum((rowSums(orgins[,"Nout",],na.rm=T)+stay)-N<1e-10)<gridsize^2){
  	warning("sum out != N-stay")
  }}
  N_new <- stay + inflow
  return(list(E=E,N=N,N_new=N_new,inflow=inflow, stay=stay,w=w,origins=orgins))
}


    	# select_o <- (orgins$i%in%(which(N_origin>0,T)[,1] -dx) & orgins$j%in%(which(N_origin>0,T)[,2] -dy)) 
    	# select_x <- (orgins$i%in%which(N_origin>0,T)[,1] & orgins$j%in%which(N_origin>0,T)[,2]) 
    	# orgins[select_o & orgins$nb==nb,]$Nout<- inflowd[inflowd>0]
# orgins[select_o & orgins$nb==nb,c("out2_i","out2_j")] <- orgins[select_x & orgins$nb==nb,c("i","j")] 
# Environment with high quality in top-left
# E <- matrix(runif(25), nrow=5)
# E[1,1] <- 5  

# # Starting abundance in center
# N <- matrix(0, nrow=5, ncol=5)
# N[3,3] <- 100  

# # Move toward high-E areas
# biased_diffuse(N, E, D=0.5, beta=1)
# rd <- rd_ar[,,1,1]
# rd <- litDif(rd)
# seeds
# difdis(rd,seeds)
# biased_diffuse(N=seeds,E=-rd,D=0.4,beta=0)$N_new -N
# #when D=1, the >>highest abundance with highest attraction will go down because all of it is going elsewhere and whats coming in is split up between others that are low
# #try different values of D and beta 

# N;E
# Nt <- biased_diffuse(N,E)
# N-Nt$N_new
# sum(N)
# sum(Nt$N_new)
# w <- exp(beta * E)
# probability of N(x) moving to N(x') is the weight of E(x') divided by the sum of neighbors weights of N(x). 
# so its normalized by proportion of weight that should add up to one 

# #   contrib_from_up    <- D * shift(N, 1, 0)  * ( w / shift(sum_w_origin, 1, 0) )
  # contrib_from_down  <- D * shift(N, -1, 0) * ( w / shift(sum_w_origin, -1, 0) )
  # contrib_from_left  <- D * shift(N, 0, 1)  * ( w / shift(sum_w_origin, 0, 1) )
  # contrib_from_right <- D * shift(N, 0, -1) * ( w / shift(sum_w_origin, 0, -1) )