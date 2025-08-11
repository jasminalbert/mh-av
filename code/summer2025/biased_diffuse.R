# Biased dispersal on 2D grid (4-neighbor). 

#N = matrix(c(0.4,0.2,0.05,0.3),ncol=2)
#E = (1-N)/3
biased_diffuse <- function(N, E, D = 0.2, beta = 1, wrap = 0, n_nbs=8) {
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
    c(-1,  0), # up
    c( 1,  0), # down
    c( 0, -1), # left
    c( 0,  1), # right
    c(-1, -1), # up-left
    c(-1,  1), # up-right
    c( 1, -1), # down-left
    c( 1,  1))  # down-right
  if (n_nbs==4){nbrs<-nbrs[1:4]}

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

  
  # Inflow from all neighbors
  inflow <- matrix(0, nr, nc)
  for (d in nbrs) {
  	d=nbrs[[i]]
    # Abundance at origin
    N_origin <- shift(N, -d[1], -d[2]);N_origin
    # Sum of weights at origin
    sumw_origin <- shift(sum_w_origin, -d[1], -d[2]);sumw_origin
    # 🔹 FIX: avoid NaNs by replacing zeros in per-origin sums
    sumw_origin[sumw_origin == 0] <- 1
    # Add contribution to current cell
    inflow <- inflow + D * N_origin * (w / sumw_origin);inflow
    i=i+1
  }

  # remaining residents that didn't disperse
  stay <- (1 - D) * N
  
  N_new <- stay + inflow
  return(list(E=E,N=N,N_new=N_new,inflow=inflow, stay=stay,w=w))
}


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