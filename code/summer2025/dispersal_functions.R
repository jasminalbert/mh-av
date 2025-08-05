# Helper function for dispersal
disperse1 <- function(mat, rate) {
  padded <- matrix(0, nrow = grid_size + 2, ncol = grid_size + 2)
  padded[2:(grid_size + 1), 2:(grid_size + 1)] <- mat
  
  neighbors <- padded[1:grid_size, 2:(grid_size + 1)] +    # up
               padded[3:(grid_size + 2), 2:(grid_size + 1)] +  # down
               padded[2:(grid_size + 1), 1:grid_size] +    # left
               padded[2:(grid_size + 1), 3:(grid_size + 2)]  # right
  
  new_mat <- (1 - rate) * mat + (rate / 4) * neighbors
  return(new_mat)
}

disperse2 <- function(mat, local_disp = 0.1, global_disp = 0.001) {
  n <- nrow(mat)
  
  # Distance-weighted 3x3 kernel (e.g., inverse-distance or Gaussian)
  kernel <- matrix(c(1/sqrt(2), 1, 1/sqrt(2),
                     1,     0, 1,
                     1/sqrt(2), 1, 1/sqrt(2)), nrow = 3, byrow = TRUE)
  kernel <- kernel / sum(kernel)  # Normalize

  # Toroidal pad function (wrap edges)
toroidal_pad <- function(mat) {
  n <- nrow(mat)
  m <- ncol(mat)
  
  # Rows
  top    <- mat[n, , drop = FALSE]
  bottom <- mat[1, , drop = FALSE]
  
  # Columns
  left   <- mat[, m, drop = FALSE]
  right  <- mat[, 1, drop = FALSE]
  
  # Corners
  topleft     <- mat[n, m, drop = FALSE]
  topright    <- mat[n, 1, drop = FALSE]
  bottomleft  <- mat[1, m, drop = FALSE]
  bottomright <- mat[1, 1, drop = FALSE]

  top_row    <- cbind(topleft, top, topright)
  middle_row <- cbind(left,    mat, right)
  bottom_row <- cbind(bottomleft, bottom, bottomright)

  padded <- rbind(top_row, middle_row, bottom_row)
  return(padded)
}



  # Apply local dispersal using convolution
  pad <- toroidal_pad(mat)
  disp_out <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      disp_out[i,j] <- sum(kernel * pad[i:(i+2), j:(j+2)])
    }
  }

  # Mix with original and global pool
  new_mat <- (1 - local_disp - global_disp) * mat +
             local_disp * disp_out +
             global_disp * mean(mat)

  return(new_mat)
}



#set.seed(12)
#test_mat <- matrix(runif(grid_size^2), nrow = grid_size, ncol = grid_size)
#test_out1 <- disperse2(test_mat)
#sum(test_mat);sum(test_out1)
#test_mat=test_out1;image(test_mat,zlim=c(0,1))
#test_out2 <- disperse1(test_mat,rate=0.1)
#sum(test_mat);sum(test_out2)
#test_mat=test_out2


disperse_gaussian <- function(mat, sigma = 1) {
  # Create 3x3 Gaussian kernel centered at (0,0)
  dx <- -1:1
  dy <- -1:1
  kernel <- outer(dx, dy, function(x, y) exp(-(x^2 + y^2) / (2 * sigma^2)))
  kernel <- kernel / sum(kernel)  # Normalize
  
  # Pad the matrix with zeros (no wraparound)
  padded <- matrix(0, nrow = nrow(mat) + 2, ncol = ncol(mat) + 2)
  padded[2:(nrow(mat)+1), 2:(ncol(mat)+1)] <- mat
  
  # Create output matrix
  result <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
  
  # Apply the kernel to each cell
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      neighborhood <- padded[i:(i+2), j:(j+2)]
      result[i, j] <- sum(kernel * neighborhood)
    }
  }
  
  return(result)
}

simpleDisperse <- function(D){
  normal <- matrix(rnorm(gridsize^2),ncol=gridsize)
  normal <- normal/sum((normal))
  id <- (1/normal*0.5)*(abs(normal))
  sum(normal+id)
}
ld <- litDif(matrix(abs(rnorm(gridsize^2)),ncol=gridsize))
ldn <- ld/(sum(abs(ld)))
sum(ldn)
seedslit=seedslit[,,5]
litdifnorm
sd1<-seedslit-0.3*litdifnorm
sd2<-seedslit-0.3*ldn
sum(seedslit)
sum(abs(sd1-seedslit))
sum(abs(sd2-seedslit))
sum(litdif)
