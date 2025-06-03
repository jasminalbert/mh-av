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
