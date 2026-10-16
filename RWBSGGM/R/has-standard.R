standardize_cp_format <- function(cp_vec, n, has_boundary = FALSE) {
  
  cp_vec <- as.integer(cp_vec)
  cp_vec <- cp_vec[!is.na(cp_vec)]
  
  if (length(cp_vec) == 0) {
    return(c(1, n))
  }
  
  if (has_boundary) {
    cp_vec[cp_vec == 0] <- 1
    cp_vec[cp_vec > n] <- n
    cp_vec <- sort(unique(cp_vec))
    return(cp_vec)
  } else {
    cp_vec <- cp_vec[cp_vec > 1 & cp_vec < n]
    cp_vec <- sort(unique(cp_vec))
    return(c(1, cp_vec, n))
  }
}


hausdorff.distance <- function(v1, v2, n) {
  p1 = length(v1)
  p2 = length(v2)
  if (p1 * p2 == 0) {
    dis = n
  } else {
    distance.mat = matrix(0, nrow = p1, ncol = p2)
    for (i in 1:p1) {
      for (j in 1:p2) {
        distance.mat[i, j] = abs(v1[i] - v2[j])
      }
    }
    dis = max(max(apply(distance.mat, 1, min)),
              max(apply(distance.mat, 2, min)))
  }
  return(dis)
}

