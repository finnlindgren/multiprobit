latent_to_wishart <- function(x, LV, df) {
  d <- nrow(LV)
  L <- diag(sqrt(qchisq(pnorm(x[seq_len(d)]), df = df - seq_len(d) + 1)))
  L[upper.tri(L)] <- x[d + seq_len(d * (d-1) / 2)]
  L <- t(L)
  LC <- LV %*% L
  list(LC = LC, L = L)
}

latent_from_wishart <- function(LC, LV, df) {
  d <- nrow(LC)
  L <- solve(LV, LC)
  x <- numeric(d * (d+1) / 2)
  x[seq_len(d)] <- qnorm(pchisq(diag(L)^2, df = df - seq_len(d) + 1))
  x[d + seq_len(d * (d-1) / 2)] <- t(L)[upper.tri(L)]
  list(x = x, L = L)
}


latent_to_nwishart <- function(x, LV, df) {
  d <- nrow(LV)
  s_vec <- c(1 / LV[1, 1], numeric(d - 1))
  L <- LC <- matrix(0, d, d)
  L[1, 1] <- LC[1, 1] <- 1
  prev_index <- cumsum(seq_len(d) - 1)
  for (k in seq_len(d - 1)) {
    L[k + 1, seq_len(k)] <- u <- x[prev_index[k] + seq_len(k)]
    lambda_vec <- (LV[k + 1, seq_len(k)] / LV[k + 1, k + 1]) %*%
      L[seq_len(k), seq_len(k), drop = FALSE]
    uu <- lambda_vec + u
    uu_sq <- sum(uu^2)
    uu_norm <- sqrt(uu_sq)
    lambda <- sum(lambda_vec^2)
    z <-
      qbeta(
        pchisq(uu_sq, df = k, ncp = lambda),
        shape1 = k / 2,
        shape2 = (df - k) / 2,
        ncp = lambda
      )
    tmp <- sqrt(z) / uu_norm
    s_vec[k + 1] <- tmp / LV[k + 1, k + 1]
    L[k + 1, k + 1] <- sqrt(1 - z) / tmp
    LC[k + 1, seq_len(k)] <- tmp * uu
    LC[k + 1, k + 1] <- sqrt(1 - z)
  }
  list(LC = LC, s = s_vec, L = L)
}

latent_from_nwishart <- function(LC, LV, df) {
  d <- nrow(LC)
  # Make sure the rows are normalised:
  LC <- rowSums(LC^2)^0.5 * LC
  s_vec <- c(1 / LV[1, 1], numeric(d - 1))
  x <- numeric(d * (d - 1) / 2)
  L <- matrix(0, d, d)
  L[1, 1] <- 1
  prev_index <- cumsum(seq_len(d) - 1)
  for (k in seq_len(d - 1)) {
    z <- 1 - LC[k + 1, k + 1]^2
    lambda_vec <-
      (LV[k + 1, seq_len(k), drop = FALSE] / LV[k + 1, k + 1]) %*%
      L[seq_len(k), seq_len(k), drop = FALSE]
    lambda <- sum(lambda_vec^2)
    uu_norm <-
      sqrt(qchisq(
        pbeta(z, shape1 = k / 2, shape2 = (df - k) / 2, ncp = lambda),
        df = k,
        ncp = lambda
      ))
    u <- (uu_norm / sqrt(z)) * LC[k + 1, seq_len(k), drop = FALSE] - lambda_vec

    # Store the latent variables
    x[prev_index[k] + seq_len(k)] <- u
    L[k + 1, seq_len(k)] <- u

    tmp <- sqrt(z) / uu_norm
    s_vec[k + 1] <- tmp / LV[k + 1, k + 1]
    L[k + 1, k + 1] <- LC[k + 1, k + 1] / tmp
  }
  list(x = x, s = s_vec, L = L)
}
