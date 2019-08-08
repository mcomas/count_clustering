ellipse = function(mu, sigma, p){
  mu = as.vector(mu)
  s = -2 * log(1 - p);
  ev = eigen(sigma * s)
  t_ = seq(0, 2 * base::pi, length.out = 500)
  a = mu + t(t(ev$vectors) * sqrt(ev$values)) %*% rbind(cos(t_), sin(t_)) 
  as.data.frame(t(a))
}
coda_ellipse = function(mu, sigma, p, B = ilr_basis(length(mu) + 1)){
  composition(ellipse(mu, sigma, p), B)
}
# x = X[1,]
# xseq = seq(-1, 5, length.out = 100)
# yseq = seq(-1, 3, length.out = 100)
data_posterior = function(xseq, yseq, x, pi, mu, sigma, B, NL = 10){
  #p0 = sum(sapply(seq_along(pi), function(i) pi[i] * dlrnm(x, mu[,i], sigma[,,i], B)))
  gH = expand.grid(h1 = xseq, h2 = yseq)
  
  gHmixt = sapply(seq_along(pi), function(i) apply(gH, 1, function(h_, i) log_join_lrnm(x, h_, mu[,i], sigma[,,i], 
                                                                                        B, constant=FALSE), i))
  m0 = max(gHmixt)
  gHmixt = gHmixt - m0
  d = rowSums(sapply(seq_along(pi), function(i) pi[i] * exp(gHmixt[,i])))
  cL = contourLines(xseq, yseq, matrix(d, nrow = length(xseq), ncol = length(yseq)), nlevels = NL)
  cL
}

dmixnorm <- function(x, Pi, Mu, S, part = 1:length(Pi), closure = TRUE) {
  # z = sample(x=1:length(Pi), size=n, prob=Pi, replace=T) rmn = matrix(0, nrow=n, ncol=nrow(Mu))
  if (is.vector(x)) {
    dmn <- 0
  } else {
    dmn <- rep(0, times = nrow(x))
  }
  for (i in part) {
    if (ncol(Mu) == 1){
      dmn <- dmn + Pi[i] * stats::dnorm(x, mean = Mu[, i], sd = sqrt(S[, , i]))
    }else{
      dmn <- dmn + Pi[i] * mvtnorm::dmvnorm(x, mean = Mu[, i], sigma = S[, , i])
    }
  }
  if(closure){
    dmn/sum(Pi[part])
  }else{
    dmn
  }
}
dmixnorm_mclust <- function(x, mclust_solution, ...) {
  func_pi <- mclust_solution$parameters$pro
  func_mean <- mclust_solution$parameters$mean
  func_sigma <- mclust_solution$parameters$variance$sigma
  dmixnorm(x, func_pi, func_mean, func_sigma, ...)
}

data_mixture = function(xseq, yseq, mixt){
  gH = expand.grid(h1 = xseq, h2 = yseq)
  d = apply(gH, 1, function(h_) dmixnorm_mclust(h_, mixt))
  cL = contourLines(xseq, yseq, matrix(d, nrow = length(xseq), ncol = length(yseq)))
  cL
}
