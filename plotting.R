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
data_posterior = function(xseq, yseq, x, pi, mu, sigma, B){
  #p0 = sum(sapply(seq_along(pi), function(i) pi[i] * dlrnm(x, mu[,i], sigma[,,i], B)))
  gH = expand.grid(h1 = xseq, h2 = yseq)
  
  gHmixt = sapply(seq_along(pi), function(i) apply(gH, 1, function(h_, i) log_join_lrnm(x, h_, mu[,i], sigma[,,i], 
                                                                                        B, constant=FALSE), i))
  m0 = max(gHmixt)
  gHmixt = gHmixt - m0
  d = rowSums(sapply(seq_along(pi), function(i) pi[i] * exp(gHmixt[,i])))
  cL = contourLines(xseq, yseq, matrix(d, nrow = length(xseq), ncol = length(yseq)))
  cL
}
