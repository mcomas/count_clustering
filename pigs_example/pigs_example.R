library(coda.base)
library(coda.count)
load('pigs_example/Pigs.rda')
load('pigs_example/Pigs.Posture.RData')

X = Pigs
B = ilr_basis(ncol(X))

pars = fit_lrnm(X, B = B, method = 'hermite', hermite.order = 3, max_iter = 1000)


Bhattacharyya = function(N1, N2){
  SIGMA = (N1$sigma + N2$sigma) / 2
  invSIGMA = solve( SIGMA )
  d = 1/8 * t(N1$mu - N2$mu) %*% invSIGMA %*% (N1$mu - N2$mu) +
    1/2 * log( det(SIGMA) / sqrt(det(N1$sigma)*det(N2$sigma)) )
  d[1]
}

dist_lrnm = function(X){
  N = nrow(X)
  D = matrix(0, nrow = N, ncol = N)
  for(i in 1:N){
    for(j in 1:N){
      D[i,j] = Bhattacharyya(lM[[i]], lM[[j]])
    }
  }
  as.dist(D)
}


M = replicate(1000, {
  M = t(apply(as.matrix(X), 1, function(x) c_rlrnm_posterior(1, x, pars$mu, pars$sigma, Binv = t(MASS::ginv(B)), r = 100)))
  as.matrix(dist(M))
})
Mm = apply(M, 1:2, mean)
hc1 = hclust(as.dist(Mm), method = 'ward.D2')

lM = lrnm_posterior_approx(as.matrix(X), pars$mu, pars$sigma, B)
D = dist_lrnm(X)
hc2 = hclust(D, method = 'ward.D2')

par(mfrow=c(1,2))
plot(hc1, xlab = 'Observations', labels = dat.pos$Treatment)
plot(hc2, xlab = 'Observations', labels = dat.pos$Treatment)
par(mfrow=c(1,1))

clusters = cutree(hc, 3)
ggtern() +
  geom_mask() +
  geom_point(data = X, aes(x=BED, y = PASSAGE, z = FEEDER, col = dat.pos$Treatment))

table(clusters, dat.pos$Treatment)

Y = as.matrix(dat.pos[,3:8])
fit_Y = fit_lrnm(Y, method = 'montecarlo', Z = Z, eps = 0.001, B = B)

fit_Y[1:2]

lM = lrnm_posterior_approx(as.matrix(Y), fit_Y$mu, fit_Y$sigma, B)
N = nrow(Y)
D = matrix(0, nrow = N, ncol = N)
for(i in 1:N){
  for(j in 1:N){
    D[i,j] = Bhattacharyya(lM[[i]], lM[[j]])
  }
}
D = as.dist(D)
hc = hclust(D, method = 'complete')
plot(hc, xlab = 'Observations')
clusters = cutree(hc, 2)

table(dat.pos$Treatment, clusters)
