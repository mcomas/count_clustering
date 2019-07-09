library(coda.base)
library(coda.count)
load('pigs_example/Pigs.rda')
load('pigs_example/Pigs.Posture.RData')

X = Pigs
B = ilr_basis(ncol(X))

Z = randtoolbox::sobol(10000, 5, normal = TRUE)
fit = fit_lrnm(X, method = 'montecarlo', Z = Z, eps = 0.001)

fit[1:2]



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

lM = lrnm_posterior_approx(as.matrix(X), fit$mu, fit$sigma, B)
D = dist_lrnm(X)
hc = hclust(D, method = 'complete')
plot(hc, xlab = 'Observations')
clusters = cutree(hc, 3)

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
