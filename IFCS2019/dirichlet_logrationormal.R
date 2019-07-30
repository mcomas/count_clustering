#devtools::install_github('mcomas/coda.base', ref = 'optimising_pb')
library(coda.base)
library(coda.count)
K = 3

load('IFCS2019/parliament_example_fitting.RData')
fit = fit_dm(X)
Xnz = t(t(X) + fit[,1])

B = pb_basis(Xnz, method = 'exact')
rownames(B) = colnames(Xnz)
colnames(B) = paste0('PB', 1:7)
B

H = coordinates(Xnz, B)

library(mclust)
m = Mclust(H)
I = 1
m$parameters$mean[,I]
m$parameters$variance$sigma[,,I]
m$z[1:5,]

L = lapply(1:m$G, function(I) lrnm_posterior_approx(X, m$parameters$mean[,I], m$parameters$variance$sigma[,,I], B))



mu = colMeans(H)
sigma = cov(H)

L = lrnm_posterior_approx(X, mu, sigma, B)
