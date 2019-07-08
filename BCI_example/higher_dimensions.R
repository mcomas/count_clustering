D = 20

set.seed(1)
X = matrix(0, ncol = D, nrow = 50)
while(min(rowSums(X)) == 0){
  IND = sample(1:ncol(BCI), D)
  X = BCI[,IND]
}

#fit = fit_lrnm(X, probs = TRUE, method = 'hermite', hermite.order = 3, max_iter = 1000, eps = 1e-5)
fit = fit_lrnm(X, probs = TRUE, method = 'montecarlo', Z = randtoolbox::sobol(1000, D-1, normal=TRUE),
               max_iter = 1000, eps = 0.001)

lM = lrnm_posterior_approx(as.matrix(X), fit$mu, fit$sigma, ilr_basis(D))
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
D = dist_lrnm(X)
hc = hclust(D, method = 'complete')
plot(hc, xlab = 'Observations')
clusters = cutree(hc, k = 3)

data(BCI.env, package = 'vegan')
ggplot() +
  geom_point(data=BCI.env, aes(x=UTM.EW, y=UTM.NS, col=factor(clusters)), size=6)

lapply(split(BCI.env, clusters), summary)
# If we take a closer look at observation 17.
# 
# ```{r}
# X[17,]
# ```
# 
# we can observe that this observation has a very small amount of counts, therefore, they behaviour is likeliky to be similar to the original gaussian modelling distribution (similar location and variability as the ones observed in the overall sample).
# 
# ```{r}
# I_ellip = ellipse(lM[[17]]$mu, lM[[17]]$sigma, 0.95)
# p + 
#   geom_text(data = H_E, aes(x = x1, y = x2, label = seq_along(x1)), alpha = 0.5, size = 4) +
#   geom_point(data = H_E[17,,drop=F], aes(x = x1, y = x2), col = 'blue') +
#   geom_path(data=I_ellip, aes(x = V1, y = V2), col = 'blue') +
#   coord_equal()
# ```
# 
# ```{r}
# data(BCI.env, package = 'vegan')
# ggplot() +
#   geom_point(data=BCI.env, aes(x=UTM.EW, y=UTM.NS, col=factor(clusters)), size=6)
# ```
# 
# ```{r}
# lapply(split(BCI.env, clusters), summary)
# ```
# 
