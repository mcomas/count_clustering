#devtools::install_github('mcomas/coda.base', ref = 'optimising_pb')
library(coda.base)
library(coda.count)
library(magrittr)
K = 3
source('plotting.R')
load('IFCS2019/parliament_example_fitting.RData')
Y = matrix(0, nrow = nrow(X), ncol = 3)
Y[,1] = X[,'catsp'] + X[,'other']
Y[,2] = X[,'jxcat'] + X[,'erc'] + X[,'cup']
Y[,3] = X[,'psc'] + X[,'cs'] + X[,'pp']
colnames(Y) = c('other', 'ind', 'esp')

fit = fit_dm(Y)
# library(zCompositions)
# Xnz = cmultRepl(X, method = 'user', t = t(matrix(prop.table(fit[,1]), nrow = ncol(X), ncol = nrow(X))), s = sum(fit))
Ynz = t(t(Y) + fit[,1])
Ynz = Ynz/rowSums(Ynz)*rowSums(Y)

Yzr = zCompositions::cmultRepl(Y)

#B = pb_basis(Xnz, method = 'exact')
B = pb_basis(Ynz, method = 'exact')
rownames(B) = colnames(Ynz)
colnames(B) = paste0('PB', 1:2)
vs <-apply(coordinates(Ynz, B), 2, var)
cumsum(vs) / sum(vs)
B

library(mclust)
B2 = B[,1:2]
B2
Hnz = coordinates(as.data.frame(Ynz), B2)
Hzr = coordinates(as.data.frame(Yzr), B2)
m0nz = Mclust(Hnz)
m0zr = Mclust(Hzr)

m0 = m0zr
H2 = Hzr
C0 = m0$classification
library(ggtern)
ggtern() +
  geom_mask() +
  geom_point(data = as.data.frame(Y), aes(x=other, y=ind, z=esp, col=factor(C0))) +
  theme_minimal() + theme(legend.position = 'none')


mA = Mclust(H2, modelNames = "VII")
#m = Mclust(H2)
summary(mA)
plot(mA, 'uncertainty')

PI = mA$parameters$pro
MU = mA$parameters$mean
SIGMA = mA$parameters$variance$sigma

set.seed(1)
library(ggplot2)
theme_set(theme_minimal() + theme(legend.position='none'))
CONF = 0.4
p1 = ggplot()
p1 = p1 +
  geom_point(data = H2, aes(x = x1, y = x2), alpha = 0.4) +
  # geom_point(data = H2, aes(x = x1, y = x2, col = factor(m$classification)), alpha = 1) +
  coord_fixed()
for(i in 1:mA$G){
  p1 = p1 + geom_path(data = ellipse(MU[,i], SIGMA[,,i], CONF),  aes(x = V1, y = V2), col = 'red')
}
p1

pt1 = ggtern()
pt1 = pt1 +
  geom_point(data = composition(H2,B2), aes(x = x1, y = x2, z = x3), size = 1, alpha = 0.4) 
  # geom_point(data = H2, aes(x = x1, y = x2, col = factor(m$classification)), alpha = 1)
for(i in 1:mA$G){
  pt1 = pt1 + geom_path(data = composition(ellipse(MU[,i], SIGMA[,,i], CONF),B2),  aes(x = x1, y = x2, z = x3), col = 'red')
}
pt1


set.seed(1)
I = 60 #413#1
H_s = c_rlrnm_mixture_posterior(200, Y[I,], PI, MU, SIGMA, B2, 1000) %>% as.data.frame()
l_post = data_posterior(xseq = seq(min(H_s[,1]), max(H_s[,1]), length.out = 100),
                        yseq = seq(min(H_s[,2]), max(H_s[,2]), length.out = 100),
                        Y[I,], PI, MU, SIGMA, B2)
p2 = p1
for(post in l_post){
  p2 = p2 + geom_path(data = as.data.frame(post),  aes(x = x, y = y), col = 'blue', alpha=0.4)
}
p2 = p2 +
  geom_point(aes(x = H2[I,'x1'], y = H2[I,'x2']), col = 'blue', size = 2)
p2

pt2 = pt1
for(post in l_post){
  pt2 = pt2 + geom_path(data = composition(as.data.frame(post)[,2:3], B2),  aes(x = x1, y = x2, z = x3), col = 'blue', alpha=0.4)
}
pt2 = pt2 +
  geom_point(aes(x = Ynz[I,1], y = Ynz[I,2], z = Ynz[I,3]), col = 'blue', size = 2)
pt2


# p2 +
#   geom_point(data=H_s, aes(x = V1, y = V2), col = 'black', alpha=0.7)+
#   geom_point(aes(x = H2[I,'x1'], y = H2[I,'x2']), col = 'blue', size = 2)
# if(I == 90){
#   p2 +
#     coord_fixed(xlim=c(-0.45,-0.4), ylim = c(-1.05,-1)) +
#     geom_point(data=H_s, aes(x = V1, y = V2), col = 'black', alpha=0.7)
# }

library(fpc)
set.seed(1)
NSIM = 10
D = matrix(0, nrow=nrow(Y), ncol=nrow(Y)) 
H = matrix(0, nrow=nrow(Y), ncol=ncol(Y)-1)
C = matrix('', ncol=NSIM, nrow=nrow(Y))
for(i in 1:NSIM){
  print(i)
  Y1 = apply(Y, 1, function(x_) c_rlrnm_mixture_posterior(1, x_, PI, MU, SIGMA, t(MASS::ginv(B2)), r = 10)) %>% t()
  H = H + Y1
  D = D + as.matrix(dist(Y1))
  m = Mclust(Y1, k = m0$G, verbose = FALSE)
  C[,i] = as.character(m$classification)
}


LABS = unique(C0)
Cv = matrix(0, ncol=NSIM, nrow=nrow(Y))
for(i in 1:NSIM){
  LABS_n = unique(C[,i])
  f = expand.grid(replicate(length(LABS_n), LABS, simplify = FALSE))
  imax = which.max(apply(f, 1, function(ind){
    names(ind) = LABS_n
    sum(ind[C[,i]] == C0)
  }))
  ind.max = f[imax,]
  names(ind.max) = LABS_n
  Cv[,i] = unlist(ind.max[C[,i]])
}


CL = apply(sapply(1:4, function(cl) rowMeans(Cv == cl)), 1, which.max)
ggtern() +
  geom_mask() +
  geom_point(data = as.data.frame(Y), aes(x=other, y=ind, z=esp, col=factor(CL))) +
  theme_minimal() + theme(legend.position = 'none')

y_1 = 100 * Y/rowSums(Y)
y_1_sub = y_1[,]
(1:nrow(Y))[50 > y_1[,2] & y_1[,2] > 40 & 60 > y_1[,3] & y_1[,3] > 50][9]

ggtern() +
  geom_mask() +
  geom_point(data = as.data.frame(Y), aes(x=other, y=ind, z=esp, col=factor(CL))) +
  geom_point(aes(x=Y[413,1],y=Y[413,2],z=Y[413,3])) +
  theme_minimal() + theme(legend.position = 'none')

Cm = matrix(0, nrow = nrow(Y), ncol = nrow(Y))
for(i in 1:nrow(Y)){
  for(j in 1:nrow(Y)){
    Cm[i,j] = mean(C[i,] == C[j,])
  }
}
hc = hclust(as.dist(1-Cm))

ggtern() +
  geom_mask() +
  geom_point(data = as.data.frame(Y), aes(x=other, y=ind, z=esp, col=factor(cutree(hc, k = 2)))) +
  theme_minimal() + theme(legend.position = 'none')


### Consensus clustering

fit = cmdscale(as.dist(D))
eigen(cov(fit))
eigen(cov(H))

hc = hclust(D)
plot(hc)

head(Y1)
p3 = ggplot() +
  geom_point(data = Y1, aes(x = V1, y = V2), alpha = 1) +
  coord_fixed()
# for(i in 1:m$G){
#   p3 = p3 + geom_path(data = ellipse(MU[,i], SIGMA[,,i], CONF),  aes(x = V1, y = V2), col = 'red')
# }
suppressMessages(p3 + coord_fixed(xlim = c(-2.5, 0.5), ylim = c(-2.5, 0.5)))

# Until this step. The approach is similar to what a replacement approach would do.
# m$parameters$mean
# m$parameters$variance$sigma
# m$z[1:5,]



# Dirichlet-multinomial is not changin our initial sample
Y[1,]
Ynz[1,]
log(Y[1,]) - log(Y[1,1])
log(Ynz[1,]) - log(Ynz[1,1])

hsample = as.data.frame(c_rlrnm_2mixture_posterior(1000, X[ISEL,], 
                                                   PI[1], MU[,1], SIGMA[,,1],
                                                   PI[2], MU[2,], SIGMA[,,2], t(MASS::ginv(B2)), r = 100))





library(ggtern)
plot_x = function(X, g = NULL, part = NULL){
  if(!is.null(part)){
    He1 = ellipse(MU[,1], PI[1] * SIGMA[,,1], 0.95)
    Xe1 = composition(He1, B2)
    He2 = ellipse(MU[,2], PI[2] * SIGMA[,,2], 0.95)
    Xe2 = composition(He2, B2)
    colnames(Xe1) = names(X[[1]])
    colnames(Xe2) = names(X[[1]])
    p = ggtern() +
      geom_path(data = Xe1,  aes_string(x = part[1], y = part[2], z = part[3]), col = 'red') +
      geom_path(data = Xe2,  aes_string(x = part[1], y = part[2], z = part[3]), col = 'red')
    if(!is.null(X)){
      for(x in X){
        if(is.null(g)){
          p0 = PI[1] * dlrnm(x, MU[,1], SIGMA[,,1], B2) + PI[2] * dlrnm(x, MU[,2], SIGMA[,,2], B2)
          gH = data.frame(expand.grid(h1 = seq(-1, 5, length.out = 100), h2 = seq(-1, 3, length.out = 100)))
          gX = composition(gH, B2)
          colnames(gX) = colnames(X)
          gH1 = exp(apply(gH,1,function(h_) log_join_lrnm(x, unlist(h_), MU[,1], SIGMA[,,1], B2)))
          gH2 = exp(apply(gH,1,function(h_) log_join_lrnm(x, unlist(h_), MU[,2], SIGMA[,,2], B2)))
          gX$density =  (PI[1] * gH1 + PI[2] * gH2)/ p0
          L = contourLines(x=seq(-1, 5, length.out = 100), y=seq(-1, 3, length.out = 100), z = matrix(gX$density, ncol = 100))
          for(l in L){
            xl = as.data.frame(composition(cbind(l$x,l$y), B2))
            colnames(xl) = names(x)
            p = p +
              geom_path(data=xl, aes_string(x = part[1], y = part[2], z = part[3]), col = 'blue')
          }
        }else{
          p = p + g
        }
      }
    }
  }else{
    p = ggplot() +
      geom_path(data = ellipse(MU[,1], PI[1] * SIGMA[,,1], 0.95),  aes(x = V1, y = V2), col = 'red') +
      geom_path(data = ellipse(MU[,2], PI[2] * SIGMA[,,2], 0.95),  aes(x = V1, y = V2), col = 'red')
    if(!is.null(X)){
      for(x in X){
        p0 = PI[1] * dlrnm(x, MU[,1], SIGMA[,,1], B2) + PI[2] * dlrnm(x, MU[,2], SIGMA[,,2], B2)
        gH = data.frame(expand.grid(h1 = seq(-1, 5, length.out = 100), h2 = seq(-1, 3, length.out = 100)))
        gH1 = exp(apply(gH,1,function(h_) log_join_lrnm(x, unlist(h_), MU[,1], SIGMA[,,1], B2)))
        gH2 = exp(apply(gH,1,function(h_) log_join_lrnm(x, unlist(h_), MU[,2], SIGMA[,,2], B2)))
        gH$density =  (PI[1] * gH1 + PI[2] * gH2)/ p0
        if(is.null(g)){
          p = p +
            geom_contour(data = gH, aes(x=h1, y=h2, z=density), bins = 20) 
        }else{
          p = p + g
        }
      }
    }
    #####
  }
  p

}


ISEL = 686
x = X[ISEL,]

plot_x(list(x)) +
  coord_cartesian(xlim = c(0, 4), ylim = c(0,3)) +
  theme_minimal()
plot_x(list(x), part = c('jxcat', 'catsp', 'cs')) +
  theme_minimal()
plot_x(list(x), part = c('pp', 'cup', 'catsp')) +
  theme_minimal()

set.seed(2)
hsample1 = as.data.frame(c_rlrnm_posterior(round(1000*PI[1]), X[ISEL,], MU[1,], SIGMA[,,1], t(MASS::ginv(B2)), r = 100))
hsample2 = as.data.frame(c_rlrnm_posterior(round(1000*PI[2]), X[ISEL,], MU[2,], SIGMA[,,2], t(MASS::ginv(B2)), r = 100))
hsample = as.data.frame(c_rlrnm_2mixture_posterior(1000, X[ISEL,], 
                                                   PI[1], MU[,1], SIGMA[,,1],
                                                   PI[2], MU[2,], SIGMA[,,2], t(MASS::ginv(B2)), r = 100))
plot_x(list(x)) + 
  geom_point(data=rbind(hsample1,hsample2), aes(x=V1, y=V2), alpha=0.7) +
  theme_minimal()
plot_x(list(x)) + 
  geom_point(data=hsample, aes(x=V1, y=V2), alpha=0.7) +
  theme_minimal()

c_rlrnm_posterior(5, X[2,], MU[1,], SIGMA[,,1], t(MASS::ginv(B2)), r = 100)




X1 = apply(X, 1, function(x_) c_rlrnm_2mixture_posterior(1, x_, 
                                                         PI[1], MU[,1], SIGMA[,,1],
                                                         PI[2], MU[2,], SIGMA[,,2], t(MASS::ginv(B2)), r = 10))
X2 = apply(X, 1, function(x_) c_rlrnm_2mixture_posterior(1, x_, 
                                                         PI[1], MU[,1], SIGMA[,,1],
                                                         PI[2], MU[2,], SIGMA[,,2], t(MASS::ginv(B2)), r = 10))
X3 = apply(X, 1, function(x_) c_rlrnm_2mixture_posterior(1, x_, 
                                                         PI[1], MU[,1], SIGMA[,,1],
                                                         PI[2], MU[2,], SIGMA[,,2], t(MASS::ginv(B2)), r = 10))
H1 = as.data.frame(t(X1))
H2 = as.data.frame(t(X2))
H3 = as.data.frame(t(X3))
Hnz = as.data.frame(coordinates(Xnz, B2))
tot = rowSums(Hnz)

ggplot() +
  geom_path(data = ellipse(MU[,1], PI[1] * SIGMA[,,1], 0.95),  aes(x = V1, y = V2), col = 'red') +
  geom_path(data = ellipse(MU[,2], PI[2] * SIGMA[,,2], 0.95),  aes(x = V1, y = V2), col = 'red') +
  geom_point(data=Hnz, aes(x=x1, y=x2, alpha=tot/max(tot)))+
  coord_cartesian(xlim = c(-3, 6), ylim = c(-1, 3)) +
  theme(legend.position='none')

ggplot() +
  geom_path(data = ellipse(MU[,1], PI[1] * SIGMA[,,1], 0.95),  aes(x = V1, y = V2), col = 'red') +
  geom_path(data = ellipse(MU[,2], PI[2] * SIGMA[,,2], 0.95),  aes(x = V1, y = V2), col = 'red') +
  geom_point(data=H1, aes(x=V1, y=V2))+
  coord_cartesian(xlim = c(-3, 6), ylim = c(-1, 3))

ggplot() +
  geom_path(data = ellipse(MU[,1], PI[1] * SIGMA[,,1], 0.95),  aes(x = V1, y = V2), col = 'red') +
  geom_path(data = ellipse(MU[,2], PI[2] * SIGMA[,,2], 0.95),  aes(x = V1, y = V2), col = 'red') +
  geom_point(data=H1, aes(x=V1, y=V2))+
  coord_cartesian(xlim = c(-3, 6), ylim = c(-1, 3))

ggplot() +
  geom_path(data = ellipse(MU[,1], PI[1] * SIGMA[,,1], 0.95),  aes(x = V1, y = V2), col = 'red') +
  geom_path(data = ellipse(MU[,2], PI[2] * SIGMA[,,2], 0.95),  aes(x = V1, y = V2), col = 'red') +
  geom_point(data=H2, aes(x=V1, y=V2))+
  coord_cartesian(xlim = c(-3, 6), ylim = c(-1, 3))

ggplot() +
  geom_path(data = ellipse(MU[,1], PI[1] * SIGMA[,,1], 0.95),  aes(x = V1, y = V2), col = 'red') +
  geom_path(data = ellipse(MU[,2], PI[2] * SIGMA[,,2], 0.95),  aes(x = V1, y = V2), col = 'red') +
  geom_point(data=H3, aes(x=V1, y=V2)) +
  coord_cartesian(xlim = c(-3, 6), ylim = c(-1, 3))

# plot_x(list(x)) + 
#   geom_path(data = ellipse(N1$mu, PI[1] * N1$sigma, 0.95),  aes(x = V1, y = V2), col = 'green') +
#   geom_path(data = ellipse(N2$mu, PI[2] * N2$sigma, 0.95),  aes(x = V1, y = V2), col = 'green')
# plot_x(list(x), geom_contour(data = mH, aes(x=h1, y=h2, z=density), bins = 20, col = 'blue'))+
#   coord_cartesian(xlim = c(0, 4), ylim = c(0,3)) +
#   theme_minimal()
# plot_x(list(x), geom_contour(data = mH, aes(x=h1, y=h2, z=density), bins = 20, col = 'blue'), ternary = TRUE) +
#   theme_minimal()


L = lapply(1:m$G, function(I) lrnm_posterior_approx(X, m$parameters$mean[,I], m$parameters$variance$sigma[,,I], B2))



mu = colMeans(H)
sigma = cov(H)

L = lrnm_posterior_approx(X, mu, sigma, B)


# p = plot_x(list(x))
# p
# N1 = lrnm_posterior_approx(matrix(x, nrow=1), MU[,1], SIGMA[,,1], B2)[[1]]
# N2 = lrnm_posterior_approx(matrix(x, nrow=1), MU[,2], SIGMA[,,2], B2)[[1]]
# Z = m$z[ISEL,]
# mH = data.frame(expand.grid(h1 = seq(-1, 5, length.out = 100), h2 = seq(-1, 3, length.out = 100)))
# mH$density = Z[1] * mvtnorm::dmvnorm(mH, N1$mu, N1$sigma) +
#   Z[2] * mvtnorm::dmvnorm(mH, N2$mu, N2$sigma)



# # head(order(rowSums(X)))
# # 3, 8, 10, 14
# Lparl[['2017']][order(rowSums(X))[I], 1]
# Lparl[['2017']][order(rowSums(X))[I], -1]
# 
# I = I + 1
# plot_x(list(X[order(rowSums(X))[I],]))
# plot_x(list(X[331,]))
# plot_x(list(X[686,]))
# plot_x(list(X[311,]))
# plot_x(list(X[220,]))
# plot_x(list(X[220,]))
# plot_x(list(X[651,]))
# 
# X[331,]
# plot_x(X[331,])