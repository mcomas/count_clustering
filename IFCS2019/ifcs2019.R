library(knitr)
library(data.table)
library(flextable)
library(ggplot2)
library(ggtern)
library(coda.base)
library(coda.count)
library(mclust)
library(magrittr)
library(fpc)
library(diceR)
source('plotting.R')
load('IFCS2019/parliament2017.RData')

Y = matrix(0, nrow = nrow(parliament), ncol = 3)
Y[,1] = with(parliament, catsp+other)
Y[,2] = with(parliament, jxcat+erc+cup) 
Y[,3] = with(parliament, cs+psc+pp)
colnames(Y) = c('other', 'ind', 'esp')


n0 = apply(Y,1,min)== 0
n1 = apply(Y,1,min)== 1
n2 = apply(Y,1,min)== 2
ggtern() +
  geom_mask() +
  geom_point(data = data.table(Y), aes(x=other, y=ind, z=esp), alpha=0.8) +
  geom_point(data = data.table(Y)[n0], aes(x=other, y=ind, z=esp), col = 'red', alpha=0.8) +
  geom_point(data = data.table(Y)[n1], aes(x=other, y=ind, z=esp), col = 'orange', alpha=0.8) +
  geom_point(data = data.table(Y)[n2], aes(x=other, y=ind, z=esp), col = 'yellow', alpha=0.8) +
  labs(title = '')

ZR_LAB = 'Zero replacement'
NZ_LAB = 'Dirichlet-multinomial smoothing'
Yzr = data.table(zCompositions::cmultRepl(Y, suppress.print = T))
fit = fit_dm(Y)
Ynz = t(t(Y) + fit[,1])
Ynz = data.table(Ynz/rowSums(Ynz)*rowSums(Y))


ggtern() +
  geom_mask() +
  geom_point(data = Yzr, aes(x=other, y=ind, z=esp), alpha=0.8) +
  geom_point(data = Yzr[n0], aes(x=other, y=ind, z=esp), col = 'red', alpha=0.8) +
  geom_point(data = Yzr[n1], aes(x=other, y=ind, z=esp), col = 'orange', alpha=0.8) +
  geom_point(data = Yzr[n2], aes(x=other, y=ind, z=esp), col = 'yellow', alpha=0.8) +
  labs(title = ZR_LAB)
ggtern() +
  geom_mask() +
  geom_point(data = Ynz, aes(x=other, y=ind, z=esp), alpha=0.8) +
  geom_point(data = Ynz[n0], aes(x=other, y=ind, z=esp), col = 'red', alpha=0.8) +
  geom_point(data = Ynz[n1], aes(x=other, y=ind, z=esp), col = 'orange', alpha=0.8) +
  geom_point(data = Ynz[n2], aes(x=other, y=ind, z=esp), col = 'yellow', alpha=0.8) +
  labs(title = NZ_LAB)

B = ilr_basis(3)[c(3,1,2),]
rownames(B) = colnames(Ynz)
colnames(B) = paste0('B', 1:2)
B2 = B[,1:2]


XLIM = c(-2, 4)
YLIM = c(0, 4.5)
Hzr = data.table(coordinates(Yzr, B2, label = 'ilr'))
p0_zr = ggplot() +
  geom_point(data = Hzr, aes(ilr1, ilr2)) +
  geom_point(data = Hzr[n0], aes(ilr1, ilr2), col = 'red') +
  geom_point(data = Hzr[n1], aes(ilr1, ilr2), col = 'orange') +
  geom_point(data = Hzr[n2], aes(ilr1, ilr2), col = 'yellow') +
  labs(title = 'Coordinates', subtitle = ZR_LAB) +
  coord_fixed(xlim = XLIM, ylim = YLIM)
p0_zr
#####
Hnz = coordinates(Ynz, B2, label = 'ilr')
p0_nz = ggplot() +
  geom_point(data = Hnz, aes(ilr1, ilr2)) +
  geom_point(data = Hnz[n0], aes(ilr1, ilr2), col = 'red') +
  geom_point(data = Hnz[n1], aes(ilr1, ilr2), col = 'orange') +
  geom_point(data = Hnz[n2], aes(ilr1, ilr2), col = 'yellow') +
  labs(title = 'Coordinates', subtitle = NZ_LAB) +
  coord_fixed(xlim = XLIM, ylim = YLIM)
p0_nz


BASIS = as.data.table(t(sign(B2)),keep.rownames=TRUE)
names(BASIS)[1]=  'Basis:'
BASIS %>%
  flextable() %>%
  colformat_int(col_keys = colnames(Y))

Mzr = Mclust(Hzr, G = 10, modelNames = "EVV")
Mnz = Mclust(Hnz, G = 10, modelNames = "EVV")
Lzr = data_mixture(xseq = seq(-2,4, length.out = 100), yseq = seq(0, 4.5, length.out = 100), Mzr) %>%
  lapply(as.data.table)
Lnz = data_mixture(xseq = seq(-2,4, length.out = 100), yseq = seq(0, 4.5, length.out = 100), Mnz) %>%
  lapply(as.data.table)

####
p1_zr = ggplot() +
  geom_point(data = Hzr, aes(ilr1, ilr2), alpha=0.2)
for(l in Lzr){
  p1_zr = p1_zr +
    geom_path(data=l, aes(x,y))
}
p1_zr = p1_zr +
  coord_fixed(xlim = XLIM, ylim = YLIM)
p1_zr + labs(title = 'Coordinates', subtitle = ZR_LAB)
#############
p1_nz = ggplot() +
  geom_point(data = Hnz, aes(ilr1, ilr2), alpha=0.2)
for(l in Lnz){
  p1_nz = p1_nz +
    geom_path(data=l, aes(x,y))
}
p1_nz = p1_nz +
  coord_fixed(xlim = XLIM, ylim = YLIM)
p1_nz + labs(title = 'Coordinates', subtitle = NZ_LAB)


####3
I = 66

set.seed(1)
#########
PI_zr = Mzr$parameters$pro
MU_zr = Mzr$parameters$mean
SIGMA_zr = Mzr$parameters$variance$sigma
H_zr_s = c_rlrnm_mixture_posterior(500, Y[I,], PI_zr, MU_zr, SIGMA_zr, B2, 1000) %>% as.data.frame()
l_post_zr = data_posterior(xseq = seq(min(H_zr_s[,1]), max(H_zr_s[,1]), length.out = 100),
                           yseq = seq(min(H_zr_s[,2]), max(H_zr_s[,2]), length.out = 100),
                           Y[I,], PI_zr, MU_zr, SIGMA_zr, B2) %>%
  lapply(as.data.table)
########
PI_nz = Mnz$parameters$pro
MU_nz = Mnz$parameters$mean
SIGMA_nz = Mnz$parameters$variance$sigma
H_nz_s = c_rlrnm_mixture_posterior(500, Y[I,], PI_nz, MU_nz, SIGMA_nz, B2, 1000) %>% as.data.frame()
l_post_nz = data_posterior(xseq = seq(min(H_nz_s[,1]), max(H_nz_s[,1]), length.out = 100),
                           yseq = seq(min(H_nz_s[,2]), max(H_nz_s[,2]), length.out = 100),
                           Y[I,], PI_nz, MU_nz, SIGMA_nz, B2) %>%
  lapply(as.data.table)

p2_zr = p1_zr
for(post in l_post_zr){
  p2_zr = p2_zr + geom_path(data = post,  aes(x = x, y = y), col = 'blue', alpha=0.4)
}
p2_zr = p2_zr +
  geom_point(data = Hzr[I,], aes(x = ilr1, y = ilr2), col = 'blue', size = 2)
# p2_zr = p2_zr +
#   geom_point(data=H_zr_s, aes(x=V1,y=V2))
p2_zr + labs(title = 'Coordinates', subtitle = ZR_LAB)
# p2_zr + geom_point(data=H_zr_s, aes(V1,V2))

p2_nz = p1_nz
for(post in l_post_nz){
  p2_nz = p2_nz + geom_path(data = post,  aes(x = x, y = y), col = 'blue', alpha=0.4)
}
p2_nz = p2_nz +
  geom_point(data = Hnz[I,], aes(x = ilr1, y = ilr2), col = 'blue', size = 2)
p2_nz + labs(title = 'Coordinates', subtitle = NZ_LAB)
p2_nz + geom_point(data=H_nz_s, aes(V1,V2))


##
Hs_zr = apply(Y, 1, function(x_) c_rlrnm_mixture_posterior(1, x_, PI_zr, MU_zr, SIGMA_zr, t(MASS::ginv(B2)), r = 100)) %>% 
  t() %>% as.data.table()
Hs_nz = apply(Y, 1, function(x_) c_rlrnm_mixture_posterior(1, x_, PI_nz, MU_nz, SIGMA_nz, t(MASS::ginv(B2)), r = 100)) %>% 
  t() %>% as.data.table()

p3_zr = ggplot() +
  geom_point(data = Hs_zr, aes(V1, V2), alpha=0.2)
for(l in Lzr){
  p3_zr = p3_zr +
    geom_path(data=l, aes(x,y))
}
p3_zr = p3_zr +
  coord_fixed(xlim = XLIM, ylim = YLIM)
p3_zr + labs(title = 'Coordinates', subtitle = ZR_LAB)

p3_nz = ggplot() +
  geom_point(data = Hs_nz, aes(V1, V2), alpha=0.2)
for(l in Lnz){
  p3_nz = p3_nz +
    geom_path(data=l, aes(x,y))
}
p3_nz = p3_nz +
  coord_fixed(xlim = XLIM, ylim = YLIM)
p3_nz + labs(title = 'Coordinates', subtitle = NZ_LAB)
