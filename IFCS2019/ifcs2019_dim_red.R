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

set.seed(1)

#Y = matrix(0, nrow = nrow(parliament), ncol = 3)
# Y[,1] = with(parliament, catsp+other)
# Y[,2] = with(parliament, jxcat+erc+cup) 
# Y[,3] = with(parliament, cs+psc+pp)
# colnames(Y) = c('other', 'ind', 'esp')
Y = as.matrix(parliament[,.(catsp,cs,cup,erc,jxcat,pp,psc)])

Yzr = data.table(zCompositions::cmultRepl(Y, suppress.print = T))
fit = fit_dm(Y)
Ynz = t(t(Y) + fit[,1])
Ynz = data.table(Ynz/rowSums(Ynz)*rowSums(Y))

B = pb_basis(Ynz, 'exact')
rownames(B) = colnames(Ynz)
colnames(B) = paste0('B', 1:ncol(B))
B2 = B[,1:2]


Hzr = data.table(coordinates(Yzr, B2, label = 'ilr'))
Hnz = coordinates(Ynz, B2, label = 'ilr')



BASIS = as.data.table(t(sign(B2)),keep.rownames=TRUE)
names(BASIS)[1]=  'Basis:'
BASIS %>%
  flextable() %>%
  colformat_int(col_keys = colnames(Y))

Mzr = Mclust(Hzr, modelNames = "EII", G=10)
Mnz = Mclust(Hnz, modelNames = "EII", G=10)
Lzr = data_mixture(xseq = seq(min(Hnz[,1])-1, max(Hnz[,1])+1, length.out = 100), 
                   yseq = seq(min(Hnz[,2])-1, max(Hnz[,2])+1, length.out = 100), Mzr) %>%
  lapply(as.data.table)
Lnz = data_mixture(xseq = seq(min(Hnz[,1])-1, max(Hnz[,1])+1, length.out = 100), 
                   yseq = seq(min(Hnz[,2])-1, max(Hnz[,2])+1, length.out = 100), Mnz) %>%
  lapply(as.data.table)

if(FALSE){
  p1_zr = ggplot() +
    geom_point(data = Hzr, aes(ilr1, ilr2), alpha=0.2)
  for(l in Lzr){
    p1_zr = p1_zr +
      geom_path(data=l, aes(x,y))
  }
  p1_zr = p1_zr +
    coord_fixed()
  p1_zr 
  #############
  p1_nz = ggplot() +
    geom_point(data = Hnz, aes(ilr1, ilr2), alpha=0.2)
  for(l in Lnz){
    p1_nz = p1_nz +
      geom_path(data=l, aes(x,y))
  }
  p1_nz = p1_nz +
    coord_fixed()
  p1_nz
}

I = 66
set.seed(1)
#########
PI_zr = Mzr$parameters$pro
MU_zr = Mzr$parameters$mean
SIGMA_zr = Mzr$parameters$variance$sigma
H_zr_s = c_rlrnm_mixture_posterior(500, Y[I,], PI_zr, MU_zr, SIGMA_zr, t(MASS::ginv(B2)), 1000) %>% as.data.frame()
l_post_zr = data_posterior(xseq = seq(min(H_zr_s[,1]), max(H_zr_s[,1]), length.out = 100),
                           yseq = seq(min(H_zr_s[,2]), max(H_zr_s[,2]), length.out = 100),
                           Y[I,], PI_zr, MU_zr, SIGMA_zr, B2) %>%
  lapply(as.data.table)
########
PI_nz = Mnz$parameters$pro
MU_nz = Mnz$parameters$mean
SIGMA_nz = Mnz$parameters$variance$sigma
H_nz_s = c_rlrnm_mixture_posterior(500, Y[I,], PI_nz, MU_nz, SIGMA_nz, t(MASS::ginv(B2)), 1000) %>% as.data.frame()
l_post_nz = data_posterior(xseq = seq(min(H_nz_s[,1]), max(H_nz_s[,1]), length.out = 100),
                           yseq = seq(min(H_nz_s[,2]), max(H_nz_s[,2]), length.out = 100),
                           Y[I,], PI_nz, MU_nz, SIGMA_nz, B2) %>%
  lapply(as.data.table)

if(FALSE){
  p2_zr = p1_zr
  for(post in l_post_zr){
    p2_zr = p2_zr + geom_path(data = post,  aes(x = x, y = y), col = 'blue', alpha=0.4)
  }
  p2_zr = p2_zr +
    geom_point(data = Hzr[I,], aes(x = ilr1, y = ilr2), col = 'blue', size = 2)
  # p2_zr = p2_zr +
  #   geom_point(data=H_zr_s, aes(x=V1,y=V2))
  p2_zr 
  
  p2_nz = p1_nz
  for(post in l_post_nz){
    p2_nz = p2_nz + geom_path(data = post,  aes(x = x, y = y), col = 'blue', alpha=0.4)
  }
  p2_nz = p2_nz +
    geom_point(data = Hnz[I,], aes(x = ilr1, y = ilr2), col = 'blue', size = 2)
  p2_nz 
}


if(FALSE){
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
    coord_fixed()
  p3_zr 
  
  p3_nz = ggplot() +
    geom_point(data = Hs_nz, aes(V1, V2), alpha=0.2)
  for(l in Lnz){
    p3_nz = p3_nz +
      geom_path(data=l, aes(x,y))
  }
  p3_nz = p3_nz +
    coord_fixed()
  p3_nz 
}

#####3 CLUSTERING
########
gen_label_zr = function(){
  Hs_zr = apply(Y, 1, function(x_) c_rlrnm_mixture_posterior(1, x_, PI_zr, MU_zr, SIGMA_zr, t(MASS::ginv(B2)), r = 100)) %>% 
    t() %>% as.data.table()
  kmeansruns(Hs_zr, krange = 1:10)$cluster
}
gen_label_nz = function(){
  Hs_nz = apply(Y, 1, function(x_) c_rlrnm_mixture_posterior(1, x_, PI_nz, MU_nz, SIGMA_nz, t(MASS::ginv(B2)), r = 100)) %>% 
    t() %>% as.data.table()
  kmeansruns(Hs_nz, krange = 1:10)$cluster
}
C_zr = replicate(100, gen_label_zr())
C_nz = replicate(100, gen_label_nz())

Cx_zr = array(0, dim = c(nrow(Y), ncol(C_zr), 1, 1), dimnames = list(1:nrow(Y), 1:ncol(C_zr), 1, max(C_zr)))
Cx_zr[,,1,1] = C_zr
C_mv_zr = majority_voting(Cx_zr)
Cx_nz = array(0, dim = c(nrow(Y), ncol(C_nz), 1, 1), dimnames = list(1:nrow(Y), 1:ncol(C_nz), 1, max(C_nz)))
Cx_nz[,,1,1] = C_nz
C_mv_nz = majority_voting(Cx_nz)

save.image('IFCS2019/ifcs2019_dim_red.RData')


