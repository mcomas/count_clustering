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

Y = matrix(0, nrow = nrow(parliament), ncol = 3)
Y[,1] = with(parliament, catsp+other)
Y[,2] = with(parliament, jxcat+erc+cup) 
Y[,3] = with(parliament, cs+psc+pp)
colnames(Y) = c('other', 'ind', 'esp')


Yzr = data.table(zCompositions::cmultRepl(Y, suppress.print = T))
fit = fit_dm(Y)
Ynz = t(t(Y) + fit[,1])
Ynz = data.table(Ynz/rowSums(Ynz)*rowSums(Y))

B = ilr_basis(3)[c(3,1,2),]
rownames(B) = colnames(Ynz)
colnames(B) = paste0('B', 1:2)
B2 = B[,1:2]


Hzr = data.table(coordinates(Yzr, B2, label = 'ilr'))
Hnz = coordinates(Ynz, B2, label = 'ilr')



BASIS = as.data.table(t(sign(B2)),keep.rownames=TRUE)
names(BASIS)[1]=  'Basis:'
BASIS %>%
  flextable() %>%
  colformat_int(col_keys = colnames(Y))

Mzr = Mclust(Hzr, G = 2:20, modelNames = "EII")
Mnz = Mclust(Hnz, G = 2:20, modelNames = "EII")
Lzr = data_mixture(xseq = seq(-2,4, length.out = 100), yseq = seq(0, 4.5, length.out = 100), Mzr) %>%
  lapply(as.data.table)
Lnz = data_mixture(xseq = seq(-2,4, length.out = 100), yseq = seq(0, 4.5, length.out = 100), Mnz) %>%
  lapply(as.data.table)

PI_zr = Mzr$parameters$pro
MU_zr = Mzr$parameters$mean
SIGMA_zr = Mzr$parameters$variance$sigma

PI_nz = Mnz$parameters$pro
MU_nz = Mnz$parameters$mean
SIGMA_nz = Mnz$parameters$variance$sigma


gen_label_zr = function(){
  Hs_zr = apply(Y, 1, function(x_) c_rlrnm_mixture_posterior(1, x_, PI_zr, MU_zr, SIGMA_zr, t(MASS::ginv(B2)), r = 10)) %>% 
    t() %>% as.data.table()
  kmeansruns(Hs_zr, krange = 1:5)$cluster
}
gen_label_nz = function(){
  Hs_nz = apply(Y, 1, function(x_) c_rlrnm_mixture_posterior(1, x_, PI_nz, MU_nz, SIGMA_nz, t(MASS::ginv(B2)), r = 10)) %>% 
    t() %>% as.data.table()
  kmeansruns(Hs_nz, krange = 1:5)$cluster
}
C_zr = replicate(100, gen_label_zr())
C_nz = replicate(100, gen_label_nz())

Cx_zr = array(0, dim = c(nrow(Y), ncol(C_zr), 1, 1), dimnames = list(1:nrow(Y), 1:ncol(C_zr), 1, max(C_zr)))
Cx_zr[,,1,1] = C_zr
C_mv_zr = majority_voting(Cx_zr, is.relabelled = FALSE)
Cx_nz = array(0, dim = c(nrow(Y), ncol(C_nz), 1, 1), dimnames = list(1:nrow(Y), 1:ncol(C_nz), 1, max(C_nz)))
Cx_nz[,,1,1] = C_nz
C_mv_nz = majority_voting(Cx_nz, is.relabelled = FALSE)

save.image('IFCS2019/ifcs2019.RData')


