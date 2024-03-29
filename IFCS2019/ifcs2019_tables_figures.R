library(knitr)
library(kableExtra)
library(ggplot2)
library(ggtern)
library(coda.base)
library(coda.count)
library(data.table)
load('IFCS2019/ifcs2019.RData')
set.seed(1)
options(kableExtra.latex.load_packages = FALSE)

parliament %>%
  head(4) %>%
  rbind(as.list(rep("$\\vdots$", 9))) %>%
  kable("latex", booktabs = T, escape = FALSE, align = 'lrrrrrrrr') %>%
  kable_styling(font_size = 7, position = "center") %>%
  write(file = 'IFCS2019/pres/table01a.tex')

cbind('mun' = parliament$mun, Y[,c(2,3,1)]) %>%
  head(4) %>%
  rbind(as.list(rep("$\\vdots$", 4))) %>%
  kable("latex", booktabs = T, escape = FALSE, align = 'lrrr') %>%
  kable_styling(font_size = 7, position = "center") %>%
  write(file = 'IFCS2019/pres/table01b.tex')

theme_set(theme_minimal() + theme(legend.position = 'none', 
                                  plot.title =  element_text(size = 14),
                                  plot.subtitle = element_text(size = 14, face = 'italic'),
                                  axis.title = element_text(size = 13)))

n0 = apply(Y,1,min)== 0
n1 = apply(Y,1,min)== 1
n2 = apply(Y,1,min)== 2
p0 = ggtern() +
  geom_mask() +
  geom_point(data = data.table(Y), aes(x=other, y=ind, z=esp), alpha=0.8) +
  geom_point(data = data.table(Y)[n0], aes(x=other, y=ind, z=esp), col = 'red', alpha=0.8) +
  # geom_point(data = data.table(Y)[n1], aes(x=other, y=ind, z=esp), col = 'orange', alpha=0.8) +
  # geom_point(data = data.table(Y)[n2], aes(x=other, y=ind, z=esp), col = 'yellow', alpha=0.8) +
  labs(title = '')
ggsave(p0, file = 'IFCS2019/pres/ternary_original.pdf', width = 7, height = 5.5)

ZR_LAB = 'Zero replacement'
NZ_LAB = 'Dirichlet-multinomial smoothing'
p1_zr = ggtern() +
  geom_mask() +
  geom_point(data = data.table(Yzr), aes(x=other, y=ind, z=esp), alpha=0.8) +
  geom_point(data = data.table(Yzr)[n0], aes(x=other, y=ind, z=esp), col = 'red', alpha=0.8) +
  # geom_point(data = data.table(Y)[n1], aes(x=other, y=ind, z=esp), col = 'orange', alpha=0.8) +
  # geom_point(data = data.table(Y)[n2], aes(x=other, y=ind, z=esp), col = 'yellow', alpha=0.8) +
  labs(title = ZR_LAB)
p1_nz = ggtern() +
  geom_mask() +
  geom_point(data = data.table(Ynz), aes(x=other, y=ind, z=esp), alpha=0.8) +
  geom_point(data = data.table(Ynz)[n0], aes(x=other, y=ind, z=esp), col = 'red', alpha=0.8) +
  # geom_point(data = data.table(Y)[n1], aes(x=other, y=ind, z=esp), col = 'orange', alpha=0.8) +
  # geom_point(data = data.table(Y)[n2], aes(x=other, y=ind, z=esp), col = 'yellow', alpha=0.8) +
  labs(title = NZ_LAB)


ggsave(p1_zr, file = 'IFCS2019/pres/ternary_zr.pdf', width = 7, height = 5.5)
ggsave(p1_nz, file = 'IFCS2019/pres/ternary_nz.pdf', width = 7, height = 5.5)

BASISm = BASIS
names(BASISm)[1] = 'Basis $\\mathcal{B}$'
BASISm[,c(1,3,4,2)] %>%
  kable("latex", booktabs = T, escape = FALSE)  %>%
  kable_styling(font_size = 7, position = "center") %>%
  write(file = 'IFCS2019/pres/basis.tex')

COORD_LABS = function(p) p + xlab(expression(esp %<->% ind)) + 
  ylab(expression(other %<->% paste("", sqrt(paste("", "esp ind"), paste("2")))))

library(latex2exp)
XLIM = c(-2, 4)
YLIM = c(0, 4.5)
Hzr = data.table(coordinates(Yzr, B2, label = 'ilr'))
p2_zr = ggplot() +
  geom_point(data = Hzr, aes(ilr1, ilr2)) +
  geom_point(data = Hzr[n0], aes(ilr1, ilr2), col = 'red') +
  # geom_point(data = Hzr[n1], aes(ilr1, ilr2), col = 'orange') +
  # geom_point(data = Hzr[n2], aes(ilr1, ilr2), col = 'yellow') +
  labs(title = 'Coordinates', subtitle = ZR_LAB) +
  coord_fixed(xlim = XLIM, ylim = YLIM)
#####
Hnz = coordinates(Ynz, B2, label = 'ilr')
p2_nz = ggplot() +
  geom_point(data = Hnz, aes(ilr1, ilr2)) +
  geom_point(data = Hnz[n0], aes(ilr1, ilr2), col = 'red') +
  # geom_point(data = Hnz[n1], aes(ilr1, ilr2), col = 'orange') +
  # geom_point(data = Hnz[n2], aes(ilr1, ilr2), col = 'yellow') +
  labs(title = 'Coordinates', subtitle = NZ_LAB) +
  coord_fixed(xlim = XLIM, ylim = YLIM)
p2_nz

ggsave(p2_zr %>% COORD_LABS(), file = 'IFCS2019/pres/coordinates_zr.pdf', width = 7, height = 5.5)
ggsave(p2_nz %>% COORD_LABS(), file = 'IFCS2019/pres/coordinates_nz.pdf', width = 7, height = 5.5)

p2_zr = ggplot() +
  geom_point(data = Hzr, aes(ilr1, ilr2)) +
  geom_point(data = Hzr[n0], aes(ilr1, ilr2), col = 'red') +
  geom_point(data = Hzr[n1], aes(ilr1, ilr2), col = 'orange') +
  geom_point(data = Hzr[n2], aes(ilr1, ilr2), col = 'yellow') +
  labs(title = 'Coordinates', subtitle = ZR_LAB) +
  coord_fixed(xlim = XLIM, ylim = YLIM)


p2_zr = ggplot() +
  geom_point(data = Hzr, aes(ilr1, ilr2)) +
  labs(title = 'Coordinates', subtitle = ZR_LAB) +
  coord_fixed(xlim = XLIM, ylim = YLIM)
#####
Hnz = coordinates(Ynz, B2, label = 'ilr')
p2_nz = ggplot() +
  geom_point(data = Hnz, aes(ilr1, ilr2)) +
  labs(title = 'Coordinates', subtitle = NZ_LAB) +
  coord_fixed(xlim = XLIM, ylim = YLIM)
p2_nz

ggsave(p2_zr%>% COORD_LABS(), file = 'IFCS2019/pres/coordinates_black_zr.pdf', width = 7, height = 5.5)
ggsave(p2_nz%>% COORD_LABS(), file = 'IFCS2019/pres/coordinates_black_nz.pdf', width = 7, height = 5.5)

Czr0 = kmeansruns(Hzr, krange = 1:10)$cluster
Cnz0 = kmeansruns(Hnz, krange = 1:10)$cluster

p1_zr = ggplot() +
  geom_point(data = Hzr, aes(ilr1, ilr2, col = factor(Czr0)), alpha=0.8)
p1_zr = p1_zr +
  coord_fixed(xlim = XLIM, ylim = YLIM)
p1_zr = p1_zr + labs(title = 'Coordinates', subtitle = ZR_LAB)

#############
p1_nz = ggplot() +
  geom_point(data = Hnz, aes(ilr1, ilr2, col = factor(Cnz0)), alpha=0.8)
p1_nz = p1_nz +
  coord_fixed(xlim = XLIM, ylim = YLIM)
p1_nz = p1_nz + labs(title = 'Coordinates', subtitle = NZ_LAB)

ggsave(p1_zr%>% COORD_LABS(), file = 'IFCS2019/pres/clustering0_zr.pdf', width = 7, height = 5.5)
ggsave(p1_nz%>% COORD_LABS(), file = 'IFCS2019/pres/clustering0_nz.pdf', width = 7, height = 5.5)

p1_zr = ggplot() +
  geom_point(data = Hzr, aes(ilr1, ilr2), alpha=0.2)
for(l in Lzr){
  p1_zr = p1_zr +
    geom_path(data=l, aes(x,y))
}
p1_zr = p1_zr +
  coord_fixed(xlim = XLIM, ylim = YLIM)
p1_zr = p1_zr + labs(title = 'Coordinates', subtitle = ZR_LAB)
#############
p1_nz = ggplot() +
  geom_point(data = Hnz, aes(ilr1, ilr2), alpha=0.2)
for(l in Lnz){
  p1_nz = p1_nz +
    geom_path(data=l, aes(x,y))
}
p1_nz = p1_nz +
  coord_fixed(xlim = XLIM, ylim = YLIM)
p1_nz = p1_nz + labs(title = 'Coordinates', subtitle = NZ_LAB)

ggsave(p1_zr%>% COORD_LABS(), file = 'IFCS2019/pres/model_zr.pdf', width = 7, height = 5.5)
ggsave(p1_nz%>% COORD_LABS(), file = 'IFCS2019/pres/model_nz.pdf', width = 7, height = 5.5)


print_posterior = function(I, LO = 100, NL = 10){
  set.seed(1)
  #########
  PI_zr = Mzr$parameters$pro
  MU_zr = Mzr$parameters$mean
  SIGMA_zr = Mzr$parameters$variance$sigma
  H_zr_s = c_rlrnm_mixture_posterior(1000, Y[I,], PI_zr, MU_zr, SIGMA_zr, B2, 1000) %>% as.data.frame()
  l_post_zr = data_posterior(xseq = seq(min(H_zr_s[,1]), max(H_zr_s[,1]), length.out = LO),
                             yseq = seq(min(H_zr_s[,2]), max(H_zr_s[,2]), length.out = LO),
                             Y[I,], PI_zr, MU_zr, SIGMA_zr, B2, NL) %>%
    lapply(as.data.table)
  ########
  PI_nz = Mnz$parameters$pro
  MU_nz = Mnz$parameters$mean
  SIGMA_nz = Mnz$parameters$variance$sigma
  H_nz_s = c_rlrnm_mixture_posterior(1000, Y[I,], PI_nz, MU_nz, SIGMA_nz, B2, 1000) %>% as.data.frame()
  l_post_nz = data_posterior(xseq = seq(min(H_nz_s[,1]), max(H_nz_s[,1]), length.out = LO),
                             yseq = seq(min(H_nz_s[,2]), max(H_nz_s[,2]), length.out = LO),
                             Y[I,], PI_nz, MU_nz, SIGMA_nz, B2, NL) %>%
    lapply(as.data.table)
  p2_zr = p1_zr
  for(post in l_post_zr){
    p2_zr = p2_zr + geom_path(data = post,  aes(x = x, y = y), col = 'blue', alpha=0.4)
  }
  p2_zr = p2_zr +
    geom_point(data = Hzr[I,], aes(x = ilr1, y = ilr2), col = 'red', size = 2)
  # p2_zr = p2_zr +
  #   geom_point(data=H_zr_s, aes(x=V1,y=V2))
  p2_zr = p2_zr + labs(title = 'Coordinates', subtitle = ZR_LAB)
  
  p2_nz = p1_nz
  for(post in l_post_nz){
    p2_nz = p2_nz + geom_path(data = post,  aes(x = x, y = y), col = 'blue', alpha=0.4)
  }
  p2_nz = p2_nz +
    geom_point(data = Hnz[I,], aes(x = ilr1, y = ilr2), col = 'red', size = 2)
  p2_nz = p2_nz + labs(title = 'Coordinates', subtitle = NZ_LAB)
  list('zr' = p2_zr, 'nz' = p2_nz)
}


# LP = parallel::mclapply(1:947, print_posterior, mc.cores = 8)
# save(LP, file = 'IFCS2019/pres/posterioris.RData')
# load('IFCS2019/pres/posterioris.RData')


I = 60
P = print_posterior(I)
ggsave(P$zr%>% COORD_LABS(), file = sprintf('IFCS2019/pres/posterior_zr_%d.pdf', I), width = 7, height = 5.5)
ggsave(P$nz%>% COORD_LABS(), file = sprintf('IFCS2019/pres/posterior_nz_%d.pdf', I), width = 7, height = 5.5)

I = 331
P = print_posterior(I)
P$zr
ggsave(P$zr%>% COORD_LABS(), file = sprintf('IFCS2019/pres/posterior_zr_%d.pdf', I), width = 7, height = 5.5)
ggsave(P$nz%>% COORD_LABS(), file = sprintf('IFCS2019/pres/posterior_nz_%d.pdf', I), width = 7, height = 5.5)

I = 90
P = print_posterior(I)
ggsave(P$zr%>% COORD_LABS(), file = sprintf('IFCS2019/pres/posterior_zr_%d.pdf', I), width = 7, height = 5.5)
ggsave(P$nz%>% COORD_LABS(), file = sprintf('IFCS2019/pres/posterior_nz_%d.pdf', I), width = 7, height = 5.5)

print_sample = function(SEED){
  set.seed(SEED)
  Hs_zr = apply(Y, 1, function(x_) c_rlrnm_mixture_posterior(1, x_, PI_zr, MU_zr, SIGMA_zr, t(MASS::ginv(B2)), r = 100)) %>% 
    t() %>% as.data.table()
  Hs_nz = apply(Y, 1, function(x_) c_rlrnm_mixture_posterior(1, x_, PI_nz, MU_nz, SIGMA_nz, t(MASS::ginv(B2)), r = 100)) %>% 
    t() %>% as.data.table()
  Hs_zr$cl = kmeansruns(Hs_zr, krange = 1:5)$cluster
  Hs_nz$cl = kmeansruns(Hs_nz, krange = 1:5)$cluster
  p3_zr = ggplot() +
    geom_point(data = Hs_zr, aes(V1, V2), alpha=0.2) +
    annotate(geom='text', x = 3.5, y = -0.1, label = sprintf("Seed %d", SEED), size = 3)
  for(l in Lzr){
    p3_zr = p3_zr +
      geom_path(data=l, aes(x,y), alpha= 0.2)
  }
  p3_zr = p3_zr +
    coord_fixed(xlim = XLIM, ylim = YLIM)
  p3_zr = p3_zr + labs(title = 'Coordinates', subtitle = ZR_LAB)
  
  p3_zr_cl = ggplot() +
    geom_point(data = Hs_zr, aes(V1, V2,col=factor(cl)), alpha=0.7) +
    annotate(geom='text', x = 3.5, y = -0.1, label = sprintf("Seed %d", SEED), size = 3)
  for(l in Lzr){
    p3_zr_cl = p3_zr_cl +
      geom_path(data=l, aes(x,y), alpha= 0.2)
  }
  p3_zr_cl = p3_zr_cl +
    coord_fixed(xlim = XLIM, ylim = YLIM) + 
    labs(title = 'Coordinates', subtitle = ZR_LAB)
  
  p3_zr = p3_zr + labs(title = 'Coordinates', subtitle = ZR_LAB)
  p3_nz = ggplot() +
    geom_point(data = Hs_nz, aes(V1, V2), alpha=0.2) +
    annotate(geom='text', x = 3.5, y = -0.1, label = sprintf("Seed %d", SEED), size = 3)
  for(l in Lnz){
    p3_nz = p3_nz +
      geom_path(data=l, aes(x,y), alpha=0.2)
  }
  p3_nz = p3_nz +
    coord_fixed(xlim = XLIM, ylim = YLIM)
  p3_nz = p3_nz + labs(title = 'Coordinates', subtitle = NZ_LAB)
  p3_nz_cl = ggplot() +
    geom_point(data = Hs_nz, aes(V1, V2,col=factor(cl)), alpha=0.7) +
    annotate(geom='text', x = 3.5, y = -0.1, label = sprintf("Seed %d", SEED), size = 3)
  for(l in Lnz){
    p3_nz_cl = p3_nz_cl +
      geom_path(data=l, aes(x,y), alpha=0.2)
  }
  p3_nz_cl = p3_nz_cl +
    coord_fixed(xlim = XLIM, ylim = YLIM) + 
    labs(title = 'Coordinates', subtitle = NZ_LAB)
  
  list('zr' = p3_zr, 'nz' = p3_nz, 'zr_cl' = p3_zr_cl, 'nz_cl' = p3_nz_cl)
}

for(SEED in 1:5){
  P = print_sample(SEED)
  ggsave(P$zr%>% COORD_LABS(), file = sprintf('IFCS2019/pres/sample_zr_%d.pdf', SEED), width = 7, height = 5.5)
  ggsave(P$nz%>% COORD_LABS(), file = sprintf('IFCS2019/pres/sample_nz_%d.pdf', SEED), width = 7, height = 5.5)
  ggsave(P$zr_cl%>% COORD_LABS(), file = sprintf('IFCS2019/pres/sample_cl_zr_%d.pdf', SEED), width = 7, height = 5.5)
  ggsave(P$nz_cl%>% COORD_LABS(), file = sprintf('IFCS2019/pres/sample_cl_nz_%d.pdf', SEED), width = 7, height = 5.5)
}

C_mv_nz = 3 - C_mv_nz

p1_zr = ggplot() +
  geom_point(data = Hzr, aes(ilr1, ilr2, col = factor(C_mv_zr)), alpha=0.8)
p1_zr = p1_zr +
  coord_fixed(xlim = XLIM, ylim = YLIM)
p1_zr = p1_zr + labs(title = 'Coordinates', subtitle = ZR_LAB)

#############
p1_nz = ggplot() +
  geom_point(data = Hnz, aes(ilr1, ilr2, col = factor(C_mv_nz)), alpha=0.8)
p1_nz = p1_nz +
  coord_fixed(xlim = XLIM, ylim = YLIM)
p1_nz = p1_nz + labs(title = 'Coordinates', subtitle = NZ_LAB)

ggsave(p1_zr%>% COORD_LABS(), file = 'IFCS2019/pres/clustering_zr.pdf', width = 7, height = 5.5)
ggsave(p1_nz%>% COORD_LABS(), file = 'IFCS2019/pres/clustering_nz.pdf', width = 7, height = 5.5)

ind_different = which(C_mv_nz != C_mv_zr)

p1_zr = ggplot() +
  geom_point(data = Hzr, aes(ilr1, ilr2, col = factor(C_mv_zr)), alpha=0.8) +
  geom_point(data = Hzr[ind_different,], aes(ilr1, ilr2), color = 'black', alpha=0.8)
p1_zr = p1_zr +
  coord_fixed(xlim = XLIM, ylim = YLIM)
p1_zr = p1_zr + labs(title = 'Coordinates', subtitle = ZR_LAB)

#############
p1_nz = ggplot() +
  geom_point(data = Hnz, aes(ilr1, ilr2, col = factor(C_mv_nz)), alpha=0.8) +
  geom_point(data = Hnz[ind_different,], aes(ilr1, ilr2), color = 'black', alpha=0.8)
p1_nz = p1_nz +
  coord_fixed(xlim = XLIM, ylim = YLIM)
p1_nz = p1_nz + labs(title = 'Coordinates', subtitle = NZ_LAB)

ggsave(p1_zr%>% COORD_LABS(), file = 'IFCS2019/pres/clustering_diff_zr.pdf', width = 7, height = 5.5)
ggsave(p1_nz%>% COORD_LABS(), file = 'IFCS2019/pres/clustering_diff_nz.pdf', width = 7, height = 5.5)

pf_zr = ggtern() +
  geom_mask() +
  geom_point(data = data.table(Y), aes(x=other, y=ind, z=esp, col = factor(C_mv_zr)), alpha=0.8) +
  labs(title = ZR_LAB)
pf_nz = ggtern() +
  geom_mask() +
  geom_point(data = data.table(Y), aes(x=other, y=ind, z=esp, col = factor(C_mv_nz)), alpha=0.8) +
  labs(title = NZ_LAB)

ggsave(pf_zr, file = 'IFCS2019/pres/clustering_tern_zr.pdf', width = 7, height = 5.5)
ggsave(pf_nz, file = 'IFCS2019/pres/clustering_tern_nz.pdf', width = 7, height = 5.5)

pf_zr = ggtern() +
  geom_mask() +
  geom_point(data = data.table(Y), aes(x=other, y=ind, z=esp, col = factor(C_mv_zr)), alpha=0.8) +
  geom_point(data = data.table(Y)[ind_different], aes(x=other, y=ind, z=esp), color = 'black', alpha=0.8) +
  labs(title = ZR_LAB)
pf_nz = ggtern() +
  geom_mask() +
  geom_point(data = data.table(Y), aes(x=other, y=ind, z=esp, col = factor(C_mv_nz)), alpha=0.8) +
  geom_point(data = data.table(Y)[ind_different], aes(x=other, y=ind, z=esp), color='black', alpha=0.8) +
  labs(title = NZ_LAB)

ggsave(pf_zr, file = 'IFCS2019/pres/clustering_tern_diff_zr.pdf', width = 7, height = 5.5)
ggsave(pf_nz, file = 'IFCS2019/pres/clustering_tern_diff_nz.pdf', width = 7, height = 5.5)

# 
# which(C_mv_nz != C_mv_zr)
# 
# I = 24
# P24 = print_posterior(I)
# P24$zr
# P24$nz
# 
# I = 516
# P516 = print_posterior(I)
# P516$zr
# P516$nz
# 
# I = 563
# P563 = print_posterior(I)
# P563$zr
# P563$nz
# 
# I = 696
# P696 = print_posterior(I)
# P696$zr
# P696$nz
# 
# I = 725
# P725 = print_posterior(I)
# P725$zr
# P725$nz
# 

# 
