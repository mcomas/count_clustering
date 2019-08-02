library(coda.count)
library(data.table)

PREBUILD = FALSE
if(PREBUILD){
  load('IFCS2019/catalan_elections.RData')
  parliament_elections = data.table(parliament_elections)
  YEAR = 2017
  years_ = unique(parliament_elections$year)
  Lparl = lapply(years_, function(YEAR){
    dcast(parliament_elections[year == YEAR], mun~party, value.var = 'votes')
  })
  names(Lparl) = years_
  
  data = Lparl[['2017']]
  X = as.matrix(data[,-1])
  B = ilr_basis(ncol(X))
  fit = fit_lrnm(X, B, method = 'hermite', hermite.order = 3, max_iter = 1000)
  save.image(file = 'IFCS2019/parliament_example_fitting.RData')
}else{
  load('IFCS2019/parliament_example_fitting.RData')
}

fit$mu
fit$sigma
