library(topicmodels)
library(tm)
library(tidytext)
library(data.table)
load('IFCS2019/parliament_example_fitting.RData')

df = melt(Lparl[['2017']], id.vars = 'mun')
PARLAMENT2017 = cast_dtm(df, mun, variable, value)

# res = CRM(PARLAMENT2017, k = 3, control = list(seed = 1, verbose = 1))
# #save.image('CTM_output.RData')

res = LDA(PARLAMENT2017, k = 3, control = list(seed = 1, verbose = 1))

# TOP = exp(res@beta)
# colnames(TOP) = names(Lparl[['2017']][,-1])
# TOP
# d = as.data.frame(res@gamma)
# head(d)
# library(ggtern)
# library(coda.base)
# ds = as.data.frame(composition(scale(coordinates(d))))
# G = apply(d, 1, which.max)
# ggtern() +
#   geom_mask() +
#   geom_point(data = ds, aes(x=x1,y=x2,z=x3,col=factor(G))) +
#   theme_minimal()

#save.image('CTM_output.RData')



