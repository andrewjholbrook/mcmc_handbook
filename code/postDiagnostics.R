setwd("~/expmDerivative/")
library(coda)

df <- read_table2("output/short.Ed.TOOBIG4GIT.log", skip = 3)
df <- df[df$state>=1000 & df$state <= 80000,]

# fixed effects
effectiveSize(df[,9:11]) #        1053.3909        1342.5106         498.3175 

# random effects
summary(effectiveSize(df[,12:1903]))


# full run
df <- read_table2("output/covid10Mar_285plusNYC.Ed.log", skip = 3)
df <- df[df$state>=800000,]

# global scale
effectiveSize(df[,1904])

# fixed effects
effectiveSize(df[,9:11])

# random effects
summary(effectiveSize(df[,12:1903]))

colMeans(df[,9:11])
coda::HPDinterval(coda::as.mcmc(df[,9:11]))

