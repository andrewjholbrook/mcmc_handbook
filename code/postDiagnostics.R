setwd("~/mcmc_handbook/")
library(coda)
library(readr)
library(ggplot2)
library(colorspace)
library(ggpubr)


df <- read_table2("output/withHmc.log", skip = 3)
df <- df[df$state>=1000,12:1903]
effs <- effectiveSize(df)

df2 <- read_table2("output/noHmc.log", skip = 3)
df2 <- df2[df2$state>=1000,12:1903]
effs2 <- effectiveSize(df2)

df3 <- data.frame(Values=effs,Algorithm="HMC")
df4 <- data.frame(Values=effs2,Algorithm="RW")

df <- rbind(df3,df4)



gg <- ggplot(df,aes(x=Values, fill=Algorithm)) +
  geom_histogram() + ggtitle("Median effective sample size") +
  #scale_fill_continuous_sequential(palette = "Heat2")+
  theme_bw()
gg




# ggsave("figures/ctmcTuning.pdf",ggarrange(gg,NULL,gg2,nrow=1,labels=c("a","","b"),
#                                           widths = c(1, -0.01, 1)),width = 10,height=4)
# system2(command = "pdfcrop",
#         args    = c("~/mcmc_handbook/figures/ctmcTuning.pdf",
#                     "~/mcmc_handbook/figures/ctmcTuning.pdf")
# )


