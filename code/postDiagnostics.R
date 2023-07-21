setwd("~/mcmc_handbook/")
library(coda)
library(readr)
library(ggplot2)
library(colorspace)
library(ggpubr)


df <- read_table2("output/withHmc.log", skip = 3)
df <- df[df$state>=1000,12:1903]
effs <- effectiveSize(df)

df2 <- read_table2("output/noHmc2.log", skip = 3)
df2 <- df2[df2$state>=1000,12:1903]
effs2 <- effectiveSize(df2)

df3 <- data.frame(Values=effs,Algorithm="HMC")
df4 <- data.frame(Values=effs2,Algorithm="RW")

df <- rbind(df3,df4)

# copy paste hours from timing files
hmcHours <- 19.033439722222223
noHmcHours <- 20.454540833333333

df$Values[df$Algorithm=="HMC"] <- df$Values[df$Algorithm=="HMC"] / 19.033439722222223
df$Values[df$Algorithm=="RW"] <- df$Values[df$Algorithm=="RW"] / 20.454540833333333

df$Algorithm[df$Algorithm=="HMC"] <- "Surrogate-trajectory\nHMC"
df$Algorithm[df$Algorithm=="RW"] <- "Random walk"


gg <- ggplot(df,aes(x=Values, fill=Algorithm)) +
  geom_histogram() + xlab(NULL) +
  scale_x_log10() +
  ggtitle("Effective sample size per hour") +
  ylab("Number of parameters") +
  scale_fill_discrete_qualitative(palette = "Harmonic") +
  theme_bw() 
gg


ggsave("figures/ctmcPerformance.pdf",gg,width = 6,height=2.4)
system2(command = "pdfcrop",
        args    = c("~/mcmc_handbook/figures/ctmcPerformance.pdf",
                    "~/mcmc_handbook/figures/ctmcPerformance.pdf")
)


