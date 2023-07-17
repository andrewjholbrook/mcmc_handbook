setwd("~/mcmc_handbook/")
library(coda)
library(readr)
library(ggplot2)
library(colorspace)
library(ggpubr)


medEssMat <- matrix(0,35,3)
minEssMat <- matrix(0,35,3)
rates     <- c(0.2,0.3,0.4,0.5,0.6,0.7,0.8)
steps     <- c(4,8,16,32,64)

r <- 1 # row
for(k in 1:7) {
  for(j in 1:5) {
    df <- read_table2(paste0("output/tuning",steps[j],rates[k],".log"), skip = 3)
    df <- df[df$state>=1000,9:1903]
    effs <- effectiveSize(df)
    medEss <- median(effs)
    minEss <- min(effs)
    medEssMat[r,] <- c(steps[j],rates[k],medEss)
    minEssMat[r,] <- c(steps[j],rates[k],minEss)
    r <- r + 1
  }
}
colnames(medEssMat) <- c("Steps","Target rate","Median")
colnames(minEssMat) <- c("Steps","Target rate","Minimum")

df <- data.frame(medEssMat)
df2 <- data.frame(minEssMat)
colnames(df) <- c("Steps","Target rate","Median")
colnames(df2) <- c("Steps","Target rate","Minimum")
df$Steps <- factor(df$Steps)
df$`Target rate` <- factor(df$`Target rate`)

df2$Steps <- factor(df2$Steps)
df2$`Target rate` <- factor(df2$`Target rate`)

gg <- ggplot(df,aes(x=`Target rate`, y=Steps, fill=Median)) +
  geom_tile(color = "white") + ylab("Leapfrog steps") + xlab("Target acceptance rate") +
  geom_text(aes(label=round(Median))) + ggtitle("Median effective sample size") +
  scale_fill_continuous_sequential(palette = "Heat2")+
  theme_bw() + theme(legend.position = "none")
gg

gg2 <- ggplot(df2,aes(x=`Target rate`, y=Steps, fill=Minimum)) +
  geom_tile() + ylab("") + xlab("Target acceptance rate") +
  geom_text(aes(label=round(Minimum))) + ggtitle("Minimum effective sample size") +
  scale_fill_continuous_sequential(palette = "Heat2")+
  theme_bw() + theme(legend.position = "none")
gg2


ggsave("figures/ctmcTuning.pdf",ggarrange(gg,NULL,gg2,nrow=1,labels=c("a","","b"),
                                          widths = c(1, -0.01, 1)),width = 10,height=4)
system2(command = "pdfcrop",
        args    = c("~/mcmc_handbook/figures/ctmcTuning.pdf",
                    "~/mcmc_handbook/figures/ctmcTuning.pdf")
)
