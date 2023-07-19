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
  geom_tile(color = "white") + ylab("Leapfrog steps") + xlab("") +
  geom_text(aes(label=round(Median))) + ggtitle("Median effective sample size") +
  scale_fill_continuous_sequential(palette = "Heat2")+
  theme_bw() + theme(legend.position = "none")
gg

gg2 <- ggplot(df2,aes(x=`Target rate`, y=Steps, fill=Minimum)) +
  geom_tile() + ylab("") + xlab("") +
  geom_text(aes(label=round(Minimum))) + ggtitle("Minimum effective sample size") +
  scale_fill_continuous_sequential(palette = "Heat2")+
  theme_bw() + theme(legend.position = "none")
gg2


# ggsave("figures/ctmcTuning.pdf",ggarrange(gg,NULL,gg2,nrow=1,labels=c("a","","b"),
#                                           widths = c(1, -0.01, 1)),width = 10,height=4)
# system2(command = "pdfcrop",
#         args    = c("~/mcmc_handbook/figures/ctmcTuning.pdf",
#                     "~/mcmc_handbook/figures/ctmcTuning.pdf")
# )

df_timing <- read_table("output/timing.txt",col_names = FALSE)
df3 <- df
colnames(df_timing) <- c("Steps","Hours")
df_timing$Hours <- df_timing$Hours * 100 / 60 # 60 minutes in an hour and 1000 its --> 100000 its

df3 <- merge(df3,df_timing)
df3$Median <- df3$Median / df3$Hours

df4 <- df2
df4 <- merge(df4,df_timing)
df4$Minimum <- df4$Minimum / df4$Hours

gg3 <- ggplot(df3,aes(x=`Target rate`, y=Steps, fill=Median)) +
  geom_tile(color = "white") + ylab("Leapfrog steps") + xlab("Target acceptance rate") +
  geom_text(aes(label=round(Median))) + ggtitle("Median effective sample size per hour") +
  scale_fill_continuous_sequential(palette = "Heat2")+
  theme_bw() + theme(legend.position = "none")
gg3

gg4 <- ggplot(df4,aes(x=`Target rate`, y=Steps, fill=Minimum)) +
  geom_tile() + ylab("") + xlab("Target acceptance rate") +
  geom_text(aes(label=round(Minimum))) + ggtitle("Minimum effective sample size per hour") +
  scale_fill_continuous_sequential(palette = "Heat2")+
  theme_bw() + theme(legend.position = "none")
gg4

ggsave("figures/ctmcTuning.pdf",ggarrange(gg,NULL,gg2,gg3,NULL,gg4,nrow=2,ncol=3,labels=c("a","","b","c","","d"),
                                          widths = c(1, -0.01, 1,1,-0.01,1)),width = 10,height=7)
system2(command = "pdfcrop",
        args    = c("~/mcmc_handbook/figures/ctmcTuning.pdf",
                    "~/mcmc_handbook/figures/ctmcTuning.pdf")
)
