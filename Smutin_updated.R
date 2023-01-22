setwd ("~/Загрузки/Telegram Desktop/statistics")

library(tidyverse)
library(ggthemes)
library(GGally)
library(Hmisc)
library(ggpubr)

#1 ----------

inb <- read.csv("./1/InbreedingWolves.csv", sep = ";")
names(inb) <- c('i','p')

#1_1.png
ggplot(data = inb, mapping = aes (x = as.character(p), y = i)) +
  geom_violin(draw_quantiles = T,
              trim = F,
              color = "cadetblue4") +
  geom_jitter(width = 0.1, color = "cadetblue4") +
  geom_smooth(mapping = aes (x = p, y = i), formula = y~x, level= 0.99, method = glm) +
  xlab("pups") +
  ylab("inbreeding coefficient") +
  theme_light()

#1_2.png
ggplot(data = inb) +
  geom_density(mapping = aes (x = i)) +
  xlab("inbreeding coefficient") +
  ylab("density") +
  theme_light()

#1_3.png
ggplot(data = inb) +
  geom_density(mapping = aes (x = p)) +
  xlab("pups") +
  ylab("density") +
  theme_light()

#correlation
cor_inb <- cor.test(y = inb$i, x = inb$p, conf.level = 0.99, method = "spearman")

#2----------

empat <- read.csv("./2/MouseEmpathy.csv", sep = ",")

m1 <- glm(data = empat, percent.stretching ~ treatment)
summary(m1)
m0 <- glm(data = empat, percent.stretching ~ percent.stretching)


pred_m1 <- data.frame("P" = m1$linear.predictors, "Tr" = empat$treatment)
pred_m0 <- data.frame("P" = m0$linear.predictors, "Tr" = empat$treatment)

#2_1.png
ggplot() +
  geom_violin(data = empat, mapping = aes (x = treatment, y = percent.stretching, color = treatment),
              draw_quantiles = T,
              trim = F) +
  geom_jitter(data = empat, mapping = aes (x = treatment, y = percent.stretching, color = treatment),
              width = 0.1) +
  geom_point(data = pred_m1, mapping = aes (x = Tr, y = P), color = "cadetblue4", shape = "_", size = 20)+
  geom_point(data = pred_m0, mapping = aes (x = Tr, y = P), color = "cadetblue1", shape = "_", size = 20)+
  xlab("treatment") +
  ylab("percent stretching") +
  theme_light()

anova(m1, m0)
AIC(m1,m0)

#3----------
prostate <-as.matrix(read.table("./3/prostate.data", sep = "\t", header = T, row.names = 1))

#взглянем на данные: 3_1-10.png

i <- 1
for (i in c(1:10)){
  if (i == 5 | i == 7 | i == 10){
    a <- ggplot() +
      geom_bar(mapping = aes (x = as.character(prostate[,i])), color = "cadetblue4") +
      xlab(colnames(prostate[i])) +
      ylab("N") +
      theme_light()
  } else {
    a <- ggplot() +
      geom_density(mapping = aes (x = prostate[,i]), color = "cadetblue4") +
      xlab(colnames(prostate[i])) +
      ylab("density") +
      theme_light()
  }
  ggsave(plot=a, str_c("3_",i,".png"))
}

#проверим на нормальность
sw_pr <- unlist(map(apply(prostate, 2, shapiro.test), "p.value"))
sw_pr > 0.01

#3_11.png
cor_p <- rcorr(prostate, type = "spearman")
ggcorr(prostate, cor_matrix = cor_p[[1]], geom = "text", mid = "cadetblue1", 
       high = "cadetblue4", low = "cadetblue4")

#3_12.png
cor_p_005 <- cor_p[[3]] < 0.05
ggcorr(prostate, cor_matrix = cor_p_005, mid = "cadetblue1", 
       high = "cadetblue4", low = "cadetblue4")

##glms
pr <- as.data.frame(prostate)
m0 <- glm(data = pr, lcavol ~ lcavol)
m1 <- glm(data = pr, lcavol ~ lcp)
m2 <- glm(data = pr, lcavol ~ lpsa)
m3 <- glm(data = pr, lcavol ~ lcp*lpsa)
m4 <- glm(data = pr, lcavol ~ lcp*lpsa*svi)
m5 <- glm(data = pr, lcavol ~ lcp*lpsa*svi*gleason)
m6 <- glm(data = pr, lcavol ~ lcp*lpsa*svi*gleason*lweight)

ANOVA_pr <- anova(m0,m1,m2,m3,m4,m5,m6)
AIC_pr <- AIC(m0,m1,m2,m3,m4,m5,m6)
c_m4 <- coef(m4)

write.table(round(ANOVA_pr,2), "ANOVA_pr.txt", sep = "\t")
write.table(round(AIC_pr,2), "AIC_pr.txt", sep = "\t")
write.table(round(c_m4,3), "c_m4.txt", sep = "\t")

#4----------

rx <- unlist(as.matrix(read.table ("./4/radsens.x", header = FALSE, quote = c(" ", "  "))))
ry <- read.table("./4/radsens.y", header = FALSE, sep = " ", fill = T)
ry <- na.omit(unlist(c(ry[1,], ry[2,])))
norm <- ry == 1

m_norm <- apply(rx[,norm], 1, mean)
sd_norm <- apply(rx[,norm], 1, sd)
n_norm <- 44

m_rad <- apply(rx[,!norm], 1, mean)
sd_rad <- apply(rx[,!norm], 1, sd)
n_rad <- 14

rx_t <- ( m_norm - m_rad ) / 
  sqrt (sd_norm/n_norm + sd_rad/n_rad)

normal_t <- rnorm(12627, mean = mean(rx_t), sd = sd(rx_t))
#4_1.png
ggplot() +
  geom_density(mapping = aes(x = normal_t), color = "cadetblue1") +
  geom_density(mapping = aes(x = rx_t), color = "cadetblue4") +
  xlab("t-test value") +
  ylab("density") +
  theme_light()

#4_2.png
ggplot() +
  geom_qq(mapping = aes(sample = rx_t), color = "cadetblue4", alpha = 0.3, inherit.aes = T) +
  geom_qq_line(mapping = aes(sample = rx_t)) +
  xlab("theoretical") +
  ylab("sample") +
  theme_light()

t.test(rx_t,normal_t)

rx_t_p <- pt(rx_t, 57)
sum(rx_t_p < 0.01/12627)

rx_t_corrected <- p.adjust(rx_t_p, "holm")
sum(rx_t_corrected < 0.01)
