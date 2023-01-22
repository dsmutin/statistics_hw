setwd ("~statistics")

library(tidyverse)
library(ggthemes)

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

anova (m1, m0)
AIC(m1,m0)

#3----------
read.csv()

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
