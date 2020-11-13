# install.packages("DirichletReg")
library(DirichletReg)
head(ArcticLake)
AL <- DR_data(ArcticLake[, 1:3])

plot(AL, cex = 0.5, a2d = list(colored = FALSE, c.grid = FALSE))
plot(rep(ArcticLake$depth, 3), as.numeric(AL), pch = 21
     , bg = rep(c("#E495A5","#86B875", "#7DB0DD"), each = 39),
     xlab = "Depth (m)", ylab = "Proportion", ylim = 0:1)

lake1 <- DirichReg(AL ~ depth, ArcticLake)
lake1




info = readRDS("data/LC_sample_info.rds")
precentLC <- DR_data(info[, 2:4]/100)
precentLC[,1] <- info[, 2]/100
precentLC[,2] <- info[, 3]/100
precentLC[,3] <- info[, 4]/100

precentLC
plot(precentLC, cex = 2, a2d = list(colored = T, c.grid = T))
plot(precentLC, cex = 0.5, a2d = list(colored = T, c.grid = F))
pre

ggplot(data= c.data) +
  geom_boxplot(aes(y= debridement.day), fill= "#CC9966", lwd=1)
  
  boxplot(x= debridement.day, data= c.data)
  
  
  
