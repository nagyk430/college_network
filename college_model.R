
library(readxl)
szocio <- read_excel("")

hist(szocio$index_ertek)
hist(szocio$dpont)
mean(szocio$index_ertek)
median(szocio$index_ertek)


#Data manipulation
szocio$nem_dummy <- ifelse(szocio$nem == "Férfi", 0, 1) #0 = férfi, 1=nő
szocio$index_bin <- ifelse(szocio$index_ertek >= 0.16, 1,0) #1 >= 0,16; 0 <- <0,16
szocio$dpont_bin <- ifelse(szocio$dpont >= 40, 1,0) #1 >= 40, 0 <- <40

szocio$szint_faktor <- as.factor(szocio$szint)
szocio$nem_faktor <- as.factor(szocio$nem)
szocio$ido_faktor <- as.factor(szocio$ido)

szocio$dpont_kat <- ifelse(szocio$dpont < 36, 1,
                           ifelse(szocio$dpont >= 36 & szocio$dpont <= 54, 2,
                                  ifelse(szocio$dpont >= 55 & szocio$dpont <= 75, 3, 4)))

hist(szocio$dpont_kat)
szocio2 <- na.omit(szocio)

#Models

summary(mod1 <- lm(index_ertek ~ ido_faktor, data=szocio))

summary(mod2 <- lm(index_ertek ~ ido_faktor + nem_faktor + szint_faktor, data=szocio))

summary(mod3 <- lm(index_ertek ~ dpont_kat, data=szocio))

#summary(mod_poisson1 <- glm(index_bin ~ bentlakas + nem_dummy, data=szocio, family=poisson))

summary(binom1 <- glm(index_bin ~ szoc_saj, data=szocio2, family=binomial))

summary(mod_poisson3 <- glm(index_bin ~ nem_dummy, data=szocio, family=binomial))

summary(mod_poisson4 <- glm(dpont_kat ~ nem_faktor + ido_faktor + szint_faktor, data=szocio, family=poisson))

summary(mod_poisson5 <- glm(index_bin ~ dpont_bin, data=szocio, family=binomial))

summary(lm(index_ertek ~ dpont + szint_faktor, data=szocio))


#Participation in programmes
summary(glm(index_bin ~ ismerk + 
              teahaz + 
              kolihet + 
              csendes + 
              Ifi1 +
              Ifi2 + 
              zold +
              angy + 
              méz + 
              diszn, data = szocio, family = binomial))

szocio$lelki <- ifelse(szocio$Ifi1 == 1 | szocio$Ifi2 == 1 | szocio$csendes == 1, 1, 0) #1 = részt vett, 0=nem vett részt
szocio$hetv <- ifelse(szocio$kolihet == 1 | szocio$zold == 1 | szocio$diszn == 1, 1, 0) #1 = részt vett, 0=nem vett részt
szocio$hetk <- ifelse(szocio$ismerk == 1 | szocio$teahaz == 1 | szocio$angy == 1 | szocio$méz == 1, 1, 0) #1 = részt vett, 0=nem vett részt


summary(glm(index_bin ~ lelki + hetv + hetk, data=szocio, family=poisson))
summary(lm(lelki ~ nem_faktor + szint_faktor + ido_faktor, data=szocio))
summary(lm(hetk ~ nem_faktor + szint_faktor + ido_faktor, data=szocio))
summary(lm(hetv ~ nem_faktor + szint_faktor + ido_faktor, data=szocio))

-----------------------------------------------------------------------------
  #####Cluster analysis
  install.packages("mvtnorm")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("fpc") 
install.packages("caret")
install.packages("mclust")
install.packages("factoextra")
install.packages("cluster")

library(readxl)
szocio <- read_excel("")

szocio$nem_dummy <- ifelse(szocio$nem == "Férfi", 0, 1) #0 = férfi, 1=nő


#After data cleaning
szocio <- na.omit(szocio)
szocio$ido <- NULL
szocio$szint <- NULL
szocio$dpont <- NULL
szocio$index_ertek <- NULL
szocio$nem <- NULL
szocio$nev <- NULL
szocio$szoc_saj <- NULL
szocio$bentlakas <- NULL


#Plot of the dendogram
d <- dist(szocio)
h <- hclust(d)
h

plot(h)

#Number of cluster
wss <- (nrow(szocio)-1)*sum(apply(szocio,2,var))
for (i in 2:15) {
  wss[i] <- sum(kmeans(szocio, centers=i)$withinss) # sum of square alapján
}
plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")

d <- dist(szocio)
h <- hclust(d)
h
plot(h)
rect.hclust(h, k=2, border= "red")

#install.packages("NbClust")
library(NbClust)
NbClust(szocio, method = 'complete', index = 'dindex')

NbClust(szocio, method = 'complete', index = 'hartigan')$Best.nc

szocio2 <- scale(szocio)
summary(szocio2)

NbClust(szocio2, method = 'complete', index = 'dindex')
NbClust(szocio2, method = 'complete', index = 'hartigan')$Best.nc #creating 4 cluster

# Hierarchical:
d <- dist(szocio,  diag = TRUE)
h <- hclust(d)
h

par(mfrow=c(1,1))

plot(h)

#Plotting 5 cluster
rect.hclust(h, k=5, border = "red")
(cn <- cutree(h, k = 5))

table(cn)

#Describe the clusters
round(aggregate(szocio, FUN = mean, by = list(cn)), 1) #mean values of the group
round(aggregate(szocio, FUN = sd, by = list(cn)), 1) #scatter values of the group
round(sapply(szocio, sd), 1)
round(apply(aggregate(szocio, FUN = mean, by = list(cn)), 2, sd), 1) #scatter of groups relative to each other

summary(aov(szocio$nem_dummy ~ cn))
sort(unlist(apply(szocio, 2, function(x) summary(aov(x ~ cn))[[1]][4][1,1])))

#K-means cluster 4
#fit <- kmeans(szocio, 4)
#aggregate(szocio,by=list(fit$cluster),FUN=mean)

#K-means cluster 5
fit <- kmeans(szocio, 5)
aggregate(szocio,by=list(fit$cluster),FUN=mean)

szocio <- data.frame(szocio, fit$cluster)

hist(szocio$fit.cluster)

#Ward Hierarchial Clustering
#install.packages("pvclust")
library(pvclust)
fit <- pvclust(t(szocio), method.hclust="ward.D2", method.dist="euclidean")
plot(fit)

pvrect(fit, alpha=.95)
fit <- kmeans(szocio, 5)

#Plotting
library(cluster)
clusplot(szocio, fit$cluster, color=TRUE, shade=TRUE, labels=2, lines=0)

fit1 <- kmeans(szocio, 5)
table(fit1$cluster, cn)

