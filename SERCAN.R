library(tidyverse)
library(GGally)
library(caTools)
library(rstatix)
library(caret)
library(patchwork)
library(ggsignif)

# DATA TIDING ----
Data <- tibble(TKS4, TKS4b, TKS4L, TKS5L, T, N, subgroup, id)
write_csv(Data, 'Sercan.csv')
Data <- read_csv('Sercan.csv')

data <- Data |> na.omit()

R <- log2(data$TKS4L)-log2(data$TKS4b)

# BOX-COX TRANSFORM 
bc1 <- BoxCoxTrans(data$TKS4)
TKS4 <- predict(bc1, data$TKS4)
bc2 <- BoxCoxTrans(data$TKS4L)
TKS4L <- predict(bc2, data$TKS4L)
bc3 <- BoxCoxTrans(data$TKS4b)
TKS4b <- predict(bc3, data$TKS4b)
bc4 <- BoxCoxTrans(data$TKS5L)
TKS5L <- predict(bc4, data$TKS5L)

# LOG2 TRANSFORM 
X <- log2(data[, 1:4])
data <- cbind(X, R, data[,5], data[, 6], data[, 7], data[, 8]) |> tibble()

# DATA ANALYSIS ----
ggpairs(data[, 1:4], 
        title = 'Correlation and density plots on NCI data',
        xlab = 'Log2(ΔCt)',
        ylab = 'Log2(ΔCt), density')

# NORMALITY ----
shapiro.test(X$TKS4[1:21])   #normal
shapiro.test(X$TKS4[22:39])  #normal
shapiro.test(X$TKS4[40:58])  #normal
shapiro.test(X$TKS4[59:75])  #normal
shapiro.test(X$TKS4[76:86])  #normal
shapiro.test(X$TKS4[87:106]) #normal

shapiro.test(X$TKS4b[1:21])   #not normal
shapiro.test(X$TKS4b[22:39])  #normal
shapiro.test(X$TKS4b[40:58])  #normal
shapiro.test(X$TKS4b[59:75])  #normal
shapiro.test(X$TKS4b[76:86])  #normal
shapiro.test(X$TKS4b[87:106]) #normal

shapiro.test(X$TKS4L[1:21])   #normal
shapiro.test(X$TKS4L[22:39])  #normal
shapiro.test(X$TKS4L[40:58])  #normal
shapiro.test(X$TKS4L[59:75])  #normal
shapiro.test(X$TKS4L[76:86])  #normal
shapiro.test(X$TKS4L[87:106]) #not normal

shapiro.test(X$TKS5L[1:21])   #normal
shapiro.test(X$TKS5L[22:39])  #normal
shapiro.test(X$TKS5L[40:58])  #normal
shapiro.test(X$TKS5L[59:75])  #not normal
shapiro.test(X$TKS5L[76:86])  #normal
shapiro.test(X$TKS5L[87:106]) #normal

shapiro.test(R[1:21])   #normal
shapiro.test(R[22:39])  #normal
shapiro.test(R[40:58])  #not normal
shapiro.test(R[59:75])  #normal
shapiro.test(R[76:86])  #not normal
shapiro.test(R[87:106]) #normal

# VARIANCES ----
bartlett.test(TKS4 ~ subgroup, data)
bartlett.test(TKS4b ~ subgroup, data)
bartlett.test(TKS4L ~ subgroup, data)
bartlett.test(TKS5L ~ subgroup, data)
bartlett.test(R ~ subgroup, data)

# NON_PARAMETRIC -----
# variances
fligner.test(TKS4 ~ subgroup, data) 
fligner.test(TKS4b ~ subgroup, data)         
fligner.test(TKS4L ~ subgroup, data)  
fligner.test(TKS5L ~ subgroup, data)
fligner.test(R ~ subgroup, data)

#means
kruskal.test(R ~ subgroup, data) #differs
kruskal.test(TKS4 ~ subgroup, data) 
kruskal.test(TKS4b ~ subgroup, data)             
kruskal.test(TKS4L ~ subgroup, data) #differs
kruskal.test(TKS5 ~ subgroup, data)  #differs

data |> wilcox_test(R ~ subgroup) #differs
data |> wilcox_test(TKS4 ~ subgroup) 
data |> wilcox_test(TKS4b ~ subgroup)  
data |> wilcox_test(TKS4L ~ subgroup) #differs
data |> wilcox_test(TKS5 ~ subgroup)  #differs

# PARAMETRIC ----------
aov <- aov(TKS4 ~ subgroup, data)
summary(aov)
plot(aov)

aov <- aov(TKS4b ~ subgroup, data)
summary(aov)
plot(aov)

# DATA VISUALIZATION ----
rd <- ggplot(data, aes(R))+
  geom_density()+
  labs(title = 'Density plot of TKS4 long:short isoform ratio', x = '2^(TKS4L:TKS4b)')

r <- ggplot(data)+aes(subgroup, y = R)+
  geom_boxplot()+
  ylim(-3.2, 13)+
  geom_signif(data = data, y_position = c(8, 9, 10, 11, 12), 
              annotations = c('**', '**', '**', '*', '*'),
              comparisons = list(c('Adjacent', 'HER2-Enriched'),
                                 c('Adjacent', 'Lum_A'),
                                 c('Adjacent', 'LumB_Her-'),
                                 c('Adjacent', 'LumB_Her+'),
                                 c('Adjacent', 'Triple_neg')))+
  labs(title = 'Comparison of TKS4 long:short isoform ratio between subgroup', y = '2^(TKS4L:TKS4b)')

t4 <- ggplot(data)+aes(subgroup)+
  geom_boxplot(aes(y = TKS4))+
  labs(title = 'TKS4', x = 'subgroup', y = 'Log2(ΔCt)')
 
t4b <- ggplot(data)+aes(subgroup)+
  geom_boxplot(aes(y = TKS4b))+
  labs(title = 'TKS4b', x = '', y = '')

t4l <- ggplot(data)+aes(subgroup, y = TKS4L)+
  geom_boxplot()+
  ylim(-8, 4)+
  geom_signif(data = data, y_position = c(2.2, 3.2), annotations = c("*", '*'),
              comparisons = list(c('Adjacent', 'LumB_Her-'), 
                                 c('Adjacent', 'Lum_A')))+
  labs(title = 'TKS4L', x = '', y = 'Log2(ΔCt)')

t5 <- ggplot(data)+aes(subgroup, y = TKS5L)+
  geom_boxplot()+
  ylim(-9.5, 1.8)+
  geom_signif(data = data, y_position = c(1, 0, -1), annotations = c("*", '*', '*'),
              comparisons = list(c('Adjacent', 'LumB_Her-'), 
                                 c('Adjacent', 'HER2-Enriched'),
                                 c('HER2-Enriched', 'Triple_neg')))+
  labs(title = 'TKS5L', x = 'subgroup', y = '')

(t4l + t4b)/(t4 + t5)

# PCA + T-SNE ------------
pca <- prcomp(data[,1:5])
tsne <- Rtsne::Rtsne(data[,1:5])

pcadata <- cbind(data[,8], pca$x)
tsnedata <- data.frame(data[, 8],
                       x = tsne$Y[,1],
                       y = tsne$Y[,2])

pc <- ggplot(pcadata, aes(col = subgroup))+ 
  geom_point(aes(x=PC1, y = PC2))
ts <- ggplot(tsnedata, aes(col = subgroup))+ 
  geom_point(aes(x=x,y=y))

ts + pc + plot_layout(guides = 'collect')

