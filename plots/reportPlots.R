library(tidyverse)
library(lmtest)
library(Rmisc)
library(R.matlab)
library(tseries)
setwd("/Users/galraz1/Developer/CoCoSci/Analysis scripts")
df1 <- data.frame(read_csv('Rdata.csv'))

# lineplot randomXbins vs disengaged
ggplot(df1, aes(binMeans, propDisengaged)) + geom_line()+
  geom_point() + geom_errorbar(data=df1, y=stdErrorDisengaged)

df2 <- data.frame(read_csv('Rdata2.csv'))

df2$quad <- df2$diffs^2

mylogit <- glm(disengaged ~ eventpos, data = df2, family = "binomial")
mylogit2 <- glm(disengaged ~ eventpos + CT, data = df2, family = "binomial")
mylogit3 <- glm(disengaged ~ diffs, data = df2, family = "binomial")

lrtest(mylogit, mylogit2)
lrtest(mylogit2, mylogit3)


tgc <- summarySE(df2, measurevar="disengaged", groupvars=c('binMeans'))

# dotplot for logistic regression
ggplot(tgc, aes(binMeans, disengaged))  + geom_point(size=3, shape=21, fill = 'black') +
  geom_errorbar(aes(ymin=disengaged-se, ymax=disengaged+se), size = 1, width=.01) +
  geom_line(color = 'black', size = 1) + 
  theme_bw() +
  ylab('P(disengage)') + xlab(expression(Delta~random(x))) +
  theme(panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust=0.5, size=26, face="bold"),
        axis.title.x = element_text(hjust=0.53, size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16)) 

ggplot(df2, aes(eventpos, diffs)) + geom_point() + 
  ylab(expression(Delta~random(x))) + xlab('position in sequence') +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust=0.5, size=26, face="bold"),
        axis.title.x = element_text(hjust=0.53, size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16)) 

ggplot(df2,  aes(eventpos, disengaged)) + geom_point(position='jitter') + 
  xlab('eventpos') + ylab('P(disengage)') +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust=0.5, size=26, face="bold"),
        axis.title.x = element_text(hjust=0.53, size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16))



# boxplot change in randomX for disengaged vs engaged 
ggplot(df2, aes(x=as.factor(disengaged), y=diffs)) + 
  geom_boxplot() + geom_point()


curves <- data.frame(read_csv('curveRData.csv'))
diffs <- data.frame(read_csv('diffRData.csv'))
sequences <- data.frame(read_csv('sequenceRData.csv'))

xlabel1 = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'

ggplot(curves, aes(x=1:29, y=curve1)) + geom_line(size = 1.4) +
 xlab('event #') + ylab(expression(random(X)/random(X[max]))) +
   ggtitle('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA') + ylim(0, 1) + theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust=0.5, size=26, face="bold"),
        axis.title.x = element_text(hjust=0.53, size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16)) 

ggplot(curves, aes(x=1:29, y=curve2)) + geom_line(size = 1.4) +
  xlab('event #') + ylab(expression(random(X)/random(X[max]))) +
  ggtitle('CBCCCAAABBABCBBACAABCBCCABACBC') + ylim(0, 1) + theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust=0.5, size=26, face="bold"),
        axis.title.x = element_text(hjust=0.53, size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16)) 


ggplot(curves, aes(x=1:29, y=curve3)) + geom_line(size = 1.4) +
  xlab('event #') + ylab(expression(random(X)/random(X[max]))) +
  ggtitle('CABCABCABCABCABCABCABCABCABCAB') +   ylim(0, 1) +
 theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust=0.5, size=26, face="bold"),
        axis.title.x = element_text(hjust=0.53, size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16)) 


avg <- as.numeric(unlist(readMat('modelHumanCorr_avg.mat')))

each <- as.numeric(unlist(readMat('modelHumanCorr_each.mat')))
each[each < 0.1] <- 0.1

# Dummy data
delta <- seq(0.25, 0.95, by=0.1)
alpha <- seq(0.05, 0.95, by=0.1)
data <- expand.grid(X=alpha, Y=delta)
data$rho <- avg

# Heatmap 
ggplot(data, aes(X, Y, fill= rho)) + 
  geom_tile() + scale_fill_gradient2(low="white", mid="blue", high="red") +
  xlim(0.05, 0.95) + ylim(0.25, 0.95) + xlab(expression(alpha)) +
  ylab(expression(delta)) + 
  theme_classic() + 
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        plot.background = element_blank(),
        plot.title = element_text(hjust=0.5, size=26, face="bold"),
        axis.title.x = element_text(vjust=6, size=20),
        axis.title.y = element_text(hjust=0.5, vjust =-5, size=20),
        axis.text.x = element_text(vjust=10, size=14, family = 'sans'),
        axis.text.y = element_text(margin = margin(l = 30), hjust = 0, size =14, family = 'sans'), 
        axis.line = element_blank(),
        axis.ticks = element_blank())


# plot for each
data <- expand.grid(X=alpha, Y=delta)
data$rho <- each

ggplot(data, aes(X, Y, fill= rho)) + 
  geom_tile() + scale_fill_gradient2(low="black", mid="blue", high="white",
                                       midpoint = 0.5) +
  xlim(0.05, 0.95) + ylim(0.25, 0.95) + xlab(expression(alpha)) +
  ylab(expression(delta)) + 
  theme_classic() + 
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        plot.background = element_blank(),
        plot.title = element_text(hjust=0.5, size=26, face="bold"),
        axis.title.x = element_text(vjust=4, size=20),
        axis.title.y = element_text(hjust=0.5, vjust =-5, size=20),
        axis.text.x = element_text(vjust=12, size=10, family = 'sans'),
        axis.text.y = element_text(margin = margin(l = 45), hjust = 0, size =12, family = 'sans'), 
        axis.line = element_blank(),
        axis.ticks = element_blank())

