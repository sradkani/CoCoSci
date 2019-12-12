library(tidyverse)
library(lmtest)
library(Rmisc)

setwd("/Users/galraz1/Developer/CoCoSci")
df1 <- data.frame(read_csv('Rdata.csv'))

# lineplot randomXbins vs disengaged
ggplot(df1, aes(binMeans, propDisengaged)) + geom_line()+
  geom_point() + geom_errorbar(data=df1, y=stdErrorDisengaged)

df2 <- data.frame(read_csv('Rdata2.csv'))

mylogit <- glm(disengaged ~ eventpos, data = df2, family = "binomial")
mylogit2 <- glm(disengaged ~ diffs + eventpos, data = df2, family = "binomial")
lrtest(mylogit, mylogit2)

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

ggplot(df2, aes(diffs, eventpos)) + geom_point() + 
  xlab(expression(Delta~random(x))) + ylab('position in sequence') +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust=0.5, size=26, face="bold"),
        axis.title.x = element_text(hjust=0.53, size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16)) 

# boxplot change in randomX for disengaged vs engaged 
ggplot(df2, aes(x=as.factor(disengaged), y=diffs)) + geom_boxplot()


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

