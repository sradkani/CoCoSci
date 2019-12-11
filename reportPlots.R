library(tidyverse)


setwd("/Users/galraz1/Developer/CoCoSci")
df1 <- data.frame(read_csv('Rdata.csv'))


# lineplot randomXbins vs disengaged
ggplot(df1, aes(binMeans, propDisengaged)) + geom_line()+ geom_point()


df2 <- data.frame(read_csv('Rdata2.csv'))

# dotplot for logistic regression
ggplot(df2, aes(diffs, disengaged))  + geom_point() + theme_classic()

# boxplot change in randomX for disengaged vs engaged 
ggplot(df2, aes(x=as.factor(disengaged), y=diffs)) + geom_boxplot()


curvedf <- data.frame(read_csv('Rdata2.csv'))


