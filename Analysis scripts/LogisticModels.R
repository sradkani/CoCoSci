library(tidyverse)
library(lmtest)
library(Rmisc)
library(R.matlab)
library(tseries)
setwd("/Users/galraz1/Developer/CoCoSci/Analysis scripts")

df_deltarandomx <- data.frame(read_csv('Rdata_deltaRandomX.csv'))
df_randomX_normalized <- data.frame(read_csv('Rdata_randomX_normalized.csv'))
df_randomX_unnormalized <- data.frame(read_csv('Rdata_randomX_unnormalized.csv'))
df_CT <- data.frame(read_csv('Rdata_CT.csv'))

df_deltarandomx$quad = df_deltarandomx$diffs^2

mylogit_delta <- glm(disengaged ~ quad + diffs + eventpos, data = df_deltarandomx, family = "binomial")

df_randomX_normalized$quad = df_randomX_normalized$randomx^2
mylogit_randomx_norm <- glm(disengaged ~  randomx + eventpos, data = df_randomX_normalized, family = "binomial")

df_randomX_unnormalized$quad = df_randomX_unnormalized$randomx^2
mylogit_randomx_unnorm <- glm(disengaged ~ quad + randomx + eventpos, data = df_randomX_unnormalized, family = "binomial")

mylogit_CT <- glm(disengaged ~  CT + eventpos, data = df_CT, family = "binomial")

summary(mylogit_delta)
summary(mylogit_randomx_norm)
summary(mylogit_randomx_unnorm)
summary(mylogit_CT)