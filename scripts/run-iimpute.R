source('./I-Impute.R')

args <- commandArgs(TRUE)
data <- read.csv(args[1],row.name=1)
iimpute(as.data.frame(t(data)), paste(args[1], '.iimpute.csv', sep=''))