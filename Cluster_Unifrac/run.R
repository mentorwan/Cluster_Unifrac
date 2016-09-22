args = commandArgs(trailingOnly=TRUE)

row <- args[1]
row = as.numeric(row)

source("./GUniFrac_single.R")
load("./otu.RData")
library(GUniFrac)
data(throat.tree)

GUniFrac_single(otu.tab.rff,throat.tree,alpha=c(0,0.5,1),row)

