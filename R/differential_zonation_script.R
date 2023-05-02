library(ggplot2)
library(DESeq2)
library(stringr)
library(foreach)

source("./R/misc_function.R")
source("./R/diffzone_function.R")

df <-  read.csv("./data/20230430_ExpSums_dl.csv", row.names = 1)
counts <- data.matrix(df)
# G <- dim(counts)[1]
# top500 <- order(rowSums(counts))[9500:G]
# 
# counts <- counts[top500, ]
# counts <- counts[1:10, 1:dim(counts)[2]]
mode(counts) <- 'integer'

dose_layer_sample <- as.data.frame(str_split_fixed(colnames(df), "_", 3))
colnames(dose_layer_sample) <- c("Dose", "Layer") 
dose_layer_sample
data <- list(countData = counts, group = dose_layer_sample$Dose, zone = dose_layer_sample$Layer)
countData = data[["countData"]]
group     = data[["group"]]
zone      = as.numeric(data[["zone"]])

diffzone_list  = diffzone(countData,group,zone)

diffzone_plot(diffzone_list, "Stat1")

