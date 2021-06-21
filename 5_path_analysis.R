# Clear
cat("\014")
rm(list = ls())

library(EnvStats)

# test
data1 <- read.csv('multi-label-emp/kingdom/k__Bacteria.csv', header = T, row.names = 1)
data2 <- read.csv('multi-label-emp/kingdom/k__Archaea.csv', header = T, row.names = 1)
data3 <- read.csv('multi-label-emp/phylum/p__Proteobacteria.csv', header = T, row.names = 1)
data4 <- read.csv('multi-label-emp/phylum/p__OP3.csv', header = T, row.names = 1)


kmers <- rownames(data1)

imp1 <- data1$TRUE.
imp2 <- data2$TRUE.
imp3 <- data3$TRUE.
imp4 <- data4$TRUE.

neg1 <- data1$FALSE.
neg2 <- data2$FALSE.

data.frame(cbind(imp1, imp2))

plot(imp1, type = 'l')
plot(imp3, type = 'l')
plot(imp4, type = 'l')

a <- summaryStats(imp1, quartiles = TRUE)