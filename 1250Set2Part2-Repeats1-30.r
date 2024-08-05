library(stringr)
library(cluster)
library(ggplot2)
library(ggdendro, lib.loc = "/mnt/Adenina/mgomez/esanchez/R/x86_64-pc-linux-gnu-library/4.2/")



ATFv1 <- read.table(file = "/mnt/Adenina/mgomez/esanchez/AUTC_Cluster/CoRaTables/OUT_ExSSs_ATFv1_1250Set2_mY_mY.txt", header = T)
ATFv2 <- read.table(file = "/mnt/Adenina/mgomez/esanchez/AUTC_Cluster/CoRaTables/OUT_ExSSs_ATFv2_1250Set2_mY_mY.txt", header = T)
BMFv1 <- read.table(file = "/mnt/Adenina/mgomez/esanchez/AUTC_Cluster/CoRaTables/OUT_ExSSs_BMFv1_1250Set2_mY_mY.txt", header = T)
BMFv2 <- read.table(file = "/mnt/Adenina/mgomez/esanchez/AUTC_Cluster/CoRaTables/OUT_ExSSs_BMFv2_1250Set2_mY_mY.txt", header = T)
BNFv1 <- read.table(file = "/mnt/Adenina/mgomez/esanchez/AUTC_Cluster/CoRaTables/OUT_ExSSs_BNFv1_1250Set2_mY_mY.txt", header = T)
BNFv2 <- read.table(file = "/mnt/Adenina/mgomez/esanchez/AUTC_Cluster/CoRaTables/OUT_ExSSs_BNFv2_1250Set2_mY_mY.txt", header = T)
FADv1 <- read.table(file = "/mnt/Adenina/mgomez/esanchez/AUTC_Cluster/CoRaTables/OUT_ExSSs_FADv1_1250Set2_mY_mY.txt", header = T)
FADv2 <- read.table(file = "/mnt/Adenina/mgomez/esanchez/AUTC_Cluster/CoRaTables/OUT_ExSSs_FADv2_1250Set2_mY_mY.txt", header = T)

ATFv1_Names <- paste0("ATFv1_", 1:nrow(ATFv1))
ATFv2_Names <- paste0("ATFv2_", 1:nrow(ATFv2))
BMFv1_Names <- paste0("BMFv1_", 1:nrow(BMFv1))
BMFv2_Names <- paste0("BMFv2_", 1:nrow(BMFv2))
BNFv1_Names <- paste0("BNFv1_", 1:nrow(BNFv1))
BNFv2_Names <- paste0("BNFv2_", 1:nrow(BNFv2))
FADv1_Names <- paste0("FADv1_", 1:nrow(FADv1))
FADv2_Names <- paste0("FADv2_", 1:nrow(FADv2))

CoRas <- as.data.frame(t(rbind(ATFv1[, 2:ncol(ATFv1)], ATFv2[, 2:ncol(ATFv2)], BMFv1[, 2:ncol(BMFv1)], BMFv2[, 2:ncol(BMFv2)], BNFv1[, 2:ncol(BNFv1)], BNFv2[, 2:ncol(BNFv2)], FADv1[, 2:ncol(FADv1)], FADv2[, 2:ncol(FADv2)])))
CoRas[CoRas > 1] <- NaN
colnames(CoRas) <- c(ATFv1_Names, ATFv2_Names, BMFv1_Names, BMFv2_Names, BNFv1_Names, BNFv2_Names, FADv1_Names, FADv2_Names)
CoRas <- CoRas[, !(colSums(is.na(CoRas)) == nrow(CoRas))]
rm(ATFv1, ATFv2, BMFv1, BMFv2, BNFv1, BNFv2, FADv1, FADv2, ATFv1_Names, ATFv2_Names, BMFv1_Names, BMFv2_Names, BNFv1_Names, BNFv2_Names, FADv1_Names, FADv2_Names)

print("I am now going to read the CalcABTC1250Set2.csv file")

CalcABTC <- read.csv("/mnt/Adenina/mgomez/esanchez/AUTC_Cluster/CalcABTC1250Set2.csv", header = T, sep = ",")
colnames(CalcABTC) <- colnames(CoRas)
rownames(CalcABTC) <- colnames(CalcABTC)

print("I have finished reading the CalcABTC1250Set2.csv file")

subsample <- function(data, percentage){
  indexes <- sort(sample(1:nrow(data), round(nrow(data)/100 * percentage), replace = FALSE))
  return(data[indexes, indexes])
}
subMaxK <- c()
subWholeK <- c()
set.seed(42)
i <- 30
print(paste0("I will start the process with the subsample of", i, "%"))
subABTC <- subsample(CalcABTC, i)
subModels <- unique(str_split_fixed(colnames(subABTC), "_", 2)[,1])
subIndexes <- sapply(subModels, function(M) grepl(paste0(M, "_"), colnames(subABTC)))
tryCatch(
  expr = {
    subModelsKCalc <- lapply(lapply(subModels, function(M1) subABTC[subIndexes[, M1], subIndexes[, M1]]), function(M2) clusGap(M2, FUN = kmeans, nstart = 15, K.max = 60, B = 100))
    subMaxK30 <- sum(sapply(subModelsKCalc, function(Model) maxSE(Model$Tab[, "gap"], Model$Tab[, "SE.sim"], method = "firstSEmax", SE.factor = 1)))
  },
  error = function(e){
    subMaxK30 <- NaN
  }
)
print(paste0("I finished the calculation for the individual models' clusters with a subsample of ", i))
tryCatch(
  expr = {
    temp <- clusGap(subABTC, FUN = kmeans, nstart = 15, K.max = subMaxK30, B = 150)
    subWholeK30<- maxSE(temp$Tab[,"gap"], temp$Tab[, "SE.sim"], method = "firstSEmax", SE.factor = 1)
  },
  error = function(e){
    subWholeK30 <- NaN
  }
)
print(paste0("I finished the calculation for the clusters of the entire subsample of ", i))

save.image(file = "/mnt/Adenina/mgomez/esanchez/AUTC_Cluster/ClusteringCalcs2-30.RData")
