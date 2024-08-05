Comparison <- function(a, b){
  if(is.nan(a) & is.nan(b)){
    return(0)
  }else if(is.nan(a) | is.nan(b)){
    return(1)
  }else{
    return(abs(a - b))
  }
}

ABTC <- function(Data){
  AreaMatrix <- matrix(0, nrow = ncol(Data), ncol = ncol(Data))
  AreaMatrix[lower.tri(AreaMatrix, diag = T)] <- unlist(sapply(c(1:(ncol(Data) - 1)), function(a) apply(sapply(Data[,(a):ncol(Data)], function(b) mapply(Comparison, Data[,a], b)), 2, sum)))
  AreaMatrix[upper.tri(AreaMatrix)] <- t(AreaMatrix)[upper.tri(AreaMatrix)]
  return(AreaMatrix)
}

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

CalcABTC <- ABTC(CoRas)
colnames(CalcABTC) <- colnames(CoRas)

write.csv(CalcABTC, "/mnt/Adenina/mgomez/esanchez/AUTC_Cluster/CalcABTC1250Set2.csv", row.names = FALSE)