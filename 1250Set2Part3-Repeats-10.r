library(stringr)
library(cluster)
library(ggplot2)
library(ggdendro, lib.loc = "/mnt/Adenina/mgomez/esanchez/R/x86_64-pc-linux-gnu-library/")
library(seriation, lib.loc = "/mnt/Adenina/mgomez/esanchez/R/x86_64-pc-linux-gnu-library/")
library(stringr)
library(dplyr)

load("/mnt/Adenina/mgomez/esanchez/AUTC_Cluster/ClusteringCalcs2-30-r10.RData")

DissimMat <- CalcABTC

pamAUTCMaxK <- pam(DissimMat, diss = T, k = subMaxK30)
ClusterMaxK <- pamAUTCMaxK$clustering
DissimMatMaxK <- DissimMat[names(sort(ClusterMaxK)), names(sort(ClusterMaxK))]
MaxKdissplot <- ggdissplot(as.dist(DissimMatMaxK), labels = sort(ClusterMaxK), method = NA, cluster_labels = T, cluster_lines = T, diag = F) + scale_fill_gradient(low = "midnightblue", high = "white", na.value = "grey75") + ggtitle("Individual Model Clustering, 30% of sample analized for clustering")

ggsave(filename = "MaxKdissplot.png", path = "/mnt/Adenina/mgomez/esanchez/AUTC_Cluster/30_Plots/r10/", plot = last_plot())

DissimMat <- CalcABTC
pamAUTCWholeK <- pam(DissimMat, diss = T, k = subWholeK30)
ClusterWholeK <- pamAUTCWholeK$clustering
DissimMatWholeK <- DissimMat[names(sort(ClusterWholeK)), names(sort(ClusterWholeK))]
WholeKdissplot <- ggdissplot(as.dist(DissimMatWholeK), labels = sort(ClusterWholeK), method = NA, cluster_labels = T, cluster_lines = T, diag = F) + scale_fill_gradient(low = "midnightblue", high = "white", na.value = "grey75") + ggtitle("AllxAll Clustering, 30% of sample analized for clustering")

ggsave(filename = "WholeKdissplot.png", path = "/mnt/Adenina/mgomez/esanchez/AUTC_Cluster/30_Plots/r10/", plot = last_plot())

clustermedoids <- CalcABTC[pamAUTCMaxK$medoids,pamAUTCMaxK$medoids]
colnames(clustermedoids) <- paste0("Cluster ", pamAUTCMaxK$clustering[pamAUTCMaxK$medoids])
rownames(clustermedoids) <- colnames(clustermedoids)
dendro <- hclust(as.dist(clustermedoids), method = "average")
MaxKDendrogram <- ggdendrogram(dendro) + ggtitle("Signature Behaviours according to Intra-Motif Clustering, 30% subsample for Clustering") + theme(plot.title = element_text(size = 32))

ggsave(filename = "MaxKDendrogram.png", path = "/mnt/Adenina/mgomez/esanchez/AUTC_Cluster/30_Plots/r10/", plot = last_plot())

clustermedoids <- CalcABTC[pamAUTCWholeK$medoids,pamAUTCWholeK$medoids]
colnames(clustermedoids) <- paste0("Cluster ", pamAUTCWholeK$clustering[pamAUTCWholeK$medoids])
rownames(clustermedoids) <- colnames(clustermedoids)
dendro <- hclust(as.dist(clustermedoids), method = "average")
WholeKDendrogram <- ggdendrogram(dendro) + ggtitle("Signature Behaviours according to AllxAll Clustering, 30% subsample for Clustering") + theme(plot.title = element_text(size = 32))

ggsave(filename = "WholeKDendrogram.png", path = "/mnt/Adenina/mgomez/esanchez/AUTC_Cluster/30_Plots/r10/", plot = last_plot())

pamAUTCMaxK_10 <- pamAUTCMaxK
pamAUTCWholeK_10 <- pamAUTCWholeK

saveRDS(pamAUTCMaxK_10, file = "/mnt/Adenina/mgomez/esanchez/AUTC_Cluster/30_Plots/r10/pamAUTCMaxK_10.RDS")
saveRDS(pamAUTCWholeK_10, file = "/mnt/Adenina/mgomez/esanchez/AUTC_Cluster/30_Plots/r10/pamAUTCWholeK_10.RDS")