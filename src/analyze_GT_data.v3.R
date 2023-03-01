rm(list=ls())

library("dplyr")
library("scales")
library("ggplot2")
library("ggdendro")
library("reshape2")
library("grid")
library("ggrepel")

args = commandArgs(trailingOnly=TRUE)

file_name <-  args[1]
out_prefix <-  args[2]

###### debug
#setwd('~/Documents/GAU/Lebensohn')
#file_name = "/Volumes/GAU/projects/Lebensohn_2021/GT_pipeline_test/test_out/Results/sub.2_ALOD4_Top_S2.gene_level.temp"
#out_prefix="test"
####


# out_prefix<-basename(out_prefix)
# out_prefix<-tools::file_path_sans_ext(out_prefix)

tdata <- read.table(file_name,sep="\t",header=T)

# for genes with multiple transcripts, take the gene with the most total insertions
tdata <- tdata[order(-tdata$InsertionsGenes_Unsorted),]
data<- tdata[!duplicated(tdata$Gene), ]


tot_control = sum(data$InsertionsGenes_Unsorted)
tot_sel_inactivating = sum(data$TotalInactivating)
tot_sel_all = sum(data$TotalAllInsertions)


data$TotalMappedControl <-matrix(tot_control, length(data$Gene),1)
data$TotalMappedSelected.Inactivating <-matrix(tot_sel_inactivating, length(data$Gene),1)
data$TotalMappedSelected.All <-matrix(tot_sel_all, length(data$Gene),1)


## compute p-values
# Original eLife P-value (uses inactivating insertions in sorted sample)
data$PValue.V1 <- matrix(0,length(data$Gene),1)
data$FDRadjustedPvalue.V1 <- matrix(0, length(data$Gene),1)

# New P-value (uses ALL gene insertions in sorted sample, not just inactivating)
data$PValue.V2 <- matrix(0,length(data$Gene),1)
data$FDRadjustedPvalue.V2 <- matrix(0, length(data$Gene),1)

for(x in 1:length(data$Gene)){
  m1 <- matrix(c(data$InsertionsGenes_Unsorted[x],tot_control-data$InsertionsGenes_Unsorted[x],data$TotalInactivating[x],tot_sel_inactivating-data$TotalInactivating[x]),2,2)
  m2 <- matrix(c(data$InsertionsGenes_Unsorted[x],tot_control-data$InsertionsGenes_Unsorted[x], data$TotalAllInsertions[x],tot_sel_all-data$TotalAllInsertions[x]),2,2)
  fs1<-fisher.test(m1, alternative="less")
  fs2<-fisher.test(m2, alternative="less")
  data$PValue.V1[x] <-fs1$p.value
  data$PValue.V2[x] <-fs2$p.value
}

# FDR correct P-value and write table

data$FDRadjustedPvalue.V1<-p.adjust(data$PValue.V1,method="BH")
data$FDRadjustedPvalue.V2<-p.adjust(data$PValue.V2,method="BH")
sortedData<-data[order(data$FDRadjustedPvalue.V2),]
out_file = paste0(out_prefix,".v3.tab")
print(paste("Writing to file", out_file))
write.table(sortedData, out_file, append =FALSE, quote=FALSE, row.names = FALSE, col.names = TRUE, sep="\t")


## make plots

# volcano v1
#vars <- c("Gene", "FDRadjustedPvalue.V1","IGTIOB", "TotalAllInsertions")
#pp <- select(data,vars)
#pp <- filter(pp,IGTIOB > 0)
#pp$nLog10FDRadjustedPvalue.V1 <- -log10(pp$FDRadjustedPvalue.V1)
#pp <- pp[order(pp$nLog10FDRadjustedPvalue.V1),]
#pp$nLog10FDRadjustedPvalue.V1[!is.finite(pp$nLog10FDRadjustedPvalue.V1)] <- -log10(.Machine$double.xmin)

#max_p1<-max(pp$nLog10FDRadjustedPvalue.V1)

#plt <- ggplot(pp, aes(x=IGTIOB, y=nLog10FDRadjustedPvalue.V1,size=rescale(TotalAllInsertions,to=c(0,1000)))) + theme_bw()
#plt <- plt + guides(size=guide_legend(title="Total Insertions"))
#plt <- plt +
#  geom_point(aes(col=nLog10FDRadjustedPvalue.V1), alpha=0.8) +
#  scale_color_gradient(low = "#0091ff", high = "#f0650e",limits=c(3,max_p1),guide = "colourbar") +
#  labs(colour = "-log10 padj") +
#  geom_text_repel(data=tail(pp, 20), aes(label=Gene)) +
#  xlab("IGTIOB") +
#  ylab("-log10 padj") +
#  xlim(0,max(abs(pp$IGTIOB))+1) +
#  ggtitle(out_prefix) + theme(plot.title = element_text(hjust = 0.5))


#out_fig<-paste0('Results/',out_prefix,"_Scatter_V1Pval.pdf")
#ggsave(out_fig,width=6.5, height=6.5, units = "in")



# volcano v2
vars <- c("Gene", "FDRadjustedPvalue.V2","IGTIOB", "TotalAllInsertions")
pp <- select(data,vars)
pp$nLog10FDRadjustedPvalue.V2 <- -log10(pp$FDRadjustedPvalue.V2)

pp <- pp[order(pp$nLog10FDRadjustedPvalue.V2),]
pp$nLog10FDRadjustedPvalue.V2[!is.finite(pp$nLog10FDRadjustedPvalue.V2)] <- -log10(.Machine$double.xmin)
max_p2<-max(pp$nLog10FDRadjustedPvalue.V2)
plt <- ggplot(pp, aes(x=IGTIOB, y=nLog10FDRadjustedPvalue.V2,size=rescale(TotalAllInsertions,to=c(0,1000)))) + theme_bw()
plt <- plt + guides(size=guide_legend(title="Total Insertions"))
plt <- plt +
  geom_point(aes(col=nLog10FDRadjustedPvalue.V2), alpha=0.8) +
  scale_color_gradient(low = "#0091ff", high = "#f0650e",limits=c(3,max_p2),guide = "colourbar") +
  labs(colour = "-log10 padj") +
  geom_text_repel(data=tail(pp, 20), aes(label=Gene)) +
  xlab("IGTIOB") +
  xlim(-max(abs(pp$IGTIOB))-1,max(abs(pp$IGTIOB))+1) +
  ggtitle(out_prefix) +
  ylab("-log10 padj") +
  theme(plot.title = element_text(hjust = 0.5))

out_fig<-paste0(out_prefix,"_Volcano_V2Pval.pdf")
ggsave(out_fig, width=6.5, height=6.5, units = "in")




# Heatmap of insertion types
#vars<-c("Gene","SenseExon","AntisenseExon","SenseIntron","AntisenseIntron")
#tmp<-filter(data, data$FDRadjustedPvalue.V2 < 0.001 & data$TotalAllInsertions > 10)
#maps <- select(tmp, vars)
#maps.scaled<- maps
#maps.scaled[,c(2:5)]<- maps[,2:5] / rowSums(maps[,2:5]) # fraction of reads

# Run clustering
#maps.matrix <- as.matrix(maps.scaled)
#rownames(maps.matrix) <- maps.scaled$Gene
#maps.dendro<- as.dendrogram(hclust(d = dist(x = maps.matrix)))

# Create dendro
#dendro.plot <- ggdendrogram(data = maps.dendro, rotate = TRUE)

# Generate heatmap
#maps.long <- melt(maps.scaled)
#maps.order <- order.dendrogram(maps.dendro)
# Order the levels according to their position in the cluster
#maps.long$Gene <- factor(x = maps.long$Gene,
#                         levels = maps.scaled$Gene[maps.order],
#                         ordered = TRUE)

#heatmap.plot <- ggplot(data = maps.long,aes(x=variable,y=Gene)) +
#  geom_tile(aes(fill = value)) +
#  scale_fill_gradient2() +
#  theme(axis.text.y = element_text(size = 6)) +
#  xlab("Insertion type") +
#  ylab("Gene (FDR Gene Pval < 0.001)") +
#  ggtitle(out_prefix) +
#  theme(plot.title = element_text(hjust = 0.5)) +
#  theme_minimal()
#out_fig = paste0(out_prefix,"_insertion_heatmap.pdf")
#ggsave(out_fig, width=8.5, height=11, units = "in")



