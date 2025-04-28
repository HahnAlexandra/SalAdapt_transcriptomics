##### Transcriptomic analysis #####
new_lib <- "~/new_R_library"
.libPaths(c(new_lib, .libPaths()))

library(tidyverse)
library(DESeq2)
library(pheatmap)
library(DEGreport)
library(UpSetR)
library(VennDiagram)
library(GO.db)
library(clusterProfiler)
library(enrichplot)
library(xlsx)
library(rtracklayer)
library(readr)
library(stringr)

setwd("~/Documents/Experiments/Transcriptomics/DESeq2")

#load files
samples <- read.csv("~/Documents/Experiments/Transcriptomics/DESeq2/samples_v3.csv")
samples$station <- factor(samples$station)
samples$time <- factor(samples$time, levels = c("t1", "t2", "t0"))#reordering is necessary to manually build matrix
samples$treatment <- factor(samples$treatment, levels = c("15", "7"))
#samples$interaction represents the collapsed interaction term

txi <- readRDS("~/Documents/Experiments/Transcriptomics/DESeq2/txi_2.rds")

##### Differential gene expression #####
#LRT model matrix ####
#manually edit the model matrices
#full model
m1 <- model.matrix(~ time * treatment * station, samples)

all.zero <- apply(m1, 2, function(x) all(x==0))
all.zero

idx <- which(all.zero)
m1 <- m1[,-idx]

#salinity
ms <- model.matrix(~time + station + time:station , samples)

#time
mt <- model.matrix(~treatment + station +  station:treatment, samples)

#population
mp <- model.matrix(~treatment + time  + time:treatment, samples)

all.zero <- apply(mp, 2, function(x) all(x==0))
idx <- which(all.zero)
mp <- mp[,-idx]

#salxtime
m2 <- model.matrix(~ time + treatment + station + time:station + station:treatment, samples)

#salxpop
m3 <- model.matrix(~ time + treatment + station + time:station + time:treatment, samples)

all.zero <- apply(m3, 2, function(x) all(x==0))
idx <- which(all.zero)
m3 <- m3[,-idx]

#salxpopxtime
m4 <- model.matrix(~time + treatment + station +
                     time:station  + time:treatment + station:treatment, samples)

all.zero <- apply(m4, 2, function(x) all(x==0))
idx <- which(all.zero)
m4 <- m4[,-idx]

#popxtime
m5 <- model.matrix(~time + treatment + station +
                     time:treatment + station:treatment, samples)

all.zero <- apply(m5, 2, function(x) all(x==0))
idx <- which(all.zero)
m5 <- m5[,-idx]


#DESeq ####
#use manually corrected model matrix as input, m1 = full model
dds  <- DESeqDataSetFromTximport(txi, colData = samples, 
                                 design = m1)

#filter out low counts, removes rows where count < 10 in more than 90% of samples
keep <- apply(counts(dds), 1, function(x) {ifelse(length(which(x > 10)) > (nrow(samples)*0.90), TRUE, FALSE)})
dds <- dds[keep,]
nrow(dds)#25567

rld <- rlog(dds, blind=FALSE)#regularized log transformation

pcaData <- plotPCA(rld, intgroup=c("station", "time", "treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

#plot time PCA
ggplot(pcaData, aes(PC1, PC2, color=station, fill = time)) +
  theme_light()+
  geom_point(size=3, shape = 21) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  scale_fill_manual(values =c("#E5EDFD", "#A7C3FB", "white"))+
  scale_color_manual(values = c("#83B3C2", "#DE3C22"))+
  stat_ellipse()

#plot treatment PCA
ggplot(pcaData, aes(PC1, PC2, color=station, fill = treatment)) +
  theme_light()+
  geom_point(size=3, shape = 21) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  scale_fill_manual(values =c("#E5EDFD", "#A7C3FB"))+
  scale_color_manual(values = c("#83B3C2", "#DE3C22"))+
  stat_ellipse()

# LRT for salinity effect ####
dds_salinity <- DESeq(dds, test = "LRT", reduced = ms)

boxplot(log10(assays(dds_salinity)[["cooks"]]), range=0, las=2)
plotDispEsts(dds_salinity)

#results
res_salinity <- results(dds_salinity, alpha = 0.05)
summary(res_salinity)
sum(res_salinity$padj < 0.05, na.rm = TRUE)#889

#subset significant genes
res_salinity_sub <- subset(res_salinity, padj < 0.05)
sig_genes_salinity <- rownames(res_salinity_sub)#889

dds_sig <- dds[sig_genes_salinity,] #subset the big set
rld <- rlog(dds_sig, blind=FALSE) #regularized log transformation

#visualize in PCA
pcaData <- plotPCA(rld, intgroup=c("station", "time", "treatment"), returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

#plot salinity PCA
ggplot(pcaData, aes(PC1, PC2, color=station, fill = time, shape = treatment)) +
  theme_light()+
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  scale_fill_manual(values =c( "#E5EDFD", "#A7C3FB", "white"))+
  scale_color_manual(values = c("#83B3C2", "#DE3C22"))+
  scale_shape_manual(values = c(24, 25))

#visualize in heatmap
dds_sig <- estimateSizeFactors(dds_sig)
norm_counts_sig <- counts(dds_sig, normalized = TRUE)

# order samples manually
desired_order <- c("1S05_1_L1", "1S05_2_L1", "1S05_3_L1","2S05_1_L5", "2S05_2_L5", "2S05_3_L5",
                   "1W05_1_L5", "1W05_2_L5", "1W05_3_L5", "2W05_1_L5", "2W05_2_L5", "2W05_3_L5",
                   "0S15_1_L5", "0S15_2_L5", "0S15_3_L5", "0W15_1_L5" ,"0W15_2_L5", "0W15_3_L1",
                   "1S15_1_L1", "1S15_2_L5","1S15_3_L5",  "1W15_1_L4", "1W15_2_L4","1W15_3_L5",
                   "2S15_1_L5", "2S15_2_L5", "2S15_3_L5", "2W15_1_L5", "2W15_2_L5", "2W15_3_L5")  

# reorder the columns of the matrix
mat_reordered <- norm_counts_sig[, desired_order]
df <- as.data.frame(colData(dds)[,c("treatment", "time", "station")])

pheatmap(mat_reordered, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames=FALSE,
         annotation= df, border_color=NA, fontsize = 10, scale="row",
         fontsize_row = 10, height=20)


# LRT for salinity x time interaction ####
dds_salxtime <- DESeq(dds, test = "LRT", reduced = m2)

#results
res_salxtime <- results(dds_salxtime, alpha = 0.05)
summary(res_salxtime)
sum(res_salxtime$padj < 0.05, na.rm = TRUE)

res_salxtime_sub <- subset(res_salxtime, padj < 0.05)
sig_genes_salxtime <- rownames(res_salxtime_sub)#105

dds_sig <- dds[sig_genes_salxtime,] #subset the big set
rld <- rlog(dds_sig, blind=FALSE) #regularized log transformation

#visualize in PCA
pcaData <- plotPCA(rld, intgroup=c("station", "time", "treatment"), returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

#plot salinity PCA
ggplot(pcaData, aes(PC1, PC2, color=station, fill = time, shape = treatment)) +
  theme_light()+
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  scale_fill_manual(values =c( "#E5EDFD", "#A7C3FB", "white"))+
  scale_color_manual(values = c("#83B3C2", "#DE3C22"))+
  scale_shape_manual(values = c(24, 25))

#visualize in heatmap
dds_sig <- estimateSizeFactors(dds_sig)
norm_counts_sig <- counts(dds_sig, normalized = TRUE)

# order samples manually
desired_order <- c("1S05_1_L1", "1S05_2_L1", "1S05_3_L1","2S05_1_L5", "2S05_2_L5", "2S05_3_L5",
                   "1W05_1_L5", "1W05_2_L5", "1W05_3_L5", "2W05_1_L5", "2W05_2_L5", "2W05_3_L5",
                   "0S15_1_L5", "0S15_2_L5", "0S15_3_L5", "0W15_1_L5" ,"0W15_2_L5", "0W15_3_L1",
                   "1S15_1_L1", "1S15_2_L5","1S15_3_L5",  "1W15_1_L4", "1W15_2_L4","1W15_3_L5",
                   "2S15_1_L5", "2S15_2_L5", "2S15_3_L5", "2W15_1_L5", "2W15_2_L5", "2W15_3_L5")  

# reorder the columns of the matrix
mat_reordered <- norm_counts_sig[, desired_order]
df <- as.data.frame(colData(dds)[,c("treatment", "time", "station")])

pheatmap(mat_reordered, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames=FALSE,
         annotation= df, border_color=NA, fontsize = 10, scale="row",
         fontsize_row = 10, height=20)

# LRT for salinity x population interaction ####
dds_salxpop <- DESeq(dds, test = "LRT", reduced = m3)

#results
res_salxpop <- results(dds_salxpop, alpha = 0.05)
summary(res_salxpop)
sum(res_salxpop$padj < 0.05, na.rm = TRUE)#2

# LRT for salinity x population x time interaction ####
dds_salxpopxtime <- DESeq(dds, test = "LRT", reduced = m4)

#results
res_salxpopxtime <- results(dds_salxpopxtime, alpha = 0.05)
summary(res_salxpopxtime)
sum(res_salxpopxtime$padj < 0.05, na.rm = TRUE)#1


## pairwise comparisons #####
#with collapsed interaction term 

#create dds object
dds_int <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ interaction)

#filtering 
keep <- apply(counts(dds_int), 1, function(x) {ifelse(length(which(x > 10)) > (nrow(samples)*0.9), TRUE, FALSE)})
dds_int <- dds_int[keep,]
nrow(dds_int)#25567(sanity check)

#run DESeq2
dds_int <- DESeq(dds_int)
dds_int_allgenes <- rownames(dds_int)

# pairwise comparisons sal x time ####
#build time sensitive gene set by comparing 15t1 and 15t2 to 15t0
#SW08 15t1 to 15t0
res_SW08_int1 <- results(dds_int, contrast = c("interaction", "1S15", "0S15"), alpha = 0.05)
summary(res_SW08_int1)
sum(res_SW08_int1$padj < 0.05, na.rm=TRUE)#90

res_SW08_int1 <- subset(res_SW08_int1, padj < 0.05)
timegenes_SW08_1 <- rownames(res_SW08_int1)

#SW08 15t2 to 15t0
res_SW08_int2 <- results(dds_int, contrast = c("interaction", "2S15", "0S15"), alpha = 0.05)
summary(res_SW08_int2)
sum(res_SW08_int2$padj < 0.05, na.rm=TRUE)#674

res_SW08_int2 <- subset(res_SW08_int2, padj < 0.05)
timegenes_SW08_2 <- rownames(res_SW08_int2)

#SW08 15t1 to 15t2
res_SW08_int3 <- results(dds_int, contrast = c("interaction", "2S15", "1S15"), alpha = 0.05)
summary(res_SW08_int3)
sum(res_SW08_int3$padj < 0.05, na.rm=TRUE)#206

res_SW08_int3 <- subset(res_SW08_int3, padj < 0.05)
timegenes_SW08_3 <- rownames(res_SW08_int3)

#WH 15t1 to 15t0
res_WH_int1 <- results(dds_int, contrast = c("interaction", "1W15", "0W15"), alpha = 0.05)
summary(res_WH_int1)
sum(res_WH_int1$padj < 0.05, na.rm=TRUE)#281

res_WH_int1 <- subset(res_WH_int1, padj < 0.05)
timegenes_WH_1 <- rownames(res_WH_int1)

#WH 15t2 to 15t0
res_WH_int2 <- results(dds_int, contrast = c("interaction", "2W15", "0W15"), alpha = 0.05)
summary(res_WH_int1)
sum(res_WH_int2$padj < 0.05, na.rm=TRUE)#341

res_WH_int2 <- subset(res_WH_int2, padj < 0.05)
timegenes_WH_2 <- rownames(res_WH_int2)

#WH 15t1 to 15t2
res_WH_int3 <- results(dds_int, contrast = c("interaction", "2W15", "1W15"), alpha = 0.05)
summary(res_WH_int3)
sum(res_WH_int3$padj < 0.05, na.rm=TRUE)#410

res_WH_int3 <- subset(res_WH_int3, padj < 0.05)
timegenes_WH_3 <- rownames(res_WH_int3)

#timegene list
timegenes_SW08 <- union(timegenes_SW08_1, timegenes_SW08_2)
timegenes_SW08 <- union(timegenes_SW08, timegenes_SW08_3)#781 genes responding between time points

timegenes_WH <- union(timegenes_WH_1, timegenes_WH_2)
timegenes_WH <- union(timegenes_WH, timegenes_WH_3)#848 genes responding between time points

#build list of salinity sensitive genes by contrasting low salinity & control per time
#SW08 5t1 to 15t1

res_SW08_sal1 <- results(dds_int, contrast = c("interaction", "1S05", "1S15"), alpha = 0.05)
summary(res_SW08_sal1)
sum(res_SW08_sal1$padj < 0.05, na.rm=TRUE)#207

res_SW08_sal1_sub<- subset(res_SW08_sal1, padj < 0.05)
salgenes1_SW08 <- rownames(res_SW08_sal1_sub)

#up and down regulated
res_SW08_sal1_up<- subset(res_SW08_sal1, padj < 0.05 & log2FoldChange > 0)
salgenes1_SW08_up <- rownames(res_SW08_sal1_up)#181

res_SW08_sal1_down<- subset(res_SW08_sal1, padj < 0.05 & log2FoldChange < 0)
salgenes1_SW08_down <- rownames(res_SW08_sal1_down)#26

#SW08 5t2 to 15t2
res_SW08_sal2 <- results(dds_int, contrast = c("interaction", "2S05", "2S15"), alpha = 0.05)
summary(res_SW08_sal2)
sum(res_SW08_sal2$padj < 0.05, na.rm=TRUE)#239

res_SW08_sal2_sub<- subset(res_SW08_sal2, padj < 0.05)
salgenes2_SW08 <- rownames(res_SW08_sal2_sub)

#up and down regulated
res_SW08_sal2_up<- subset(res_SW08_sal2, padj < 0.05 & log2FoldChange > 0)
salgenes2_SW08_up <- rownames(res_SW08_sal2_up)#168

res_SW08_sal2_down<- subset(res_SW08_sal2, padj < 0.05 & log2FoldChange < 0)
salgenes2_SW08_down <- rownames(res_SW08_sal2_down)#71

#WH 5t1 to 15t1
res_WH_sal1 <- results(dds_int, contrast = c("interaction", "1W05", "1W15"), alpha = 0.05)
summary(res_WH_sal1)
sum(res_WH_sal1$padj < 0.05, na.rm=TRUE)#165

res_WH_sal1_sub <- subset(res_WH_sal1, padj < 0.05)
salgenes1_WH <- rownames(res_WH_sal1_sub)

#up and down regulated
res_WH_sal1_up <- subset(res_WH_sal1, padj < 0.05 & log2FoldChange > 0)
salgenes1_WH_up <- rownames(res_WH_sal1_up)#108

res_WH_sal1_down <- subset(res_WH_sal1, padj < 0.05 & log2FoldChange < 0)
salgenes1_WH_down <- rownames(res_WH_sal1_down)#57

#WH 5t2 to 15t2
res_WH_sal2 <- results(dds_int, contrast = c("interaction", "2W05", "2W15"), alpha = 0.05)
summary(res_WH_sal2)
sum(res_WH_sal2$padj < 0.05, na.rm=TRUE)#370

res_WH_sal2_sub <- subset(res_WH_sal2, padj < 0.05)
salgenes2_WH <- rownames(res_WH_sal2_sub)

#up and down regulated
res_WH_sal2_up <- subset(res_WH_sal2, padj < 0.05 & log2FoldChange > 0)
salgenes2_WH_up <- rownames(res_WH_sal2_up)#118

res_WH_sal2_down <- subset(res_WH_sal2, padj < 0.05 & log2FoldChange < 0)
salgenes2_WH_down <- rownames(res_WH_sal2_down)#252

#subtract time sensitive genes
salgenes1_SW08 <- setdiff(salgenes1_SW08, timegenes_SW08)#186
salgenes2_SW08 <- setdiff(salgenes2_SW08, timegenes_SW08)#201

salgenes1_WH <- setdiff(salgenes1_WH, timegenes_WH)#144
salgenes2_WH <- setdiff(salgenes2_WH, timegenes_WH)#299

#combine into 1 set with both time points
salgenes_SW08 <- union(salgenes1_SW08, salgenes2_SW08)#341
salgenes_WH <- union(salgenes1_WH, salgenes2_WH)#415

salgenes_inter <- intersect(salgenes_SW08, salgenes_WH)#154
salgenes_combined <- union(salgenes_SW08, salgenes_WH)#602
#of combined genes 516 overlap with LRT salinity genes
salgenes_LRT_pair <- intersect(salgenes_combined, sig_genes_salinity)#516

allgenes_dds_gene_names <- gff_df %>%
  filter(TRINITY %in% dds_int_allgenes)%>%
  distinct(TRINITY, .keep_all = TRUE)%>%
  filter(Gene != "---NA---")

#unique for populations
salgenes_SW08_u <- setdiff(salgenes_SW08, salgenes_inter)#187
salgenes_SW08_u_names <- gff_df %>% 
  filter(TRINITY %in% salgenes_SW08_u)

salgenes_WH_u <- setdiff(salgenes_WH, salgenes_inter)#261
salgenes_WH_u_names <- gff_df %>% 
  filter(TRINITY %in% salgenes_WH_u)

#shared and unique for time points
salgenes_inter_t1 <- intersect(salgenes1_SW08, salgenes1_WH)#79
salgenes_inter_t2 <- intersect(salgenes2_SW08, salgenes2_WH)#91

salgenes_SW08_u_t1 <- setdiff(salgenes1_SW08, salgenes_inter)#99
salgenes_SW08_u_t2 <- setdiff(salgenes2_SW08, salgenes_inter)#98

salgenes_WH_u_t1 <- setdiff(salgenes1_WH, salgenes_inter)#63
salgenes_WH_u_t2 <- setdiff(salgenes2_WH, salgenes_inter)#202

#add log2fold changes
extract_lfc <- function(res_obj, col_name) {
  df <- data.frame(
    TRINITY = rownames(res_obj),  # Extract row names as a column
    log2FoldChange = res_obj@listData$log2FoldChange  # Extract log2FoldChange values
  )
  setNames(df, c("TRINITY", col_name))  # Rename column dynamically
}

#extract log2fold changes
lfc_SW08_t1 <- extract_lfc(res_SW08_sal1, "lfc_SW08_t1")#use not sub-setted results object to get also non-sig lfc
lfc_SW08_t2 <- extract_lfc(res_SW08_sal2, "lfc_SW08_t2")
lfc_WH_t1 <- extract_lfc(res_WH_sal1, "lfc_WH_t1")
lfc_WH_t2 <- extract_lfc(res_WH_sal2, "lfc_WH_t2")

#load annotation file
gff_file <- "blast2go_gff_export_20170613_0815.gff.gz"
gff <- import(gff_file)

#convert to data frame
gff_df <- as.data.frame(gff) %>%
  dplyr::filter(type == "CDS") %>%
  dplyr::select(seqnames, start, end, strand, ID, Description, Gene) %>%
  dplyr::mutate(
    transcript_id = sapply(strsplit(as.character(ID), " "), `[`, 1),
    gene_id = sapply(strsplit(as.character(Gene), ";"), `[`, 1)
  )

gff_df <- gff_df %>%
  mutate(TRINITY = sub("_i[0-9]+$", "", ID))

#extract gene names
salgenes_inter_names <- gff_df %>%
  filter(TRINITY %in% salgenes_inter)

salgenes_combined_names <- gff_df %>%
  filter(TRINITY %in% salgenes_combined)

salgenes_LRT_pair_names <- gff_df %>%
  filter(TRINITY %in% salgenes_LRT_pair)

salgenes_LRT_names <- gff_df %>%
  filter(TRINITY %in% sig_genes_salinity)

salgenes_SW08_u_names <- gff_df %>% 
  filter(TRINITY %in% salgenes_SW08_u)

salgenes_WH_u_names <- gff_df %>% 
  filter(TRINITY %in% salgenes_WH_u)

salgenes_inter_t1_names <- gff_df %>% 
  filter(TRINITY %in% salgenes_inter_t1)

salgenes_inter_t2_names <- gff_df %>% 
  filter(TRINITY %in% salgenes_inter_t2)

salgenes_SW08_u_t1_names <- gff_df %>% 
  filter(TRINITY %in% salgenes_SW08_u_t1)

salgenes_SW08_u_t2_names <- gff_df %>% 
  filter(TRINITY %in% salgenes_SW08_u_t2)

salgenes_WH_u_t1_names <- gff_df %>% 
  filter(TRINITY %in% salgenes_WH_u_t1)

salgenes_WH_u_t2_names <- gff_df %>% 
  filter(TRINITY %in% salgenes_WH_u_t2)

#get one big table that lists all functions and gene sets

salgenes_list <- list(
  LRT = salgenes_LRT_names,
  pairwise = salgenes_combined_names,
  pairwise_shared = salgenes_inter_names,
  pairwise_SW08 = salgenes_SW08_u_names,
  pairwise_WH = salgenes_WH_u_names,
  shared_t1 = salgenes_inter_t1_names,
  shared_t2 = salgenes_inter_t2_names,
  SW08_t1 = salgenes_SW08_u_t1_names,
  SW08_t2 = salgenes_SW08_u_t2_names,
  WH_t1 = salgenes_WH_u_t1_names,
  WH_t2 = salgenes_WH_u_t2_names
)

#bind all data sets, with row for origin
salgenes_list_combined <- bind_rows(
  lapply(names(salgenes_list), function(name) {
    salgenes_list[[name]] %>% mutate(Origin = name)
  })
)

#filter out genes with "---NA---" as description
salgenes_list_filtered <- salgenes_list_combined %>%
  filter(Description != "---NA---")

#create presence/absence and filter TRINITY duplicates
salgenes_list_combined <- salgenes_list_filtered %>%
  mutate(Present = 1)%>%
  pivot_wider(names_from = Origin, values_from = Present, values_fill = 0)%>%
  dplyr::select(TRINITY, Gene, Description,
         LRT, pairwise, pairwise_shared, pairwise_SW08, pairwise_WH,
         shared_t1, shared_t2, SW08_t1, SW08_t2, WH_t1, WH_t2)%>%
  distinct(TRINITY, .keep_all = TRUE)#add Gene to avoid losing Genes with same TRINITY

salgenes_list_combined <- salgenes_list_combined %>%
  left_join(lfc_SW08_t1, by = "TRINITY")%>%
  left_join(lfc_SW08_t2, by = "TRINITY")%>%
  left_join(lfc_WH_t1, by = "TRINITY")%>%
  left_join(lfc_WH_t2, by = "TRINITY")

write.csv(salgenes_list_combined, "Salinity_genes_table_lfc-uniq.csv", row.names=FALSE, quote=TRUE)

# visualize pairwise comparison sets ####
#visualize with heatmaps

#Baltic Sea
dds_sig_int <- dds_int[salgenes_SW08,]

#extract normalized counts
norm_counts_sig <- counts(dds_sig_int, normalized = TRUE)

#subset sample meta data
norm_counts_sig_SW08 <- norm_counts_sig[, grep("S", colnames(norm_counts_sig))]

# order samples manually
desired_order_SW08 <- c( "1S05_1_L1", "1S05_2_L1", "1S05_3_L1","2S05_1_L5", "2S05_2_L5", "2S05_3_L5",
                         "0S15_1_L5", "0S15_2_L5", "0S15_3_L5",
                         "1S15_1_L1", "1S15_2_L5","1S15_3_L5",  "2S15_1_L5", "2S15_2_L5", "2S15_3_L5")

# reorder the columns of the matrix
mat_reordered <- norm_counts_sig[, desired_order_SW08]

df <- as.data.frame(colData(dds_int)[,c("treatment", "time", "station")])

# plot heat map
pheatmap(mat_reordered, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames=FALSE,
         annotation= df, border_color=NA, fontsize = 10, scale="row",
         fontsize_row = 10, height=20)

#North Sea
dds_sig_int <- dds_int[salgenes_WH,]

#extract normalize counts
norm_counts_sig <- counts(dds_sig_int, normalized = TRUE)

#subset sample meta data
norm_counts_sig_WH <- norm_counts_sig[, grep("W", colnames(norm_counts_sig))]

# order samples manually
desired_order_WH<- c("1W05_1_L5", "1W05_2_L5", "1W05_3_L5", "2W05_1_L5", "2W05_2_L5", "2W05_3_L5" ,
                     "0W15_1_L5" ,"0W15_2_L5", "0W15_3_L1",
                     "1W15_1_L4", "1W15_2_L4","1W15_3_L5", "2W15_1_L5", "2W15_2_L5", "2W15_3_L5")  


# reorder the columns of the matrix
mat_reordered <- norm_counts_sig_WH[, desired_order_WH]

df <- as.data.frame(colData(dds_int)[,c("treatment", "time", "station")])

# plot heat map
pheatmap(mat_reordered, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames=FALSE,
         annotation= df, border_color=NA, fontsize = 10, scale="row",
         fontsize_row = 10, height=20)


#visualize with venn diagram
length(salgenes_SW08)#341
length(salgenes_WH)#415
intersect(salgenes_SW08, salgenes_WH)#154 overlap

grid.newpage()
draw.pairwise.venn(area1= 415,area2 = 341, cross.area = 154,
                   category=c("North","Baltic"),fill=c("Red","Blue"))

#visualize with upSet plot
gene_sets <- list(
  Baltic_t1 = salgenes1_SW08,
  Baltic_t2 = salgenes2_SW08,
  North_t1 = salgenes1_WH,
  North_t2 = salgenes2_WH
)


# create a matrix with genes as rows and sets as columns
all_genes <- unique(unlist(gene_sets))#gets unique genes from sets
presence_absence_matrix <- sapply(gene_sets, function(genes) all_genes %in% genes)
rownames(presence_absence_matrix) <- all_genes

presence_absence_matrix <- as.data.frame(presence_absence_matrix)

#plot
upset(fromList(gene_sets), order.by = "freq",nsets = 4, keep.order = TRUE)

#heat map for shared genes
dds_sig <- dds_int[salgenes_inter,]#subset the big set
norm_counts_sig <- counts(dds_sig, normalized = TRUE)

# order samples manually
desired_order <- c("1S05_1_L1", "1S05_2_L1", "1S05_3_L1","2S05_1_L5", "2S05_2_L5", "2S05_3_L5",
                   "1W05_1_L5", "1W05_2_L5", "1W05_3_L5", "2W05_1_L5", "2W05_2_L5", "2W05_3_L5",
                   "0S15_1_L5", "0S15_2_L5", "0S15_3_L5", "0W15_1_L5" ,"0W15_2_L5", "0W15_3_L1",
                   "1S15_1_L1", "1S15_2_L5","1S15_3_L5",  "1W15_1_L4", "1W15_2_L4","1W15_3_L5",
                   "2S15_1_L5", "2S15_2_L5", "2S15_3_L5", "2W15_1_L5", "2W15_2_L5", "2W15_3_L5")  

# reorder the columns of the matrix
mat_reordered <- norm_counts_sig[, desired_order]
df <- as.data.frame(colData(dds_int)[,c("treatment", "time", "station")])

pheatmap(mat_reordered, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames=FALSE,
         annotation= df, border_color=NA, fontsize = 10, scale="row",
         fontsize_row = 10, height=20)

#heat map for unique SW08 genes
dds_sig <- dds_int[salgenes_SW08_u,]#subset the big set

rld <- rlog(dds_sig, blind=TRUE)

#does it make sense to visualize all with the significant genes from one station?
norm_counts_sig <- counts(dds_sig, normalized = TRUE)

# order samples manually
desired_order <- c("1S05_1_L1", "1S05_2_L1", "1S05_3_L1","2S05_1_L5", "2S05_2_L5", "2S05_3_L5",
                   "1W05_1_L5", "1W05_2_L5", "1W05_3_L5", "2W05_1_L5", "2W05_2_L5", "2W05_3_L5",
                   "0S15_1_L5", "0S15_2_L5", "0S15_3_L5", "0W15_1_L5" ,"0W15_2_L5", "0W15_3_L1",
                   "1S15_1_L1", "1S15_2_L5","1S15_3_L5",  "1W15_1_L4", "1W15_2_L4","1W15_3_L5",
                   "2S15_1_L5", "2S15_2_L5", "2S15_3_L5", "2W15_1_L5", "2W15_2_L5", "2W15_3_L5")  

# reorder the columns of the matrix
mat_reordered <- norm_counts_sig[, desired_order]
df <- as.data.frame(colData(dds_int)[,c("treatment", "time", "station")])

pheatmap(mat_reordered, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames=FALSE,
         annotation= df, border_color=NA, fontsize = 10, scale="row",
         fontsize_row = 10, height=20)

#heat map for unique WH genes
#visualize in heatmap
dds_sig <- dds_int[salgenes_WH_u,]#subset the big set

rld <- rlog(dds_sig, blind=TRUE)

#does it make sense to visualize all with the significant genes from one station?
norm_counts_sig <- counts(dds_sig, normalized = TRUE)

# order samples manually
desired_order <- c("1S05_1_L1", "1S05_2_L1", "1S05_3_L1","2S05_1_L5", "2S05_2_L5", "2S05_3_L5",
                   "1W05_1_L5", "1W05_2_L5", "1W05_3_L5", "2W05_1_L5", "2W05_2_L5", "2W05_3_L5",
                   "0S15_1_L5", "0S15_2_L5", "0S15_3_L5", "0W15_1_L5" ,"0W15_2_L5", "0W15_3_L1",
                   "1S15_1_L1", "1S15_2_L5","1S15_3_L5",  "1W15_1_L4", "1W15_2_L4","1W15_3_L5",
                   "2S15_1_L5", "2S15_2_L5", "2S15_3_L5", "2W15_1_L5", "2W15_2_L5", "2W15_3_L5")  

# reorder the columns of the matrix
mat_reordered <- norm_counts_sig[, desired_order]
df <- as.data.frame(colData(dds_int)[,c("treatment", "time", "station")])

pheatmap(mat_reordered, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames=FALSE,
         annotation= df, border_color=NA, fontsize = 10, scale="row",
         fontsize_row = 10, height=20)

##### Functional enrichment #####
#extract genes remaining after filtering
filtered_genes <- rownames(dds_int)
gene2GO <- read.delim2("atonsa_trinotate_go_annotations.txt", header=FALSE)
colnames(gene2GO) <- c("gene", "GO")

#build gene universe from all genes
gene2GO <- gene2GO %>% 
  filter(gene %in% filtered_genes)

gene2GO$GO <- gsub(",", " ", gene2GO$GO)

gene2GO<- gene2GO %>%
  mutate(GO = strsplit(GO, " "))

#clusterProfiler needs unnested format
gene2GO_long <- gene2GO %>%
  unnest(GO)

GO2gene_long <- gene2GO_long %>%
  dplyr::select(gene, GO) %>%
  dplyr::select(GO, gene)


GO_ids <- unique(GO2gene_long$GO)
GO2term <- go2term(GO_ids)#associates each GO ID with term

# use GO.db database to split the output into different ontologies
go_table <- as_tibble( toTable(GOTERM), .name_repair = "universal" )
go_table$go_id <- go_table$go_id...1

go_table <- go_table %>%
  dplyr::select(go_id, Term, Ontology) %>%
  distinct()

GO_overview <- go_table %>%
  left_join(GO2gene_long, by = c("go_id" = "GO"))

# enrichment for salinity ####
res_salinity_ordered <- res_salinity[order(res_salinity$log2FoldChange, decreasing = TRUE),]
res_salinity_ordered <- na.omit(res_salinity_ordered)

list_salinity_gsea <- res_salinity_ordered$stat
#log2FoldChange
names(list_salinity_gsea) <- rownames(res_salinity_ordered)
list_salinity_gsea <- sort(list_salinity_gsea, decreasing = TRUE)

res_salinity <- subset(res_salinity, padj < 0.05)
sig_genes_salinity <- rownames(res_salinity)

#GSEA
set.seed(123)
gsea_result <- GSEA(gene=list_salinity_gsea,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    TERM2GENE = GO2gene_long,
                    TERM2NAME = GO2term)

gsea_result_tib <- as_tibble(gsea_result)

gsea_result_tib_list <- left_join(gsea_result_tib, go_table, by = c("ID" = "go_id"))
#9 GOs
GSEA_salinity <- paste0("GSEA_salinity.xlsx")
write.xlsx(file = GSEA_salinity, gsea_result_tib_list)

dotplot(gsea_result)

#ORA
ora_result <- enricher(gene = sig_genes_salinity, 
                       TERM2GENE = GO2gene_long,
                       TERM2NAME = GO2term,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH"
                       )

ora_result_tib <- as_tibble(ora_result)

ora_result_tib_list <- left_join(ora_result_tib, go_table, by = c("ID" = "go_id"))
ORA_salinity <- paste0("ORA_salinity.xlsx")
write.xlsx(file = ORA_salinity, ora_result_tib_list)
#82 GOs (23 after filtering out counts lower than 10)

dotplot(ora_result, showCategory = 20)

ID_ora <- ora_result_tib$Description
ID_gsea <- gsea_result_tib$Description

intersect(ID_ora, ID_gsea)#8 shared

# enrichment for salinity x time ####
res_salxtime_ordered <- res_salxtime[order(res_salxtime$log2FoldChange, decreasing = TRUE),]
res_salxtime_ordered <- na.omit(res_salxtime_ordered)

list_salxtime_gsea <- res_salxtime_ordered$stat
#log2FoldChange
names(list_salxtime_gsea) <- rownames(res_salxtime_ordered)
list_salxtime_gsea <- sort(list_salxtime_gsea, decreasing = TRUE)

#GSEA
set.seed(123)
gsea_result <- GSEA(gene=list_salxtime_gsea,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    TERM2GENE = GO2gene_long,
                    TERM2NAME = GO2term)

gsea_result_tib <- as_tibble(gsea_result)

gsea_result_tib_list <- left_join(gsea_result_tib, go_table, by = c("ID" = "go_id"))
#31 one GOs
GSEA_salinityxtime <- paste0("GSEA_salinityxtime.xlsx")
write.xlsx(file = GSEA_salinityxtime, gsea_result_tib_list)

dotplot(gsea_result)

#ORA
ora_result <- enricher(gene = sig_genes_salxtime, 
                       TERM2GENE = GO2gene_long,
                       TERM2NAME = GO2term,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH")

ora_result_tib <- as_tibble(ora_result)

ora_result_tib_list <- left_join(ora_result_tib, go_table, by = c("ID" = "go_id"))
#44 GOs
ORA_salinityxtime <- paste0("ORA_salinityxtime.xlsx")
write.xlsx(file = ORA_salinityxtime, ora_result_tib_list)

ID_ora <- ora_result_tib$ID
ID_gsea <- gsea_result_tib$ID

intersect(ID_ora, ID_gsea)#18 shared

# enrichment for salinity x population ####
res_salxpop_ordered <- res_salxpop[order(res_salxpop$log2FoldChange, decreasing = TRUE),]
res_salxpop_ordered <- na.omit(res_salxpop_ordered)

list_salxpop_gsea <- res_salxpop_ordered$stat
#log2FoldChange
names(list_salxpop_gsea) <- rownames(res_salxpop_ordered)
list_salxpop_gsea <- sort(list_salxpop_gsea, decreasing = TRUE)

res_salxpop <- subset(res_salxpop, padj < 0.05)
sig_genes_salxpop <- rownames(res_salxpop)

#GSEA
set.seed(123)
gsea_result <- GSEA(gene=list_salxpop_gsea,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    TERM2GENE = GO2gene_long,
                    TERM2NAME = GO2term)

#no enrichment 

#ORA
# no enrichment done for 2 significant genes

# enrichment for salinity x population x time ####
res_salxpopxtime_ordered <- res_salxpopxtime[order(res_salxpopxtime$log2FoldChange, decreasing = TRUE),]
res_salxpopxtime_ordered <- na.omit(res_salxpopxtime_ordered)

list_salxpopxtime_gsea <- res_salxpopxtime_ordered$stat
#log2FoldChange
names(list_salxpopxtime_gsea) <- rownames(res_salxpopxtime_ordered)
list_salxpopxtime_gsea <- sort(list_salxpopxtime_gsea, decreasing = TRUE)

res_salxpopxtime <- subset(res_salxpopxtime, padj < 0.05)
sig_genes_salxpopxtime <- rownames(res_salxpopxtime)

#GSEA
set.seed(123)
gsea_result <- GSEA(gene=list_salxpopxtime_gsea,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    TERM2GENE = GO2gene_long,
                    TERM2NAME = GO2term)
# no enrichment

#ORA
# no enrichment done for 1 significant gene

# enrichment for pairwise comparisons ####
# combined analysis with compare cluster for unique and shared genes per time
#unique for WH t1
WH_u_1 <- as.data.frame(salgenes_WH_u_t1)
WH_u_1$gene <- WH_u_1$salgenes_WH_u_t1
WH_u_1$station <- "WH"
WH_u_1$time <- "t1"

#unique for WH t2
WH_u_2 <- as.data.frame(salgenes_WH_u_t2)
WH_u_2$gene <- WH_u_2$salgenes_WH_u_t2
WH_u_2$station <- "WH"
WH_u_2$time <- "t2"

#unique for SW08 t1
SW08_u_1 <- as.data.frame(salgenes_SW08_u_t1)
SW08_u_1$gene <- SW08_u_1$salgenes_SW08_u_t1
SW08_u_1$station <- "SW08"
SW08_u_1$time <- "t1"

#unique for SW08 t2
SW08_u_2 <- as.data.frame(salgenes_SW08_u_t2)
SW08_u_2$gene <- SW08_u_2$salgenes_SW08_u_t2
SW08_u_2$station <- "SW08"
SW08_u_2$time <- "t2"

#shared at t1
WH_SW08_inter_1 <- as.data.frame(salgenes_inter_t1)
WH_SW08_inter_1$gene <- WH_SW08_inter_1$salgenes_inter_t1
WH_SW08_inter_1$station <- "both"
WH_SW08_inter_1$time <- "t1"

#shared at t2
WH_SW08_inter_2 <- as.data.frame(salgenes_inter_t2)
WH_SW08_inter_2$gene <- WH_SW08_inter_2$salgenes_inter_t2
WH_SW08_inter_2$station <- "both"
WH_SW08_inter_2$time <- "t2"


gene_list_t <- rbind(SW08_u_1[, c("gene", "station", "time")], 
                     SW08_u_2[, c("gene", "station", "time")], 
                     WH_u_1[, c("gene", "station", "time")],
                     WH_u_2[, c("gene", "station", "time")],
                     WH_SW08_inter_1[, c("gene", "station", "time")],
                     WH_SW08_inter_2[, c("gene", "station", "time")])

enrich_cluster <- compareCluster(gene ~ station + time,
                                 data = gene_list_t, 
                                 fun = enricher,
                                 TERM2GENE = GO2gene_long,
                                 TERM2NAME = GO2term,
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = "BH")

dev.new()
dotplot(enrich_cluster, showCategory = 100)


#subset for ontologies
#BP
enrich_cluster_tib <- as_tibble(enrich_cluster)

enrich_cluster_tib <- enrich_cluster_tib %>%
  left_join(go_table, by = c("ID" = "go_id"))
  
GO_BP <- enrich_cluster_tib[enrich_cluster_tib$Ontology == "BP", "ID", drop = TRUE]

enrich_cluster_BP <- enrich_cluster %>%
  filter(ID %in% GO_BP)

enrich_cluster_BP_tib <- as_tibble(enrich_cluster_BP)

dev.new()
dotplot(enrich_cluster_BP, showCategory = 45)+
  geom_point(aes(color=p.adjust)) + 
  scale_color_gradientn(colours=c("#DA3529","#FC9D64" ,"#FEE598", "#D7EDF4","#78A7CE", "#4575B4"),
                        guide=guide_colorbar(reverse=TRUE, order=1))+
  theme(strip.background =element_rect(fill="white"))+
  facet_wrap(~station, scales = "free_x")+
  theme(strip.text = element_text(colour = 'black', size = 12))+
  theme(axis.text.y = element_text(size = 7))

#MF
GO_MF <- enrich_cluster_tib[enrich_cluster_tib$Ontology == "MF", "ID", drop = TRUE]

enrich_cluster_MF <- enrich_cluster %>%
  filter(ID %in% GO_MF)

dev.new()
dotplot(enrich_cluster_MF, showCategory = 50)+
  geom_point(aes(color=p.adjust)) + 
  scale_color_gradientn(colours=c("#DA3529","#FC9D64" ,"#FEE598", "#D7EDF4","#78A7CE", "#4575B4"),
                        guide=guide_colorbar(reverse=TRUE, order=1))+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black', size = 12))+
  theme(axis.text.y = element_text(size = 7))

#save GO tables for each ontology

#BP
GO_BP <- enrich_cluster_tib[enrich_cluster_tib$Ontology == "BP", "ID", drop = TRUE]
enrich_cluster_BP <- enrich_cluster %>%
  filter(ID %in% GO_BP)
enrich_cluster_tib_BP <- as_tibble(enrich_cluster_BP)

#MF
GO_MF <- enrich_cluster_tib[enrich_cluster_tib$Ontology == "MF", "ID", drop = TRUE]
enrich_cluster_MF <- enrich_cluster %>%
  filter(ID %in% GO_MF)
enrich_cluster_tib_MF <- as_tibble(enrich_cluster_MF)

#CC
GO_CC <- enrich_cluster_tib[enrich_cluster_tib$Ontology == "CC", "ID", drop = TRUE]
enrich_cluster_CC <- enrich_cluster %>%
  filter(ID %in% GO_CC)
enrich_cluster_tib_CC <- as_tibble(enrich_cluster_CC)

write.csv(enrich_cluster_tib_BP, "GO_BP_cluster_pop_time.csv", row.names=FALSE, quote=TRUE)
write.csv(enrich_cluster_tib_MF, "GO_MF_cluster_pop_time.csv", row.names=FALSE, quote=TRUE)
write.csv(enrich_cluster_tib_CC, "GO_CC_cluster_pop_time.csv", row.names=FALSE, quote=TRUE)

#transcript expression
#normalize counts
norm_counts <- counts(dds_int, normalized = TRUE)

#visualize counts for all significant shared genes
norm_counts_df <- as.data.frame(norm_counts)

norm_counts_inter <- norm_counts_df[rownames(norm_counts_df) %in% salgenes_inter,]
norm_counts_inter<- norm_counts_inter %>% rownames_to_column(var = "Gene")

#turn counts into column 
norm_counts_inter_long <- norm_counts_inter %>%
  pivot_longer(
    cols = -Gene,                  # All columns except 'Gene'
    names_to = "Sample",           # New column for sample names
    values_to = "Counts"           # New column for counts
  
  )

#mutate sample name into extra time, station, treatment and replicate column
norm_counts_inter_long <- norm_counts_inter_long %>%
  mutate(
    Time = substr(Sample, 1, 1),                    # Extract Time (1st character)
    Station = substr(Sample, 2, 2),                 # Extract Station (2nd character)
    Treatment = substr(Sample, 3, 4),               # Extract Treatment (3rd and 4th characters)
    Replicate = as.integer(substr(Sample, 6, 6))    # Extract Replicate (6th character)
  )

norm_counts_inter_long$Treatment <- as.factor(norm_counts_inter_long$Treatment)
norm_counts_inter_long$Time <- as.factor(norm_counts_inter_long$Time)
norm_counts_inter_long$Station <- as.factor(norm_counts_inter_long$Station)

#plot normalized but untransformed counts
ggplot(norm_counts_inter_long, aes(x = Sample, y = Counts)) +
  geom_point() +                    
  geom_line(aes(group = Gene))+  
  theme_light()+
  facet_wrap(~Station , scales = "free_x")

ggplot(norm_counts_inter_long, aes(x = Time, y = log(Counts), color = Treatment)) +
  geom_boxplot(outlier.colour = NA)+
  facet_wrap(~Station)+
  theme_light()

##### DEGReport - pattern analysis #####
#log transformed data
dds_sig <- dds_int[salgenes_inter,]#subset the big set (154 shared genes)
dds_sig <- dds_int[salgenes_combined,]#subset the big set (602 combined genes WH and SW08)

rld <- rlog(dds_sig, blind=TRUE)

#visualize clusters
rld_mat <- assay(rld)

meta_dds <- samples
meta_dds$time <- factor(meta_dds$time, levels = c("t0", "t1", "t2"))
rownames(meta_dds) <- meta_dds$sample
meta_dds$station_treatment <- paste(meta_dds$station, meta_dds$treatment, sep = "_")

# use the `degPatterns` function show gene clusters across sample groups
#combine station and treatment to prevent overwriting of station variable
clusters_sal <- degPatterns(rld_mat, metadata = meta_dds, minc = 5, 
                            time = "time", col= "station_treatment")

clusters_sal_df <- clusters_sal$plot$data
clusters_sal_df$expression <- clusters_sal_df$value

#subset clusters to plot
clusters_sal_1 <- clusters_sal_df %>%
  filter(cluster == 1)#very high at t1, elevated t2

genes_cluster_sal_1 <- unique(clusters_sal_1$genes)
genes_cluster_1 <- gff_df %>% 
  filter(TRINITY %in% genes_cluster_sal_1)

clusters_sal_2 <- clusters_sal_df %>%
  filter(cluster == 2)#elevated at t1, very high at t2

genes_cluster_sal_2 <- unique(clusters_sal_2$genes)
genes_cluster_2 <- gff_df %>% 
  filter(TRINITY %in% genes_cluster_sal_2)

clusters_sal_13 <- clusters_sal_df %>%
  filter(cluster == 13)#overall lower for 7 PSU WH

genes_cluster_sal_13 <- unique(clusters_sal_13$genes)
genes_cluster_13 <- gff_df %>% 
  filter(TRINITY %in% genes_cluster_sal_13)

clusters_sal_14 <- clusters_sal_df %>%
  filter(cluster == 14)#similar for both

genes_cluster_sal_14 <- unique(clusters_sal_14$genes)
genes_cluster_14 <- gff_df %>% 
  filter(TRINITY %in% genes_cluster_sal_14)

clusters_sal_7 <- clusters_sal_df %>%
  filter(cluster == 7)#similar for both

genes_cluster_sal_7 <- unique(clusters_sal_7$genes)
genes_cluster_7 <- gff_df %>% 
  filter(TRINITY %in% genes_cluster_sal_7)

clusters_sal_8 <- clusters_sal_df %>%
  filter(cluster == 8)#overall higher for 7 PSU SW08

genes_cluster_sal_8 <- unique(clusters_sal_8$genes)
genes_cluster_8 <- gff_df %>% 
  filter(TRINITY %in% genes_cluster_sal_8)

clusters_sal_6 <- clusters_sal_df %>%
  filter(cluster == 6)#lower t2, more for WH

genes_cluster_sal_6 <- unique(clusters_sal_6$genes)
genes_cluster_6 <- gff_df %>% 
  filter(TRINITY %in% genes_cluster_sal_6)


clusters_sal_3 <- clusters_sal_df %>%
  filter(cluster == 3)

clusters_sal_4 <- clusters_sal_df %>%
  filter(cluster == 4)

clusters_sal_7 <- clusters_sal_df %>%
  filter(cluster == 7)

clusters_sal_16 <- clusters_sal_df %>%
  filter(cluster == 16)

clusters_sal_14 <- clusters_sal_df %>%
  filter(cluster == 14)

#plot all identified clusters
dev.new()
degPlotCluster(clusters_sal_df, time = "time", color ="station_treatment", min_genes = 10,
               process = TRUE, points = FALSE, boxes = FALSE, smooth = FALSE,
               lines = TRUE, facet = FALSE, cluster_column = "cluster")+
  facet_wrap(~title)+
  theme_light()+
  geom_boxplot(position = position_dodge(width = 0.1), alpha = 0.2) +
  geom_point( position = position_jitterdodge(dodge.width = 0.1, jitter.width = 0.5), alpha = 0.2) +
  scale_color_manual(values = c("#579fb5", "#2894b5","#e88f80" ,"#c4331b"))+
  scale_fill_manual(values = c("#579fb5", "#2894b5","#e88f80" ,"#c4331b"))



#plot only one cluster
degPlotCluster(clusters_sal_2, time = "time", color ="station_treatment", min_genes = 10,#change cluster number to see others
               process = TRUE, points = FALSE, boxes = FALSE, smooth = FALSE,
               lines = TRUE, facet = FALSE, cluster_column = "cluster")+
  theme_light()+
  geom_boxplot(position = position_dodge(width = 0.1), alpha = 0.2) +
  geom_point( position = position_jitterdodge(dodge.width = 0.1, jitter.width = 0.5), alpha = 0.2) +
  scale_color_manual(values = c("#98BBC4", "#4f7786","#BC9B9A" ,"#714b4a"))+
  scale_fill_manual(values = c("#98BBC4", "#4f7786","#BC9B9A" ,"#714b4a"))+
  ylim(c(-2.5, 2.5))

"#c4331b",  "#71ABBD"
"#BAD7E0","#e88f80"

#please note that some plots generated here have been modified in inkscape,
#final appearance may differ
