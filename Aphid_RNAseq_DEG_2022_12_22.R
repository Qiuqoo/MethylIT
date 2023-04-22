#============================================================================= #

#

# ================= Gene Expression of c17 fruit development ================= #

#

#============================================================================= #
library(data.table)
library( rtracklayer )
library( DESeq2 )
library( edgeR )

# ============================================================================= #
#
# ================================= /data2/users/pjw/fruit_development_RNAseq/ ===============================
#
# ============================================================================= #

# GENES <- import("/data/Work/Memory/Arabidopsis_thaliana.TAIR10.43.gff3.gz")
# GENES <- import("/data/Work/TAIR10_gff3/Arabidopsis_thaliana.TAIR10.37.gff3.gz")
#GENES <- import("/data1/users/yxd/Cs_Annotation/Cucumis_sativus.gene_with_exon.withName.gtf")
#head (GENES)
#seqlevels(GENES, pruning.mode = "coarse") <- paste(chr1:7)
# GENES <- GENES[, c("type", "gene_biotype", "gene_id" ,"gene_name")]
# 
# GENES <- GENES[ which(GENES$type == "gene") ]
Cs_C17_Genes = import("/data1/users/yxd/Cs_Annotation/Cucumis_sativus.gene_with_exon.withName.gff")
GENES = Cs_C17_Genes[ Cs_C17_Genes$type == "gene", c( "ID", "Name" ) ]
seqlevels(GENES, pruning.mode = "coarse") <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7")
# Cs_C17_TEs = import("/data1/users/yxd/Cs_Annotation/Cucumis_sativus.all_without_trf.repeat.gff")
# TEs = Cs_C17_TEs[ , c("type", "ID", "Target" ) ]
# seqlevels(TEs, pruning.mode = "coarse") <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7")
# GENES <- GENES[ which(GENES$biotype == "protein_coding") ]
# ============================ Expressed genes =============================== #
# sampleFiles <- list.files( pattern = ".geneCounts.formatted.for.DESeq.txt.gz" )
# sampleName = sub("_QC.geneCounts.formatted.for.DESeq.txt.gz","", sampleFiles )
setwd("/data4/users/Qiuqoo/Aphid_RNAseq_C17/DEGs/QoRTs_out")
# setwd("/data4/yxd/QoRTs_out/QoRTs_out")

sampleFiles <- list.files(pattern = ".QC.geneCounts.formatted.for.DESeq.txt.gz" )
sampleName <- sub( ".QC.geneCounts.formatted.for.DESeq.txt.gz" ,"", sampleFiles )
#sampleName <- gsub("[-]", "_", sampleName)

sampleCondition <- factor( c( "T06h", "T06h", "T06h",
                              "T24h", "T24h", "T24h",
                              "T72h", "T72h", "T72h",
                              "WT06h", "WT06h", "WT06h",
                              "WT24h", "WT24h", "WT24h",
                              "WT72h", "WT72h", "WT72h"),
                           levels = c("T06h","T24h", "T72h", "WT06h", "WT24h", "WT72h"))
sampleTable <- data.frame( sampleName = sampleName,
                           fileName = sampleFiles,
                           condition = sampleCondition )
condition = sampleCondition

expr_data <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                        directory = "/data4/users/Qiuqoo/Aphid_RNAseq_C17/DEGs/QoRTs_out",
                                        design= ~ condition)
# rownames(expr_data) <- gsub("gene-","",rownames(expr_data))

#View(expr_data&dds))
# idx <- match(rownames(expr_data), TEs$ID)
# sum(is.na(idx)) # 24029
#idx1 <- which(is.na(idx))
#rownames(expr_data)[idx1]
idx <- match(GENES$ID, rownames(expr_data))
sum(is.na(idx)) # 0
expr_data1 <- expr_data[idx]
# #sum(is.na(idx)) # 0
# length(TEs)
# #320780
# TEs <- TEs[idx]
# te <- data.frame(TEs)[, -4]
# head(counts(expr_data))
# dim(counts(expr_data))
#23494
exp_count <- DGEList(counts = counts(expr_data1), group = sampleCondition,
                     genes = data.frame(genes = GENES,length = width(GENES)))
exp_count$counts
write.csv(exp_count,file = "/data4/users/Qiuqoo/Aphid_RNAseq_C17/DEGs/Aphid_DEGs_counts_2022_12_22.csv")

#####  Determine which genes have sufficiently large counts to be retained in a statistical analysis.
keep <- filterByExpr(exp_count)

exp_count <- exp_count[keep, ]
# Counts per Million
exp_cpm <- cpm(exp_count)
head(exp_cpm)

# Reads per Kilobase per Million (RPKM)
exp_rpkm <- rpkm(exp_count)
head(exp_cpm)
write.csv(exp_rpkm,"/data4/users/Qiuqoo/Aphid_RNAseq_C17/DEGs/Aphid_DEGs_exp_rpkm_2022_12_22.csv")

# fragments per kilobase per million mapped fragments
rowRanges(expr_data1) <- split(GENES, as.factor(GENES$ID))
expr_data_fpkm <- fpkm(expr_data1)
write.csv(expr_data_fpkm,"/data4/users/Qiuqoo/Aphid_RNAseq_C17/DEGs/Aphid_DEGs_expr_data_fpkm_2022_12_22.csv")
# head(exp_rpkm)
# Maximizes the negative binomial likelihood to give the estimate of the common, 
# trended and tagwise dispersions across all tags.
exp_count <- estimateDisp(exp_count)
