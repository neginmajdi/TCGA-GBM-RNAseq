library("TCGAbiolinks")
setwd("C:/Users/mitteam/Desktop/GBM/R_scripts")

# get a list of all GDC projects 
all_GDC_projects<-getGDCprojects()
View(all_GDC_projects)

# Get a summary of all GBM samples in TCGA
gbm_data <- getSampleFilesSummary("TCGA-GBM")
View(gbm_data)

# Prepare a query to access TCGA-GBM transcriptome data (STAR - Counts)
query_TCGA<-GDCquery(project = "TCGA-GBM",
                     data.category = "Transcriptome Profiling",
                     workflow.type = "STAR - Counts")

# Fetch and save the results 
query_result<- getResults(query_TCGA)
write.table(query_result, file = "all samples of GBM.txt", quote = F, sep = "\t")

# Separate tumor and normal samples from the query results
sample_list<- getResults(query_TCGA, cols="cases")

# Extract tumor ("TP") and normal ("NT") sample barcodes
sample_t <- TCGAquery_SampleTypes(barcode = sample_list, typesample = "TP")
sample_n <- TCGAquery_SampleTypes(barcode = sample_list, typesample = "NT")
length(sample_t)    #372 tumor samples
length(sample_n)    #5 normal samples


# Select first 5 tumor and 5 normal samples for comparison
sample_t_5<- sample_t[1:5]
sample_n_5<- sample_n[1:5]
selected_samples<- c(sample_n_5, sample_t_5)

View(selected_samples)
write.table(selected_samples,"selected_samples.txt", quote = F, sep = "\t", row.names = F, col.names = F)


# Download all transcriptome data for the selected samples
query_TCGA_2<- GDCquery(project = "TCGA-GBM",
                        data.category = "Transcriptome Profiling",
                        workflow.type = "STAR - Counts",
                        barcode = selected_samples,
                        data.type = "Gene Expression Quantification") 

GDCdownload(query_TCGA_2)


# Prepare the data for analysis 
data_prepare<- GDCprepare(
  query_TCGA_2, save=TRUE, save.filename="data.rda"
)

# Load SummarizedExperiment library to access assay function
library(SummarizedExperiment)

# Extract expression matrix from the prepared data
data<- assay(data_prepare)

ls()

write.table(data, "data_1.txt",quote = FALSE, sep = "\t")


library(readxl)

# Read the manually prepared annotation file
# (The original TSV was extracted from the TCGA ZIP, converted to XLS manually, 
# headers adjusted, and only necessary columns kept)
annotation<- read_excel("annotation_file.xls", col_names =TRUE)

annotation<- annotation[,1:3]
View(annotation)

data<- read.table("data_1.txt")
View(data)


annotation<- as.data.frame(annotation)
class(annotation)

rownames(annotation)<- annotation[,1]
View(annotation)
annotation<- annotation[,-1]

View(annotation)
write.table(annotation, "annotation.txt", quote=F, sep = "\t")

annotation[1,1]


# Remove rows with all zeros (genes with no expression) to clean the dataset
data_filtered<- data[which(rowSums(data)!=0),]
View(data_filtered)
nrow(data_filtered)

write.table(data_filtered, "data_filtered.txt", sep = "\t", quote = F)

# Extract gene names and gene types from the annotation file
gene_name<- annotation[rownames(data_filtered), 1]

gene_type<- annotation[rownames(data_filtered),2]


# Combine the filtered data with gene information
data_2<- data_filtered
data_2<- cbind(data_2,gene_name)
data_2<- cbind(data_2,gene_type)
data_2<- as.data.frame(data_2)

# Summarize the number of genes by RNA type
rna_table<-table(data_2$gene_type)
rna_table<- as.data.frame(rna_table)
View(rna_table)
colnames(rna_table)<- c("type", "number")
View(rna_table)

# Subset only protein-coding genes for downstream analysis
mrna_df<- subset(data_2, data_2$gene_type=="protein_coding")

write.table(mrna_df, "mrna.txt", quote = F, sep = "\t")

# Plot the distribution of RNA types
library(ggplot2)

plot <- ggplot(rna_table, aes(x = type, y = number)) +
  geom_bar(stat = "identity", fill = "yellow") +
  labs(title = "Barplot",
       x = "RNA type",
       y = "Number") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave("data_rna_info.png", plot, width = 15, height = 7)



selected_samples<- read.table("selected_samples.txt")
View(selected_samples)

# Create a vector indicating sample type
sample_type<- c(rep("Normal", 5), rep("Tumor",5))


# Combine sample names and types into metadata
meta_data<- data.frame(selected_samples,sample_type)
View(meta_data)

colnames(meta_data)<- c("sample_name","sample_type")

write.table(meta_data,"meta_data.txt", quote = F, sep = "\t", row.names = F)

meta_data_2<- read.table("meta_data.txt", sep="\t",header = T)
View(meta_data_2)

mrna_df<- read.table("mrna.txt")
mrna_df<- as.data.frame(mrna_df)
View(mrna_df)

mrna_raw<- mrna_df[,-11:-12]
View(mrna_raw)

# Convert all columns to numeric
for (i in 1: ncol(mrna_raw)){
  mrna_raw[,i]<- as.numeric(mrna_raw[,i])
}

all(colnames(mrna_raw) == rownames(meta_data))

View(mrna_raw)

class(mrna_raw)

write.table(mrna_raw, "mrna_raw.txt", quote = F, sep = "\t" )

nrow(mrna_raw)
View(meta_data)


selected_samples<- read.table("selected_samples.txt")
selected_samples<- as.data.frame (selected_samples)
meta_data<- read.table("meta_data.txt", sep = "\t", header = T) 
meta_data<- as.data.frame(meta_data)

mrna_raw<- read.table("mrna_raw.txt")
mrna_raw<- as.data.frame(mrna_raw)

nrow(meta_data)
nrow(mrna_raw)

library(DESeq2)

# Create DESeq2 dataset using raw counts and metadata
dd1<- DESeqDataSetFromMatrix(
  countData=mrna_raw,
  colData=meta_data,
  design=~sample_type   
)

# Perform rlog transformation 
data_rlog<- rlogTransformation(d1)
data_rlog_df<- data.frame(assay(data_rlog))

# Estimate size factors 
d2<- estimateSizeFactors(d1)

# Run DESeq2 differential expression analysis
d3<- DESeq(d2)

# Extract results with FDR-adjusted p-values
deseq_result<- results(d3,pAdjustMethod = "fdr")
final_result<- data.frame(deseq_result)
write.table(final_result,"DEG.txt", quote = F, sep="\t")
saveRDS(final_result, "DEG.rds")

deg <- readRDS("C:/Users/mitteam/Desktop/GBM/Data/DEG.rds")

View(deg)


#Annotate DEGs with gene names

final_result<- deg
class(final_result)
annot_df<- read.table("annotation.txt", sep = "\t")

class(annot_df)
colnames(annot_df)<- c("gene_name","Bio_type")


# Function to extract gene name from annotation dataframe
annot_fx <- function(x) {
  k <- annot_df[x, 1]
  return(k)
}

annot_df[1,]

gene_name<- annot_fx(rownames(final_result))
length(gene_name)== nrow(final_result)

# Combine gene names with DESeq2 results
final_result_annot<- cbind(gene_name)
final_result_annot<- cbind(final_result_annot, final_result[1:6])

rownames(final_result_annot)<- rownames(final_result)
colnames(final_result_annot)[1]<- "symbol"
View(final_result_annot)

final_result_annot<- final_result_annot[order(final_result_annot$padj, decreasing=F), ]
View(final_result_annot)

write.table(final_result_annot, "final_result_annot.txt", quote=F, sep="\t")


#heatmap plot
library(pheatmap)
library(ggplot2)

# Use rlog-transformed data for visualization
normdata<- data_rlog_df

# Select top 100 DEGs based on adjusted p-value
top_genes_count <- data_rlog_df[rownames(final_result_annot)[1:100], ]

# Assign gene symbols as rownames
gene_names_2<- annot_fx(rownames(top_genes_count))
rownames(top_genes_count)<- gene_names_2
View(top_genes_count)

#rename columns to reflect sample type
colnames(top_genes_count)<- c("normal 1", "normal 2","normal 3", "normal 4","normal 5",
                              "tumor 1", "tumor 2", "tumor 3", "tumor 4","tumor 5")
View(top_genes_count)



# Select top 60 DEGs for heatmap visualization
top_60_genes <- top_genes_count[1:60, ]

# Generate heatmap for top 60 DEGs
library(RColorBrewer)

heat_plot <- pheatmap(top_60_genes,
                      color = colorRampPalette(c("green", "black", "red"))(256),
                      scale = "row",
                      cellwidth = 6,
                      cellheight = 6,   
                      fontsize_row = 6,
                      fontsize_col = 6,
                      main = "Top 60 DEGs",
                      show_rownames = TRUE)  



png("heatmap_plot.png",width = 1500, height = 15)
heat_plot
dev.off()

#volcano plot
library(EnhancedVolcano)
final_result_annot<- read.table("final_result_annot.txt")
final_result_annot<- as.data.frame(final_result_annot)
colnames(final_result_annot)

# Select top 10 genes to highlight in the plot
top_genes<- final_result_annot[1:10,1]

my_plot <- EnhancedVolcano(
  final_result_annot,
  title = "Gene volcano plot",
  lab = final_result_annot$symbol,
  x = "log2FoldChange",
  y = "pvalue",
  selectLab = top_genes,  
  xlab = bquote(~log[2]~"fold change"),
  pCutoff = 0.01,
  pointSize = 1.5,  
  labSize = 4.0,
  labCol = "black",
  labFace = "bold",
  boxedLabels = TRUE,
  colAlpha = 0.5,  
  legendPosition = "right",
  legendLabSize = 14,
  legendIconSize = 4.0,
  drawConnectors = TRUE,  
  colConnectors = "black"
)

ggsave("volcano_plot.png",plot=my_plot, width = 12, height = 10)


#pca plot 
norm_data<- log2(top_genes_count+1)
pc<- prcomp(norm_data)
pcr<- data.frame(pc$rotation[,1:3],meta_data$sample_type)
group<- meta_data$sample_type

library(ggplot2)
pca_plot<- ggplot(pcr,aes(PC1, PC2,color=group))+
  geom_point(size=3)
View(pca_plot)

ggsave("pca_plot.png", plot = pca_plot, width = 7, height = 6)

