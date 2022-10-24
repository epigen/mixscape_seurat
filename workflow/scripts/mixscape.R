#### load libraries & utility function 
library(Seurat)
library(ggplot2)

# source utility functions
# source("workflow/scripts/utils.R")
snakemake@source("./utils.R")

# inputs
filtered_object_path <- snakemake@input[[1]] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/scrnaseq_processing_seurat/condition_24h_cytokines/FILTERED_object.rds"

# outputs
mixscape_object_path <- snakemake@output[["mixscape_object"]] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/mixscape_seurat/condition_24h_cytokines/MIXSCAPE_ALL_object.rds" 
mixscape_plot_path <- snakemake@output[["mixscape_plot"]] 

# parameters
assay <- snakemake@params[["assay"]] #"SCT" #"RNA"
variable_features_only <- snakemake@params[["variable_features_only"]] 
calcPerturbSig_params <- snakemake@params[["CalcPerturbSig_params"]]
runMixscape_params <- snakemake@params[["RunMixscape_params"]]
grna_split_symbol <- snakemake@params[["grna_split_symbol"]]

if (calcPerturbSig_params[["split_by_col"]]==''){
    perturbSig_split_by <- NULL
}else{
    perturbSig_split_by <- calcPerturbSig_params[["split_by_col"]]
}

if (runMixscape_params[["split_by_col"]]==''){
    mixscape_split_by <- NULL
}else{
    mixscape_split_by <- runMixscape_params[["split_by_col"]]
}

result_dir <- dirname(mixscape_object_path)
# make directories if not exist
if (!dir.exists(result_dir)){
        dir.create(result_dir, recursive = TRUE)
    }

### load filtered data
data <- readRDS(file = file.path(filtered_object_path))
DefaultAssay(object = data) <- assay

### Prepare assay for dimensionality reduction: Normalize data, find variable features and scale data.

# check if data is normalized, if not normalize the data
if (all(GetAssayData(object = data, slot = "counts", assay = assay) == GetAssayData(object = data, slot = "data", assay = assay))){
    data <- NormalizeData(object = data)
}

# set features for the analysis
if (variable_features_only==1){
    # check if data has variable features, if not determine them
    if (length(VariableFeatures(object = data))==0){
        data <- FindVariableFeatures(object = data)
    }
    features <- VariableFeatures(object = data)
}else{
    features <- rownames(data)
}

# check if data is scaled, if not scale data
if (dim(GetAssayData(object = data, slot = "scale.data", assay = assay))[1]==0){
    data <- ScaleData(object = data, features = features)
}

# Run Principle Component Analysis (PCA) to reduce the dimensionality of the data.
data <- RunPCA(object = data, features = features)

### Mixscape Analysis

# Calculate perturbation signature (PRTB).
# https://satijalab.org/seurat/reference/calcperturbsig
data <- CalcPerturbSig(
  object = data,
  assay = assay,
  features = features,
  slot = "data",
  gd.class = calcPerturbSig_params[["gene_col"]],
  nt.cell.class = calcPerturbSig_params[["nt_term"]],
  split.by = perturbSig_split_by,
  num.neighbors = calcPerturbSig_params[["n_neighbors"]],
  reduction = "pca",
  ndims = calcPerturbSig_params[["ndims"]],
  new.assay.name = "PRTB",
  verbose = TRUE
)

# Run mixscape.
# https://satijalab.org/seurat/reference/runmixscape
data <- RunMixscape(
  object = data,
  assay = "PRTB",
  slot = "scale.data",
  labels = calcPerturbSig_params[["gene_col"]],
  nt.class.name = calcPerturbSig_params[["nt_term"]],
  new.class.name = "mixscape_class",
  min.de.genes = runMixscape_params[["min_de_genes"]],
  min.cells = runMixscape_params[["min_cells"]],
  de.assay = "RNA",
  logfc.threshold = runMixscape_params[["lfc_th"]],
  iter.num = 10,
  verbose = FALSE,
  split.by = mixscape_split_by,
  fine.mode = as.logical(runMixscape_params[["mixscape_fine_mode"]]),
  fine.mode.labels = runMixscape_params[["grna_col"]],
  prtb.type = runMixscape_params[["prtb_type"]]
)


### Plot Mixscape Statistics

# plot specs
width <- 15
height <- 2*ceiling((dim(unique(data[[calcPerturbSig_params[["gene_col"]]]]))[1]-1)/10)

# Calculate percentage of KO cells for all target gene classes.
stat_table <- table(data$mixscape_class.global, unlist(data[[calcPerturbSig_params[["grna_col"]]]]))
df <- prop.table(stat_table,2)

df2 <- reshape2::melt(df)
df2$Var2 <- as.character(df2$Var2)
test <- df2[which(df2$Var1 == "KO"),]
test <- test[order(test$value, decreasing = T),]
new.levels <- test$Var2
df2$Var2 <- factor(df2$Var2, levels = new.levels )
df2$Var1 <- factor(df2$Var1, levels = c(calcPerturbSig_params[["nt_term"]], "NP", "KO"))
df2$gene <- sapply(as.character(df2$Var2), function(x) head(strsplit(x, split = grna_split_symbol)[[1]],1))
df2$guide_number <- sapply(as.character(df2$Var2), function(x) tail(strsplit(x, split = grna_split_symbol)[[1]],1)) 
# remove NT
df3 <- df2[-c(which(df2$gene == calcPerturbSig_params[["nt_term"]])),]
df3$Var1 <- factor(df3$Var1, levels = c("NP", "KO"))
df3 <- df3[!is.na(df3$Var1),]
                           
p1 <- ggplot(df3, aes(x = guide_number, y = value*100, fill= Var1)) +
  geom_bar(stat= "identity") +
  theme_classic()+
  scale_fill_manual(values = c("grey79","coral1")) + 
  ylab("% of cells") +
  xlab("sgRNA")

p1 <- p1 + theme(axis.text.x = element_text(size = 10, hjust = 1), 
           axis.text.y = element_text(size = 10), 
           axis.title = element_text(size = 8), 
           strip.text = element_text(size=8, face = "bold")) + 
  facet_wrap(vars(gene),ncol = 10, scales = "free") +
  labs(fill = "mixscape class") +theme(legend.title = element_text(size = 12),
          legend.text = element_text(size = 12))

# p1
                           
ggsave_new(filename = "MIXSCAPE_ALL_stats", 
           results_path=dirname(mixscape_plot_path), 
           plot=p1, 
           width=width, 
           height=height)


### save results
# save seurat object and metadata
save_seurat_object(seurat_obj=data, 
                   result_dir=dirname(mixscape_object_path), 
                   prefix="MIXSCAPE_ALL_")
 
# save mixscape statistics
write.csv(t(stat_table), file=file.path(result_dir, paste0('MIXSCAPE_ALL_stats.csv')), row.names=TRUE)
                           
# save matrix of PRTB values
slot <- "data"                           
write.csv(GetAssayData(object = data, slot = slot, assay = "PRTB"), file=file.path(result_dir, paste0('MIXSCAPE_ALL_PRTB_',slot,'.csv')), row.names=TRUE)
  
                           
# slots counts and scale.data are emtpy -> only save data slot
# for (slot in c("counts","data","scale.data")){
#     write.csv(GetAssayData(object = data, slot = slot, assay = "PRTB"), file=file.path(result_dir, paste0('MIXSCAPE_ALL_PRTB_',slot,'.csv')), row.names=TRUE)
# }