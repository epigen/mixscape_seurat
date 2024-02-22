#### load libraries & utility function 
library("Seurat")
library("ggplot2")

# source utility functions
# source("workflow/scripts/utils.R")
snakemake@source("./utils.R")

# inputs
filtered_object_path <- snakemake@input[[1]]

# outputs
mixscape_object_path <- snakemake@output[["mixscape_object"]]
stat_plots_path <- snakemake@output[["stat_plots"]]
prtb_data_path <- snakemake@output[["prtb_data"]]
mixscape_stats_path <- snakemake@output[["mixscape_stats"]]

# parameters
assay <- snakemake@config[["assay"]]
variable_features_only <- snakemake@config[["variable_features_only"]] 
calcPerturbSig_params <- snakemake@config[["CalcPerturbSig"]]
runMixscape_params <- snakemake@config[["RunMixscape"]]
grna_split_symbol <- snakemake@config[["grna_split_symbol"]]

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
print("Calculating perturbation signature (PRTB)")
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
print("Running Mixscape")
data <- RunMixscape(
  object = data,
  assay = "PRTB",
  slot = "scale.data",
  labels = calcPerturbSig_params[["gene_col"]],
  nt.class.name = calcPerturbSig_params[["nt_term"]],
  new.class.name = "mixscape_class",
  min.de.genes = runMixscape_params[["min_de_genes"]],
  min.cells = runMixscape_params[["min_cells"]],
  de.assay = "RNA", # TODO: why not configured assay?
  logfc.threshold = runMixscape_params[["lfc_th"]],
  iter.num = 10,
  verbose = TRUE,
  split.by = mixscape_split_by,
  fine.mode = as.logical(runMixscape_params[["mixscape_fine_mode"]]),
  fine.mode.labels = runMixscape_params[["grna_col"]],
  prtb.type = runMixscape_params[["prtb_type"]]
)


### Plot Mixscape Statistics
print("Plot Mixscape statistics")
# plot specs
# width <- 15
# height <- 2*ceiling((dim(unique(data[[calcPerturbSig_params[["gene_col"]]]]))[1]-1)/10)

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
                           
# p1 <- ggplot(df3, aes(x = guide_number, y = value*100, fill= Var1)) +
#   geom_bar(stat= "identity") +
#   theme_classic()+
#   scale_fill_manual(values = c("grey79","coral1")) + 
#   ylab("% of cells") +
#   xlab("sgRNA")

# p1 <- p1 + theme(axis.text.x = element_text(size = 10, hjust = 1), 
#            axis.text.y = element_text(size = 10), 
#            axis.title = element_text(size = 8), 
#            strip.text = element_text(size=8, face = "bold")) + 
#   facet_wrap(vars(gene),ncol = 10, scales = "free") +
#   labs(fill = "mixscape class") +theme(legend.title = element_text(size = 12),
#           legend.text = element_text(size = 12))

# p1

for (gene in unique(df3$gene)){
   tmp_plot <- ggplot(df3[df3$gene==gene,], aes(x = guide_number, y = value*100, fill= Var1)) +
    geom_bar(stat= "identity") +
    theme_classic()+
    scale_fill_manual(values = c("grey79","coral1")) +
    ylab("% of cells") +
    xlab("sgRNA") +
    ggtitle(gene) +
    theme(axis.text.x = element_text(size = 10), 
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 10),
          strip.text = element_text(size = 10, face = "bold"),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
         plot.title = element_text(hjust = 0.5)) +
    labs(fill = "mixscape class")
    
    ggsave_new(filename = gene, 
           results_path=dirname(stat_plots_path), 
           plot=tmp_plot, 
           width=4, 
           height=4)
}
                           
# options(repr.plot.width=4, repr.plot.height=4)

### save results
print("Saving results")
# save seurat object and metadata
save_seurat_object(seurat_obj=data, 
                   result_dir=dirname(mixscape_object_path), 
                   prefix="ALL_")
 
# save mixscape statistics
fwrite(as.data.frame(t(stat_table)), file=file.path(mixscape_stats_path), row.names=TRUE)
                           
# save matrix of PRTB values; slots counts and scale.data are emtpy -> only save data slot
fwrite(as.data.frame(GetAssayData(object = data, slot = "data", assay = "PRTB")), file=file.path(prtb_data_path), row.names=TRUE)
