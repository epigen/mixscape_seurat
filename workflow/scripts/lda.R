#### load libraries & utility function 
library(Seurat)
library(ggplot2)
library(scales)

# source utility functions
# source("workflow/scripts/utils.R")
snakemake@source("./utils.R")

# inputs
mixscape_object_path <- snakemake@input[["mixscape_object"]] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/mixscape_seurat/condition_24h_cytokines/MIXSCAPE_ALL_object.rds" 

# outputs
lda_object_path <- snakemake@output[["lda_object"]] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/mixscape_seurat/condition_24h_cytokines/MIXSCAPE_FILTERED_LDA_object.rds" 
lda_plot_path <- snakemake@output[["lda_plot"]] 

# parameters
assay <- snakemake@params[["assay"]] #"SCT" #"RNA"
calcPerturbSig_params <- snakemake@params[["CalcPerturbSig_params"]]
runMixscape_params <- snakemake@params[["RunMixscape_params"]]
mixscapeLDA_params <- snakemake@params[["MixscapeLDA_params"]]

result_dir <- dirname(lda_object_path)
# make directories if not exist
if (!dir.exists(result_dir)){
        dir.create(result_dir, recursive = TRUE)
    }

### load mixscape data
data <- readRDS(file = file.path(mixscape_object_path))
DefaultAssay(object = data) <- assay

# Remove non-perturbed cells
Idents(data) <- "mixscape_class.global"
sub <- subset(data, idents = c(runMixscape_params[["prtb_type"]], calcPerturbSig_params[["nt_term"]]))

### perform Linear Discriminant Analysis (LDA)
# run LDA to reduce the dimensionality of the data
# https://satijalab.org/seurat/reference/mixscapelda
sub <- MixscapeLDA(
  object = sub,
  assay = assay,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.key = "LDA_",
  seed = 42,
  pc.assay = "PRTB",
  labels = calcPerturbSig_params[["gene_col"]],
  nt.label = calcPerturbSig_params[["nt_term"]],
  npcs = mixscapeLDA_params[["npcs"]],
  verbose = TRUE,
  logfc.threshold = runMixscape_params[["lfc_th"]]
)

lda_data <- Embeddings(object = sub, reduction = "lda")

### Visualize results

# Use LDA results to run UMAP and visualize cells in 2-D
# https://satijalab.org/seurat/reference/runumap
sub <- RunUMAP(
  object = sub,
  dims = 1:ncol(lda_data),#(length(unique(sub$mixscape_class))-1),
  reduction = 'lda',
  reduction.key = 'ldaumap',
  reduction.name = 'ldaumap')

# plot UMAP
width <- 10
height <- 10


# if only 3 classes remain, then LDA projection is already 2D, no UMAP necessary
# if ((length(unique(sub$mixscape_class))-1)==2){
if (ncol(lda_data)==2){
    reduction <- 'lda'
    x_label <- "LDA 1"
    y_label <- "LDA 2"
}else{
    reduction <- 'ldaumap'
    x_label <- "UMAP 1"
    y_label <- "UMAP 2"
}

# Visualize UMAP clustering results.
Idents(sub) <- "mixscape_class"
sub$mixscape_class <- as.factor(sub$mixscape_class)

# Set colors for each perturbation.
col = setNames(object = hue_pal()(length(unique(sub$mixscape_class))),nm = unique(sub$mixscape_class))
col[[calcPerturbSig_params[["nt_term"]]]] <- "#D3D3D3"

p <- DimPlot(object = sub, 
             reduction = reduction, 
             repel = T, 
             label.size = 4, 
             label = T, 
             cols = col,
             pt.size=0.1,
             label.box=T) 

p2 <- p+ 
  scale_color_manual(values=col, drop=FALSE) + 
  ylab(y_label) +
  xlab(x_label) +
  custom_theme + NoLegend()

# p2

ggsave_new(filename = "MIXSCAPE_FILTERED_LDA_UMAP", 
           results_path=dirname(lda_plot_path), 
           plot=p2, 
           width=width, 
           height=height)


### save results
# save seurat object and metadata
save_seurat_object(seurat_obj=sub, 
                   result_dir=dirname(lda_object_path), 
                   prefix="MIXSCAPE_FILTERED_LDA_")

# save matrix of LDA values
write.csv(lda_data, file=file.path(result_dir, paste0('MIXSCAPE_FILTERED_LDA_data.csv')), row.names=TRUE)
                           
# save matrix of PRTB values
slot <- "data"                           
write.csv(GetAssayData(object = sub, slot = slot, assay = "PRTB"), file=file.path(result_dir, paste0('MIXSCAPE_FILTERED_PRTB_',slot,'.csv')), row.names=TRUE)
