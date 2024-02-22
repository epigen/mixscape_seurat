#### load libraries & utility function 
library("Seurat")
library("ggplot2")
library("patchwork")

# source utility functions
# source("workflow/scripts/utils.R")
snakemake@source("./utils.R")

# inputs
mixscape_object_path <- snakemake@input[["mixscape_object"]]

# outputs
prtb_score_plot_path <- snakemake@output[["prtb_score_plots"]]
post_prob_plot_path <- snakemake@output[["post_prob_plots"]]
ab_expr_plot_path <- snakemake@output[["ab_expr_plots"]]

# parameters
calcPerturbSig_params <- snakemake@config[["CalcPerturbSig"]]
runMixscape_params <- snakemake@config[["RunMixscape"]]
antibody_capture_flag <- snakemake@config[["Antibody_Capture"]]

# make plot directories
dir.create(post_prob_plot_path, recursive = TRUE)
dir.create(prtb_score_plot_path, recursive = TRUE)
if (antibody_capture_flag != "") dir.create(ab_expr_plot_path, recursive = TRUE)

### load mixscape data
data <- readRDS(file = file.path(mixscape_object_path))
Idents(data) <- "mixscape_class"

### Visualize Mixscape analysis results

# plot specifications
width <- 5
height <- 3

# Explore the perturbation scores of cells.
for (split_by in c("default", calcPerturbSig_params[["split_by_col"]], runMixscape_params[["split_by_col"]])){
    if (split_by==''){
        next
    }
    
    if (split_by=="default"){
        split_by <- NULL
    }

    perturb_scores <- list()

    for (gene in unique(unlist((data[[calcPerturbSig_params[["gene_col"]]]])))){
        
        # skip the gene if NT, only KO or NP cells exist, or one group (KO,NP) have only one cell
        if (gene==calcPerturbSig_params[["nt_term"]] | length(unique(data$mixscape_class.global[data[[calcPerturbSig_params[["gene_col"]]]]==gene]))==1 | sum(table(data$mixscape_class.global[data[[calcPerturbSig_params[["gene_col"]]]]==gene])==1)>0){
            next
        }

        # plot perturbation scores for each target gene
        tmp_perturb_scores <- PlotPerturbScore(
          object = data,
          target.gene.class = calcPerturbSig_params[["gene_col"]],
          target.gene.ident = gene,
          mixscape.class = "mixscape_class",
          col = "orange2",
          split.by = split_by,
          before.mixscape = FALSE,
          prtb.type = runMixscape_params[["prtb_type"]]
        ) + custom_theme + theme(legend.title = element_blank())
        
        ggsave_new(filename = paste0(gene,"_",split_by), 
               results_path=prtb_score_plot_path, 
               plot=tmp_perturb_scores, 
               width=width, 
               height=height)
    }
}


# Inspect the posterior probability values in NP and KO cells.
for (split_by in c("default", calcPerturbSig_params[["split_by_col"]], runMixscape_params[["split_by_col"]])){
    if (split_by==''){
        next
    }
    
    if (split_by=="default"){
        split_by <- NULL
    }
    post_probs <- list()

    for (gene in unique(unlist((data[[calcPerturbSig_params[["gene_col"]]]])))){

        if (gene==calcPerturbSig_params[["nt_term"]]){
            next
        }

        # posterior probability values for each target gene
        tmp_p <- VlnPlot(
            object = data, 
            features = "mixscape_class_p_ko", 
            idents = c(calcPerturbSig_params[["nt_term"]], unique(data[["mixscape_class"]][data[[calcPerturbSig_params[["gene_col"]]]]==gene])),
            split.by=split_by
        ) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.text = element_text(size = 8),
          legend.text = element_text(size=8),
          plot.title = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank())

        if (is.null(split_by)){
            tmp_p <- tmp_p + NoLegend()
        }else{
            tmp_p <- tmp_p 
        }
        
        ggsave_new(filename = paste0(gene,"_",split_by), 
               results_path=post_prob_plot_path,
               plot=tmp_p, 
               width=width, 
               height=height)
        
    }
}

# check if antibody_capture_flag is present
if (antibody_capture_flag != ""){
    
    # check if data is normalized, if not normalize the data
    if (all(GetAssayData(object = data, slot = "counts", assay = antibody_capture_flag) == GetAssayData(object = data, slot = "data", assay = antibody_capture_flag))){
        data <- NormalizeData(data, normalization.method = "CLR", margin = 2, assay = antibody_capture_flag)
    }

    # plot surface protein expression measurements split by mixscape class
    DefaultAssay(object = data) <- antibody_capture_flag
    Idents(object = data) <- calcPerturbSig_params[["gene_col"]]
    Idents(object = data) <- factor(x = Idents(data), levels = unique(c(calcPerturbSig_params[["nt_term"]],sort(levels(data)))))

    features <- rownames(GetAssayData(data, slot = "data", assay = antibody_capture_flag))
    
    width <- 0.5 * length(unique(unlist(data[[calcPerturbSig_params[["gene_col"]]]])))

    for (feature in features){
        tmp_violin_plot <- VlnPlot(object = data,
            features = feature,
            cols = c("coral3","grey79","grey39"),
            pt.size = 0,
            idents = NULL,
            sort = FALSE,
            assay = antibody_capture_flag,
            group.by = NULL,
            split.by = "mixscape_class.global",
            adjust = 1,
            y.max = NULL,
            same.y.lims = FALSE,
            log = FALSE,
            ncol = NULL,
            slot = "data",
            split.plot = FALSE,
            stack = FALSE,
            combine = FALSE,
            fill.by = "feature",
            flip = FALSE
           )[[1]] + theme(axis.text = element_text(size = 8),
                     legend.text = element_text(size=8),
                     axis.title.y = element_text(size=8),
                     axis.title.x = element_blank())
        
        ggsave_new(filename = feature, 
               results_path = ab_expr_plot_path,
               plot = tmp_violin_plot, 
               width = width, 
               height = height)
    }
}
