#### load libraries & utility function 
library(Seurat)
library(ggplot2)
library(patchwork)

# source utility functions
# source("workflow/scripts/utils.R")
snakemake@source("./utils.R")

# inputs
mixscape_object_path <- snakemake@input[["mixscape_object"]] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/mixscape_seurat/condition_24h_cytokines/MIXSCAPE_ALL_object.rds" 

# outputs
post_prob_plot_path <- snakemake@output[["post_prob_plot"]] #"/nobackup/lab_bock/projects/macroIC/results/AKsmall/mixscape_seurat/condition_24h_cytokines/plots/MIXSCAPE_ALL_PerturbScores_.png" 

# parameters
calcPerturbSig_params <- snakemake@params[["CalcPerturbSig_params"]]
runMixscape_params <- snakemake@params[["RunMixscape_params"]]
antibody_capture_flag <- snakemake@params[["antibody_capture_flag"]]

# for testing
# calcPerturbSig_params <- list()
# calcPerturbSig_params[["split_by_col"]] <- "batch"
# calcPerturbSig_params[["gene_col"]] <- "KOcall"
# calcPerturbSig_params[["nt_term"]] <- "NonTargeting"
# runMixscape_params <- list()
# runMixscape_params[["split_by_col"]] <- "condition"
# runMixscape_params[["prtb_type"]] <- "KO"
# antibody_capture_flag <- "AB"

### load mixscape data
data <- readRDS(file = file.path(mixscape_object_path))
Idents(data) <- "mixscape_class"

### Visualize Mixscape analysis results

# plot specifications
n_col <- 5
width <- 5
height <- 3
width_panel <- n_col * width

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
        perturb_scores[[gene]] <- PlotPerturbScore(
          object = data,
          target.gene.class = calcPerturbSig_params[["gene_col"]],
          target.gene.ident = gene,
          mixscape.class = "mixscape_class",
          col = "orange2",
          split.by = split_by,
          before.mixscape = FALSE,
          prtb.type = runMixscape_params[["prtb_type"]]
        ) + custom_theme + theme(legend.title = element_blank())
    }

    width_panel <- n_col * width
    height_panel <- ceiling(length(perturb_scores)/n_col) * height

    perturb_scores_panel <- wrap_plots(perturb_scores, ncol = n_col)

    # save plot
#     options(repr.plot.width=width_panel, repr.plot.height=height_panel)
#     print(perturb_scores_panel)
    ggsave_new(filename = paste0("MIXSCAPE_ALL_PerturbScores_",split_by), 
               results_path=dirname(post_prob_plot_path), 
               plot=perturb_scores_panel, 
               width=width_panel, 
               height=height_panel)
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
            post_probs[[gene]] <- tmp_p + NoLegend()
        }else{
            post_probs[[gene]] <- tmp_p 
        }
        
    }
    height_panel <- ceiling(length(post_probs)/n_col) * height

    post_probs_panel <- wrap_plots(post_probs, ncol = n_col)

    # save plot
#     options(repr.plot.width=width_panel, repr.plot.height=height_panel)
#     print(post_probs_panel)
    ggsave_new(filename = paste0("MIXSCAPE_ALL_PosteriorProbabilities_",split_by), 
               results_path=dirname(post_prob_plot_path), 
               plot=post_probs_panel, 
               width=width_panel, 
               height=height_panel)
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

    violin_plots <- VlnPlot(object = data,
            features = features,
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
           )

    # adapt theme for each plot          
    violin_plots <- lapply(violin_plots, function(tmp_plot) tmp_plot + theme(axis.text = element_text(size = 8),
                                                                             legend.text = element_text(size=8),
                                                                             axis.title.y = element_text(size=8),
                                                                             axis.title.x = element_blank()))             

    violin_plots_panel <- wrap_plots(violin_plots, ncol = 1)

    height_panel <- height * length(features)
    width_panel <- 0.5 * length(unique(unlist(data[[calcPerturbSig_params[["gene_col"]]]])))

    # save plot
#     options(repr.plot.width=width_panel, repr.plot.height=height_panel)
#     print(violin_plots_panel)
    ggsave_new(filename = paste0("MIXSCAPE_ALL_",antibody_capture_flag,"_expression"), 
               results_path=dirname(post_prob_plot_path), 
               plot=violin_plots_panel, 
               width=width_panel, 
               height=height_panel)
}