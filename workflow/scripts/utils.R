### Utility Functions
library("data.table")

# save Seurat object results
save_seurat_object <- function (seurat_obj, result_dir, prefix){
    
    # save seurat object
    saveRDS(seurat_obj, file=file.path(result_dir, paste0(prefix,"object.rds")))
    
    # save metadata
    fwrite(as.data.frame(seurat_obj[[]]), file=file.path(result_dir, paste0(prefix,"metadata.csv")), row.names=TRUE)
    
    # save stats
    stats <- paste0("cells: ",ncol(seurat_obj),"\nfeatures: ",nrow(seurat_obj))
    write(stats, file=file.path(result_dir, paste0(prefix,"stats.txt")))
}


# extended ggsave
ggsave_new <- function(filename, results_path, plot, width=5, height=5){
    
    # failsafe for large plots
    width <- min(100, width)
    height <- min(100, height)
    
    for (format in c('png')){
        ggsave(
          filename = paste0(filename,'.',format),
          plot = plot,
          device = format,
          path = file.path(results_path),
          scale = 1,
          dpi = 300,
            width = width,
            height = height,
          limitsize = FALSE,
            units="in"
        )
    }
}


# custom theme settings
font <- "Arial"
size <- 8

custom_theme <- theme_bw(
        base_size=size,
        base_family = font
    ) %+replace% 
    
    theme(
      #text elements
      plot.title = element_text(             #title
                   family = font,            #set font family
                   size = size,                #set font size
                   face = 'bold',            #bold typeface
                   hjust = 0.5,                #center align
                   vjust = 2),               #raise slightly
      
      plot.subtitle = element_text(          #subtitle
                   family = font,            #font family
                   size = size),               #font size
      
      plot.caption = element_text(           #caption
                   family = font,            #font family
                   size = size,                 #font size
                   hjust = 0.5),               #center align
      
      axis.title = element_text(             #axis titles
                   family = font,            #font family
                   size = size),               #font size
      
      axis.text = element_text(              #axis text
                   family = font,            #axis famuly
                   size = size),                #font size
    )