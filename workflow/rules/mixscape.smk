
# perform Mixscape analysis
rule mixscape:
    input:
        get_sample_paths
    output:
        mixscape_object = os.path.join(result_path,'{sample}','ALL_object.rds'),
        metadata = os.path.join(result_path,'{sample}','ALL_metadata.csv'),
        prtb_data = os.path.join(result_path,'{sample}','ALL_PRTB_data.csv'),
        mixscape_stats = os.path.join(result_path,'{sample}','mixscape_stats.csv'),
        stat_plots = report(directory(os.path.join(result_path,'{sample}','plots','stats')),
                            patterns=["{ko}.png"],
                          caption="../report/mixscape.rst", 
                          category="{}_{}".format(config["project_name"], module_name),
                          subcategory="{sample}",
                            labels={
                                "name": "{ko}",
                                "type": "Statistic",
                                "misc": "PNG",
                                  }),
    params:
        utils_path = workflow.source_path("../scripts/utils.R"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: 8*config.get("threads", 1)
    conda:
        "../envs/seurat_mixscape.yaml"
    log:
        os.path.join("logs","rules","mixscape_{sample}.log"),
    script:
        "../scripts/mixscape.R"

# perform LDA on perturbed subset
rule lda:
    input:
        mixscape_object = os.path.join(result_path,'{sample}','ALL_object.rds'),
    output:
        lda_object = os.path.join(result_path,'{sample}','FILTERED_object.rds'),
        filtered_metadata = os.path.join(result_path,'{sample}','FILTERED_metadata.csv'),
        lda_data = os.path.join(result_path,'{sample}','LDA_data.csv'),
        filtered_assay_data = os.path.join(result_path,'{sample}','FILTERED_{}_data.csv'.format(config["assay"])),
        filtered_prtb_data = os.path.join(result_path,'{sample}','FILTERED_PRTB_data.csv'),
        lda_plot = report(os.path.join(result_path,'{sample}','plots','LDA_UMAP.png'), 
                          caption="../report/lda_umap.rst", 
                          category="{}_{}".format(config["project_name"], module_name),
                          subcategory="{sample}",
                          labels={
                              "name": "LDA",
                              "type": "UMAP",
                              "misc": "PNG",
                                  }),
    params:
        utils_path = workflow.source_path("../scripts/utils.R"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat_lda.yaml"
    log:
        os.path.join("logs","rules","lda_{sample}.log"),
    script:
        "../scripts/lda.R"