
# perform Mixscape analysis
rule mixscape:
    input:
        get_sample_paths
    output:
        mixscape_object = os.path.join(result_path,'{sample}','MIXSCAPE_ALL_object.rds'),
        mixscape_data = os.path.join(result_path,'{sample}','MIXSCAPE_ALL_PRTB_data.csv'),
        mixscape_metadata = os.path.join(result_path,'{sample}','MIXSCAPE_ALL_metadata.csv'),
        mixscape_plot = report(os.path.join(result_path,'{sample}','plots','MIXSCAPE_ALL_stats.png'), 
                          caption="../report/mixscape.rst", 
                          category="{}_mixscape_seurat".format(config["project_name"]), 
                          subcategory="{sample}"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","mixscape_{sample}.log"),
    params:
        partition=config.get("partition"),
        assay = config["assay"],
        variable_features_only = config["variable_features_only"],
        CalcPerturbSig_params = config["CalcPerturbSig"],
        RunMixscape_params = config["RunMixscape"],
        grna_split_symbol = config["grna_split_symbol"],
    script:
        "../scripts/mixscape.R"
      
# perform LDA on perturbed subset
rule lda:
    input:
        mixscape_object = os.path.join(result_path,'{sample}','MIXSCAPE_ALL_object.rds'),
    output:
        lda_object = os.path.join(result_path,'{sample}','MIXSCAPE_FILTERED_LDA_object.rds'),
        lda_data = os.path.join(result_path,'{sample}','MIXSCAPE_FILTERED_LDA_data.csv'),
        filtered_prtb_data = os.path.join(result_path,'{sample}','MIXSCAPE_FILTERED_PRTB_data.csv'),
        filtered_metadata = os.path.join(result_path,'{sample}','MIXSCAPE_FILTERED_LDA_metadata.csv'),
        lda_plot = report(os.path.join(result_path,'{sample}','plots','MIXSCAPE_FILTERED_LDA_UMAP.png'), 
                          caption="../report/lda_umap.rst", 
                          category="{}_mixscape_seurat".format(config["project_name"]), 
                          subcategory="{sample}"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","lda_{sample}.log"),
    params:
        partition=config.get("partition"),
        assay = config["assay"],
        CalcPerturbSig_params = config["CalcPerturbSig"],
        RunMixscape_params = config["RunMixscape"],
        MixscapeLDA_params = config["MixscapeLDA"],
    script:
        "../scripts/lda.R"