
# visualize Mixscape analysis results
rule visualize:
    input:
        mixscape_object = os.path.join(result_path,'{sample}','MIXSCAPE_ALL_object.rds'),
    output:
        prtb_score_plot = report(expand(os.path.join(result_path,'{{sample}}','plots','MIXSCAPE_ALL_PerturbScores_{split_by}.png'), split_by=['',config["CalcPerturbSig"]["split_by_col"], config["RunMixscape"]["split_by_col"]]), 
                          caption="../report/PerturbScores.rst", 
                          category="{}_mixscape_seurat".format(config["project_name"]), 
                          subcategory="{sample}"),
        post_prob_plot = report(expand(os.path.join(result_path,'{{sample}}','plots','MIXSCAPE_ALL_PosteriorProbabilities_{split_by}.png'), split_by=['',config["CalcPerturbSig"]["split_by_col"], config["RunMixscape"]["split_by_col"]]), 
                          caption="../report/PosteriorProbabilities.rst", 
                          category="{}_mixscape_seurat".format(config["project_name"]), 
                          subcategory="{sample}"),
        ab_expr_plot = report(os.path.join(result_path,'{sample}','plots','MIXSCAPE_ALL_{}_expression.png'.format(config["Antibody_Capture"])),
                          caption="../report/AntibodyExpression.rst", 
                          category="{}_mixscape_seurat".format(config["project_name"]), 
                          subcategory="{sample}") if config["Antibody_Capture"]!="" else "",
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat.yaml"
    log:
        os.path.join("logs","rules","visualize_{sample}.log"),
    params:
        partition=config.get("partition"),
        CalcPerturbSig_params = config["CalcPerturbSig"],
        RunMixscape_params = config["RunMixscape"],
        antibody_capture_flag = config["Antibody_Capture"],
    script:
        "../scripts/visualize.R"