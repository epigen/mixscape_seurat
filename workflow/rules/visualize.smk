
# visualize Mixscape analysis results
rule visualize:
    input:
        mixscape_object = os.path.join(result_path,'{sample}','ALL_object.rds'),
    output:
        prtb_score_plots = report(directory(os.path.join(result_path,'{sample}','plots','PerturbScore')), 
                          patterns=["{ko}.png"],
                                  caption="../report/PerturbScores.rst", 
                          category="{}_{}".format(config["project_name"], module_name),
                          subcategory="{sample}",
                          labels={
                              "name": "{ko}",
                              "type": "Perturb Score",
                              "misc": "PNG",
                                  }),
        post_prob_plots = report(directory(os.path.join(result_path,'{sample}','plots','PosteriorProbability')), 
                          patterns=["{ko}.png"],
                                 caption="../report/PosteriorProbabilities.rst", 
                          category="{}_{}".format(config["project_name"], module_name),
                          subcategory="{sample}",
                          labels={
                              "name": "{ko}",
                              "type": "Posterior Probability",
                              "misc": "Violin",
                                  }),
        ab_expr_plots = report(directory(os.path.join(result_path,'{sample}','plots','{}_expression'.format(config["Antibody_Capture"]))),
                          patterns=["{antibody}.png"],
                               caption="../report/AntibodyExpression.rst", 
                          category="{}_{}".format(config["project_name"], module_name),
                          subcategory="{sample}",
                          labels={
                              "name": "{antibody}",
                              "type": "Expression",
                              "misc": "Violin",
                                  }) if config["Antibody_Capture"]!="" else [],
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: config.get("threads", 1)
    conda:
        "../envs/seurat_mixscape.yaml"
    log:
        os.path.join("logs","rules","visualize_{sample}.log"),
    script:
        "../scripts/visualize.R"
