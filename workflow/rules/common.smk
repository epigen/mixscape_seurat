##### utility functions #####

def get_sample_paths(wildcards):
    return annot.loc[wildcards.sample,'data']
