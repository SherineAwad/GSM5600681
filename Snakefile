rule all:
    input:
        expand("noDoublets.h5ad")


rule preprocess:
    input:
        samples=config['samples_file']
    params:
        myObject=config["myObject"]
    output:
        expand("{myObject}.h5ad", myObject=config["myObject"])
    shell:
        "python preprocess.py {params} {input}"

rule doublets: 
     input: 
        expand("{myObject}.h5ad", myObject=config["myObject"])
     params: 
        threshold = config['doublet_threshold'] 
     output: 
        expand("noDoublets.h5ad", myObject=config["myObject"])
     shell: 
        "python detect_doublets.py {input} {params}" 

rule analyse:
    input:
        expand("noDoublets.h5ad", myObject=config["myObject"])
    output:
        expand("analysed.h5ad", myObject=config["myObject"])
    shell:
        "python analyse.py {input}"


rule cluster:
    input:
        expand("analysed.h5ad", myObject=config["myObject"])
    output:
        expand("clustered.h5ad", myObject=config["myObject"])
    shell:
        "python cluster.py {input}"


rule annotate:
    input:
        expand("clustered.h5ad", myObject=config["myObject"])
    params:
        annotations=config["annotation_file"]
    output:
        expand("annotated.h5ad", myObject=config["myObject"])
    shell:
        "python annotate.py {input} {params}"

