rule all:
    input:
        expand("{myObject}.h5ad", myObject=config["myObject"]),
        expand("analysed.h5ad"), 
        expand("clustered.h5ad"),
        expand("annotated.h5ad")


rule preprocess:
    input:
        samples=config['samples_file']
    params:
        myObject=config["myObject"]
    output:
        expand("{myObject}.h5ad", myObject=config["myObject"])
    shell:
        "python preprocess.py {params} {input}"

rule analyse:
    input:
        expand("{myObject}.h5ad", myObject=config["myObject"]), 
    output:
        expand("analysed.h5ad")
    shell:
        "python analyse.py {input}"


rule cluster:
    input:
        expand("analysed.h5ad")
    output:
        expand("clustered.h5ad")
    shell:
        "python cluster.py {input}"


rule annotate:
    input:
        expand("clustered.h5ad")
    params:
        annotations=config["annotation_file"]
    output:
        expand("annotated.h5ad")
    shell:
        "python annotate.py {input} {params}"

