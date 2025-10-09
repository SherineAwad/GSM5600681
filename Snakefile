rule all:
    input:
        expand(
            "annotated_clustered_analysed_{myObject}.h5ad",
            myObject=config["myObject"]
        )


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
        expand("{myObject}.h5ad", myObject=config["myObject"])
    output:
        expand("analysed_{myObject}.h5ad", myObject=config["myObject"])
    shell:
        "python analyse.py {input}"


rule cluster:
    input:
        expand("analysed_{myObject}.h5ad", myObject=config["myObject"])
    output:
        expand("clustered_analysed_{myObject}.h5ad", myObject=config["myObject"])
    shell:
        "python cluster.py {input}"


rule annotate:
    input:
        expand("clustered_analysed_{myObject}.h5ad", myObject=config["myObject"])
    params:
        annotations=config["annotation_file"]
    output:
        expand("annotated_clustered_analysed_{myObject}.h5ad", myObject=config["myObject"])
    shell:
        "python annotate.py {input} {params}"

