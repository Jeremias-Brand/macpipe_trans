# This is a template script to run trinity
# Using a JSON file as input will be better
# Sample is a wildcard equivalent to .+ 

rule trim_and_trinity:
    input:
        forward="/home/jeremias/dev/macpipe_trans/data/{sample}_R1_subsample.fastq.gz",
        reverse="/home/jeremias/dev/macpipe_trans/data/{sample}_R2_subsample.fastq.gz"
    output:
        assembly="{sample}.trinity.Trinity.fasta"
        #trinity_name="{sample}.trinity"
    params: adapters="/home/jeremias/soft/Trimmomatic-0.36/adapters/Truseq_barcodes.fa"
    threads: 36 # threads only works if --cores is set to the actual number of cores when running the snakemake
    message: "Executing somecommand with {threads} threads on the following files {input}."
    shell:  "/home/jeremias/soft/trinityrnaseq-2.2.0/Trinity --trimmomatic --CPU {threads} --output {wildcards.sample}'.trinity' --left {input.forward}  --right {input.reverse} \
--quality_trimming_params 'ILLUMINACLIP:{params.adapters}:2:40:15 LEADING:2 TRAILING:2 MINLEN:25'"
