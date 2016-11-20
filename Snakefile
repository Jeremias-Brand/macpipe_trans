rule trim_and_trinity:
    input:
        forward="/home/jeremias/dev/macpipe_trans/data/{sample}_R1_subsample.fastq.gz",
        reverse="/home/jeremias/dev/macpipe_trans/data/{sample}_R2_subsample.fastq.gz"
    output:
        assembly="{sample}.trinity.Trinity.fasta"
        #trinity_name="{sample}.trinity"
    params: adapters="/home/jeremias/soft/Trimmomatic-0.36/adapters/Truseq_barcodes.fa"
    threads: 36
    message: "Executing somecommand with {threads} threads on the following files {input}."
    shell:  "/home/jeremias/soft/trinityrnaseq-2.2.0/Trinity --trimmomatic --CPU {threads} --output {wildcards.sample}'.trinity' --left {input.forward}  --right {input.reverse} \
--quality_trimming_params 'ILLUMINACLIP:{params.adapters}:2:40:15 LEADING:2 TRAILING:2 MINLEN:25'"
