# This is a template script to run trinity
# Using a JSON file as input will be better
# Sample is a wildcard equivalent to .+

rule trim_and_trinity:
    input:
        forward="/home/jeremias/dev/macpipe_trans/data/{sample}_R1_subsample.fastq.gz",
        reverse="/home/jeremias/dev/macpipe_trans/data/{sample}_R2_subsample.fastq.gz"
    output:
        "{sample}.Trinity.fasta"
        #trinity_name="{sample}.trinity"
    params: adapters="/home/jeremias/soft/Trimmomatic-0.36/adapters/Truseq_barcodes.fa"
    threads: 36 # threads only works if --cores is set to the actual number of cores when running the snakemake
    message: "Executing somecommand with {threads} threads on the following files {sample}."
    # We also want a logfile
    log: "logs/{sample}.trinity.log"
    shell:  "/home/jeremias/soft/trinityrnaseq-2.2.0/Trinity --seqType fq --trimmomatic --CPU {threads} --max_memory 240G \
            --output {wildcards.sample}'.trinity' \
            --left {input.forward}  \
            --right {input.reverse} \
            --quality_trimming_params 'ILLUMINACLIP:{params.adapters}:2:40:15 LEADING:2 TRAILING:2 MINLEN:25' && \
            mv {wildcards.sample}.trinity/Trinity.fasta {wildcards.sample}.Trinity.fasta"

rule run_Transdecoder:
    input:
        "{sample}.Trinity.fasta"
    output:
        "{sample}.Trinity.fasta.transdecoder_dir/longest_orfs.pep"
    shell:
        "TransDecoder.LongOrfs -t {wildcards.sample}.Trinity.fasta"


rule creat_gene_mapping:
    input:
        "{sample}.Trinity.fasta"
    output:
        "{sample}.Trinity.fasta.gene_trans_map"
    shell:
        "/home/jeremias/soft/trinityrnaseq-2.2.0/util/support_scripts/get_Trinity_gene_to_trans_map.pl {wildcards.sample}.Trinity.fasta \
        > {wildcards.sample}.Trinity.fasta.gene_trans_map"


rule build_db:
    input:
        geneMap="{sample}.Trinity.fasta.gene_trans_map",
        assembly="{sample}.Trinity.fasta",
        transPep="{sample}.Trinity.fasta.transdecoder_dir/longest_orfs.pep"
    output:
        "{sample}.Trinotate.sqlite"
    shell:
        """
        cp /home/jeremias/soft/Trinotate-3.0.1/Trinotate.sqlite ./{wildcards.sample}.Trinotate.sqlite && \
        Trinotate {wildcards.sample}.Trinotate.sqlite init --gene_trans_map {input.geneMap} \
        --transcript_fasta {input.assembly}  --transdecoder_pep {input.transPep}
        """

rule run_trinotate:
    input:
        sqlDb="{sample}.Trinotate.sqlite",
        geneMap="{sample}.Trinity.fasta.gene_trans_map",
        assembly="{sample}.Trinity.fasta"
    output:
        "{sample}.Trinity.fasta.transdecoder.gff3",
        "{sample}.Trinity.fasta.rnammer.gff"
        #"chkpt.BLASTX_SPROT_PEP.ok",
        #"chkpt.BLASTX_SPROT_TRANS.ok",
        #"chkpt.PFAM.ok",
        #"chkpt.RNAMMER.ok",
        #"chkpt.SIGNALP.ok",
        #"chkpt.TMHMM.ok",
        #"chkpt.TRANSDECODER_LONGORF.ok",
        #"chkpt.TRANSDECODER_PREDICT.ok"
    threads: 36
    shell:
        """
        autoTrinotate.pl --CPU {threads} --conf ~/soft/Trinotate-3.0.1/config.trinotate.txt \
        --Trinotate_sqlite {input.sqlDb} \
        --transcripts {input.assembly} \
        --gene_to_trans_map {input.geneMap}
        """

# rule annotate_db:
#     input:

#     output:
#         "{sample}trinotate_annotation_report.xls"
#     shell:
#         """
#         Trinotate Trinotate.sqlite LOAD_swissprot_blastp swissprot.blastp.outfmt6
#         Trinotate Trinotate.sqlite LOAD_swissprot_blastx swissprot.blastx.outfmt6
#         Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
#         Trinotate Trinotate.sqlite LOAD_tmhmm tmhmm.out
#         Trinotate Trinotate.sqlite LOAD_signalp signalp.out
#         Trinotate Trinotate.sqlite report -E 0.0001 > {wildcards.sample}trinotate_annotation_report.xls
#         """