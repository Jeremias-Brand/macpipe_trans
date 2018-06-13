"""
Author: J. Brand
Affiliation: Unibas
Aim: Workflow for transcriptome assembly, annotation and assessment
Date: 20. Nov. 2016
Run: snakemake   -s Snakefile
Latest modification:
  - add more flexible input
  - add reporting
  - add checks for gzip
  - integrate gzip of intermediate results
  - add a block to run if fail or succeed
  - add evironment file for bioconda
  - plotting of busco automated
  - move transcriptomes to better places
  - add trinity version to software.versions file
"""
# import statements
from snakemake.utils import report
from os.path import join
import os
import time


# Uses a yaml file as input
configfile: "config.yaml"


# Preparation------------------------------------------------------------------
#TIMESTAMP = time.strftime("%Y%m%d")
if config["TIMESTAMP"]:
    TIMESTAMP = config["TIMESTAMP"]
else:
    TIMESTAMP = time.strftime("%Y%m%d")
# check if the necessary dirs exist and if not creates them
def save_mkdir( dirs ):
	for d in dirs:
		if not os.path.isdir(d):
			os.mkdir(d)
			print("Creating directory: " + d)

def check_trailing_slash( path ):
    if path[-1] != "/":
        path = path + "/"
    return path

dirs = ["logs", "logs/trinity", "logs/transrate",
        "logs/rcorrector", "trinity", "transrate",
        "transrate/summary", "busco", "logs/busco",
        "expression", "expression/salmon", "expression/RSEM",
        "transdecoder", "logs/transdecoder", "mapping", "logs/trimmomatic",
        "logs/pfam"]
save_mkdir(dirs)


# Globals ---------------------------------------------------------------------

# define global things here

# Full path to a folder that holds all of your FASTQ files.
FASTQ_DIR = check_trailing_slash(config["fastqdir"])
# loads all of the path into variable
TRANSRATE_DIR = check_trailing_slash(config["transrateDir"])
HOME_DIR = check_trailing_slash(config["homedir"])
VOUCHER_DIR = check_trailing_slash(config["voucher_dir"])


# A Snakemake regular expression matching the forward mate FASTQ files.
# we get a list with the names of all the files in the working directory
SAMPLES, = glob_wildcards(join(FASTQ_DIR, '{sample}_R1.fq.gz'))

# Patterns for the 1st mate and the 2nd mate using the 'sample' wildcard.
PATTERN_R1 = '{sample}_R1.fq.gz'
PATTERN_R2 = '{sample}_R2.fq.gz'



# Rules -----------------------------------------------------------------------



rule all:
    input:
        "report_" + TIMESTAMP + ".html"

# This is a template script to run trinity
# Using a JSON file as input will be better
# Sample is a wildcard equivalent to .+
# rule All: # this top rule is here to define the goal of the script
#     input:
#         "{sample}trinotate_annotation_report.xls"
#     log:"snakemake.report.log"

# TODO software version printout and link to report

rule run_rcorrector:
    input:
        forward = join(FASTQ_DIR, PATTERN_R1),
        reverse = join(FASTQ_DIR, PATTERN_R2)
    output:
        FASTQ_DIR + "{sample}_R1.cor.fq.gz",
        FASTQ_DIR + "{sample}_R2.cor.fq.gz"
    log:
    	"logs/rcorrector/{sample}.log"
    benchmark:
        "benchmarks/{sample}.run_rcorrector.txt"
    params:
        outdir = FASTQ_DIR
    threads: 14
    shell:
        """
        ~/perl5/perlbrew/perls/5.16.2t/bin/perl ~/bin/run_rcorrector.pl -k 31 -t {threads} \
        -od {params.outdir}  -1 {input.forward}  -2 {input.reverse} &> {log}
        """

rule trimmomatic:
    input:
        forward  = FASTQ_DIR + "{sample}_R1.cor.fq.gz",
        reverse  = FASTQ_DIR + "{sample}_R2.cor.fq.gz"
    output:
        forward_trimmed  = FASTQ_DIR + "{sample}_R1_cor_trim.fq.gz",
        reverse_trimmed  = FASTQ_DIR + "{sample}_R2_cor_trim.fq.gz"
    log: 
        "logs/trimmomatic/{sample}.log"
    benchmark:
        "benchmarks/{sample}.trimmomatic.txt"
    params:
        sample_name = FASTQ_DIR + "{sample}",
        adapters = config["adapters_fasta"],
        trimmomatic_string = config["trimmomatic_string"]
    threads: 14
    shell:
        # this step also renames the files
        """
        trimmomatic PE -threads {threads}\
        {input.forward} {input.reverse} \
        {params.sample_name}_R1_cor_trim.fq.gz {params.sample_name}_R1_cor_unpaired.fq.gz \
        {params.sample_name}_R2_cor_trim.fq.gz {params.sample_name}_R2_cor_unpaired.fq.gz \
        ILLUMINACLIP:{params.adapters}{params.trimmomatic_string} &> {log}
        """

rule remove_rRNA:
    input:
        forward  = FASTQ_DIR + "{sample}_R1_cor_trim.fq.gz",
        reverse  = FASTQ_DIR + "{sample}_R2_cor_trim.fq.gz"
    output:
        # here we want the misses because the hits are most likely rRNA
        frRNA  = FASTQ_DIR + "{sample}_rRNA.1.gz",
        rrRNA  = FASTQ_DIR + "{sample}_rRNA.2.gz",
        fclean  = FASTQ_DIR + "{sample}_clean.1.gz",
        rclean  = FASTQ_DIR + "{sample}_clean.2.gz"
    log: "logs/bowtie_{sample}.log"
    benchmark:
        "benchmarks/{sample}.remove_rRNA.txt"
    params:
        rRNAdb = config["rRNAdb"],
        statfile = "bowtie_stats",
        sample_name = FASTQ_DIR + "{sample}",
        err = "logs/bowtie_{sample}.sam"
    threads: 14
    shell:
        """
        (bowtie2 --very-sensitive-local --phred33 -q  \
        -x {params.rRNAdb} -1 {input.forward}  -2 {input.reverse} \
        --threads {threads} --met-file {params.statfile} \
        -S {params.err} \
        --al-conc-gz {params.sample_name}_rRNA.gz --un-conc-gz {params.sample_name}_clean.gz) 2> {log}
        """

rule rename_remove_rRNA:
    input:
        frRNA  = FASTQ_DIR + "{sample}_rRNA.1.gz",
        rrRNA  = FASTQ_DIR + "{sample}_rRNA.2.gz",
        fclean  = FASTQ_DIR + "{sample}_clean.1.gz",
        rclean  = FASTQ_DIR + "{sample}_clean.2.gz"
    output:
        # here we want the misses because the hits are most likely rRNA
        forward_clean  = FASTQ_DIR + "{sample}_R1_cor_trim_clean.fq.gz",
        reverse_clean  = FASTQ_DIR + "{sample}_R2_cor_trim_clean.fq.gz",
        forward_rRNA  = FASTQ_DIR + "{sample}_R1_cor_trim_rRNA.fq.gz",
        reverse_rRNA  = FASTQ_DIR + "{sample}_R2_cor_trim_rRNA.fq.gz"
    log: "logs/bowtie"
    shell:
        """
        mv {input.frRNA}   {output.forward_rRNA}
        mv {input.rrRNA}   {output.reverse_rRNA}
        mv {input.fclean}  {output.forward_clean}
        mv {input.rclean}  {output.reverse_clean}
        """

rule fastqc:
    input:
        f  = FASTQ_DIR + "{sample}_R1.fq.gz",
        r  = FASTQ_DIR + "{sample}_R2.fq.gz",
        f_trim  = FASTQ_DIR + "{sample}_R1_cor_trim.fq.gz",
        r_trim  = FASTQ_DIR + "{sample}_R2_cor_trim.fq.gz",
        f_clean  = FASTQ_DIR + "{sample}_R1_cor_trim_clean.fq.gz",
        r_clean  = FASTQ_DIR + "{sample}_R2_cor_trim_clean.fq.gz"
    output:
        f_qc         = "QC/{sample}_R1_fastqc.html",
        r_qc         = "QC/{sample}_R2_fastqc.html",
        f_trim       = "QC/{sample}_R1_cor_trim_fastqc.html",
        r_trim       = "QC/{sample}_R2_cor_trim_fastqc.html",
        f_trimmed_qc = "QC/{sample}_R1_cor_trim_clean_fastqc.html",
        r_trimmed_qc = "QC/{sample}_R2_cor_trim_clean_fastqc.html"
    log:
        "logs/fastqc/{sample}_fastqc.log"
    benchmark:
        "benchmarks/{sample}.fastqc.txt"
    params:
    threads: 6
    shell:
        """
        fastqc {input.f} {input.r} \
        {input.f_trim} {input.r_trim} \
        {input.f_clean} {input.r_clean} \
        --outdir QC/  -t {threads} &> /dev/null
        for f in {output};do
            echo $f
            file=${{f%%.html}}
            z=$file".zip"
            unzip -qq -d ./QC/ $z
            scripts/collect_fastqcData.sh $file
        done
        """

# this gives just a rough estimate of the mapping rate
rule map_snap:
    input:
        f_clean  = FASTQ_DIR + "{sample}_R1_cor_trim_clean.fq.gz",
        r_clean  = FASTQ_DIR + "{sample}_R2_cor_trim_clean.fq.gz"
    output:
        expand("mapping/{{sample}}_{T}_snap_to_ref.txt", T=TIMESTAMP)
    benchmark:
        "benchmarks/{sample}.map_snap.txt"
    threads: 14
    params:
        map_yn = config["map_yn"],
        reference_snap = config["reference_snap"]
    run:
        if str(params.map_yn) == "y":
            shell("snap-aligner paired  {params.reference_snap}  {input.f_clean} {input.r_clean} -I -d 25 -t {threads} 1> {output} 2> /dev/null")
        else:
            with open(str(output), "a") as f:
                f.write("Mapping to reference was not performed.\n")


rule trinity_on_rRNA:
    input:
        forward  = FASTQ_DIR + "{sample}_R1_cor_trim_rRNA.fq.gz",
        reverse  = FASTQ_DIR + "{sample}_R2_cor_trim_rRNA.fq.gz"
    output:
        assembly = expand("rRNA/{{sample}}_{T}_rRNA.Trinity.fasta", T=TIMESTAMP),
        log      = expand("logs/rRNA/{{sample}}_{T}_trinity.shell.log", T=TIMESTAMP)
    log:
        expand("logs/rRNA/{{sample}}_{T}.log", T=TIMESTAMP)
    benchmark:
        "benchmarks/{sample}.trinity_on_rRNA.txt"
    params:
        outdir   =  FASTQ_DIR,
        trinity  =  config["trinity"],
        stranded_string =  config["stranded_string"]
    threads: 14
    run:
        try:
            shell("""{params.trinity} --seqType fq \
        --CPU {threads} --max_memory 150G \
        --output 'trinity/'{wildcards.sample}_rRNA_{TIMESTAMP}'.trinity' \
        --left   {input.forward}  \
        --right  {input.reverse} \
        --normalize_reads  --normalize_max_read_cov 30 \
        {params.stranded_string} > ./logs/trinity/{wildcards.sample}_rRNA_{TIMESTAMP}_trinity.shell.log 2> {log} && \
        cp trinity/{wildcards.sample}_rRNA_{TIMESTAMP}.trinity/Trinity.fasta {output.assembly} && \
        sed -i 's/>TRINITY/>{wildcards.sample}_rRNA_{TIMESTAMP}/' {output.assembly} && \ 
        sed -E 's/^(>[^ ]+) ([^ ]+) .+$/\1 \2/' {output.assembly}
        find  ./trinity/{wildcards.sample}_rRNA_{TIMESTAMP}.trinity/ -type f -not -name 'Trinity.fasta' -not -name 'Trinity.fasta.gene_trans_map' -print0 | xargs -0 rm
        cd trinity/{wildcards.sample}_rRNA_{TIMESTAMP}.trinity/
        rm -rdf chrysalis insilico_read_normalization read_partitions""")
        except:
        # We need to convert the variables to a sting because otehrwise they are of type <class 'snakemake.io.Namedlist'>
            print("Trinity run on rRNA did not return an assembly!\n This could be because there were not enough input reads\n")
            open(str(output.assembly), 'a').close()
            with open(str(output.log), "a") as log_f:
                log_f.write("Trinity run did not return an assembly.\n")


rule trinity_on_trimmed:
    input:
        forward  = FASTQ_DIR + "{sample}_R1_cor_trim_clean.fq.gz",
        reverse  = FASTQ_DIR + "{sample}_R2_cor_trim_clean.fq.gz"
    output:
        assembly = expand("trinity/{{sample}}_{T}.Trinity.fasta", T=TIMESTAMP),
        log      = expand("logs/trinity/{{sample}}_{T}.shell.log", T=TIMESTAMP)
    log:
        expand("logs/trinity/{{sample}}_{T}.log", T=TIMESTAMP)
    benchmark:
        "benchmarks/{sample}.trinity_on_trimmed.txt"
    params:
        outdir   =  FASTQ_DIR,
        stranded_string =  config["stranded_string"],
        trinity  =  config["trinity"]
    threads: 14  # threads only works if --cores is set to the actual number of cores when running the snakemake
    # message: expand("Executing with {threads} threads on the following files {sample}.", sample=SAMPLES)
    # We also want a logfile
    # log: expand("logs/{sample}.trinity.log", sample=SAMPLES)
    shell:
        """
        {params.trinity} --seqType fq \
        --CPU {threads} --max_memory 150G \
        --output 'trinity/'{wildcards.sample}_{TIMESTAMP}'.trinity' \
        --left   {input.forward}  \
        --right  {input.reverse} \
        --normalize_reads  --normalize_max_read_cov 30 \
        {params.stranded_string} > ./logs/trinity/{wildcards.sample}_{TIMESTAMP}.shell.log 2> {log} && \
        cp trinity/{wildcards.sample}_{TIMESTAMP}.trinity/Trinity.fasta {output.assembly} && \
        sed -i 's/>TRINITY/>{wildcards.sample}_{TIMESTAMP}/' {output.assembly} 
        find  ./trinity/{wildcards.sample}_{TIMESTAMP}.trinity/ -type f -not -name 'Trinity.fasta' -not -name 'Trinity.fasta.gene_trans_map' -print0 | xargs -0 rm
        cd trinity/{wildcards.sample}_{TIMESTAMP}.trinity/
        rm -rdf chrysalis insilico_read_normalization read_partitions
        """

rule basic_assembly_stats:
    input:
        assembly = expand("trinity/{{sample}}_{T}.Trinity.fasta", T=TIMESTAMP),
        forward  = FASTQ_DIR + "{sample}_R1_cor_trim_clean.fq.gz",
        reverse  = FASTQ_DIR + "{sample}_R2_cor_trim_clean.fq.gz"
    output:
        basic = expand("assembly_stats/{{sample}}_{T}_N50.txt", T=TIMESTAMP),
        salmon_quant = expand("expression/salmon/{{sample}}_{T}/quant.sf", T=TIMESTAMP),
        # RSEM_quant = expand("expression/RSEM/{{sample}}_{T}/RSEM.isoforms.results", T=TIMESTAMP),
        ExN50  = expand("assembly_stats/{{sample}}_{T}_ExN50.txt", T=TIMESTAMP),
        out_salmon  = expand("expression/salmon/{{sample}}_{T}/", T=TIMESTAMP)
        # out_RSEM  = expand("expression/RSEM/{{sample}}_{T}/", T=TIMESTAMP)
    threads: 14
    params:
        trinity_home = config["trinity_home"]
    log:
        expand("logs/basic_assembly_stats/{{sample}}_{T}.log", T=TIMESTAMP)
    benchmark:
        "benchmarks/{sample}.basic_assembly_stats.txt"
    
    shell:
        """
        {params.trinity_home}'util/TrinityStats.pl' {input.assembly} > {output.basic}
        
        {params.trinity_home}'util/align_and_estimate_abundance.pl' --thread_count {threads} --prep_reference \
        --transcripts {input.assembly} --seqType fq --left {input.forward} --right {input.reverse} \
        --est_method salmon --aln_method bowtie2 --trinity_mode  --output_dir {output.out_salmon} &>> {log}


        {params.trinity_home}'util/misc/contig_ExN50_statistic.pl' {output.salmon_quant} {input.assembly} | tee {output.ExN50} &> /dev/null
        
        R/plot_ExN50_statistic.Rscript {output.ExN50} 
        """
        
#        {params.trinity_home}'util/align_and_estimate_abundance.pl' --thread_count {threads} --prep_reference \
#        --transcripts {input.assembly} --seqType fq --left {input.forward} --right {input.reverse} \
#        --est_method RSEM --aln_method bowtie2 --trinity_mode --output_dir {output.out_RSEM} &>> {log}


rule run_transdecoder_longORFs:
    input:
        assembly = expand("trinity/{{sample}}_{T}.Trinity.fasta", T=TIMESTAMP)
    output:
        longest_pep = expand("transdecoder/{{sample}}_{T}.Trinity.fasta.transdecoder_dir/longest_orfs.pep", T=TIMESTAMP),
        longest_cds = expand("transdecoder/{{sample}}_{T}.Trinity.fasta.transdecoder_dir/longest_orfs.cds", T=TIMESTAMP),
        length_pep  = expand("transdecoder/{{sample}}_{T}_longest_orfs.pep.length", T=TIMESTAMP),
        length_cds  = expand("transdecoder/{{sample}}_{T}_longest_orfs.cds.length", T=TIMESTAMP)
    threads: 1
    params:
        T = TIMESTAMP,
        H = HOME_DIR
    benchmark:
        "benchmarks/{sample}.run_transdecoder.txt"
    log:
        expand("logs/transdecoder/{{sample}}_{T}_longORFs.log", T=TIMESTAMP)
    shell:
        """
        cd transdecoder
        pyenv which TransDecoder.LongOrfs
        TransDecoder.LongOrfs -m 100 -t {params.H}{input.assembly} &> ../{log} 
        python {params.H}scripts/seq_length.py {params.H}{output.longest_pep} > {params.H}{output.length_pep} 
        python {params.H}scripts/seq_length.py {params.H}{output.longest_cds} > {params.H}{output.length_cds}
        """



rule run_pfam:
    input:
        longest_pep = expand("transdecoder/{{sample}}_{T}.Trinity.fasta.transdecoder_dir/longest_orfs.pep", T=TIMESTAMP),
        longest_cds = expand("transdecoder/{{sample}}_{T}.Trinity.fasta.transdecoder_dir/longest_orfs.cds", T=TIMESTAMP),
        length_pep  = expand("transdecoder/{{sample}}_{T}_longest_orfs.pep.length", T=TIMESTAMP),
        length_cds  = expand("transdecoder/{{sample}}_{T}_longest_orfs.cds.length", T=TIMESTAMP)
    output:
        pfam = expand("transdecoder/{{sample}}_{T}.pfam.domtblout", T=TIMESTAMP)
    threads: 14
    params:
        T = TIMESTAMP,
        H = HOME_DIR,
        pfam = config["pfam"]
    benchmark:
        "benchmarks/{sample}.run_pfam.txt"
    log:
        expand("logs/pfam/{{sample}}_{T}.pfam.log", T=TIMESTAMP)
    shell:
        """
        hmmscan --cpu {threads} --domtblout {output.pfam} {params.pfam} {input.longest_pep} &> {log}
        """


rule run_transdecoder_predict:
    input:
        assembly = expand("trinity/{{sample}}_{T}.Trinity.fasta", T=TIMESTAMP),
        pfam = expand("transdecoder/{{sample}}_{T}.pfam.domtblout", T=TIMESTAMP)
    output:
        pep = expand("transdecoder/{{sample}}_{T}.transdecoder.pep", T=TIMESTAMP),
        cds = expand("transdecoder/{{sample}}_{T}.transdecoder.cds", T=TIMESTAMP),
        pep_length = expand("transdecoder/{{sample}}_{T}.transdecoder.pep.length", T=TIMESTAMP),
        cds_length = expand("transdecoder/{{sample}}_{T}.transdecoder.cds.length", T=TIMESTAMP)
    threads: 1
    params:
        T = TIMESTAMP,
        H = HOME_DIR
    benchmark:
        "benchmarks/{sample}.run_transdecoder.txt"
    log:
        expand("logs/transdecoder/{{sample}}_{T}_predict.log", T=TIMESTAMP)
    shell:
        """
        cd transdecoder
        TransDecoder.Predict  -t {params.H}{input.assembly} --retain_pfam_hits {params.H}{input.pfam} --single_best_only  &> ../{log} 
        cd ..
        mv transdecoder/{wildcards.sample}_{params.T}.Trinity.fasta.transdecoder.pep {output.pep}
        mv transdecoder/{wildcards.sample}_{params.T}.Trinity.fasta.transdecoder.cds {output.cds}
        python {params.H}scripts/seq_length.py {output.pep} > {output.pep_length} 
        python {params.H}scripts/seq_length.py {output.cds} > {output.cds_length}
        """
    


rule run_transrate:
    input:
        assembly =  expand("trinity/{{sample}}_{T}.Trinity.fasta", T=TIMESTAMP),
        forward  =  FASTQ_DIR + "{sample}_R1.fq.gz",
        reverse  =  FASTQ_DIR + "{sample}_R2.fq.gz"
    output:
        expand("transrate/{{sample}}_{T}.transrate/{{sample}}_{T}.Trinity/good.{{sample}}_{T}.Trinity.fasta", T=TIMESTAMP),
        expand("transrate/{{sample}}_{T}.transrate/assemblies.csv", T=TIMESTAMP),
        expand("transrate/{{sample}}_{T}.transrate/{{sample}}_{T}.Trinity/contigs.csv", T=TIMESTAMP)
    log:
        expand("logs/transrate/{{sample}}_{T}.log", T=TIMESTAMP)
    benchmark:
        "benchmarks/{sample}.transrate.txt"
    params:
        transDir  =  TRANSRATE_DIR,
        dataDir   =  FASTQ_DIR,
        homeDir   =  HOME_DIR,
        transrate =  config["transrate"],
        transrate_cores = config["transrate_cores"]
    threads: 7
    shell:
        """
        {params.transrate} --assembly {params.homeDir}{input.assembly} \
        --output {params.transDir}{wildcards.sample}_{TIMESTAMP}.transrate \
        --threads {params.transrate_cores} \
        --left {input.forward} \
        --right {input.reverse} > logs/transrate/{wildcards.sample}.transrate.log 2> {log}
        """


rule gather_transrate:
    # move all transrate results to same folder for easier summary later
    # First I used the same wildcard input output names for this rule as I did
    # for the gather busco rule but that somehow caused it not to be plotted in the
    # dag. More of a cosmetic issue.
    input:
        a = expand("transrate/{{sample}}_{T}.transrate/assemblies.csv", T=TIMESTAMP),
        c  = expand("transrate/{{sample}}_{T}.transrate/{{sample}}_{T}.Trinity/contigs.csv", T=TIMESTAMP)
    output:
        a = expand("transrate/summary/{{sample}}_{T}.assemblies.csv", T=TIMESTAMP),
        c = expand("transrate/summary/{{sample}}_{T}.contigs.csv", T=TIMESTAMP)
    params:
        T = TIMESTAMP,
        dir = expand("transrate/{{sample}}_{T}.transrate/{{sample}}_{T}.Trinity/", T=TIMESTAMP)
    shell:
        """
        cp {input.a} {output.a}
        cp {input.c} {output.c}
        rm -rdf {params.dir}*.bam {params.dir}single_component_bad
        rm -rdf {params.dir}{wildcards.sample}_{params.T}.Trinity
        """



# Here the assemblies are being summarized and we need to expanda list
rule summarize_transrate_assembly:
    input:
        expand("transrate/summary/{sample}_{T}.assemblies.csv",sample=SAMPLES, T=TIMESTAMP)
    output:
        "transrate/summary/run_" + TIMESTAMP + "_assemblies.csv"
    run:
        assemblies = input
        #generate first line we need in outfile:
        with open(output[0], "w") as outfile:
            with open(assemblies[0], "r") as infile:
                lines = infile.readlines()
                outfile.write(lines[0])

        # now we open the outfile in appending mode to add all the results
        with open(output[0], "a") as outfile:
            for assembly in assemblies:
                with open(assembly, "r") as infile:
                    lines = infile.readlines()
                    outfile.write(lines[1])


rule fish_4_genes:
    input:
        assemblies = expand("trinity/{sample}_{T}.Trinity.fasta", sample=SAMPLES, T=TIMESTAMP),
        assembly_dir = "trinity/",
        voucher_dir = VOUCHER_DIR,
        vouchers    = config["vouchers"],
        one = config["ref_28S"],
        four = config["ref_4genes"]
    output:
        one = "blast_results/blastn_summary_" + TIMESTAMP + "_fish_28S",
        four ="blast_results/blastn_summary_" + TIMESTAMP + "_fish_4genes"
    params:
        timestamp = TIMESTAMP
    threads: 12
    shell:
        """
        python3 scripts/fish_genes.py \
            -a {input.assembly_dir} \
            -q {input.voucher_dir} -c {input.vouchers} \
            -v {input.one} -T {threads} \
            -o {output.one} &&  python3 scripts/fish_genes.py \
            -a {input.assembly_dir} \
            -q {input.voucher_dir} -c {input.vouchers} \
            -v {input.four} -T {threads} -o {output.four}

        """

# transrate plotting
# TODO change locations of file to something sensible
# ownimgg directory?
rule plot_transrate_contigs:
    input:
        expand("transrate/summary/{{sample}}_{T}.contigs.csv", T=TIMESTAMP)
    output:
        expand("QC/{{sample}}_{T}_contig_dist.png", T=TIMESTAMP)
    script:
        "R/contigs.plot.R"


# add busco rule
rule busco:
    input:
        assembly =   expand("trinity/{{sample}}_{T}.Trinity.fasta", T=TIMESTAMP)
    output:
        out      =   expand(HOME_DIR + "busco/run_{{sample}}_{T}.busco/short_summary_{{sample}}_{T}.busco.txt", T=TIMESTAMP)
    params:
        py3      =   config["py3"],
        BUSCO    =   config["BUSCO"],
        BuscoLib =   config["BuscoLib"],
        homedir  =   HOME_DIR
    threads: 14
    log:
        expand("logs/busco/{{sample}}_{T}.log", T=TIMESTAMP)
    benchmark:
        "benchmarks/{sample}.busco.txt"
    shell:
        # cd {params.outdir}
        """
        {params.py3}  {params.BUSCO} -f \
                -i {params.homedir}{input.assembly} \
                -o {wildcards.sample}_{TIMESTAMP}.busco              \
                -l  {params.BuscoLib}    \
                -m  tran -c {threads} &>> {log} && \
        mv  run_{wildcards.sample}_{TIMESTAMP}.busco  busco/

        """


rule busco_on_transrate:
    # run busco on transrate optimized assembly
    input:
        assembly =   expand("transrate/{{sample}}_{T}.transrate/{{sample}}_{T}.Trinity/good.{{sample}}_{T}.Trinity.fasta", T=TIMESTAMP)
    output:
        out      =   expand(HOME_DIR + "busco/run_{{sample}}_{T}.good.busco/short_summary_{{sample}}_{T}.good.busco.txt", T=TIMESTAMP),
        outass   =   expand("trinity/{{sample}}_{T}_transrate_optimized.Trinity.fasta", T=TIMESTAMP)
    params:
        outdir   =   HOME_DIR + "busco/",
        py3      =   config["py3"],
        BUSCO    =   config["BUSCO"],
        BuscoLib =   config["BuscoLib"],
        homedir  =   HOME_DIR
    threads: 14
    log:
        expand("logs/busco/{{sample}}_{T}.good.log", T=TIMESTAMP)
    # I move the assessed transrate assembly over to the trinity folder
    benchmark:
        "benchmarks/{sample}.busco_on_transrate.txt"
    shell:
        # cd {params.outdir}
        """
        {params.py3}  {params.BUSCO} -f \
                -i {params.homedir}{input.assembly} \
                -o {wildcards.sample}_{TIMESTAMP}.good.busco              \
                -l  {params.BuscoLib}    \
                -m  tran -c {threads} &>> {log} && \
        mv  run_{wildcards.sample}_{TIMESTAMP}.good.busco  busco/ && \
        mv {input.assembly} {output.outass}
        """


# TODO rule decision on assembly
rule gather_busco:
    # move all transrate results to same folder for easier summary later
    # we do not need the expand statement here because plot busco will just call this function many times
    input:
        normal   = expand(HOME_DIR + "busco/run_{{sample}}_{T}.busco/short_summary_{{sample}}_{T}.busco.txt", T=TIMESTAMP),
        good     = expand(HOME_DIR + "busco/run_{{sample}}_{T}.good.busco/short_summary_{{sample}}_{T}.good.busco.txt", T=TIMESTAMP)
    output:
        normal   = expand("busco/summary/short_summary_{{sample}}_{T}.busco.txt", T=TIMESTAMP),
        good     = expand("busco/summary/short_summary_{{sample}}_{T}.good.busco.txt", T=TIMESTAMP)
    shell:
        """
        cp {input.normal} {output.normal}
        cp {input.good}   {output.good}

        """


rule plot_busco:
    input:
        expand("transrate/summary/{sample}_{T}.contigs.csv", T=TIMESTAMP, sample=SAMPLES),
        transrate_summary = "transrate/summary/run_" + TIMESTAMP + "_assemblies.csv",
        normal            = expand("busco/summary/short_summary_{sample}_{T}.busco.txt", sample=SAMPLES, T=TIMESTAMP),
        good              = expand("busco/summary/short_summary_{sample}_{T}.good.busco.txt", sample=SAMPLES, T=TIMESTAMP)
    output:
        busco_pdf         = "busco/summary/run_" + TIMESTAMP + "_busco.pdf",
        busco_trans_pdf   = "busco/summary/run_" + TIMESTAMP + "_busco_with_transrate.pdf"
    params:
        timestamp = TIMESTAMP,
        busco_folder      = "busco/summary/"
    script:
        "R/plot_busco_summary.R"


# TODO rename assembly headers
# TODO
# add plotting for BUSCO
# The default plotting does not work and we have to add
# options(bitmapType='cairo')
# at the beginning of the script. -.-
rule versions:
    output:
        log         = "logs/run_" + TIMESTAMP + "_software.versions"
    params:
        trinity     = config["trinity"],
        transrate   = config["transrate"],
        BUSCO       = config["BUSCO"],
        time        = TIMESTAMP
    shell:
        """
        scripts/versions.sh {params.trinity} {params.transrate} {params.BUSCO} > {output.log}
        """

rule collect_all_stats:
    input:
        busco               = expand(HOME_DIR + "busco/run_{sample}_{T}.busco/short_summary_{sample}_{T}.busco.txt", sample=SAMPLES, T=TIMESTAMP),
        transrate_summary   = "transrate/summary/run_" + TIMESTAMP + "_assemblies.csv",
        assembly_plot       = expand("QC/{sample}_{T}_contig_dist.png", sample=SAMPLES, T=TIMESTAMP),
        forward_qc          = expand("QC/{sample}_R1_fastqc.html", sample=SAMPLES),
        reverse_qc          = expand("QC/{sample}_R2_fastqc.html", sample=SAMPLES),
        forward_trimmed_qc  = expand("QC/{sample}_R1_cor_trim_clean_fastqc.html", sample=SAMPLES),
        reverse_trimmed_qc  = expand("QC/{sample}_R2_cor_trim_clean_fastqc.html", sample=SAMPLES),
        assembly            = expand("rRNA/{sample}_{T}_rRNA.Trinity.fasta", sample=SAMPLES, T=TIMESTAMP),
        basic               = expand("assembly_stats/{sample}_{T}_N50.txt", sample=SAMPLES, T=TIMESTAMP),
        salmon_quant        = expand("expression/salmon/{sample}_{T}/quant.sf", sample=SAMPLES, T=TIMESTAMP),
        ExN50               = expand("assembly_stats/{sample}_{T}_ExN50.txt", sample=SAMPLES, T=TIMESTAMP),
        pep                 = expand("transdecoder/{sample}_{T}.transdecoder.pep", sample=SAMPLES, T=TIMESTAMP),
        cds                 = expand("transdecoder/{sample}_{T}.transdecoder.cds", sample=SAMPLES, T=TIMESTAMP),
        pep_length          = expand("transdecoder/{sample}_{T}.transdecoder.pep.length", sample=SAMPLES, T=TIMESTAMP),
        cds_length          = expand("transdecoder/{sample}_{T}.transdecoder.cds.length", sample=SAMPLES, T=TIMESTAMP),
        mapping             = expand("mapping/{sample}_{T}_snap_to_ref.txt", sample=SAMPLES, T=TIMESTAMP)
    output:
        a = "assembly_stats/All_stats_" + TIMESTAMP + ".txt",
        qc = "QC/All_qc.txt"
    params:
        s = SAMPLES,
        T = TIMESTAMP,
        script = "scripts/collect_macpipe_stats.sh"
    shell:
        """
        set -euo pipefail
        for i in {params.s}; do

            {params.script} "$i"_{params.T} > assembly_stats/"$i"_{params.T}_all_stats.txt
        done
        
        # take first element to have the header line
        arr=({params.s})
        head -n 1 assembly_stats/${{arr[0]}}_{params.T}_all_stats.txt > {output.a}
        for s in $( find assembly_stats/ -name "*_all_stats.txt" );do
            tail -n 1 $s >> {output.a}
        done
        
        echo -e "No_reads\tPercent_deduplicated\tGC_content" > {output.qc}
        cat QC/*_stats.txt >> {output.qc} 
        """

# TODO report file should contain:
# plot of coverage
# busco and transrate stats
# software versions
# timestamp
rule report:
    input:
        all_stats           = "assembly_stats/All_stats_" + TIMESTAMP + ".txt",
        qc                  = "QC/All_qc.txt",
        busco_pdf           = "busco/summary/run_" + TIMESTAMP + "_busco.pdf",
        busco_trans_pdf     = "busco/summary/run_" + TIMESTAMP + "_busco_with_transrate.pdf",
        transrate_summary   = "transrate/summary/run_" + TIMESTAMP + "_assemblies.csv",
        assembly_plot       = expand("QC/{sample}_{T}_contig_dist.png", sample=SAMPLES, T=TIMESTAMP),
        forward_qc          = expand("QC/{sample}_R1_fastqc.html", sample=SAMPLES),
        reverse_qc          = expand("QC/{sample}_R2_fastqc.html", sample=SAMPLES),
        forward_trimmed_qc  = expand("QC/{sample}_R1_cor_trim_clean_fastqc.html", sample=SAMPLES),
        reverse_trimmed_qc  = expand("QC/{sample}_R2_cor_trim_clean_fastqc.html", sample=SAMPLES),
        versions            = "logs/run_" + TIMESTAMP + "_software.versions",
        fish_4genes         = "blast_results/blastn_summary_" + TIMESTAMP + "_fish_4genes",
        fish_28S            = "blast_results/blastn_summary_" + TIMESTAMP + "_fish_28S",
        assembly            = expand("rRNA/{sample}_{T}_rRNA.Trinity.fasta", sample=SAMPLES, T=TIMESTAMP),
        basic               = expand("assembly_stats/{sample}_{T}_N50.txt", sample=SAMPLES, T=TIMESTAMP),
        salmon_quant        = expand("expression/salmon/{sample}_{T}/quant.sf", sample=SAMPLES, T=TIMESTAMP),
        ExN50               = expand("assembly_stats/{sample}_{T}_ExN50.txt", sample=SAMPLES, T=TIMESTAMP),
        pep                 = expand("transdecoder/{sample}_{T}.transdecoder.pep", sample=SAMPLES, T=TIMESTAMP),
        cds                 = expand("transdecoder/{sample}_{T}.transdecoder.cds", sample=SAMPLES, T=TIMESTAMP),
        mapping             = expand("mapping/{sample}_{T}_snap_to_ref.txt", sample=SAMPLES, T=TIMESTAMP)
    output:
        report_name         = "report_" + TIMESTAMP + ".html"
    run:
        from snakemake.utils import report
        run_name = TIMESTAMP
        report("""
        Results of the macpipe_trans run for transcriptome assembly and assessment.
        Run timestamp/ID: {run_name}
        """, output[0], **input)


rule clean_up:
    input:
    output:
        "cleanup"
    run:
        shell("""
            while true; do
                read -p "Do you want to delete all output of the snakemake pipeline?\n WARNING THIS CAN NOT BE REVERSED!" yn
                case $yn in
                    [Yy]* ) rm -rdf data/*cor* trinity/ transrate/ expression/ QC/ benchmarks/ transdecoder/ mapping/ busco/ assembly_stats/ logs/ rRNA/ blast_results/ && touch cleanup ; break;;
                    [Nn]* ) exit;;
                    * ) echo "Please answer yes or no.";;
                esac
            done
                """)



# Finishing up --------------------------------------------------------------
onsuccess:
    print("Workflow finished, no error")
    shell("mail -s 'Snakemake finished with no error' jeremias.brand@unibas.ch < {log}")
onerror:
    print("An error occurred with the snakemake run")
    # here it would be good to include timstamping
    shell("mail -s 'Error in Snakemake run' jeremias.brand@unibas.ch < {log}")

# TODO add separate snakefile for annotation
