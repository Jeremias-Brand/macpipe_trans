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
  - add a plock to run if fail or succeed

"""
# import statements
from os.path import join
import os
import time



# Uses a yaml file as input
configfile: "config.yaml"

# Preparation------------------------------------------------------------------
TIMESTAMP = time.strftime("%Y%m%d")

# check if the necessary dirs exist and if not creates them
def save_mkdir( dirs ):
	for d in dirs:
		if not os.path.isdir(d):
			os.mkdir(d)
			print("Creating directory: " + d)


dirs = ["logs", "logs/trinity", "logs/transrate",
        "logs/rcorrector", "trinity", "transrate",
        "busco", "logs/busco"]
save_mkdir(dirs)

# Globals ---------------------------------------------------------------------

# define global things here

# Full path to a folder that holds all of your FASTQ files.
FASTQ_DIR = config["fastqdir"]

# A Snakemake regular expression matching the forward mate FASTQ files.
SAMPLES, = glob_wildcards(join(FASTQ_DIR, '{sample}_R1.fq'))

# Patterns for the 1st mate and the 2nd mate using the 'sample' wildcard.
PATTERN_R1 = '{sample}_R1.fq'
PATTERN_R2 = '{sample}_R2.fq'

# Rules -----------------------------------------------------------------------



rule all:
    input:
        "trans.sentinel",
        "report.html"

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
        FASTQ_DIR + "{sample}_R1.cor.fq",
        FASTQ_DIR + "{sample}_R2.cor.fq"
    log:
    	"logs/rcorrector/{sample}.log"
    params:
        outdir = config["fastqdir"]
    shell:
        """
        ~/perl5/perlbrew/perls/5.16.2t/bin/perl ~/bin/run_rcorrector.pl -k 31 -t 30 \
        -od {params.outdir}  -1 {input.forward}  -2 {input.reverse} &> {log}
        """


rule trim_and_trinity:
    input:
        forward = FASTQ_DIR + "{sample}_R1.cor.fq",
        reverse = FASTQ_DIR + "{sample}_R2.cor.fq"
    output:
        "{sample}.Trinity.fasta",
        "logs/trinity/{sample}.shell.log"
    log:
        "logs/trinity/{sample}.log"
    params:
        adapters =  config["adapters_fasta"],
        outdir   =  config["fastqdir"]
    threads: 28  # threads only works if --cores is set to the actual number of cores when running the snakemake
    # message: expand("Executing with {threads} threads on the following files {sample}.", sample=SAMPLES)
    # We also want a logfile
    # log: expand("logs/{sample}.trinity.log", sample=SAMPLES)
    shell:  
        """
        /home/jeremias/soft/trinityrnaseq-2.2.0/Trinity --seqType fq \
        --trimmomatic --CPU {threads} --max_memory 150G \
        --output 'trinity/'{wildcards.sample}'.trinity' \
        --left   {input.forward}  \
        --right  {input.reverse} \
        --normalize_reads  --normalize_max_read_cov 30 \
        --quality_trimming_params 'ILLUMINACLIP:{params.adapters}:2:40:15 LEADING:2 TRAILING:2 MINLEN:25' > ./logs/trinity/{wildcards.sample}.shell.log 2> {log} && \
        mv trinity/{wildcards.sample}.trinity/Trinity.fasta {wildcards.sample}.Trinity.fasta
        """

rule run_transrate:
    input:
        assembly =  "{sample}.Trinity.fasta",
        forward  =  FASTQ_DIR + "{sample}_R1.cor.fq",
        reverse  =  FASTQ_DIR + "{sample}_R2.cor.fq"
    output:
        "logs/transrate/{sample}.transrate.log",
        "transrate/{sample}.transrate/{sample}.Trinity/good.{sample}.Trinity.fasta"
    log:
    	"logs/transrate/{sample}.log"
    params:
        transDir =  config["transrateDir"],
        dataDir  =  config["fastqdir"],
        homeDir  =  config["homedir"]
    threads: 28
    shell:
        """
        transrate --assembly {params.homeDir}{input.assembly} \
        --output {params.transDir}{wildcards.sample}.transrate \
        --threads {threads} \
        --left {input.forward} \
        --right {input.reverse} > logs/transrate/{wildcards.sample}.transrate.log 2> {log}
        """

rule gather_transrate:
    input:
        assembly = expand("transrate/{samples}.transrate/assemblies.csv", samples=SAMPLES)
        contigs  = expand("transrate/{samples}.transrate/{samples}.Trinity/contigs.csv", samples=SAMPLES)
    output:
        assembly = 
    shell:
        """
        cp {
# TODO
# transrate plotting

rule cat_transrate:
    input:
        expand("logs/transrate/{sample}.transrate.log", sample = SAMPLES)
    output:
        "trans.sentinel"
    shell:
        "touch trans.sentinel"

# add busco rule
rule busco:
    input:
        assembly =   "{sample}.Trinity.fasta"
    output:
        out      =   config["homedir"] + "busco/run_{sample}.busco/short_summary_{sample}.busco.txt"
    params:
        py3      =   config["py3"],
        BUSCO    =   config["BUSCO"],   
        BuscoLib =   config["BuscoLib"],
        homedir  =   config["homedir"],
        folder   =   "run_{sample}.busco"
    threads: 14
    log:
        "logs/busco/{sample}.log"

    shell:
        # cd {params.outdir}
        """
        {params.py3}  {params.BUSCO} -f \
                -i {params.homedir}{input.assembly} \
                -o {wildcards.sample}.busco              \
                -l  {params.BuscoLib}    \
                -m  tran -c {threads} &>> {log} && \
        mv  {params.folder}  busco/

        """


rule busco_on_transrate:
    # run busco on transrate optimized assembly
    input:
        assembly =   "transrate/{sample}.transrate/{sample}.Trinity/good.{sample}.Trinity.fasta"
    output:
        out      =   config["homedir"] + "busco/run_{sample}.good.busco/short_summary_{sample}.good.busco.txt"
    params:
        outdir   =   config["homedir"] + "busco/",
        py3      =   config["py3"],
        BUSCO    =   config["BUSCO"],   
        BuscoLib =   config["BuscoLib"],
        homedir  =   config["homedir"],   
        folder   =   "run_{sample}.good.busco"
    threads: 14
    log:
        "logs/busco/{sample}.good.log"

    shell:
        # cd {params.outdir}
        """
        {params.py3}  {params.BUSCO} -f \
                -i {params.homedir}{input.assembly} \
                -o {wildcards.sample}.good.busco              \
                -l  {params.BuscoLib}    \
                -m  tran -c {threads} &>> {log} && \
        mv  {params.folder}  busco/

        """
# TODO rule decision on assembly

# TODO rename assembly headers

# TODO 


# TODO report file should contain:
# plot of coverage
# busco and transrate stats
# software versions
# timestamp
rule report:
    input:
        expand("logs/transrate/{sample}.transrate.log", sample=SAMPLES),
        expand( config["homedir"] + "busco/run_{sample}.busco/short_summary_{sample}.busco.txt", sample=SAMPLES),
        expand( config["homedir"] + "busco/run_{sample}.good.busco/short_summary_{sample}.good.busco.txt", sample=SAMPLES)
    output:
        "report.html"
    run:
        from snakemake.utils import report
        with open(input[0]) as trans.logs:
            print(trans.logs)

        report("""
        An example variant calling workflow
        ===================================

        Reads were mapped to the Yeast
        reference genome and variants were called jointly with
        SAMtools/BCFtools.

        This resulted in {n_calls} variants (see Table T1_).
        """, output[0], T1=input[0])


# Finishing up --------------------------------------------------------------
# onsuccess:
#     print("Workflow finished, no error")

# onerror:
#     print("An error occurred with the snakemake run")
#     # here it would be good to include timstamping
#     # TODO config for email
#     shell:
#     	"mail -s "an error occurred" jeremias.brand@unibas.ch < {log}"

# TODO add separate snakefile for annotation
