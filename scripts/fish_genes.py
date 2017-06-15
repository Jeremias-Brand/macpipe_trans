# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 12:06:18 2017

@author: jeremias
"""
import subprocess
import argparse
import glob
import os


def make_nt_blastdb(infile, outname):
    """
    infile: path to nucleotide fasta
    outname: name of the blast db that is created
    """
    # check if the db allready exists
    if not os.path.isfile(outname + ".nhr"):
        subprocess.check_call(["makeblastdb", "-in", infile,
                               "-dbtype", "nucl", "-out", str(outname)],
                              shell=False)


def blastn(db, query, outfile, evalue="1e-40", threads=12, max_target_seqs="40"):
    """
    blasts query against db and prints the output in the specified format
    column 3 is the name of the sequence hit.
    """
    threads = str(threads)
    if os.path.isfile(outfile):
        print("blastn error, " + outfile + " allready exists! \n skipping blast")
    else:
        tabular_fields = ("6",
                          "qseqid", "qlen", "sseqid", "slen", "frames", "pident",
                          "nident", "length", "mismatch", "gapopen", "qstart",
                          "qend", "sstart", "send", "evalue", "bitscore", "stitle")
        tabular_fields_str = ' '.join(tabular_fields)
        cmd = ["blastn", "-evalue", evalue, "-num_threads", threads, "-db", db,
               "-query", query, "-max_target_seqs", max_target_seqs, "-out",
               outfile, "-outfmt", tabular_fields_str]
        subprocess.check_call(cmd)


def extract_seq_fasta(fasta, name, outfile):
    name = name.strip()
    sed_cmd = "/" + str(name) + "/,/>/p"
    if os.path.isfile(fasta):
        with open(outfile, "a") as f:
            p1 = subprocess.Popen(["sed", "-n", sed_cmd, fasta],
                                  stdout=subprocess.PIPE)
            p2 = subprocess.Popen(["head", "-n", "-1"],
                                  stdin=p1.stdout, stdout=f)
            p1.stdout.close()
    else:
        print("Could not find: " + fasta)


def print_blasthit_names_to_file(blastoutput):
    outfile = blastoutput + "_contig_names"
    with open(outfile, "w") as f:
        subprocess.check_call(("awk", '{print $3}', blastoutput), stdout=f)


def check_hits_against_db(hit, db, outfile, evalue="1e-40", threads="12"):
    from Bio import SeqIO
    # TODO *import* sequences
    hits_fasta = SeqIO.parse(hit, "fasta")
    with open(outfile, "w") as f:
        f.write("""Results from performing blastn of hits with
         INSERT_BAIT_NAME, against""" + db + """.\n
         following we see the name of the sequence and then the top 4 blast hits.\n""")
        for seq in hits_fasta:
            f.writelines("\n          " + seq.id + "\n")
            with open("tmp_seq", "w") as tmp_seq:
                SeqIO.write(seq, tmp_seq, "fasta")
            blastn(db, "tmp_seq", "tmp_file", evalue=evalue, threads=threads, max_target_seqs="4")
            with open("tmp_file", "r") as tmp_f:
                for l in tmp_f.readlines():
                    f.writelines(l)
            os.remove("tmp_file")
        # blastn(barcode_db, seq, outfile, evalue="1e-40", threads="12")
        # print_blasthit_names_to_file(outfile)

        # # we need a function that formats 
        # transcriptome resulted in len(outfile) hits to query
        # <order by length>
        # seq name
        #     blastn results

def check_multiple_hits_against_db(hits, db, outfile, evalue="1e-40", threads="12"):
    from Bio import SeqIO
    if(os.path.isfile(outfile)):
        print(str(outfile) + " allready exists. Aborting check_multiple_hits_against_db.")
    else:   
        with open(outfile, "w") as f:
            f.write("""Results from performing blastn of hits with:\n""")
            for name in hits:
                f.write(name.rsplit("/", 1)[1] + "\n")

            f.write("""\n against:\n"""
              + db + 
              """\n
               Shown is the name of the sequence followed by the top 4 blast hits.
               \n""")
            if(os.path.isfile("tmp_file")):
                os.remove("tmp_file")
            if(os.path.isfile("tmp_seq")):
                os.remove("tmp_seq")
            for hit in hits:
                seqs_in_hit = SeqIO.parse(hit, "fasta")
                f.writelines("\nWe found " + str(len(list(SeqIO.parse(hit, "fasta")))) + " contigs that matched "+ hit +"\n")
                for seq in seqs_in_hit:
                    f.writelines("##" + seq.id + " len= " + str(len(seq.seq)) + "\n")
                    with open("tmp_seq", "w") as tmp_seq:
                        SeqIO.write(seq, tmp_seq, "fasta")
                    try:
                        blastn(db, "tmp_seq", "tmp_file", evalue=evalue, threads=threads, max_target_seqs="4")
                    except subprocess.CalledProcessError:
                        with open("tmp_file", "w") as tmp_file:
                            tmp_file.writelines("NO SEQUENCE FOUND\n")
                        print(str(seq.id) + " is empty!")
                    with open("tmp_file", "r") as tmp_f:
                        for l in tmp_f.readlines():
                            f.writelines(l)
                    os.remove("tmp_file")
                    os.remove("tmp_seq")

def main():
    """
    User input
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--assemblies", type=str, required=True,
                        help="directory containing assemblies to fish. assemblies must end with .fasta")

    parser.add_argument("-q", "--query", type=str, required=True,
                        help="""path to directory with multiple baits\n""")
    parser.add_argument("-c", "--config", type=str, required=True,
                        help="""path to config file containing one bait per line, used with -q""")
    parser.add_argument("-T", "--threads", type=int,
                        help="number of threads available for blastn")
    parser.add_argument("-t", "--timestamp", type=str, required=True,
                        help="timestamp that will be printed in the outfile")
    args = parser.parse_args()
    """
    Checking what arguments have been supplied and supplying default if needed
    """
    TIMESTAMP = args.timestamp
    if(args.threads):
        THREADS = args.threads
    else:
        print("No threads argument given, setting threads to 4.")
        THREADS = "4"
    if(args.config):
        query_dir = args.query
        with open(args.config) as f:
            queries = [line.rstrip('\n') for line in f]
    else:
        query_dir = args.query.rsplit("/", 1)[0]
        queries = list(args.query.rsplit("/", 1)[1])

    assembly_dir = args.assemblies
    # collect all assemblies in the data folder
    assemblies_abs = glob.glob(assembly_dir + "*Trinity.fasta")
    assemblies = [x.rsplit("/", 1)[1] for x in assemblies_abs]


    """
    Iterates through the assemblies and for each blasts with the content of query.
    The result is then extracted and stored in a separate file.
    """
    try:
        os.mkdir("blast_results")
    except OSError:
        print("Directory blast results exits - not created new")

    for assembly in assemblies:
        # create blast dbs
        # assembly_path = os.path.join(assembly_dir, assembly)
        assembly_path = os.path.join(assembly_dir, assembly)
        make_nt_blastdb(assembly_path, assembly_path)
        for query in queries:
            query_path = os.path.join(query_dir, query)
            # simplify the query and transcriptome names
            if assembly.endswith("Trinity.fasta"):
                short_ass = assembly[:-14]
            else:
                short_ass = assembly
            if query.endswith(".fasta"):
                short_query = query[:-6]
            else: 
                short_query = query
            result_name = "blast_results/" + short_ass + "_blast_" + short_query
            if os.path.isfile(result_name):
                print(result_name + "allready exists.")
            else:
                # replace query here with the name in the query seq
                blastn(assembly_path, query_path, result_name, threads=THREADS)
                print_blasthit_names_to_file(result_name)
                with open(result_name + "_contig_names", "r") as inf:
                    for line in inf.readlines():
                        print(line)
                        extract_seq_fasta(assembly_path, line,
                                          result_name+".fasta")
                os.remove(result_name + "_contig_names")

    """
    After we have blasted the targets against the bait we then blast the hits
    against a small database to see what species they likely belong too.
    """

    hits = glob.glob("blast_results/*.fasta")
    check_multiple_hits_against_db(hits, "dbs/vouchers/Macrostomidae_28S.fasta" ,
     "blast_results/blastn_summary_" + TIMESTAMP + "_fish_28S", threads=THREADS)
    check_multiple_hits_against_db(hits, "dbs/vouchers/Macrostomidae_4genes_db.fasta" ,
     "blast_results/blastn_summary_" + TIMESTAMP + "_fish_4genes", threads=THREADS)


if __name__ == '__main__':
    main()
