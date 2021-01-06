
# Author: Andrew Valesano
# Purpose: Combine the consensus files and coverage data.

# Usage: python ivar_CombineAndExport.py --run-info test --min-length 100

# ======================= Import modules ======================

import argparse
import glob
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ========================= Main =============================

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--run-info', action="store", dest="runinfo")
    parser.add_argument('--min-length', action="store", dest="min", type = int)
    args = parser.parse_args()

    ### Edit and combine the consensus files into one fasta ###
    os.system("rm data/ivar_output/consensus/*.qual.txt")
    consensus_files = "data/ivar_output/consensus/*.fa"
    all_fasta = list()
    tmp_filename = "data/ivar_output/" + args.runinfo + ".consensus.fasta"
    final_filename = "data/ivar_output/all.consensus.fasta"

    for file in glob.glob(consensus_files):

        print("Working on file: " + file)
        filename_only = file.split("/")[3]
        id = filename_only.split(".")[0]
        for record in SeqIO.parse(file, "fasta"):
            record.id = id
            length = sum(map(lambda x: record.seq.count(x), ["a", "t", "g", "c", "A", "T", "G", "C"]))
            if(length >= args.min):
                all_fasta.append(record)
            else:
                print("Sequence " + id + " was too short! Excluding from final file.")

    with open(tmp_filename, 'w') as full_fasta:
        SeqIO.write(all_fasta, full_fasta, "fasta")

    os.system("sed '/^>/ s/ .*//' " + tmp_filename + " > " + final_filename)
    os.system("rm " + tmp_filename)

    ### Combine coverage files ###
    cov_files = "data/ivar_output/coverage/*.csv"
    cov_filename = "data/ivar_output/coverage.csv"

    df_list = []
    for filename in glob.glob(cov_files):
        filename_only = filename.split("/")[3]
        id = filename_only.split(".")[0]
        df = pd.read_csv(filename, sep = '\t')
        df.columns = ["chr", "pos", "cov"]
        df['ID'] = id
        df_list.append(df)
    df_full = pd.concat(df_list, axis=0, ignore_index=True)
    df_full['run'] = args.runinfo
    df_full.to_csv(cov_filename, index = False)


if __name__ == "__main__":
    main()
