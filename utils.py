from os import path, system, mkdir
import glob
import argparse

import numpy as np
import pandas as pd

from fasta_from_ena import to_excel


def EzClerment(fp, output, script="EzClermont.py"):
    OUT = path.basename(fp).replace(".fasta", ".txt")
    system(f'python "{script}" "{fp}" 1> "{output}{OUT}" 2>&1')


def run_phylo(fp, script="EzClermont.py"):
    out = path.join(fp, "phylo")
    if not path.isdir(out):
        mkdir(out)
    for f in glob.glob(path.join(fp, "*.fasta")):
        EzClerment(fp=f, output=out, script=script)


def fetch_phylogroup(ena : pd.DataFrame, fp):

    def phylo(x):
        try:
            p = open(path.join(fp, f"{x}.txt")).readlines()[-1]
            result = p.split()[1]
            return result
        except FileNotFoundError:
            return np.NaN

    ena["phylogroup"] = ena.index.map(phylo)
    return ena


def fetch_serotype(ena : pd.DataFrame, fp):

    def ectype(x):
        PATH = path.join(fp, f"{x}.tabular")
        try:
            result = open(PATH).readlines()[-1].split()
        except FileNotFoundError:
            return np.NaN
        
        Otype = result[1]
        Htype = result[2]

        response = " ".join([Otype, Htype])

        return response

    ena["serotype"] = ena.index.map(ectype)
    return ena


def fetch_sequence_type(ena: pd.DataFrame, fp: str):
    mlst = pd.read_table(fp, index_col=0, names=["index","ST"], usecols=[0, 2])
    mlst.index = mlst.index.str.replace(".fasta", "")
    ena = ena.merge(mlst, left_index=True, right_index=True)
    return ena


def main(fp, phylo, sero, st, script=None):

    ena_fp = path.join(fp, "summary.xlsx")
    ena = pd.read_excel(ena_fp, index_col=0)

    if phylo:
        run_phylo(fp, script=script)
        out = path.join(fp, "phylo")
        ena = fetch_phylogroup(ena, out)
        to_excel(ena, fp)

    if sero:
        ena = fetch_serotype(ena, sero)
        to_excel(ena, fp)

    if st:
        fetch_sequence_type(ena, st)
        to_excel(ena, fp)



def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('-i', '--input', type=str)
    parser.add_argument('-s', '--serotype', type=str, default=False)
    parser.add_argument('-p', '--phylo', type=bool, default=False)
    parser.add_argument('-t', '--st', type=str, default=False)
    parser.add_argument('--script', type=str, default="EzClermont.py")

    args = parser.parse_args()
    main(args.input, args.phylo, args.serotype, args.st,script=args.script)


if __name__ == "__main__":
    parse_arguments()


