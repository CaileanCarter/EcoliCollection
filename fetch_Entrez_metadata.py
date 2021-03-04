from Bio import Entrez
import numpy as np
import pandas as pd
from os import path

from fasta_from_ena import get_FASTA, unzip_gz


f"""
fetch_Entrez_metadata.py

This script is a continuation from fasta_from_ena.py and requires the 'ena' output from it.
ena is a DataFrame containing results from ENA. This collection of scripts gathers extra information.
"""

def fetch_id(term):
    """Given an accession or assembly name, fetch the assembly ID"""
    handle = Entrez.esearch(db="assembly", term=term)
    record = Entrez.read(handle)
    idlist = record["IdList"]

    if len(idlist) > 1:
        print(f"{term} returned more than 1 result")
    return idlist[0]


def fetch_summary(id_num):
    """Given an assembly ID, fetch the document summary from Entrez"""
    esummary_handle = Entrez.esummary(db="assembly", id=id_num, report="full")
    esummary_record = Entrez.read(esummary_handle, validate=False)
    summary = esummary_record['DocumentSummarySet']['DocumentSummary'][0]

    return summary


def fetch_parts(summary):
    """Picks out elements of the document summary from Entrez"""
    AssemblyAccession = summary["AssemblyAccession"]
    AssemblyStatus = summary["AssemblyStatus"]
    WGS = summary["WGS"]
    Coverage = summary["Coverage"]
    SubmissionDate = summary["SubmissionDate"]
    LastUpdateDate = summary["LastUpdateDate"]

    return AssemblyAccession, AssemblyStatus, WGS, Coverage, SubmissionDate, LastUpdateDate


def fetch_all(term):
    id_num = fetch_id(term)
    summary = fetch_summary(id_num)
    return fetch_parts(summary)


def fetch_sequence(ena : pd.DataFrame, DIR : str):
    """"Fetch missing sequences from NCBI"""
    for index in ena[ena["FASTA"].isna()].index:
        id_num = fetch_id(index)
        summary = fetch_summary(id_num)
        ftp = summary["FtpPath_RefSeq"]
        

        unzipped = path.join(DIR, f"{index}.fasta")
        zipped = unzipped + ".gz"

        label = path.basename(ftp)
        if not ftp:
            continue
        link = ftp +"/"+ label+'_genomic.fna.gz'
        ena.at[index, "FASTA"] = link

        get_FASTA(url=link, fp=zipped)
        unzip_gz(zipped, unzipped)

    return ena


def main(ena : pd.DataFrame) -> pd.DataFrame:
    """Fetch metadata for Ecoli ENA collection"""
    dump = {"AssemblyAccession" : [],
        "AssemblyStatus" : [],
        "WGS" : [],
        "Coverage" : [],
        "SubmissionDate" : [],
        "LastUpdateDate" : []}

    for index in ena.index:
        results = fetch_all(index)
        for key, result in zip(dump.keys(), results):
            dump[key].append(result if result else np.NaN)
        
    result = pd.DataFrame(dump, index=ena.index)
    ena = ena.merge(result, left_index=True, right_index=True)
    return ena

