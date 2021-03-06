[![Total alerts](https://img.shields.io/lgtm/alerts/g/CaileanCarter/EcoliCollection.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/CaileanCarter/EcoliCollection/alerts/)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/CaileanCarter/EcoliCollection.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/CaileanCarter/EcoliCollection/context:python)

# EcoliCollection
Scripts used for handling and analysing my E. coli genome assembly collection.

The workflow starts by pulling genome assemblies from ENA using a saved search result from the site.

## Tutorial

1. Enter a search term into ENA search
2. Select <b>Assembly</b> results
3. Download ENA records: XML
4. Create an output directory
5. Run fasta_from_ena.py with XML file as input and directory created as output
6. A summary Excel file is created from search result
7. Wait for FASTA files to be downloaded and unzipped. This can take a while...
8. Run fetch_Entrez_metadata.py with directory as positional argument (you can include Entrez API and your account email to speed up the run)

For serotyping:
1. Use the Ectyper tool on Galaxy
2. Download results into a folder
3. Run utils.py with `-s` or `--sero` flag to input path for directory containing serotpe data. `--input` should be path for working directory

For phylotyping:
1. Run utils.py with `--input` flag for working directory and use `-p` or `--phylo` flag.
2. Optionally, you can specify location of EzClermont.py script using `--script`, but script is included in package.

## Conda support
1. Run `setup_conda` script in command-line 
or
1. Run:
```
$ conda env create --name ecoli-collection --file environment.yaml
$ conda activate ecoli-collection
```

## Snakemake support
1. Have XML file in EcoliCollection directory
2. Rename XML file to `ena_genome_assembly.xml`
3. Run `snakemake --cores 1`

---

This is working progress:<br>
- [x] Tidy scripts
- [x] Add docstrings
- [x] Make changes to EzClermont
- [ ] Finish README
- [x] Include tutorial
- [x] Trialed and tested
- [x] CLI support
- [x] Snakemake support
- [x] conda support
- [ ] Have Excel spreadsheet by put in parent directory instead of FASTA directory