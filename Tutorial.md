# Tutorial

1. Enter your search term into ENA search
2. Select <b>Assembly</b> results
3. Download ENA records: XML

4. Create an output directory
5. Run fasta_from_ena.py with XML file as input and directory created as output
6. A summary Excel file is created from search result
7. Wait for FASTA files to be downloaded and unzipped. This can take a while...
8. Run fetch_Entrez_metadata with directory as positional argument (you can include Entrez API and your account email to speed up the run)


For serotyping:
1. Use the Ectyper tool on Galaxy
2. Download results into a folder
3. Run utils.py with `-s` or `--sero` flag to input path for directory containing serotpe data. `--input` should be path for working directory

For phylotyping:
1. Run utils.py with `--input` flag for working directory and use `-p` or `--phylo` flag.
2. Optionally, you can specify location of EzClermont.py script using `--script`, but script is included in package.
