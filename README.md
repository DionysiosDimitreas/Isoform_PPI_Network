# Isoform PPI Network

## Building an isoform-aware CellChat database

1. Prepare 3did domain interaction rulebook

```
# Filter out PFAM domain-domain interaction information
less databases/3did/v2022_01/3did_flat.gz | grep "^#=ID" | cut -f4,5 | perl -ane '$F[0]=~ s/.*(PF\d+).*/$1/; $F[1]=~ s/.*(PF\d+).*/$1/; print "$F[0]\t$F[1]\n$F[1]\t$F[0]\n"' | sort -u | gzip > results/3did_pfamInteractions_1804.tsv.gz
```

2. Prepare Ensembl fasta files for PFAM domain identification

```
# Create directory structure for Ensembl protein FASTA files
mkdir -p fasta/ensembl_v75/human

# Create ENSEMBL protein FASTA files
zless databases/ensembl/v75/Homo_sapiens.GRCh37.75.pep.all.fa.gz | perl -ne 'if(/^>(ENSP\d+) .*/){$i=$1; $h{$i}=$_; next} else{$h{$i}.=$_} END{foreach $i(sort keys %h){open(O, ">fasta/ensembl_v75/human/$i.fasta"); print O $h{$i}; close(O) }}'
```

3. Run PfamScan for domain identification
This repository includes a script to automate running `hmmscan` on a batch of FASTA files and merging the results.

  a) Prerequisites:

You must have [HMMER](http://hmmer.org/) (which includes `hmmscan` and `hmmpress`) installed and available in your system's PATH.

  b) Run the Script

First, make the script executable (you only need to do this once):
```bash
chmod +x scripts/run_pfam_scan.sh
```

Now, run the script by providing three arguments:
i.  Your input FASTA directory
ii.  The path to your `Pfam-A.hmm` database
iii.  A directory name for the results

Example:
```bash
./scripts/run_pfam_scan.sh fasta/ensembl_v75/human/ databases/Pfam/Pfam-A.hmm PfamScan_results/
```

The PfamScan results will be merged in one file, and cleaned.

```
# Merge results
find PfamScan_results/ -name "*.domtblout" -print0 | xargs -0 cat | gzip > ensp_pfam_v115.domtblout.txt.gz

# Clean merged results
python3 scripts/clean_pfam.py
```

4. Get Ligand-Receptor gene pairs from CellChatDB in TSV format

```
Rscript scripts/extract_cellchat_pairs.R
```

5. Gene Identifier Map

This project uses `biomart_export.tsv` as a translation table to map various gene identifiers. This file was generated using **Ensembl BioMart v115**.

The table includes the following columns:
* `ENSG`: Ensembl Gene ID
* `ENST`: Ensembl Transcript ID
* `ENSP`: Ensembl Protein ID
* `HGNC symbol`: Official gene name

6. Create PPI network

```
python3 build_isoform_network.py
```
