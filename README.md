# FLANDERS
Functional Landscape And Neighbor Determining gEnomic Region Search

The FLANDERS pipeline enables high throughput examination of neighboring regions of a query gene cluster and its homologs.  FLANDERS utilizes ncbi-acc-download (https://github.com/kblin/ncbi-acc-download) to obtain sequences, antiSMASH (https://docs.antismash.secondarymetabolites.org/install/) to identify known biosynthetic gene clusters, signalP (https://services.healthtech.dtu.dk/cgi-bin/sw_request) to predict signal peptides, and Multigene BLAST (dependcy included in this repository) to identify conserved putative biosynthetic gene clusters of the neighboring regions of a query gene cluster.

# INSTALLATION

In addition to its .py main file, FLANDERS requires a series of dependencies:

1. The Multigene BLAST directory available at https://drive.google.com/drive/folders/1j_Sdj281jqMjUGwsSlVRCKpJUHz51CbD?usp=sharing is to be downloaded and extracted into the same directory as flanders_mgbloop.py
2. Install antiSMASH according to instructions at https://docs.antismash.secondarymetabolites.org/install/ (for use with FLANDERS, Bioconda installation is recommended)
3. Within conda environment created for antiSMASH, install signalP according to instructions at https://services.healthtech.dtu.dk/cgi-bin/sw_request
4. Within conda environment created for antiSMASH, install ncbi-acc-download:
    $pip install ncbi-acc-download
5. Test installation using test_cluster_list.txt file included in repository


# USAGE

FLANDERS accepts a list of genome and gene accession numbers in .csv (.txt) format:

  (query 1) genome acc, gene acc, gene acc, gene acc, gene acc
  (query 2) genome acc, gene acc, gene acc, gene acc, gene acc
  etc.

$flanders_mgbloop.py -i <input file> -o <output parent directory> -r <range of genes upstream and downstream to analyze> -a <perform antismash analysis> -c <threads for antismash computation> -t <taxon for antismash (bacteria or fungi)> -s <perform signalP analysis>
