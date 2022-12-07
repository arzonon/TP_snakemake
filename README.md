# Cloud_Snakemake

TP Snakemake
The project is done in my master degree for the Parallel Calculations module. The objective is to re-analyse ATAC-seq data with a compute cluster. The analysis was on a B3 mouse cell line. The project is to make a workflow for the ATAC-seq analysis with snakemake.

Experimental Design

Dataset collection
Datasets were collected from European Nucleotide Archive (ENA) https://www.ebi.ac.uk/ena on 2021-07-26 by Nadia Goué with enaDataGet tool https://www.ebi.ac.uk/about/news/service-news/new-tools-download-data-ena included in inhouse script ena_array_atac.slurm

Cells collected for 0 and 24 hours post-treatment with tamoxifen.
3 biological replicates of ~50,000 cells


Raw dataset

SRR4785152      50k-Rep1-0h-sample.0h        GSM2367179      0.7G
SRR4785153      50k-Rep2-0h-sample.0h        GSM2367180      0.7G
SRR4785154      50k-Rep3-0h-sample.0h        GSM2367181      0.7G
SRR4785341      50k-24h-R1-sample.24h.2      GSM2367368      0.6G
SRR4785342      50k-24h-R2-sample.24h.2      GSM2367369      0.7G
SRR4785343      50k-24h-R3-sample.24h.2      GSM2367370      0.6G


Workflow
 I - Data pre-processing 

Data quality control with  fastqc  tool
Cleaning of reads with  trimmomatic  tool
Data quality control with  fastqc  tool

 II - Mapping 

Alignment with  bowtie2  tool
Alignment cleaning with  picard  tool

 III - Data mining  with  deeptools  tool
 IV - Identification of DNA accessibility sites  with  macs2  tool
 V - Identification of unique and common accessibility sites between cell stages  with  bedtools  tool

Run the program
To launch the workflow, execute the Snakefile script in a conda environment with the command:

conda activate snakemake
snakemake --cores all --use-conda --snakefile scripts/Snakefile


Pay attention to the paths in the different scripts before launching the workflow

Tools

snakemake : https://snakemake.readthedocs.io/en/stable/

fastqc : https://manpages.ubuntu.com/manpages/focal/man1/fastqc.1.html

trimmomatic : http://www.usadellab.org/cms/index.php?page=trimmomatic

bowtie2 : https://bowtie-bio.sourceforge.net/bowtie2/index.shtml

picard : https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates

deeptools : https://deeptools.readthedocs.io/en/develop/

macs2 : https://pypi.org/project/MACS2

bedtools : https://github.com/hbctraining/Intro-to-ChIPseq/blob/master/lessons/data_visualization_with_bedtools.md

IGV : http://software.broadinstitute.org/software/igv/



Reference
Gomez-Cabrero et al. (2019). STATegra, a comprehensive multi-omics dataset of B-cell differentiation in mouse. Scientific Data, 6(1). https://doi.org/10.1038/s41597-019-0202-7.

Authors

Aristides ZONON aristides.zonon@etu.uca.fr

## Directory structure

├── .gitignore
├── README.md
├──Snakefile
├── config
│   ├── config.yaml
│   └── env
│        ├── qc.yaml
│        ├── trim.yaml
│        ├── bowtie2.yaml
│        ├── deeptools.yaml
│        ├── picard.yaml
│        └── macs2.yaml
├── data
│   └── data
│        ├── subset
│        
│               
└── results

 How to run the workflow
To run the snakemake workflow use the command : snakemake --cores all --use-conda --snakefile scripts/Snakefile