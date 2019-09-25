# Forked from Hoohm/dropSeqPipe for TataLab

pipeline for Crop-Seq

## Installation

please make sure dropSeqPip is installed. And make sure you installation is like this

```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
git clone https://github.com/jianhong/dropSeqPipe.git
cd dropSeqPipe
conda create -n dropSeqPipe
conda install -n dropSeqPipe -c bioconda -c conda-forge snakemake
```

## How to run

The files in your WORKING_DIR
```
/path/to/your/WORKING_DIR/
| -- RAW_DATA/
| -- -- sample_name1_R1.fastq.gz
| -- -- sample_name1_R2.fastq.gz
| -- -- sample_name2_R1.fastq.gz
| -- -- sample_name2_R2.fastq.gz
| samples.csv
| gRNAseq.fa
```

- fastq files must be put in RAW_DATA folder.
- gRNAseq.fa must be available.
-	samples.csv must be available.
-	samples.csv template:

```
samples,expected_cells,read_lengths,batch
sample_name1,500,100,Batch1
sample_name2,500,100,Batch2
```

sample of gRNAseq.fa
```
> gRNA1
CACCGTGGGCGGAGACCGTCCTAAT
> gRNA2
CACCGCAGGTCTTCTCGCTACCGA
> gRNA3
CACCGTAGGACGGTCTCCGCCCACC
```

expected_cells is the amount of cells you expect from your sample.

read_length is the read length of the mRNA (Read2). This is necessary for STAR index generation

batch is the batch of your sample. If you are added new samples to the same experiment, this is typically a good place to add the main batch.

run cropseq.sh

```
path/to/cropseq.sh -e your@email.addr -u your.username -d /path/to/dropSeqPipe -i /path/to/igenome/UCSC
```

To download igenome: https://support.illumina.com/sequencing/sequencing_software/igenome.html
