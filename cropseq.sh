#!/bin/bash
usage="$(basename "$0") [-h] [-seupabmrd] -- CropSeq_pipeline

where:
    -e  email address.
    -u  your name.
    -d  the directory of dropSeqPipe.
    -i  the directory of ignome. for example /work/jo117/igenome/UCSC

    -h  show this help text
    -s  species; default mm10
    -p  5-prime-smart-adapter; default: AATGATACGGCGACCACCGAGATCTACACGCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC
    -a  adapters; default: NexteraPE-PE.fa; Could be NexteraPE-PE.fa, TruSeq2-PE.fa, TruSeq2-SE.fa, TruSeq3-PE-2.fa, TruSeq3-PE.fa, TruSeq3-SE.fa or any customized fasta file
    -b  strand strategy. default: SENSE. could be SENSE, ANTISENSE and BOTH.
    -m  minimum UMI/Gene pair needed to be counted as one. default: 1
    -l  the filename of your whitelist fi you have one. Well plate base protocols often have one.
	-r  bsub or not; default yes

example:
	$(basename "$0") -e jianhong.ou@duke.edu -u jianhong -d /work/jo117/TataLab/DropSeq/dropSeqPipe -i /work/jo117/igenome/UCSC

Note:
	fastq files must be put in RAW_DATA folder.
	gRNA.gtf and gRNA.fa must be available.
	samples.csv must be available.
	samples.csv template:
samples,expected_cells,read_lengths,batch
sample_name1,500,100,Batch1
sample_name2,500,100,Batch2

expected_cells is the amount of cells you expect from your sample.
read_length is the read length of the mRNA (Read2). This is necessary for STAR index generation
batch is the batch of your sample. If you are added new samples to the same experiment, this is typically a good place to add the main batch.


/path/to/your/WORKING_DIR/
| -- RAW_DATA/
| -- -- sample_name1_R1.fastq.gz
| -- -- sample_name1_R2.fastq.gz
| -- -- sample_name2_R1.fastq.gz
| -- -- sample_name2_R2.fastq.gz
| samples.csv
| gRNAseq.fa

sample of gRNAseq.fa
> gRNA1
CACCGTGGGCGGAGACCGTCCTAAT
> gRNA2
CACCGCAGGTCTTCTCGCTACCGA
> gRNA3
CACCGTAGGACGGTCTCCGCCCACC

please make sure dropSeqPip is installed. And make sure you installation is like this

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
git clone https://github.com/jianhong/dropSeqPipe.git
cd dropSeqPipe
conda create -n dropSeqPipe
conda install -n dropSeqPipe -c bioconda -c conda-forge snakemake

"


email=""
person=""
d=`date +%m-%d-%Y-%H%M%S`
pA="AATGATACGGCGACCACCGAGATCTACACGCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC"
adapters="NexteraPE-PE.fa"
strand="SENSE"
minNumUMI=1
dsp=""
species=mm10
igenome=""
run="yes"
whitelist=""
upsgRNA="gagggcctatttcccatgattccttcatatttgcatatacgatacaaggctgttagagagataattagaattaatttgactgtaaacacaaagatattagtacaaaatacgtgacgtagaaagtaataatttcttgggtagtttgcagttttaaaattatgttttaaaatggactatcatatgcttaccgtaacttgaaagtatttcgatttcttggctttatatatcttgtggaaaggacgaaa"
dwsgRNA="gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgctttttt"

while getopts ':hs:e:u:p:a:b:m:d:r:i:l:' option; do
	case "$option" in
		h) echo "$usage"
		   exit
		   ;;
		s) species=$OPTARG
		   ;;
		e) email=$OPTARG
		   ;;
		u) person=$OPTARG
		   ;;
		p) pA=$OPTARG
		   ;;
		a) adapters=$OPTARG
		   ;;
		b) strand=$OPTARG
		   ;;
		m) minNumUMI=$OPTARG
		   ;;
		d) dsp=$OPTARG
		   ;;
		r) run=$OPTARG
		   ;;
		i) igenome=$OPTARG
		   ;;
		l) whitelist=$OPTARG
		   ;;
		:) printf "missing argument for -%s\n" "$OPTARG" >&2
       	   echo "$usage" >&2
           exit 1
           ;;
		\?) printf "illegal option: -%s\n" "$OPTARG" >&2
           echo "$usage" >&2
           exit 1
           ;;
	esac
done



if [ "$email" == "" ] || [ "$person" == "" ] || [ "$dsp" == "" ] || [ "$igenome" == "" ]; then
	printf "\-e\-u\-d\-i argument must be set\n\n" >&2
	echo "$usage" >&2
	exit
fi

if [ ! -f "gRNAseq.fa" ]; then
	printf "gRNAseq.fa do not exist\n\n" >&2
	exit
fi

if [ -f "$dsp/templates/$adapters" ]; then
	adapters="$dsp/templates/$adapters"
fi

# check Suerat
Rscript -e "packageVersion('Seurat')"

igenome=$igenome/$species

if [ ! -f "${igenome}/Sequence/WholeGenomeFasta/genome.fa" ] || [ ! -f "${igenome}/Annotation/Genes/ensembl.gtf" ]; then
	printf "${igenome}/Sequence/WholeGenomeFasta/genome.fa or ${igenome}/Annotation/Genes/ensembl.gtf do not exist\n\n" >&2
	printf "Please download igenome from https://support.illumina.com/sequencing/sequencing_software/igenome.html" >&2
	exit
fi

userID=`whoami`
tmpFolder="/work/$userID/tmp"
mkdir -p $tmpFolder

#git clone https://github.com/Hoohm/dropSeqPipe.git

## prepare samples.csv
#samples,expected_cells,read_lengths,batch
#sample_name1,500,100,Batch1
#sample_name2,500,100,Batch2

#mkdir -p RAW_DATA
## put fastq files into fastq folder
pd=${PWD}
## prepare genome and gtf files
## the gtf files must contain exon in column 3

## prepare gRNA.gtf and gRNA.fa
rm -f gRNA.gtf
rm -f gRNA.fa

seq=""
newheader=""
while IFS= read -r line
do
 if [ "${line:0:1}" == ">" ]; then
   header=$newheader
   newheader=${line#">"}
   newheader=${newheader//[[:blank:]]/}
   if [ "$header" != "" ]; then
     echo ">chr_$header" >> gRNA.fa
     echo "$upsgRNA$seq$dwsgRNA" >> gRNA.fa
     end=$(( 247 + ${#seq}))
     echo -e "chr_${header}\tspikein\tgene\t246\t${end}\t.\t+\t.\tgene_id \"${header}\"; gene_name \"${header}\"; gene_source \"spikein\"; gene_biotype \"protein_coding\"" >> gRNA.gtf
     echo -e "chr_${header}\tspikein\ttranscript\t246\t${end}\t.\t+\t.\tgene_id \"${header}\"; gene_name \"${header}\"; gene_source \"spikein\"; gene_biotype \"protein_coding\"; transcript_id \"${header}-transcript\"; transcript_name \"${header}-transcript\"; transcript_source \"spikein\"; transcript_biotype \"protein_coding\"" >> gRNA.gtf
     seq=""
   fi
 else
   seq="$seq$line"
 fi
done < "gRNAseq.fa"
header=$newheader
echo ">chr_$header" >> gRNA.fa
echo "$upsgRNA$seq$dwsgRNA" >> gRNA.fa
echo -e "chr_${header}\tspikein\tgene\t246\t${end}\t.\t+\t.\tgene_id \"${header}\"; gene_name \"${header}\"; gene_source \"spikein\"; gene_biotype \"protein_coding\"" >> gRNA.gtf
echo -e "chr_${header}\tspikein\ttranscript\t246\t${end}\t.\t+\t.\tgene_id \"${header}\"; gene_name \"${header}\"; gene_source \"spikein\"; gene_biotype \"protein_coding\"; transcript_id \"${header}-transcript\"; transcript_name \"${header}-transcript\"; transcript_source \"spikein\"; transcript_biotype \"protein_coding\"" >> gRNA.gtf
seq=""

mkdir -p ${species}_gRNA_${d}_1

grep "gene_name" ${igenome}/Annotation/Genes/ensembl.gtf > ${species}_gRNA_${d}_1/ensembl.tmp.gtf
cat ${species}_gRNA_${d}_1/ensembl.tmp.gtf gRNA.gtf > ${species}_gRNA_${d}_1/annotation.gtf
rm ${species}_gRNA_${d}_1/ensembl.tmp.gtf

cat ${igenome}/Sequence/WholeGenomeFasta/genome.fa gRNA.fa > ${species}_gRNA_${d}_1/genome.fa

## prepare config.yaml
cat <<EOT >> config.yaml
CONTACT:
    email: "$email"
    person: "$person"
LOCAL:
    temp-directory: "$tmpFolder"
    memory: "40g"
    raw_data: "RAW_DATA"
    results: "RESULTS"
META:
    species:
        ${species}_gRNA:
            build: "${d}"
            release: 1
    ratio: 0.2
    reference-directory: "$pd"
    gtf_biotypes: "$dsp/templates/gtf_biotypes.yaml"
FILTER:
    barcode-whitelist: "$whitelist"
    5-prime-smart-adapter: "$pA"
    cell-barcode:
        start: 1
        end: 12
    UMI-barcode:
        start: 13
        end: 20
    cutadapt:
        adapters-file: "$adapters"
        R1:
            quality-filter: 20
            maximum-Ns: 0
            extra-params: ''
        R2:
            quality-filter: 20
            minimum-adapters-overlap: 6
            minimum-length: 15
            extra-params: ''
MAPPING:
    STAR:
        genomeChrBinNbits: 18
        outFilterMismatchNmax: 10
        outFilterMismatchNoverLmax: 0.3
        outFilterMismatchNoverReadLmax: 1
        outFilterMatchNmin: 0
        outFilterMatchNminOverLread: 0.66
        outFilterScoreMinOverLread: 0.66
EXTRACTION:
    LOCUS:
        - CODING
        - UTR
    strand-strategy: $strand
    UMI-edit-distance: 1
    minimum-counts-per-UMI: $minNumUMI

EOT

pt=`which python`
pt=`echo "${pt/bin\/python/lib}"`

cat <<EOT > dropseqpipe.sh
#!/bin/bash
#SBATCH -J dropSeqPipe #jobname
#SBATCH --array=1 # job array, can be 1,3,5,7 or 1-7:2 '%' will limit the number of job
#SBATCH -o log/dropSeqPipe.out.%A_%a.txt # %A == $SLURM_ARRAY_JOB_ID #log folder must be there
#SBATCH -e log/dropSeqPipe.err.%A_%a.txt # %a == $SLURM_ARRAY_TASK_ID
#SBATCH --mem-per-cpu=7G #5G
#SBATCH -c 8 # cpus-per-task

# module load STAR/2.5.3a
export PATH=$PATH:/opt/apps/rhel7/STAR-2.5.4/STAR/bin/Linux_x86_64/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$pt

i=$(($SLURM_ARRAY_TASK_ID-1))

cd $dsp

source activate dropSeqPipe
snakemake --use-conda --cores \$SLURM_CPUS_PER_TASK --directory $pd
conda deactivate


EOT

mkdir -p log

if [ "$run" == "yes" ]; then
sbatch dropseqpipe.sh
fi

