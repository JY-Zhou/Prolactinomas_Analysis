#! /bin/bash
#$ -S /bin/bash
#$ -j yes
#$ -cwd
#$ -q week.q
#$ -N exec_salmon_index
#$ -l p=8,h_vmem=50G,mem_free=50G

source /data/home/jyzhou/.bashrc

set -x

ANNOTATION_GFF=./data/GCF_000001405.25_GRCh37.p13_genomic.gff
REFGENOME_FASTA=./data/GCF_000001405.25_GRCh37.p13_genomic.fna

rm -rf ./index
mkdir ./index

salmon index \
    -t ./data/GCF_000001405.25_GRCh37.p13_rna.fna \
    -i ./index
