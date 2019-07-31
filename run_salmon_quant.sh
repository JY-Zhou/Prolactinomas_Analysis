#! /bin/bash
#$ -S /bin/bash
#$ -j yes
#$ -cwd
#$ -q week.q
#$ -N exec_salmon_quant
#$ -l p=8,h_vmem=50G,mem_free=50G
#$ -l h=n001|n002|n003|n004|n005|n006|n008|n009|n010|n011|n012
#$ -t 1-22

source /data/home/jyzhou/.bashrc

set -x

LIST=./samples.list
RAW_READS=./data/raw_reads

IFS=$'\n' filelist=($(cat $LIST))
filename=${filelist[$((SGE_TASK_ID-1))]}
echo '>>> Now processing ' ${filename}
rm -rf PA_${filename}
mkdir PA_${filename}

salmon quant -p 8 -l A --gcBias \
    -i ./index \
    -1 ${RAW_READS}/${filename}_1.clean.fq.gz \
    -2 ${RAW_READS}/${filename}_2.clean.fq.gz \
    -o PA_${filename}

echo '^^^ Finished! ' ${filename}
