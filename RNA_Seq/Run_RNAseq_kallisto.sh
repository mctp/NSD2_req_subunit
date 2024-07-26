#!/usr/bin/env bash

#With thanks to James George 
#directory where raw fastq files are stored
rawFileDirectory_original="C:\files\raw" #Windows drive format
rawFileDirectory="${rawFileDirectory_original//\\//}"
rawFileDirectory="/${rawFileDirectory#*/}"
mainFileDirectory=${rawFileDirectory%/*}
organism="human" 

#id file number range
SI_file_range={1..5}


#Kallisto input files
#mouse
m_index="/mus_musculus/m38.99/grcm38.99.idx"
m_gtf="/mus_musculus/m38.99/grcm38.99.clean.gtf"
m_chromosomes="/mus_musculus/m38.99/mouse_chromosomes.tsv"
#human
h_index="/homo_sapiens/hg38_alignment_input_files/hg38.transcriptome.idx"
h_gtf="/homo_sapiens/hg38_alignment_input_files/Homo_sapiens.GRCh38.96.gtf"
h_chromosomes="/homo_sapiens/hg38_alignment_input_files/hg38.chrom.sizes.txt"


# # Begin alignment process
#print out conditions entered
echo "organism: ${organism}"
echo "file_range: ${SI_file_range}"

#1. Merge reads from different lanes
#Enter directory with raw fastqs.  New merged fastqs will be deposited in main file directory
echo "(Step 1 of 3) Combine the libraries distributed across different wells"

cd ${rawFileDirectory}
for file in *_SI_*_*_*_1.fq.gz; do
	new_name=${file#*_}
	new_name=${new_name%_*_*_1.fq.gz}
	cat "$file">>../${new_name}_R1.fq.gz
done

for file in *_SI_*_*_*_2.fq.gz; do
	new_name=${file#*_}
	new_name=${new_name%_*_*_2.fq.gz}
	cat "$file">>../${new_name}_R2.fq.gz
done

# # Enter main file directory
cd ${mainFileDirectory}

#2. Check fastq quality
echo "(Step 2 of 3) Start fastqc analysis"
mkdir fastqcResults
fastqc -t 10 *.fq.gz --outdir=./fastqcResults

cd fastqcResults
unzip '*.zip' '*/summary.txt'
echo "!!!!!!  FASTQC BASE SEQUENCE QUALITY SCREEN  !!!!!"
grep $'Per base sequence quality' */*.txt
grep $'PASS/tPer base sequence quality' */*.txt | wc -l
numberoftotalfiles=$(find *.zip | wc -l)
numberofpassedfiles=$(grep $'PASS\tPer base sequence quality' */*.txt | wc -l)


# Note: -eq is used for numerical comparisons, while = and == are used for string comparisons
if [ ${numberofpassedfiles} -eq ${numberoftotalfiles} ]; then echo 'ALL FILES PASSED!'; else echo 'ERROR: One or more files failed quality check!'; fi

if [ ${numberofpassedfiles} < ${numberoftotalfiles} ]; then
    echo 'Exiting script due to error condition' >&2
    exit 1
fi

#3. Align reads to genome using kallisto (psuedoalignment)
echo "(Step 3 of 3) Begin alignment with kallisto"

# Enter main file directory
cd ${mainFileDirectory}


#change conda environment to access kallisto package
source /miniconda3/bin/activate python36

if [ ${organism} == "mouse" ];
then
	echo "Aligning reads to mouse genome!"
	for i in $(eval echo ${SI_file_range})
	do
		kallisto quant -i ${m_index} -o SI_$i -b 100 --rf-stranded --threads=5 --genomebam --gtf ${m_gtf} --chromosomes ${m_chromosomes} SI_${i}_R1.fq.gz SI_${i}_R2.fq.gz
	done
fi



if [ ${organism} == "human" ]
then
	echo "Aligning reads to human genome!"
	for i in $(eval echo ${SI_file_range})
	do
		nohup kallisto quant -i ${h_index} -o SI_$i -b 100 --rf-stranded --threads=5 --gtf ${h_gtf} --chromosomes ${h_chromosomes} SI_${i}_R1.fq.gz SI_${i}_R2.fq.gz > $i"_log.out" 2>&1 &
	done
fi

echo "Alignment Completed!"

