#! /bin/bash
echo $1

##Runs inside a docker container, /data2 is your reference folder, /data is your working folder with fastqs as input. docker run chipimage:latest /bin/bash 
##Trim data before running!##
# for i in *.fastq.gz; do java -jar trimmomatic-0.39.jar SE -threads 30 $i $i"Out.fq.gz" ILLUMINACLIP:Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 MINLEN:30; done
												

prefix1=$1
#prefix2=basename $PWD 
input1=$prefix1"_1P.fq.gz"
input2=$prefix1"_2P.fq.gz"
echo $prefix1 "started. Using "$input1 "and" $input2 "as input" $(date) >> Log.log
#bwa mem -T0 -t22 /data2/genome.fa $input1 -o $prefix1"_aligned.sam"
bwa mem -5SP -T0 -t5 /data2/genome.fa $input1 $input2 -o $prefix1"_aligned.sam" 
samtools sort $prefix1"_aligned.sam" -o $prefix1"sorted_aligned.sam"
samtools view -hq 20 $prefix1"sorted_aligned.sam" -o $prefix1"filtered_sorted_aligned.sam"
echo "bwa finished..." $(date) >> Log.log;
java -jar /picard/build/libs/picard.jar MarkDuplicates -INPUT $prefix1"filtered_sorted_aligned.sam" -OUTPUT $prefix1"_aligned_PCRDupes.bam" -ASSUME_SORTED true -METRICS_FILE $prefix1"Aligned_Sorted_PCRDupes.txt" -VALIDATION_STRINGENCY SILENT ;
echo "picard finished..." $(date) >> Log.log;
samtools index $prefix1"_aligned_PCRDupes.bam"
samtools view -b -h -F 0x900 $prefix1"_aligned_PCRDupes.bam" | bedtools bamtobed -i stdin > $prefix1".primary.aln.bed"
#comment out one depending on sample running 
macs2 callpeak -t $prefix1".primary.aln.bed" -n $prefix1".macs2" -B
macs2 callpeak -t $prefix1".primary.aln.bed" -n $prefix1"_NoModel.macs2" --nomodel -B
echo "mac2 finished..." $(date) >> Log.log;
bedtools intersect -v -a $prefix1".macs2_peaks.narrowPeak" -b /data2/FixBlacklist_ENCFF356LFX.bed > $prefix1"NarrowPeakNoBL.bed"
wigToBigWig /data/$prefix1".macs2_treat_pileup.bdg" /data2/hg38.genome /data/$prefix1".bw" -clip
wigToBigWig /data/$prefix1"_NoModel.macs2_treat_pileup.bdg" /data2/hg38.genome /data/$prefix1"_NoModel.bw" -clip
echo "starting flagstat" $(date) >> Log.log;
samtools view -Sb $prefix1"_aligned.sam" -o $prefix1"sorted_aligned.sam"
samtools flagstat $prefix1"sorted_aligned.sam" >> $prefix1"_flagstat.txt"
echo "DONE DONE DONE" $(date) >> Log.log;
done

