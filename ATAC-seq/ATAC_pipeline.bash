# the fastp script 
#!/bin/bash
INDIR=...
DESTDIR=...
mkdir -p $DESTDIR/fastp_trimseq
mkdir -p $DESTDIR/fastp_json
# list of folder names of original files
list1="..."
# list of files names for new files after
list2="..."
array1=($list1)
array2=($list2)
count=${#array1[@]}
for i in `seq 1 $count`
do
mkdir -p $DESTDIR/fastp_trimseq/"${array2[$i-1]}"
fastp -j $DESTDIR/fastp_json/"${array2[$i-1]}".fastp.json -w 30  -i $INDIR/"${array1[$i-1]}"/"${array1[$i-1]}"_R1.fastq.gz -I $INDIR/"${array1[$i-1]}"/"${array1[$i-1]}"_R2.fastq.gz -o $DESTDIR/fastp_trimseq/"${array2[$i-1]}"/"${array2[$i-1]}"_R1_trim.fastq.gz -O $DESTDIR/fastp_trimseq/"${array2[$i-1]}"/"${array2[$i-1]}"_R2_trim.fastq.gz
done
######################################################################################################
########################################################################################################
# the  script for mapping 
#!/bin/bash
INDIR=...
DESTDIR=...
REF=.../hg38 #Ref genome
mkdir -p $DESTDIR
source /home/miabdelh/miniconda3/etc/profile.d/conda.sh
conda activate bowtie2
GB20=$(echo "20*1024^3" | bc)
echo "Start bowtie"
list2="..."
array2=($list2)
count=${#array2[@]}
for i in `seq 1 $count`
do
echo "${array2[$i-1]}"
bowtie2 -x $REF -1 $INDIR/"${array2[$i-1]}"/"${array2[$i-1]}"_R1_trim*.gz -2 $INDIR/"${array2[$i-1]}"/"${array2[$i-1]}"_R2_trim*.gz -p60 | samtools view -Su - | samtools sort - -m $GB20 -@ 8 -o  $DESTDIR/"${array2[$i-1]}".sorted.bam 
samtools index -@ 24 $DESTDIR/"${array2[$i-1]}".sorted.bam
echo "done"
done&
wait
flagstat_out=/projects/imb-pkbphil/moab/Project_Galigniana-ATAC2-2023/Thermo_reads/Mapped_reads_trim/flagstat_out
mkdir -p $flagstat_out
echo "samtools flagstat"
for x in $DESTDIR/*.sorted.bam ; do
 qname2=`basename $x `
 echo "$qname2"
 samtools flagstat -@ 12 $x > $flagstat_out/$qname2.flagstats.txt
done
########################################################################################################
###########################################################################################################
# script for filtering bam files
#### Remove the mitochondrial reads after alignment 
INDIR=...
DESTDIR=...
removeChrom=./removeChrom.py
mkdir -p $DESTDIR
GB20=$(echo "20*1024^3" | bc)
source activate bowtie2
for i in $INDIR/*.bam ; do
name=$(basename -s .sorted.bam $i)
echo "$name"
samtools view -h  $i  |  python $removeChrom - - chrM | samtools view -Su - | samtools sort - -m $GB20 -@ 4 -o  $DESTDIR/"$name"_mtremoved_sorted.bam
samtools index -@ 24 $DESTDIR/"$name"_mtremoved_sorted.bam
echo "done"
done&
wait
#### get flagstat
mkdir -p $DESTDIR/flagstat_out
echo "samtools flagstat"
for x in $DESTDIR/*_mtremoved_sorted.bam ; do
 qname2=`basename $x `
 echo "$qname2"
 samtools flagstat -@ 12 $x > $DESTDIR/flagstat_out/$qname2.flagstats.txt
done
echo "done"
##########################################################################################
############################################################################################
#### remove PCR duplicates
srcd_output=...
INDIR=...
tmp_dir=$srcd_output/tmp # tmp directory for picard
mkdir -p $tmp_dir
echo "MarkDuplicates"
for i in $INDIR/*_mtremoved_sorted.bam ; do
 COUNTER=$((COUNTER+1)) 
 qname=$(basename -s _mtremoved_sorted.bam $i)
 echo -e "["$COUNTER"]"'\t'$qname 
  picard MarkDuplicates I=$i O=$srcd_output/$qname.Mkdup.bam M=$srcd_output/$qname_Mkdup_metrics.txt ASSUME_SORTED=false REMOVE_DUPLICATES=true CREATE_INDEX=true TMP_DIR=$tmp_dir 
  echo "done"
done&
wait
#### get flagstat
echo "flagstat MarkDuplicates "
mkdir -p $srcd_output/flagstat_out
echo "samtools flagstat"
for x in $srcd_output/*Mkdup.bam ; do
 qname2=`basename $x `
 echo "$qname2"
 samtools flagstat -@ 12 $x > $srcd_output/flagstat_out/$qname2.flagstats.txt
done
echo "done"
####################################################################
#####################################################################
##### remove Non-unique alignments
echo "remove Non-unique alignments"
DESTDIR=...
INDIR=...
mkdir -p $DESTDIR
for i in $INDIR/*.Mkdup.bam ; do
 name=$(basename -s .Mkdup.bam $i)
 echo "$name"
 samtools view -b  -q 10 $i | samtools sort - -m $GB20 -@ 4 -o  $DESTDIR/"$name"_filteredfinal_sorted.bam
 samtools index -@ 24 $DESTDIR/"$name"_filteredfinal_sorted.bam
echo "done"
done&
wait
#### get flagstat
echo "flagstat remove Non-unique alignments"
mkdir -p $DESTDIR/flagstat_out
for x in $DESTDIR/*_filteredfinal_sorted.bam ; do
 qname2=`basename $x `
 echo "$qname2"
 samtools flagstat -@ 12 $x > $DESTDIR/flagstat_out/$qname2.flagstats.txt
done
echo "done"
########################################################################################################################
##############################################################################################################
# the script used for peakcalling 
export TMPDIR=/data/Moab/Project_Osterbye/tmp
source /home/miabdelh/miniconda3/etc/profile.d/conda.sh
INDIR=...
DESTDIR=...
GB20=$(echo "20*1024^3" | bc)
mkdir -p $DESTDIR/peaks
mkdir -p $DESTDIR/reads_bed
mkdir -p $DESTDIR/pq_bed
mkdir -p $INDIR/sorted_queryname
# sort man by queryname
echo "sort man by queryname"
conda activate bowtie2
for i in $INDIR/*.bam ; do
name=$(basename -s _filteredfinal_sorted.bam $i)
echo "$name"
samtools sort $i -n -m $GB20 -@ 8 -O bam -o $INDIR/sorted_queryname/"$name"_queryname_sorted.bam
echo "done"
done&
wait
# to start Genrich
conda activate base
echo "start Genrich"
for i in $INDIR/sorted_queryname/*.bam ; do
name=$(basename -s _queryname_sorted.bam $i)
echo "$name"
Genrich -t $i  -o $DESTDIR/peaks/"$name".bed -j -y -v -b $DESTDIR/reads_bed/"$name".bed -f $DESTDIR/pq_bed/"$name".bed
echo "done"
done

