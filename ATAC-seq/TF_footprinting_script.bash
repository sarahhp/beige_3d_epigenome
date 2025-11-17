# the script for TF footprinting script
#!/bin/bash
INDIR_peak=... # merged_peaks  
INDIR_bam=... 
INDIR_genome=./hg38.fa
INDIR_blacklist=./hg38.blacklist.ENCFF419RSJ.bed
DESTDIR_bam=...
DESTDIR_tobias=...
source /home/miabdelh/miniconda3/etc/profile.d/conda.sh
export TMPDIR=/data/Moab/tmp
GB20=$(echo "20*1024^3" | bc)
mkdir -p $DESTDIR_bam
mkdir -p $DESTDIR_tobias
conda activate bowtie2
# merge bam files 
list1="hASCD13A_bD15 hASCD13A_wD15"
array1=($list1)
count=${#array1[@]}
echo "samtools merge"
for i in `seq 1 $count`
do
 echo "${array1[$i-1]}"
 samtools merge -@ 60 $DESTDIR_bam/"${array1[$i-1]}"_merge.bam $INDIR_bam/"${array1[$i-1]}"_rep1_*_sorted.bam $INDIR_bam/"${array1[$i-1]}"_rep2_*_sorted.bam
 echo "done"
done&
wait
echo "sort_index"
for x in $DESTDIR_bam/*.bam ; do
 name=$(basename -s .bam $x)
 echo "$name"
 samtools sort -@ 8 -m $GB20 $x > $DESTDIR_bam/"$name".sorted.bam
 samtools index -@ 60 $DESTDIR_bam/"$name".sorted.bam
echo "done"
done
# TOBIAS ATACorrect
echo " TOBIAS ATACorrect"
conda activate tobias
echo "hASCD13A_bD15"
mkdir -p $DESTDIR_tobias/bD15
TOBIAS ATACorrect --bam $DESTDIR_bam/hASCD13A_bD15_merge.bam --genome $INDIR_genome --peaks $INDIR_peak/hASCD13A_bD15_Genrich_merged_peaks.bed --blacklist $INDIR_blacklist --outdir $DESTDIR_tobias/bD15  --cores 48
echo "done"
echo "hASCD13A_wD15"
mkdir -p $DESTDIR_tobias/wD15
TOBIAS ATACorrect --bam $DESTDIR_bam/hASCD13A_wD15_merge.bam --genome $INDIR_genome --peaks $INDIR_peak/hASCD13A_wD15_Genrich_merged_peaks.bed --blacklist $INDIR_blacklist --outdir $DESTDIR_tobias/wD15  --cores 48
echo "done"
#################################################################################
#################################################################################
# the second step FootprintScores
INDIR_signal=... 
INDIR_regions=ATAC_peaks.bed
DESTDIR=...
mkdir -p $DESTDIR
source /home/miabdelh/miniconda3/etc/profile.d/conda.sh
export TMPDIR=/data/Moab/Project_Osterbye/tmp
# TOBIAS FootprintScores
echo " TOBIAS FootprintScores"
conda activate tobias
echo "hASCD13A_bD15"
mkdir -p $DESTDIR/bD15
TOBIAS FootprintScores --signal $INDIR_signal/bD15/hASCD13A_bD15_merge_corrected.bw --regions $INDIR_regions --output $DESTDIR/bD15/hASCD13A_bD15_merge_footprints.bw --cores 48
echo "done"
echo "hASCD13A_wD15"
mkdir -p $DESTDIR/wD15
TOBIAS FootprintScores --signal $INDIR_signal/wD15/hASCD13A_wD15_merge_corrected.bw --regions $INDIR_regions --output $DESTDIR/wD15/hASCD13A_wD15_merge_footprints.bw --cores 48
echo "done"
###########################################################################################
#############################################################################################
# the Third step : BINDetect
INDIR_signal=...
INDIR_motifs=./JASPAR2022_CORE_non-redundant_pfms_jaspar.txt
INDIR_genome=./hg38.fa
INDIR_peaks=...
DESTDIR=...
mkdir -p $DESTDIR
source /home/miabdelh/miniconda3/etc/profile.d/conda.sh
export TMPDIR=/data/Moab/tmp
# TOBIAS BINDetect 
echo " TOBIAS BINDetect "
conda activate tobias
TOBIAS BINDetect --motifs $INDIR_motifs --signals $INDIR_signal/bD15/hASCD13A_bD15_merge_footprints.bw $INDIR_signal/wD15/hASCD13A_wD15_merge_footprints.bw  --genome $INDIR_genome --peaks $INDIR_peaks --outdir $DESTDIR --cond_names bD15 wD15 --cores 48
echo "done"
