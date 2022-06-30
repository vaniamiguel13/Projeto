#! /usr/bin/bash
#conda install -c bioconda samtools
#conda install -c bioconda bwa
#conda install -c bioconda picard
#conda install -c bioconda gatk4
#conda install -c bioconda snpeff
#conda install -c conda-forge r
#conda install r-gsalib

echo 'Insert Path of Reference Genome'
read RG

echo 'Insert Path of Fastq_1'
read FQ1

echo 'Insert Path of Fastq_2'
read FQ2

echo 'Insert wanted file name (no extension)'
read NAME

bwa index $RG

bwa mem -K 100000000 -Y -R "@RG\tID:$NAME\tLB:$NAME\tPL:illumina\tPM:hiseq\tSM:$NAME" $RG $FQ1 $FQ2 > $NAME.sam

_s='.sam'
SAMP=$NAME$_s

gatk SortSam -I $SAMP -O $NAME.sorted.sam -SO coordinate

_ss='.sorted.sam'
SorSAM=$NAME$_ss

gatk MarkDuplicates -I $SorSAM -M $NAME.dedup_metric.txt -O $NAME.sorted_dedup.bam

_sdb='.sorted_dedup.bam'
BAMP=$NAME$_sdb

samtools index $BAMP 

gatk ValidateSamFile -R $RG -I $BAMP

samtools depth -a $BAMP > $NAME.depth.txt

gatk CollectAlignmentSummaryMetrics -R $RG -I $BAMP -O $NAME.alignment_metrics.txt 

gatk CollectInsertSizeMetrics -I $BAMP -O $NAME.insert_metrics2.txt -H $NAME.insert_size_histogram.pdf

samtools faidx $RG

ref=${RG%.*}

gatk CreateSequenceDictionary -R $RG -O $ref.dict

gatk HaplotypeCaller -R $RG -I $BAMP -O $NAME._raw_variants.vcf

rvcf1='._raw_variants.vcf'
RV1=$NAME$rvcf1

gatk SelectVariants -R $RG -V $RV1 --select-type-to-include SNP -O $NAME.raw_snps.vcf
rs1='.raw_snps.vcf'
SVCF=$NAME$rs1

gatk SelectVariants -R $RG -V $RV1 --select-type-to-include INDEL -O $NAME.raw_indel.vcf
ri1='.raw_indel.vcf'
IVCF=$NAME$ri1

gatk VariantFiltration -R $RG -V $SVCF -O $NAME.filtered_snps.vcf --filter-name "QD_filter" -filter "QD < 2.0" --filter-name "FS_filter" -filter "FS > 60.0" --filter-name "MQ_filter" -filter "MQ < 40.0" --filter-name "SOR_filter" -filter "SOR > 4.0" --filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" --filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
fsv='.filtered_snps.vcf'
SFVCF=$NAME$fsv

gatk VariantFiltration -R $RG -V $IVCF -O $NAME.filtered_indels.vcf --filter-name "QD_filter" -filter "QD < 2.0" --filter-name "FS_filter" -filter "FS > 200.0" --filter-name "SOR_filter" -filter "SOR > 10.0"
fiv='.filtered_indels.vcf'
IFVCF=$NAME$fiv

gatk SelectVariants --exclude-filtered -V $SFVCF -O $NAME.bqsr_snps.vcf
gatk SelectVariants --exclude-filtered -V $IFVCF -O $NAME.bqsr_indels.vcf

bqsrs='.bqsr_snps.vcf'
bqsri='.bqsr_indels.vcf'
BQSRS=$NAME$bqsrs
BQSRI=$NAME$bqsri

gatk BaseRecalibrator -R $RG -I $BAMP --known-sites $BQSRS --known-sites $BQSRI -O $NAME.recal_data.table
rt='.recal_data.table'
RT=$NAME$rt

gatk ApplyBQSR -R $RG -I $BAMP -bqsr $RT -O $NAME.recal_reads.bam
rrb='.recal_reads.bam'
RBAM=$NAME$rrb

gatk BaseRecalibrator -R $RG -I $RBAM --known-sites $BQSRS --known-sites $BQSRI -O $NAME.post_recal_data.table
prt='.post_recal_data.table'
PRT=$NAME$prt

gatk AnalyzeCovariates -before $RT -after $PRT -plots $NAME.recalibration_plots.pdf

gatk HaplotypeCaller -R $RG -I $RBAM -O $NAME.raw_variants_recal.vcf
rvr='.raw_variants_recal.vcf'
RVRVCF=$NAME$rvr

gatk SelectVariants -R $RG -V $RVRVCF --select-type-to-include SNP -O $NAME.raw_snps_recal.vcf
gatk SelectVariants -R $RG -V $RVRVCF --select-type-to-include INDEL -O $NAME.raw_indels_recal.vcf

rsrv='.raw_snps_recal.vcf'
rirv='.raw_indels_recal.vcf'
RSRVCF=$NAME$rsrv
RIRVCF=$NAME$rirv

gatk VariantFiltration -R $RG -V $RSRVCF -O $NAME.filtered_snps_final.vcf --filter-name "QD_filter" -filter "QD < 2.0" --filter-name "FS_filter" -filter "FS > 60.0" --filter-name "MQ_filter" -filter "MQ < 40.0" --filter-name "SOR_filter" -filter "SOR > 4.0" --filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" --filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

gatk VariantFiltration -R $RG -V $RIRVCF -O $NAME.filtered_indels_final.vcf --filter-name "QD_filter" -filter "QD < 2.0" --filter-name "FS_filter" -filter "FS > 200.0" --filter-name "SOR_filter" -filter "SOR > 10.0"

fsf='.filtered_snps_final.vcf'
NFSF=$NAME$fsf
fif='.filtered_indels_final.vcf'
NFIF=$NAME$fif

snpEff -v Plasmodium_falciparum -s SNP_summary.html $NFSF > $NAME.filtered_snps_final.ann.vcf
snpEff -v Plasmodium_falciparum -s Indel_summary.html $NFIF > $NAME.filtered_indels_final.ann.vcf
