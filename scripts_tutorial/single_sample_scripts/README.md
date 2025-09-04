## Single samples full worked example (worked example)

This is a single sample full protocol steps: from download reference genome and paired-ended reads, quality control of reads and alignments, and trimming, mapping, sort. dedupping, and realignment of indels in all aligned reads, until a variant calling  procedure from the realigned reads.

It is worth mentioning that *mapping* the reads and perform a *variant calling* procedure using a **single sample** will introduce a lot of bias from a perspective of a Population Genetics analysis, since one individual cannot be a full representative organism of a entire species specifically when its genetic variances, in this case based on SNPS, will be studied. Therefore, more samples should be included in this type of studies. There is not a maximum limit, but a thumbs up rule for a minimum number of samples could be 10 samples. This requirement, however, cannot be achieve due a several reasons (technical, biological, geographical, economical, etc.), so we understand that this is an ideal case and sometimes it cannot be possible.   

### Protcol

1. Create working directory.

```bash
# current location: $HOME // ~
cd $HOME
mkdir book_chapter_SNPs
cd book_chapter_SNPs;
mkdir programs snps_identification && cd snps_identification && mkdir reference_genome raw_reads mapped_unmapped alignments variant_calling;
cd $HOME;

# * WORKDIR STRUCTURE * 
tree book_chapter_SNPs/
# book_chapter_SNPs/
# ├── programs
# └── snps_identification
#     ├── alignments
#     ├── mapped_unmapped
#     ├── raw_reads
#     ├── reference_genome
#     └── variant_calling

# working directory
workdir="$HOME/book_chapter_SNPs/";

# define other variables:
mdkir tmpDIR;
tmpdir="$HOME/tmpDIR";

raw_reads="$workdir/raw_reads";
mapped_unmapped="$workdir/mapped_unmapped";
BAM_alignments="$workdir/alignments";
variants_DIR="$workdir/variant_calling";

java22="$workdir/programs/java22/bin/java";                     # for pircard, java
java8="$workdir/book_chapter_SNPs/programs/java8/bin/java";     # for gatk
trimmomatic="$workdir/programs/Trimmomatic-0.39/trimmomatic-0.39.jar";
picard="$workdir/programs/picard.jar";
gatk_v381="$workdir/programs/gatk_v381/GenomeAnalysisTK.jar";

cd $workdir/snps_identification/raw_reads;
```

2. Download a single sample and transform it into FASTQ file.

```bash
prefetch SRR11006725;
wait;
sleep 2;

fasterq-dump –progress --split 2 --verbose --threads 10 SRR11006725.sra
wait;
sleep 2;

# compress
ls file.*.fastq | parallel 'gzip {} > {}.gz'; 
```

3. Prepare the **Reference Genome**.

This step can be done while FASTQ files are created. This `reference genome` step has been fully described in in the script section. Just [follow the same procedure presented in this file](https://github.com/naborlozada/BookChapter_SNP_identification_tutorial/blob/main/scripts_tutorial/scripts_description/02_reference_genome.md).


4. Trimming adapters from alinged reads.


```bash
# Use latest JAVA version for trimmomatic

time $java22 -Djava.io.tmpdir=$tmpdir -jar $trimmomatic PE -threads 32 -phred33 \
            -trimlog $sortedDIR/SRR11006725.R1_R2.trimmomatic.logfile.txt \
            $raw_reads/SRR11006725.file_1.fastq.gz \
            $raw_reads/SRR11006725.file_2.fastq.gz \
            $mapped_unmapped/SRR11006725.R1.paired.trimclip.fastq.gz \
            $mapped_unmapped/SRR11006725.R1.unpaired.trimclip.fastq.gz \
            $mapped_unmapped/SRR11006725.R2.paired.trimclip.fastq.gz \
            $mapped_unmapped/SRR11006725.R2.unpaired.trimclip.fastq.gz \
            ILLUMINACLIP:$workdir/programs/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
wait;

```

5. Alignment of reads against reference genome sequence.

```bash

time bwa mem -t 32 \
            $reference_genome \
            $mapped_unmapped/SRR11006725.R1.paired.trimclip.fastq.gz \
            $mapped_unmapped/SRR11006725.R2.paired.trimclip.fastq.gz \
            -R '@RG\tID:'SRR11006725'\tSM:'SRR11006725'\tPL:illumina' | samtools view -b --threads 32 - > $BAM_alignments/SRR11006725.paired.aln_raw.bam
wait;

```

6. Sort aligned reads in the BAM file. 

```bash
time $java22 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=32 -Xmx60G -jar $picard SortSam \
            INPUT=$BAM_alignments/SRR11006725.paired.aln_raw.bam \
            OUTPUT=$BAM_alignments/SRR11006725.paired.aln_sort.bam \
            SORT_ORDER=coordinate \
            TMP_DIR=$tmpdir
wait;
```

7. Deduplicate BAM file.

```bash
time $java22 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=32 -Xmx60G -jar $picard MarkDuplicates \
            INPUT=$BAM_alignments/SRR11006725.paired.aln_sort.bam \
            OUTPUT=$BAM_alignments/SRR11006725.paired.dedup.bam \
            METRICS_FILE=$BAM_alignments/SRR11006725.paired.dedup.metrics.txt \
            TMP_DIR=$tmpdir
wait;
```

8. Make indexes for the `deduped` BAM file.

```bash
# create indexes:
# bwa index
bwa index $BAM_alignments/SRR11006725.paired.dedup.bam
wait;

# samtools index
samtools index -@ 20 $BAM_alignments/SRR11006725.paired.dedup.bam
wait;

```

9. Check alignment file integrity: verify for errors in the alignment.

```bash
time $java22 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=32 -Xmx60G -jar $picard ValidateSamFile \
            I=$BAM_alignments/SRR11006725.paired.dedup.bam \
            O=$BAM_alignments/SRR11006725.paired.dedup.ValidateSamFile.txt \
            MODE=SUMMARY \
            TMP_DIR=$tmpdir
wait;
```

10. Correction in the alignment: indels re-alignment.

This procedure involves two steps: (a) indentify all low indel regions, and (b) realign those regions. For the first step we use the and second steps we use the functions `RealignerTargetCreator` and  `IndelRealigner` from GATK v.3.8.1, respectively. 

```bash
# First step (multi-thread):
time $java8 -XX:ParallelGCThreads=32 -Xmx40G -jar $gatk_v381 \
    -T RealignerTargetCreator \
    -R $reference_genome \
    -I $BAM_alignments/SRR11006725.paired.dedup.bam \
    -known $knownINDELs \           # ---> optional step if you have a 'golden indels dataset'.
    -o $BAM_alignments/SRR11006725.paired.dedup.indels.intervals \
    --num_threads 32

wait;

# Second step (single-thread):
time java8 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=1 -Xmx5G -jar $gatk_v381 \
    -T IndelRealigner \
    -R $reference_genome \
    -I $BAM_alignments/SRR11006725.paired.dedup.bam \
    -known $knownINDELs \           # ---> optional step if you have a 'golden indels dataset'.
    -targetIntervals $BAM_alignments/SRR11006725.paired.dedup.indels.intervals \
    -o $BAM_alignments/SRR11006725.paired.dedup.indels.raln.bam

wait;

# make an index of the latest BAM file.
time $samtools index -@ 20 $BAM_alignments/SRR11006725.paired.dedup.indels.raln.bam
```

11. Fix mate information.

```bash
time $java22 -Djava.io.tmpdir=$tmpdir -Xmx40G -XX:+UseParallelGC -XX:ParallelGCThreads=16 -jar $picard FixMateInformation \
    INPUT=$BAM_alignments/SRR11006725.paired.dedup.indels.raln.bam \
    OUTPUT=$BAM_alignments/SRR11006725.paired.dedup.indels.raln.fixed_mate.bam \
    ADD_MATE_CIGAR=true \
    TMP_DIR=$tmpdir

wait;

# make index
time $samtools index -@ 16 $BAM_alignments/SRR11006725.paired.dedup.indels.raln.fixed_mate.bam
wait;
```

12. Variant calling using BCFTOOLS.


```bash

# [1] call variants
bcftools mpileup -f $reference_genome  $BAM_alignments/SRR11006725.paired.dedup.indels.raln.bam | bcftools call -mv -Oz -o $variants_DIR/SRR11006725.aln.srtd.dedup.raln.fm.vcf.gz;
wait;
sleep 2;

# [2] Suggestion (recommended): Remove SNPs with a close proximity to indes of 10 bp
bcftools filter --threads 3 --SnpGap 10 -Oz  -o $variants_DIR/SRR11006725.closeprox_rmvd.vcf.gz  $variants_DIR/SRR11006725.aln.srtd.dedup.raln.fm.vcf.gz
wait;


# [3] Extract only SNP variants:
bcftools view -v snps  $variants_DIR/SRR11006725.closeprox_rmvd.vcf.gz  -Oz -o $variants_DIR/SRR11006725.aln.srtd.dedup.raln.fm.snps.raw.vcf.gz
wait;
sleep 2;


# filter SNPs by retaining only SNPs with a minor allele frequency of 0.01 AND mapping quality bigger than 20:
bcftools view -q 0.01:minor  $variants_DIR/SRR11006725.aln.srtd.dedup.raln.fm.snps.raw.vcf.gz | bcftools filter -i '%QUAL>20' -Oz -o $variants_DIR/SRR11006725.aln.srtd.dedup.raln.fm.snps.fltrd.vcf.gz
wait;
sleep 2;


# Extract only biallelic SNPs
bcftools -m2 -M2  -Oz -o $variants_DIR/SRR11006725.aln.srtd.dedup.raln.fm.snps.fltrd.biallelic.vcf.gz  $variants_DIR/SRR11006725.aln.srtd.dedup.raln.fm.snps.fltrd.vcf.gz
wait;
sleep 2;

# Count number of SNPs: first, index the VCF file (option: index -t), then count then (option: index -n). This will create a file with the same name of the VCF file, but with a file extension *vcf.gz.tbi
bcftools index -t $variants_DIR/SRR11006725.aln.dedup.raln_indels.fm.snps.fltrd.biallelic.vcf.gz
wait;
sleep 2;

# count SNPs:
bcftools index -n $variants_DIR/SRR11006725.aln.dedup.raln_indels.fm.snps.fltrd.biallelic.vcf.gz
```

### Quality samples: fastq and bam files

Here, this step can be done in parallel after the step 2 (compressed fastq.gz files) and 10 (dedupped alignment). 

A) Quality of the reads using `FASTQC`.

```bash
cd $workdir;
mkdir quality_control && cd quality_control;

fastqc $mapped_unmapped/SRR11006725.R1.paired.trimclip.fastq.gz  $mapped_unmapped/SRR11006725.R2.paired.trimclip.fastq.gz  -o  $workdir/quality_control;
```
 and involve two anlaysis

B) Statistics for aligned reads `samtools`.

```bash
samtools flagstats $BAM_alignments/SRR11006725.paired.dedup.bam > $BAM_alignments/SRR11006725.paired.dedup.samtools_stats.txt;
```

c) Merge all quality data and statistics into one single graphical report using `MULTIQC`. After running this command, you will open the HTML report.

```bash
# assuming your current workdir is $workdir/quality_control
multiqc .
```
