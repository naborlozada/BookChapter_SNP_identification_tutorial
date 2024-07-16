#!/usr/bin/bash


# variables:
# $tmpdir     : location of a directory where temporal files will be created. 
# $trimmomatic : location of the binary program 'trimmomatic'.
# $java22      : location of the binary java file (version 22).
# $picard      : location of Picard program.  

# $reference_DIR   : location of the reference genome (fasta) and associated index files.
# $alignments_DIR  : location of BAM alignment files.
# $reads_DIR       : location of the *.FASTQ.gz files (raw reads).
# $paired_unpaired : location of the paired/unpaired reads files (mate1, R1; mate2, R2).

# Single FASTQ pair files (mate 1 and mate 2). Replace these names for your samples' names or defined in a variable, as described in the book chapter.
#   sample1.R1.fastq.gz
#   sample1.R2.fastq.gz



# trimming
# -----------------------------------------------------------------------------
# Define adapters directory
adapters_DIR="$trimmomatic/adapters";

# trim adapters
time $java22 -Djava.io.tmpdir=$tmpdir -jar $trimmomatic PE \
    -threads 32 -phred33 \
    -trimlog $paired_unpaired/sample1.R1_R2.trimmomatic.logfile.txt \
    $reads_DIR/sample1.R1.fastq.gz \
    $reads_DIR/sample1.R2.fastq.gz \
    $paired_unpaired/sample1.R1.paired.trimclip.fastq.gz \
    $paired_unpaired/sample1.R1.unpaired.trimclip.fastq.gz \
    $paired_unpaired/sample1.R2.paired.trimclip.fastq.gz \
    $paired_unpaired/sample1.R2.unpaired.trimclip.fastq.gz \
    ILLUMINACLIP:$adapters_DIR/TruSeq3-PE.fa:2:30:10
wait;
sleep 2;


# make alignment
# -----------------------------------------------------------------------------
bwa mem -t 20 \
    $reference_DIR/reference_parsed.fasta \
    $paired_unpaired/sample1.R1.paired.trimclip.fastq.gz \
    $paired_unpaired/sample1.R2.paired.trimclip.fastq.gz \
    -R '@RG\tID:'sample1'\tSM:'sample1'\tPL:illumina' |\
    samtools view -b --threads 20 - > $alignments_DIR/sample1.aln.raw.bam
wait;
sleep 2;


# sort alignment
# -----------------------------------------------------------------------------
$java22 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=20 \
    -Xmx60G -jar $picard SortSam \
    INPUT=$alignments_DIR/sample1.aln.raw.bam \
    OUTPUT=$alignments_DIR/sample1.aln.sorted.bam \
    SORT_ORDER=coordinate \
    TMP_DIR=$tmpdir
wait;
sleep 2;


# Dedupping alignment
# -----------------------------------------------------------------------------
$java22 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=32 \
    -Xmx60G -jar $picard MarkDuplicates \
    INPUT=$alignments_DIR/sample1.aln.sorted.bam \
    OUTPUT=$alignments_DIR/sample1.aln.dedup.bam \
    METRICS_FILE=$alignments_DIR/sample1.aln.dedup.metrics.txt \
    TMP_DIR=$tmpdir
wait;
sleep 2;


# compress the *.txt file:
gzip $alignments_DIR/sample1.aln.dedup.metrics.txt;
wait;
sleep 2;


# Verify errors in deduped alignment
# -----------------------------------------------------------------------------
$java22 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=20 \
    -Xmx60G -jar $picard ValidateSamFile \
    INPUT=$alignments_DIR/sample1.aln.dedup.bam \
    OUTPUT=$alignments_DIR/sample1.aln.dedup.ValidateSam.txt \
    MODE=SUMMARY \
    TMP_DIR=$tmpdir
wait;
sleep 2;

