#!/usr/bin/bash


# variables:
# $tmpdir     : location of a directory where temporal files will be created. 
# $java8      : location of the binary java file (version 8)
# $gatk381    : location of GATK (version 3.8.1)  

# $reference_DIR  : location of the reference genome (fasta) and associated index files.
# $alignments_DIR : location of BAM alignment files.



# Recalibration of the alignment based on Indel realigments
# -----------------------------------------------------------------------------

# Dedup BAM (script 3) as arget file: $alignments_DIR/sample1.aln.dedup.bam

# identify conflictive signal regions (run it using several threads):
$java8 -XX:ParallelGCThreads=31 -Xmx40G -jar $gatk381 \
    -T RealignerTargetCreator \
    -R $reference_DIR/reference_parsed.fasta \
    -I $alignments_DIR/sample1.aln.dedup.bam \
    -o $alignments_DIR/sample1.aln.dedup.indels.intervals \
    --num_threads 15
wait;
sleep 2;


# Make correction on those regions:
# -----------------------------------------------------------------------------
$ $java8 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=1 -Xmx40G -jar $gatk381 \
    -T IndelRealigner \
    -R $reference_DIR/reference_parsed.fasta \
    -I $alignments_DIR/sample1.aln.dedup.bam \
    -known indels.db.vcf.gz \                   # this line is optional
    -targetIntervals $alignments_DIR/sample1.aln.dedup.indels.intervals \
    -o $alignments_DIR/sample1.aln.dedup.raln_indels.bam
wait;
sleep 2;


# First check: Verify errors after the indel realignment
# -----------------------------------------------------------------------------
$java22 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=20 \
    -Xmx60G -jar $picard ValidateSamFile \
    INPUT=$alignments_DIR/sample1.aln.dedup.raln_indels.bam \
    OUTPUT=$alignments_DIR/sample1.aln.dedup.raln_indels.ValidateSam.txt \
    MODE=SUMMARY \
    TMP_DIR=$tmpdir
wait;
sleep 2;


# Make correction if errors were found in the realignment
# -----------------------------------------------------------------------------
$java22 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=20 \
    -Xmx60G -jar $picard FixMateInformation \
    INPUT=$alignments_DIR/sample1.aln.dedup.raln_indels.bam \
    OUTPUT=$alignments_DIR/sample1.aln.dedup.raln_indels.fm.bam \
    ADD_MATE_CIGAR=true
wait;
sleep 2;


# Second check: Verify errors after the indel realignment
# -----------------------------------------------------------------------------
$java22 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=20 \
    -Xmx60G -jar $picard ValidateSamFile \
    INPUT=$alignments_DIR/sample1.aln.dedup.raln_indels.fm.bam \
    OUTPUT=$alignments_DIR/sample1.aln.dedup.raln_indels.fm.ValidateSam.txt \
    MODE=SUMMARY \
    TMP_DIR=$tmpdir
wait;
sleep 2;

# The content of the second check file should be "No errors found.". 
