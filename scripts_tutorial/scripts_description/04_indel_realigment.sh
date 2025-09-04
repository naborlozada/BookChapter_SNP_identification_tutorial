#!/usr/bin/bash



# // Correction of aligments for each single sample in a single Directory //
# ------------------------------------------------------------------------------------------------------------------------------------------------------- #
# Bash script to process ONE or MORE alignments in a LOOP, that is, a single run per paired reads sample.
# This script is particularly useful when there is only 1 computer with not much computational power or Server with low numer of threads accesible.
#
# NOTE 1: MULTI-THREADS (or multithreading): 
# Overall, this script operate with 1 thread/core use. However, the first step uses multi-threads: "RealignerTargetCreator" (--num_threads NUMBER)
# Where NUMBER is the number of threads. To control this just add desire number of threads to be used. Before running this script, please, check the number 
# of threads available.
# Additionally, all java process here use a multi-threads option named "Garbage Collector (GC)" (-XX:-UseParallelGC), and basically it handles the startup memory
# of each process. To control the number of threads of this process change the number in the option: -XX:ParallelGCThreads=31
# More information: https://sematext.com/java-garbage-collection-tuning/ 
# 
# NOTE 2: Set properly all directory full path names for programs, scripts, and files used to run this bash script.
#
# Run this script as follows (better under a screen linux session):
#       nohup bash 04_indel_realignment.sh 2>> 04_indel_realignment.loop_job.stderr.log &> 04_indel_realignment.loop_job.stderr.log &
# ------------------------------------------------------------------------------------------------------------------------------------------------------- #

# variables:
# $tmpdir     : location of a directory where temporal files will be created. 
# $java8      : location of the binary java file (version 8)
# $gatk381    : location of GATK (version 3.8.1)  

# $reference_DIR  : location of the reference genome (fasta) and associated index files.
# $alignments_DIR : location of BAM alignment files.



# // paths & programs //
# --------------------------------------------------
# working directory
workdir="$HOME/book_chapter_SNPs";

# programs
java8="$workdir/programs/java8/bin/java8.jar";            # write full path directory location of java version 8 for GATK 3.8 
java22="$workdir/programs/java22/bin/java22.jar";         # write full path directory location of java version 22 for Picard 
gatk381="$workdir/programs/gatk_v38/bin/gatk_v38.jar";    # write full path directory location of GATK 3.8

# infiles
reference_genome=$workdir/reference_genome/Aedes-aegypti-LVP_AGWG_CHROMOSOMES.AaegL5_2.fasta;
bams_dir=$workdir/alignments/*.paired.dedup_reads.bam;
bams_path=$workdir/alignments/;
tmpdir=$workdir/TMP_DIR/;




echo
echo
echo
echo JOB_STARTED: `date`
echo
echo `pwd`
echo
echo
echo
echo "// [MAIN] PROCESSING REALINGMENT CORRECTION PER SAMPLE //"
echo "####################################################################################################################################################################################"

for BAMFILE in $bams_dir
do
        echo
        echo "[1] // Recalibration of the alignment based on Indel realigments //"
        echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
        FILENAME=$(basename ${BAMFILE} .bam);
        
        # identify conflictive signal regions (run it using several threads):
        echo "$java8 -XX:ParallelGCThreads=31 -Xmx40G -jar $gatk381 -T RealignerTargetCreator -R $reference_genome -I $BAMFILE -o $workdir/$bams_path/${FILENAME}.indels.intervals --num_threads 15";
        
        # main
        $java8 -XX:ParallelGCThreads=31 -Xmx40G -jar $gatk381 \
            -T RealignerTargetCreator \
            -R $reference_genome \
            -I $BAMFILE \
            -o $workdir/$bams_path/${FILENAME}.indels.intervals \
            --num_threads 15
        wait;
        sleep 2;
        
        
        echo "[2] // Make correction on those regions //"
        echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
        echo "$java8 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=1 -Xmx40G -jar $gatk381 -T IndelRealigner -R $reference_genome -I $BAMFILE -known indels.db.vcf.gz -targetIntervals $workdir/$bams_path/${FILENAME}.indels.intervals -o $workdir/$bams_path/${FILENAME}.raln_indels.bam";
        
        # main:        
        $ $java8 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=1 -Xmx40G -jar $gatk381 \
            -T IndelRealigner \
            -R $reference_genome \
            -I $BAMFILE \
            -known indels.db.vcf.gz \                   # this line is optional
            -targetIntervals $workdir/$bams_path/${FILENAME}.indels.intervals \
            -o $workdir/$bams_path/${FILENAME}.raln_indels.bam
        wait;
        sleep 2;
        
        
        echo "[3] // First check: Verify errors after the indel realignment //"
        echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
        echo "$java22 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=20  -Xmx60G -jar $picard ValidateSamFile INPUT=$workdir/$bams_path/${FILENAME}.raln_indels.bam OUTPUT=$workdir/$bams_path/${FILENAME}.raln_indels.ValidateSam.txt MODE=SUMMARY TMP_DIR=$tmpdir";
        
        # main:
        $java22 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=20  -Xmx60G -jar $picard ValidateSamFile \
            INPUT=$workdir/$bams_path/${FILENAME}.raln_indels.bam \
            OUTPUT=$workdir/$bams_path/${FILENAME}.raln_indels.ValidateSam.txt \
            MODE=SUMMARY \
            TMP_DIR=$tmpdir
        wait;
        sleep 2;
        
        
        echo "[4] // Make correction if errors were found in the realignment //"
        echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
        echo "$java22 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=20  -Xmx60G -jar $picard FixMateInformation  INPUT=$workdir/$bams_path/${FILENAME}.raln_indels.bam  OUTPUT=$workdir/$bams_path/${FILENAME}.raln_indels.fm.bam  ADD_MATE_CIGAR=true";
        
        # main
        $java22 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=20 \
            -Xmx60G -jar $picard FixMateInformation \
            INPUT=$workdir/$bams_path/${FILENAME}.raln_indels.bam \
            OUTPUT=$workdir/$bams_path/${FILENAME}.raln_indels.fm.bam \
            ADD_MATE_CIGAR=true
        wait;
        sleep 2;
        
        
        echo "[5] // Second check: Verify errors after the indel realignment //"
        echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
        echo "$java22 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=20 -Xmx60G -jar $picard ValidateSamFile INPUT=$workdir/$bams_path/${FILENAME}.raln_indels.fm.bam OUTPUT=$workdir/$bams_path/${FILENAME}.raln_indels.fm.ValidateSam.txt MODE=SUMMARY TMP_DIR=$tmpdir";        
        
        # main:
        $java22 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=20 \
            -Xmx60G -jar $picard ValidateSamFile \
            INPUT=$workdir/$bams_path/${FILENAME}.raln_indels.fm.bam \
            OUTPUT=$workdir/$bams_path/${FILENAME}.raln_indels.fm.ValidateSam.txt \
            MODE=SUMMARY \
            TMP_DIR=$tmpdir
        wait;
        sleep 2;
        
        # NOTE: The content of the 'second check' file should be "No errors found.". 


        echo
        endtime=`date +"%s"`;
        duration_sec=$((endtime - starttime));
        duration_min=$(( $duration_sec / 60 ));
        duration_hrs=$( echo $duration_sec / 3600 | bc -l );    # add precision count
        wait;
        sleep 2;


        echo
        echo
        echo JOB_ENDED: `date`
        echo
        echo
        echo "STAT:startTime:$starttime";
        echo "STAT:doneTime:$endtime";
        echo "STAT:runtime_sec:$duration_sec";
        echo "STAT:runtime_min:$duration_min";
        echo "STAT:runtime_hrs:$duration_hrs";
        echo
        echo

        endtime=0;
        duration_sec=0;
        duration_min=0;
        duration_hrs=0;
done

echo
echo
echo "####################################################################################################################################################################################"
echo
echo
echo JOB_ENDED: `date`
echo
echo
