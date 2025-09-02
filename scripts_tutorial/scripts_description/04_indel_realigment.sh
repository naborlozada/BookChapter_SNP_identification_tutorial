#!/usr/bin/bash



# // Correction of aligments for each single sample in a single Directory //

# Bash script to process one or several alignments in a LOOP, that is, a single run per paired reads sample.
# This is particularly useful when there is only 1 computer and/or 1 server several threads that can be used.
# It is important to mention that the first two steps are the main processes that run in parallel, while the rest does run with one thread.
# IMPORTANT: Do the proper changes in the directory full path names for programs/scripts ('programs' section) used to run this bash script,  
#            and also for the infiles, outfiles (if required), and output directory names. 
#
# Run this script as follow (better under a screen linux session):
#       nohup bash 04_indel_realignment.sh 2>> 04_indel_realignment.loop_job.stderr.log &> 04_indel_realignment.loop_job.stderr.log &


# variables:
# $tmpdir     : location of a directory where temporal files will be created. 
# $java8      : location of the binary java file (version 8)
# $gatk381    : location of GATK (version 3.8.1)  

# $reference_DIR  : location of the reference genome (fasta) and associated index files.
# $alignments_DIR : location of BAM alignment files.



# // paths & programs //
# --------------------------------------------------
# working directory
workdir=/home/username/book_snps_identification

# programs
java8=/programs/java8/bin/java8.jar            # write full path directory location of java version 8 for GATK 3.8 
java22=/programs/java22/bin/java22.jar         # write full path directory location of java version 22 for Picard 
gatk381=/programs/gatk_v38/bin/gatk_v38.jar    # write full path directory location of GATK 3.8

# infiles
reference_genome=/home/username/book_snps_identification/reference_genome/Aedes-aegypti-LVP_AGWG_CHROMOSOMES.AaegL5_2.fasta;
bams_dir=/home/username/alignments/*.paired.dedup_reads.bam;
bams_path=/home/username/alignments/;
tmpdir=/home/username/TMP_DIR/;




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

for BAMFILE in $reads_dir
do
        echo
        echo "[1] // Recalibration of the alignment based on Indel realigments //"
        echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
       FILENAME=$(basename ${BAMFILE} .bam);
        
        # identify conflictive signal regions (run it using several threads):
        echo "$java8 -XX:ParallelGCThreads=31 -Xmx40G -jar $gatk381 -T RealignerTargetCreator -R $reference_genome -I $workdir/$bams_path/$BAMFILE -o $workdir/$bams_path/${FILENAME}.indels.intervals --num_threads 15";
        
        # main
        $java8 -XX:ParallelGCThreads=31 -Xmx40G -jar $gatk381 \
            -T RealignerTargetCreator \
            -R $reference_genome \
            -I $workdir/$bams_path/$BAMFILE \
            -o $workdir/$bams_path/${FILENAME}.indels.intervals \
            --num_threads 15
        wait;
        sleep 2;
        
        
        echo "[2] // Make correction on those regions //"
        echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
        echo "$java8 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=1 -Xmx40G -jar $gatk381 -T IndelRealigner -R $reference_genome -I $workdir/$bams_path/$BAMFILE -known indels.db.vcf.gz -targetIntervals $workdir/$bams_path/${FILENAME}.indels.intervals -o $workdir/$bams_path/${FILENAME}.raln_indels.bam";
        
        # main:        
        $ $java8 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=1 -Xmx40G -jar $gatk381 \
            -T IndelRealigner \
            -R $reference_genome \
            -I $workdir/$bams_path/$BAMFILE \
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
        
        # The content of the second check file should be "No errors found.". 


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
