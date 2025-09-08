#!/usr/bin/bash


# author: Alejandro Nabor Lozada-Chávez
# questions: nabor.lozada@gmail.com
#
# Script: correction of aligments for each single sample present in a single directory



# workding directory structure:
# 
# processing_samples/
# ├── programs
# └── snps_identification
#     ├── TMP_DIR
#     ├── alignments
#     ├── paired_mapped
#     ├── paired_unmapped
#     ├── quality_control
#     ├── raw_reads
#     ├── reference_genome
#     └── variant_calling




# Record the start time in seconds since the epoch
START_TIME=$(date +%s);



# // paths & programs //
# --------------------------------------------------
# working directory
workdir="$HOME/processing_samples/snps_identification";

# programs
java8=$workdir/...      # write full path directory location of java version 8 for GATK 3.8 
java22=$workdir/...     # write full path directory location of java version 22 for Picard 
gatk381=$workdir/...    # write full path directory location of GATK 3.8

# infiles
reference_genome=$workdir/reference_genome/my_fasta_seq.fasta;
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

for BAMFILE in $reads_dir
do
        echo
        echo "[1] // Recalibration of the alignment based on Indel realigments //"
        echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
       FILENAME=$(basename ${BAMFILE} .bam);
        
        # identify conflictive signal regions (run it using several threads):
        echo "$java8 -XX:ParallelGCThreads=31 -Xmx40G -jar $gatk381 -T RealignerTargetCreator -R $reference_genome -I $bams_path/$BAMFILE -o $bams_path/${FILENAME}.indels.intervals --num_threads 15";
        
        # main
        $java8 -XX:ParallelGCThreads=31 -Xmx40G -jar $gatk381 \
            -T RealignerTargetCreator \
            -R $reference_genome \
            -I $bams_path/$BAMFILE \
            -o $bams_path/${FILENAME}.indels.intervals \
            --num_threads 15
        wait;
        sleep 2;
        
        
        echo "[2] // Make correction on those regions //"
        echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
        echo "$java8 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=1 -Xmx40G -jar $gatk381 -T IndelRealigner -R $reference_genome -I $bams_path/$BAMFILE -known indels.db.vcf.gz -targetIntervals $bams_path/${FILENAME}.indels.intervals -o $bams_path/${FILENAME}.raln_indels.bam";
        
        # main:        
        $ $java8 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=1 -Xmx40G -jar $gatk381 \
            -T IndelRealigner \
            -R $reference_genome \
            -I $bams_path/$BAMFILE \
            -known indels.db.vcf.gz \                   # this line is optional
            -targetIntervals $bams_path/${FILENAME}.indels.intervals \
            -o $bams_path/${FILENAME}.raln_indels.bam
        wait;
        sleep 2;
        
        
        echo "[3] // First check: Verify errors after the indel realignment //"
        echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
        echo "$java22 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=20  -Xmx60G -jar $picard ValidateSamFile INPUT=$bams_path/${FILENAME}.raln_indels.bam OUTPUT=$bams_path/${FILENAME}.raln_indels.ValidateSam.txt MODE=SUMMARY TMP_DIR=$tmpdir";
        
        # main:
        $java22 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=20  -Xmx60G -jar $picard ValidateSamFile \
            INPUT=$bams_path/${FILENAME}.raln_indels.bam \
            OUTPUT=$bams_path/${FILENAME}.raln_indels.ValidateSam.txt \
            MODE=SUMMARY \
            TMP_DIR=$tmpdir
        wait;
        sleep 2;
        
        
        echo "[4] // Make correction if errors were found in the realignment //"
        echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
        echo "$java22 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=20  -Xmx60G -jar $picard FixMateInformation  INPUT=$bams_path/${FILENAME}.raln_indels.bam  OUTPUT=$bams_path/${FILENAME}.raln_indels.fm.bam  ADD_MATE_CIGAR=true";
        
        # main
        $java22 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=20 \
            -Xmx60G -jar $picard FixMateInformation \
            INPUT=$bams_path/${FILENAME}.raln_indels.bam \
            OUTPUT=$bams_path/${FILENAME}.raln_indels.fm.bam \
            ADD_MATE_CIGAR=true
        wait;
        sleep 2;
        
        
        echo "[5] // Second check: Verify errors after the indel realignment //"
        echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
        echo "$java22 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=20 -Xmx60G -jar $picard ValidateSamFile INPUT=$bams_path/${FILENAME}.raln_indels.fm.bam OUTPUT=$bams_path/${FILENAME}.raln_indels.fm.ValidateSam.txt MODE=SUMMARY TMP_DIR=$tmpdir";        
        
        # main:
        $java22 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=20 \
            -Xmx60G -jar $picard ValidateSamFile \
            INPUT=$bams_path/${FILENAME}.raln_indels.fm.bam \
            OUTPUT=$bams_path/${FILENAME}.raln_indels.fm.ValidateSam.txt \
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



        END_TIME=$(date +%s)
        
        # Calculate the elapsed time in seconds
        TOTAL_TIME_SECONDS=$((END_TIME - START_TIME))
        # Convert elapsed seconds to minutes
        TOTAL_TIME_MINUTES=$((TOTAL_TIME_SECONDS / 60))
        # Calculate remaining seconds for precise display
        REMAINING_SECONDS=$((TOTAL_TIME_SECONDS % 60))
       # Convert elapsed seconds to hours
        TOTAL_TIME_HOURS=$((TOTAL_TIME_SECONDS / 3600 ))
        wait;
        sleep 2;



        echo
        echo
        echo JOB_COMPLETED: `date`;
        echo
        echo
        echo "Start Time: $(date -d @$START_TIME)";
        echo "End Time: $(date -d @$END_TIME)";
        echo "Total Time: ${TOTAL_TIME_SECONDS} seconds (total elapsed: $REMAINING_SECONDS)";
        echo "Total Time: ${TOTAL_TIME_MINUTES} minutes";
        echo "Total Time: ${TOTAL_TIME} hours";
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
