#!/usr/bin/bash


# ------------------------------------------------ #
#      *** Alejandro Nabor Lozada-Chávez ***       #
# Questions: alejandro.chavez@unipv.it             #
#            nabor.lozada@gmail.com                #
# ------------------------------------------------ #

# Script: correction of aligments for each single sample present in a single directory



#SBATCH --account=
#SBATCH --partition=
#SBATCH --cpus-per-task=20
#SBATCH --mem=
#SBATCH --ntasks=1
#SBATCH --job-name=
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --time=
#SBATCH --error=quality_control_fastqc_n_samtools.report.%A_%a.stderr.log
#SBATCH --output=quality_control_fastqc_n_samtools.report.%A_%a.stderr.log




echo
echo
echo "            **** SBATCH JOBS SUBMISSION wih SLURM system ****                 "
echo 
echo

echo "------------------------------------------------------------------------------"
echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "JOB submission name            = $SLURM_JOB_NAME"
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"
echo "User account name              = $SLURM_JOB_ACCOUNT"
echo "Partition name                 = $SLURM_JOB_PARTITION"
echo "Export environment DIR         = $SLURM_EXPORT_ENV"
echo "Task PID                       = $SLURM_TASK_PID $SLURM_JOB_ID"
echo "------------------------------------------------------------------------------"
echo
echo
echo



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
bam_files=$workdir/alignments/*.paired.dedup_reads.bam;
bams_path=$workdir/alignments/;
tmpdir=$workdir/TMP_DIR/;

# infiles
BAM_FILES=($bam_files)
BAM_FILE=${BAM_FILES[$SLURM_ARRAY_TASK_ID]}
FILENAME="$(basename ${FILE1} .bam)"




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

echo
echo "[1] // Recalibration of the alignment based on Indel realigments //"
echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";

# identify conflictive signal regions (run it using several threads):
echo "$java8 -XX:ParallelGCThreads=31 -Xmx40G -jar $gatk381 -T RealignerTargetCreator -R $reference_genome -I $bams_path/$BAM_FILE -o $bams_path/${FILENAME}.indels.intervals --num_threads 1";

# main
$java8 -XX:ParallelGCThreads=31 -Xmx40G -jar $gatk381 \
    -T RealignerTargetCreator \
    -R $reference_genome \
    -I $bams_path/$BAM_FILE \
    -o $bams_path/${FILENAME}.indels.intervals \
    --num_threads 1
wait;
sleep 2;



echo "[2] // Make correction on those regions //"
echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
echo "$java8 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=1 -Xmx40G -jar $gatk381 -T IndelRealigner -R $reference_genome -I $bams_path/$BAM_FILE -known indels.db.vcf.gz -targetIntervals $bams_path/${FILENAME}.indels.intervals -o $bams_path/${FILENAME}.raln_indels.bam";

# main:        
$ $java8 -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=1 -Xmx40G -jar $gatk381 \
    -T IndelRealigner \
    -R $reference_genome \
    -I $bams_path/$BAM_FILE \
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
 Convert elapsed seconds to hours
TOTAL_TIME_HOURS=$((TOTAL_TIME_SECONDS / 3600 ))
wait;
sleep 2;


echo
echo
echo "####################################################################################################################################################################################"
echo
echo
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
