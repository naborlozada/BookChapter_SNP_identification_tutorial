#!/usr/bin/bash


# ------------------------------------------------ #
#      *** Alejandro Nabor Lozada-Chávez ***       #
# Questions: alejandro.chavez@unipv.it             #
#            nabor.lozada@gmail.com                #
# ------------------------------------------------ #

# Script: performs a trimming and alignment of raw reads on a reference genome, while 
#         at the end sort and dedup all aligned reads present in a single directory


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




# working directory
workdir="$HOME/processing_samples/snps_identification";



# // paths & programs //
# --------------------------------------------------
# programs
trimmomatic=$workdir/programs/Trimmomatic-0.39/trimmomatic-0.39.jar
trimmomatic_adapters=$workdir/programs/Trimmomatic-0.39/adapters

SAMTOOLS=$workdir/programs/samtools/samtools
BCFTOOLS=$workdir/programs/bcftools/bcftools
PICARD=$workdir/programs/picard.jar
BWA=$workdir/programs/bwa/bwa

# infiles
reference_genome=$workdir/reference_genome/my_reference_genomes.fasta;
RAW_READS_R1_FILES=$workdir/raw_reads/*1.fastq.gz;
reads_path=$workdir/raw_reads;
tmpdir=$workdir/TMP_DIR/;

# outfiles
paired_dir_path=$workdir/paired_mapped;
unpaired_dir_path=$workdir/paired_unmapped;
alignments=$workdir/alignments;


# infiles
FILES1=($RAW_READS_R1_FILES)
FILE1=${FILES1[$SLURM_ARRAY_TASK_ID]}
FILE2=$(basename ${FILE1} 1.fastq.gz)_2.fastq.gz

newfile="$(basename $FILE1 1.fastq.gz)"




echo
echo
echo
echo JOB_STARTED: `date`
echo
echo `pwd`
echo
echo
echo
echo "##############################################################################################################################"
echo
echo
echo "READS FILES R1 R2: $FILE1 $FILE2 ..."
echo 
# -----------------------------------------------------------------------------
echo "[1] // Trimmomatic: Trimming reads adapters //"
echo "command_line: java -Djava.io.tmpdir=$tmpdir -jar $trimmomatic PE -threads 32 -phred33  -trimlog $paired_dir_path/${newfile}.R1_R2.trimmomatic.logfile.txt  $FILE1  $reads_path/$FILE2  $paired_dir_path/${newfile}.R1.paired.trimclip.fastq.gz  $unpaired_dir_path/${newfile}.R1.unpaired.trimclip.fastq.gz  $paired_dir_path/${newfile}.R2.paired.trimclip.fastq.gz $unpaired_dir_path/${newfile}.R2.unpaired.trimclip.fastq.gz ILLUMINACLIP:$trimmomatic_adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
echo

time java -Djava.io.tmpdir=$tmpdir -jar $trimmomatic PE -threads 32 -phred33 \
            -trimlog $paired_dir_path/${newfile}.R1_R2.trimmomatic.logfile.txt \
            $FILE1 \
            $reads_path/$FILE2 \
            $paired_dir_path/${newfile}.R1.paired.trimclip.fastq.gz \
            $unpaired_dir_path/${newfile}.R1.unpaired.trimclip.fastq.gz \
            $paired_dir_path/${newfile}.R2.paired.trimclip.fastq.gz \
            $unpaired_dir_path/${newfile}.R2.unpaired.trimclip.fastq.gz \
            ILLUMINACLIP:$trimmomatic_adapters/TruSeq3-PE.fa:2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
wait;
sleep 2;



# -----------------------------------------------------------------------------
echo
echo "[2] // BWA: Alignment of reads to  he reference genome // "
echo "command_line: $BWA mem -t 32  $reference_genome  $paired_dir_path/${newfile}.R1.paired.trimclip.fastq.gz  $paired_dir_path/${newfile}.R2.paired.trimclip.fastq.gz  -R '@RG\tID:'$newfile'\tSM:'$newfile'\tPL:illumina' | $SAMTOOLS view -b --threads 32 - > $alignments/${newfile}.paired.aligned_raw_reads.bam "
echo

time $BWA mem -t 32 \
            $reference_genome \
            $paired_dir_path/${newfile}.R1.paired.trimclip.fastq.gz \
            $paired_dir_path/${newfile}.R2.paired.trimclip.fastq.gz \
            -R '@RG\tID:'$newfile'\tSM:'$newfile'\tPL:illumina' | $SAMTOOLS view -b --threads 32 - > $alignments/${newfile}.paired.aligned_raw_reads.bam
wait;
sleep 2;



# -----------------------------------------------------------------------------
echo
echo "[3] // SortSam: Sort alignment // "
echo "command_line: java -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=32 -Xmx60G -jar $PICARD SortSam  INPUT=$alignments/${newfile}.paired.aligned_raw_reads.bam  OUTPUT=$alignments/${newfile}.paired.aligned_reads_sorted.bam  SORT_ORDER=coordinate  TMP_DIR=$tmpdir "
echo

time java -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=32 -Xmx60G -jar $PICARD SortSam \
            INPUT=$alignments/${newfile}.paired.aligned_raw_reads.bam \
            OUTPUT=$alignments/${newfile}.paired.aligned_reads_sorted.bam \
            SORT_ORDER=coordinate \
            TMP_DIR=$tmpdir
wait;
sleep 2;



# -----------------------------------------------------------------------------
echo
echo "[4] // MarkDuplicates //"
echo "command_line:  java -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=32 -Xmx60G -jar $PICARD MarkDuplicates  INPUT=$alignments/${newfile}.paired.aligned_reads_sorted.bam  OUTPUT=$dedupDIR/${newfile}.paired.dedup_reads.bam  METRICS_FILE=$dedupDIR/${newfile}.paired.dedup_reads.metrics.txt  TMP_DIR=$tmpdir "
echo

time java -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=32 -Xmx60G -jar $PICARD MarkDuplicates \
            INPUT=$alignments/${newfile}.paired.aligned_reads_sorted.bam \
            OUTPUT=$alignments/${newfile}.paired.dedup_reads.bam \
            METRICS_FILE=$alignments/${newfile}.paired.dedup_reads.metrics.txt \
            TMP_DIR=$tmpdir
wait;
sleep 2;



# -----------------------------------------------------------------------------
echo
echo "[5] // Make BWA-based index... //"
echo "command_line: $BWA index $alignments/${newfile}.paired.dedup_reads.bam "
echo

$BWA index $alignments/${newfile}.paired.dedup_reads.bam
wait;
sleep 2;



# -----------------------------------------------------------------------------
echo
echo "[6] // Make samtools-based index... //"
echo "command_line: $SAMTOOLS index -@ 64 $alignments/${newfile}.paired.dedup_reads.bam"
echo

$SAMTOOLS index -@ 64 $alignments/${newfile}.paired.dedup_reads.bam
wait;
sleep 2;



# -----------------------------------------------------------------------------
echo
echo "[7] // Validate SAM/BAM error in alignment... //"
echo "command_line: java -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=32 -Xmx60G -jar $PICARD ValidateSamFile  I=$dedupDIR/${newfile}.paired.dedup_reads.bam  O=$dedupDIR/${newfile}.paired.dedup_reads.ValidateSamFile.txt  MODE=SUMMARY  TMP_DIR=$tmpdir "
echo

time java -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=32 -Xmx60G -jar $PICARD ValidateSamFile \
            I=$dedupDIR/${newfile}.paired.dedup_reads.bam \
            O=$dedupDIR/${newfile}.paired.dedup_reads.ValidateSamFile.txt \
            MODE=SUMMARY \
            TMP_DIR=$tmpdir
wait;
sleep 2;



# -----------------------------------------------------------------------------
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
echo "##############################################################################################################################"
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
echo
