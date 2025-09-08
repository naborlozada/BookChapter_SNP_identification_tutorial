
#!/usr/bin/bash


# author: Alejandro Nabor Lozada-Chávez
# questions: nabor.lozada@gmail.com
#
# Script: Series of command lines used to perform QUALITY CONTROL analyses for FASTQC raw reads and BAM files 
#         (aligned reads, either dedupped or realigned), which are present in a single directory. The script will run all analysis
#         taking advantage of the power of the `GNU parallel` program.




START_TIME=$(date +%s);

# A single directory can contain all raw reads and BAM files or they can be separated. They have to be declare/defined in the
# variables below. FASTQ extension files are sometimes different (e.g. *.fq, *.fastq, *.fq.gz, *.fastq.gz), therefore, it needs to be
# defined in the command lines below for FASTQC and SAMtools programs.
RAW_READS=/path_to/fastq_files/
BAM_FILES=/path_to/bam_files/

# Define the number of threads to be used to run in parallel.
N_THREADS=20;

# assuming...
workdir="$HOME/processing_samples/snps_identification";

cd $workdir;
mkdir quality_control && cd quality_control;


echo
echo
echo
echo JOB_STARTED: `date`
echo
echo
echo
echo "// [MAIN] QUALITY CONTROL ANALYSIS: FASTQC & SAMTOOLS //"
echo "####################################################################################################################################################################################"
echo
echo
echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
echo "RUNING 'FASTQC' in PARALLEL for ALL SAMPLES..."; 
echo
ls $RAW_READS/*.fastq.gz | parallel --jobs $N_THREADS 'fastqc {} -o $workdir/quality_control';
wait;
sleep 2;


echo
echo
echo
echo "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
echo "RUNING 'SAMTOOLS' in PARALLEL for ALL SAMPLES..."; 
echo
ls $BAM_FILES/*.raln_indels.fm.bam | parallel --jobs  'samtools flagstats {} > $workdir/quality_control/{/.}.samtools_stats.txt';
wait;
sleep 2;

# assuming your current workdir is $workdir/quality_control
cd  $workdir/quality_control/;
multiqc .

# Done!

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
echo "####################################################################################################################################################################################"
echo
echo
echo "Start Time: $(date -d @$START_TIME)";
echo "End Time: $(date -d @$END_TIME)";
echo "Total Time: ${TOTAL_TIME_SECONDS} seconds (total elapsed: $REMAINING_SECONDS)";
echo "Total Time: ${TOTAL_TIME_MINUTES} minutes";
echo "Total Time: ${TOTAL_TIME} hours";
echo
echo
