#!/usr/bin/bash


# ------------------------------------------------ #
#      *** Alejandro Nabor Lozada-Chávez ***       #
# Questions: alejandro.chavez@unipv.it             #
#            nabor.lozada@gmail.com                #
# ------------------------------------------------ #
#
# Script: performs a trimming and alignment of raw reads on a reference genome, while at the end sort and dedup all aligned reads present in a single directory.




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
reads_dir=$workdir/raw_reads/*1.fastq.gz;
reads_path=$workdir/raw_reads;
tmpdir=$workdir/TMP_DIR/;

# outfiles
paired_dir_path=$workdir/paired_mapped;
unpaired_dir_path=$workdir/paired_unmapped;
alignments=$workdir/alignments;




echo
echo
echo
echo JOB_STARTED: `date`
echo
echo `pwd`
echo
echo
echo
echo "// [MAIN] PRE-PROCESSING SAMPLES //"
echo "####################################################################################################################################################################################"

for R1 in $reads_dir
do
        echo
        echo
        echo "[1] // Trimmomatic_process //"
        echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
        R2=${R1//_1.fastq.gz/_2.fastq.gz};
        trimclip=${R1##*/};
        trimclip_logfile=`echo "$trimclip" | sed "s/\_1.fastq.gz/.R1_R2.trimmomatic.logfile.txt/"`;
        #echo "pair-ended files: $R1 ---- $R2"

        R1paired=`echo "$R1" | sed "s/\.fastq.gz/_paired.trimclip.fastq/"`;
        R1unpaired=`echo "$R1" | sed "s/\.fastq.gz/_unpaired.trimclip.fastq/"`;
        R2paired=`echo "$R2" | sed "s/\.fastq.gz/_paired.trimclip.fastq/"`;
        R2unpaired=`echo "$R2" | sed "s/\.fastq.gz/_unpaired.trimclip.fastq/"`;
        
        R1_paired=${R1paired##*/};
        R2_paired=${R2paired##*/};
        R1_unpaired=${R1unpaired##*/};
        R2_unpaired=${R2unpaired##*/};
        
        # show the full command line used in this step, so we can check what options were used in this process when we look at the standard output file. 
        echo ">>>conmmand_line.trimmomatic: java -jar $trimmomatic PE -threads 30 -phred33 -trimlog $paired_dir_path/$trimclip_logfile $R1 $R2  $paired_dir_path/$R1_paired  $unpaired_dir_path/$R1_unpaired $paired_dir_path/$R2_paired  $unpaired_dir_path/$R2_unpaired  ILLUMINACLIP:$trimmomatic_adapters/TruSeq3-PE.fa:2:30:10  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36";
        
        # main command execution:
        time java -jar $trimmomatic PE -threads 30 -phred33 \
            -trimlog $paired_dir_path/$trimclip_logfile \
            $R1 \
            $R2 \
            $paired_dir_path/$R1_paired \
            $unpaired_dir_path/$R1_unpaired \
            $paired_dir_path/$R2_paired \
            $unpaired_dir_path/$R2_unpaired \
            ILLUMINACLIP:$trimmomatic_adapters/TruSeq3-PE.fa:2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        wait;
        sleep 2;


        echo
        echo
        echo "[2] // Mapping_reads //"
        echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
        NEWTAG=${R1##*/};
        ID_TAG=${NEWTAG//_1.fastq.gz/};
        BAM_R1R2_paired_outfile=`echo "$R1_paired" | sed "s/\_1_paired.trimclip.fastq/.paired.aligned_raw_reads.bam/"`;
        echo ">>>conmmand_line.BWA: $BWA mem -t 30  $reference_genome  $paired_dir_path/$R1_paired  $paired_dir_path/$R2_paired  -R '@RG\tID:$ID_TAG\tSM:$ID_TAG\tPL:illumina' | $SAMTOOLS view -b -@ 30 - >  $alignments/$BAM_R1R2_paired_outfile";

        # Paired-ended samples (trimmed quallity and clippled)
        time $BWA mem -t 30 \
            $reference_genome \
            $paired_dir_path/$R1_paired \
            $paired_dir_path/$R2_paired \
            -R '@RG\tID:'$ID_TAG'\tSM:'$ID_TAG'\tPL:illumina' | $SAMTOOLS view -b -@ 30 - > $alignments/$BAM_R1R2_paired_outfile
        wait;
        sleep 2;


        echo
        echo
        echo "[3] // Sort bam reads file: SortSam //";
        echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
        BAM_SORTED_outfile=`echo "$BAM_R1R2_paired_outfile" | sed "s/\.paired.aligned_raw_reads.bam/.paired.aligned_reads_sorted.bam/"`;
        echo ">>>conmmand_line.SortSam: java -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=20 -Xmx60G -jar $PICARD SortSam  INPUT=$alignments/$BAM_R1R2_paired_outfile  OUTPUT=$alignments/$BAM_SORTED_outfile  SORT_ORDER=coordinate  TMP_DIR=$tmpdir";

        # [SortSam] sort bam file reads (by coordinates, default & mandatory required from GATK)
        time java -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=20 -Xmx60G -jar $PICARD SortSam \
            INPUT=$alignments/$BAM_R1R2_paired_outfile \
            OUTPUT=$alignments/$BAM_SORTED_outfile \
            SORT_ORDER=coordinate \
            TMP_DIR=$tmpdir
        wait;
        sleep 2;


        echo
        echo        
        echo "[4] // MarkDuplicates:_Mark_duplicates //";
        echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
        DUPLICATES_outfile=`echo "$BAM_SORTED_outfile" | sed "s/\.paired.aligned_reads_sorted.bam/.paired.dedup_reads.bam/"`;
        DUPLICATES_metrics=`echo "$BAM_SORTED_outfile" | sed "s/\.paired.aligned_reads_sorted.bam/.paired.dedup_reads.metrics.txt/"`;
        echo ">>>conmmand_line.MarkDuplicates: java -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=32 -Xmx60G -jar $PICARD MarkDuplicates  INPUT=$alignments/$BAM_SORTED_outfile  OUTPUT=$alignments/$DUPLICATES_outfile  METRICS_FILE=$alignments/$DUPLICATES_metrics  TMP_DIR=$tmpdir";

        # MarkDuplicates
        time java -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=32 -Xmx60G -jar $PICARD MarkDuplicates \
            INPUT=$alignments/$BAM_SORTED_outfile \
            OUTPUT=$alignments/$DUPLICATES_outfile \
            METRICS_FILE=$alignments/$DUPLICATES_metrics \
            TMP_DIR=$tmpdir
        wait;
        sleep 2;


        echo
        echo
        echo " [5] // Make an index files for the sorted bam files //";
        echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
        echo ">>>conmmand_line.BWA: $BWA index  $alignments/$DUPLICATES_outfile";
        $BWA index  $alignments/$DUPLICATES_outfile;
        echo
        echo
        echo ">>>conmmand_line.SAMTOOLS: $SAMTOOLS index  $alignments/$DUPLICATES_outfile";
        $SAMTOOLS index $alignments/$DUPLICATES_outfile;
        wait;
        sleep 2;


        echo
        echo
        echo "[6] // ValidateSamFile:_check_errors_in_SAM_file //"
        echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
        BAM_DEDUP_ValidateSamFile=`echo "$DUPLICATES_outfile" | sed "s/\.paired.dedup_reads.bam/.paired.dedup_reads.ValidateSamFile.txt/"`
        
        echo ">>>conmmand_line.ValidateSamFile: java  -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=20 -Xmx60G -jar $PICARD ValidateSamFile I=$alignments/$DUPLICATES_outfile  O=$alignments/$BAM_DEDUP_ValidateSamFile   MODE=SUMMARY   TMP_DIR=$tmpdir";

        # validate the format and other feature that can cause errors during the next analysis.
        time java -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=20 -Xmx50G -jar $PICARD ValidateSamFile \
            I=$alignments/$DUPLICATES_outfile \
            O=$alignments/$BAM_DEDUP_ValidateSamFile \
            MODE=SUMMARY \
            TMP_DIR=$tmpdir
        wait;
        sleep 2;


    echo
        echo
        echo "[7] // Compress logfile from trimmomatic //"
        echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
        #trimclip_logfile_GZIP=`echo "$trimclip_logfile" | sed "s/.R1_R2.trimmomatic.logfile.txt/.R1_R2.trimmomatic.logfile.txt.gz/"`;

        echo ">>>conmmand_line.compress_logfile: gzip $trimclip_logfile ";

        time gzip $trimclip_logfile;
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
echo ALL_JOBS_COMPLETED: `date`
echo
echo
echo


