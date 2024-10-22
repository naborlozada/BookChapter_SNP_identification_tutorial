#!/usr/bin/bash

# // Process all samples in a single Directory //

# Bash script to process several alignments in a LOOP. That is, a single run per sample. This is particularly useful when there is only 1 computer and/or 1 server several threads that can be used.
# It is important to mention that the first two steps are the main processes that run in parallel, while the rest does run with one thread.
# IMPORTANT: Do the proper changes in the directory full path names for programs/scripts ('programs' section) used to run this bash script, and also for the infiles, outfiles (if required), and output directory names. 
#
# Run this script as follow:
#       nohup bash process_all_samples.loop_job.sh 2>> process_all_samples.loop_job.stderr.log &> process_all_samples.loop_job.stderr.log &




# // paths & programs //
# --------------------------------------------------
# programs
trimmomatic=/home/username/programs/bioinformatics/Trimmomatic-0.38/trimmomatic-0.38.jar
SAMTOOLS=/home/username/programs/bioinformatics/samtools/samtools
PICARD=/home/username/programs/bioinformatics/picard.jar
BWA=/home/username/programs/bioinformatics/bwa/bwa

# infiles
reference_genome=Aedes-aegypti-LVP_AGWG_CHROMOSOMES.AaegL5_2.fasta;
reference_dir=/home/username/book_variant_calling/steps/2_reference_genome;
reads_dir=/home/username/raw_sequencing_data/Aedes_aegypti/population_name/*1.fastq.gz;
reads_path=/home/username/raw_sequencing_data/Aedes_aegypti/population_name;
tmpdir=/home/username/TMP_DIR/;

# outfiles
paired_dir_path=/home/username/Aedes-aegypti/population_name/processing/outfiles/paired;
unpaired_dir_path=/home/username/Aedes-aegypti/population_name/processing/outfiles/unpaired;
outfiles=/home/username/Aedes-aegypti/population_name/processing/outfiles/main;




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
        R2=${R1//_1.fastq.gz/_2.fastq.gz}
        trimclip=${R1##*/}
        trimclip_logfile=`echo "$trimclip" | sed "s/\_1.fastq.gz/.R1_R2.trimmomatic.logfile.txt/"`
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
        echo ">>>conmmand_line.trimmomatic: java -jar $trimmomatic PE -threads 30 -phred33 -trimlog $outfiles/$trimclip_logfile $R1 $R2  $paired_dir_path/$R1_paired  $unpaired_dir_path/$R1_unpaired $paired_dir_path/$R2_paired  $unpaired_dir_path/$R2_unpaired  ILLUMINACLIP:/home/username/programs/bioinformatics/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36";
        
        # main command execution:
        time java -jar $trimmomatic PE -threads 30 -phred33 \
            -trimlog $outfiles/$trimclip_logfile \
            $R1 \
            $R2 \
            $paired_dir_path/$R1_paired \
            $unpaired_dir_path/$R1_unpaired \
            $paired_dir_path/$R2_paired \
            $unpaired_dir_path/$R2_unpaired \
            ILLUMINACLIP:/home/username/programs/bioinformatics/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 \
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
        echo ">>>conmmand_line.BWA: $BWA mem -t 30  $reference_dir/$reference_genome  $paired_dir_path/$R1_paired  $paired_dir_path/$R2_paired  -R '@RG\tID:$ID_TAG\tSM:$ID_TAG\tPL:illumina' | $SAMTOOLS view -b -@ 30 - >  $outfiles/$BAM_R1R2_paired_outfile";

        # Paired-ended samples (trimmed quallity and clippled)
        time $BWA mem -t 30 \
            $reference_dir/$reference_genome \
            $paired_dir_path/$R1_paired \
            $paired_dir_path/$R2_paired \
            -R '@RG\tID:'$ID_TAG'\tSM:'$ID_TAG'\tPL:illumina' | $SAMTOOLS view -b -@ 30 - > $outfiles/$BAM_R1R2_paired_outfile
        wait;
        sleep 2;


        echo
        echo
        echo "[3] // Sort bam reads file: SortSam //";
        echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
        BAM_SORTED_outfile=`echo "$BAM_R1R2_paired_outfile" | sed "s/\.paired.aligned_raw_reads.bam/.paired.aligned_reads_sorted.bam/"`;
        echo ">>>conmmand_line.SortSam: java -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=20 -Xmx60G -jar $PICARD SortSam  INPUT=$outfiles/$BAM_R1R2_paired_outfile  OUTPUT=$outfiles/$BAM_SORTED_outfile  SORT_ORDER=coordinate  TMP_DIR=$tmpdir";

        # [SortSam] sort bam file reads (by coordinates, default & mandatory required from GATK)
        time java -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=20 -Xmx60G -jar $PICARD SortSam \
            INPUT=$outfiles/$BAM_R1R2_paired_outfile \
            OUTPUT=$outfiles/$BAM_SORTED_outfile \
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
        echo ">>>conmmand_line.MarkDuplicates: java -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=32 -Xmx60G -jar $PICARD MarkDuplicates  INPUT=$outfiles/$BAM_SORTED_outfile  OUTPUT=$outfiles/$DUPLICATES_outfile  METRICS_FILE=$outfiles/$DUPLICATES_metrics  TMP_DIR=$tmpdir";

        # MarkDuplicates
        time java -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=32 -Xmx60G -jar $PICARD MarkDuplicates \
            INPUT=$outfiles/$BAM_SORTED_outfile \
            OUTPUT=$outfiles/$DUPLICATES_outfile \
            METRICS_FILE=$outfiles/$DUPLICATES_metrics \
            TMP_DIR=$tmpdir
        wait;
        sleep 2;


        echo
        echo
        echo " [5] // Make an index files for the sorted bam files //";
        echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
        echo ">>>conmmand_line.BWA: $BWA index  $outfiles/$DUPLICATES_outfile";
        $BWA index  $outfiles/$DUPLICATES_outfile;
        echo
        echo
        echo ">>>conmmand_line.SAMTOOLS: $SAMTOOLS index  $outfiles/$DUPLICATES_outfile";
        $SAMTOOLS index $outfiles/$DUPLICATES_outfile;
        wait;
        sleep 2;


        echo
        echo
        echo "[6] // ValidateSamFile:_check_errors_in_SAM_file //"
        echo "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------";
        BAM_DEDUP_ValidateSamFile=`echo "$DUPLICATES_outfile" | sed "s/\.paired.dedup_reads.bam/.paired.dedup_reads.ValidateSamFile.txt/"`
        
        echo ">>>conmmand_line.ValidateSamFile: java  -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=20 -Xmx60G -jar $PICARD ValidateSamFile I=$outfiles/$DUPLICATES_outfile  O=$outfiles/$BAM_DEDUP_ValidateSamFile   MODE=SUMMARY   TMP_DIR=$tmpdir";

        # validate the format and other feature that can cause errors during the next analysis.
        time java -Djava.io.tmpdir=$tmpdir -XX:ParallelGCThreads=20 -Xmx50G -jar $PICARD ValidateSamFile \
            I=$outfiles/$DUPLICATES_outfile \
            O=$outfiles/$BAM_DEDUP_ValidateSamFile \
            MODE=SUMMARY \
            TMP_DIR=$tmpdir
        wait;
        sleep 2;

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
```

