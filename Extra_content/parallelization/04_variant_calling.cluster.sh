#!/usr/bin/bash


# ------------------------------------------------ #
#      *** Alejandro Nabor Lozada-Chávez ***       #
# Questions: alejandro.chavez@unipv.it             #
#            nabor.lozada@gmail.com                #
# ------------------------------------------------ #

# Script: Variant calling procedure based on a list of BAM files using BCFTOOLS for variant calling and filtering 
#         present in a single directory. BCFtools will produce a VCF file filtered in a compressed format (GZIP).



#SBATCH --account=
#SBATCH --partition=
#SBATCH --cpus-per-task=1
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
bcftools=$workdir/programs/bcftools/bin/bcftools;

# infiles
reference_genome=$workdir/reference_genome/my_fasta_seq.fasta;
bams_dir=$workdir/alignments/*.raln_indels.fm.bam;
bams_path=$workdir/alignments/;
variants_DIR=$workdir/variant_calling;
tmpdir=$workdir/TMP_DIR/;





# ------------------------------------------------------------------------------------------------
#                                 /// MAIN STEPS ///
# ------------------------------------------------------------------------------------------------
echo
echo
echo
echo JOB_STARTED: `date`
echo
echo `pwd`
echo
echo
echo
echo "// [MAIN] VARIANT CALLING VCF AND HARD FILTERING  //"
echo "####################################################################################################################################################################################"
echo
echo
# list all BAM files of interest (e.g. from a specific location o region) in a text file.
ls $bams_dir > $variants_DIR/all_BAM_files_list.txt;


# call variants
bcftools mpileup -f ${reference_genome}.fasta  -b $variants_DIR/all_BAM_files_list.txt | bcftools call -mv -Oz -o $variants_DIR/mexico.tapachula.aln.dedup.raln_indels.fm.vcf.gz;
wait;
sleep 2;

# Suggestion (recommended): Remove SNPs with a close proximity to indes of 10 bp:
bcftools filter --threads 3 --SnpGap 10 -Oz -o $variants_DIR/mySNPs.closeprox_rmvd.vcf.gz  $variants_DIR/mexico.tapachula.aln.dedup.raln_indels.fm.vcf.gz
wait;


# Extract only SNP variants:
bcftools view -v snps  $variants_DIR/mySNPs.closeprox_rmvd.vcf.gz  -Oz -o $variants_DIR/mexico.tapachula.aln.dedup.raln_indels.fm.snps.vcf.gz
wait;
sleep 2;


# filter SNPs by retaining only SNPs with a minor allele frequency of 0.01 && mapping quality bigger than 20:
bcftools view -q 0.01:minor  $variants_DIR/mexico.tapachula.aln.dedup.raln_indels.fm.snps.vcf.gz | bcftools filter -i '%QUAL>20' -Oz -o $variants_DIR/mexico.tapachula.aln.dedup.raln_indels.fm.snps.fltrd.vcf.gz
wait;
sleep 2;


# Extract only biallelic SNPs
bcftools -m2 -M2  -Oz -o $variants_DIR/mexico.tapachula.aln.dedup.raln_indels.fm.snps.fltrd.biallelic.vcf.gz  $variants_DIR/mexico.tapachula.aln.dedup.raln_indels.fm.snps.fltrd.vcf.gz
wait;
sleep 2;

# Count number of SNPs: first, index the VCF file (option: index -t), then count then (option: index -n). This will create a file with the same name of the VCF file, but with a file extension *vcf.gz.tbi
bcftools index -t $variants_DIR/mexico.tapachula.aln.dedup.raln_indels.fm.snps.fltrd.biallelic.vcf.gz
wait;
sleep 2;

# count SNPs:
bcftools index -n $variants_DIR/mexico.tapachula.aln.dedup.raln_indels.fm.snps.fltrd.biallelic.vcf.gz

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
