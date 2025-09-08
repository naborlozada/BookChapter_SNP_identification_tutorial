#!/usr/bin/bash


# ------------------------------------------------ #
#      *** Alejandro Nabor Lozada-Chávez ***       #
# Questions: alejandro.chavez@unipv.it             #
#            nabor.lozada@gmail.com                #
# ------------------------------------------------ #

# Script: Variant calling procedure based on a list of BAM files using BCFTOOLS for variant calling and filtering present in a single directory.
#         BCFtools will produce a VCF file filtered in a compressed format (GZIP).




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
