#!/usr/bin/bash


# This script will call variants (SNPs and INDELs) for a set of samples in a population. In our 4 examples, 2 pairs belong to two different populations, Thies (Senegal)  and Tapachula (Mexico).
# Although these are two examples of two populations, the number of samples per population should be expected to be higher, so we can have a better quality call (identification) of each variant. However, this is an ideal situation and we know that sometimes this is no feasible.

# First, we have to create two TXT files with the names of the BAM files of each population separately. Next, use them to make the calls.

vi mexico.tapachula.txt  # (SRR11196649, SRR11196650)
vi senegal.thies.txt     # (SRR11006725, SRR11006726)

# Second, set variables for infiles:
# reference_DIR : location of the reference genome (fasta) and associated index files.
# variants_DIR  : location where the file with the variants (variant calling file, VCF) will be saved.

reference_genome=/home/username/book_snps_identification/reference_genome/Aedes-aegypti-LVP_AGWG_CHROMOSOMES.AaegL5_2.parsed;
variants_DIR=/home/username/variant_calling;


# ------------------------------------------------------------------------------------------------
#                                 /// MAIN STEPS ///
# ------------------------------------------------------------------------------------------------

# call variants
bcftools mpileup -f ${reference_genome}.fasta  -b $variants_DIR/mexico.tapachula.txt | bcftools call -mv -Oz -o $variants_DIR/mexico.tapachula.aln.dedup.raln_indels.fm.vcf.gz;
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
# This procedure has to be repeated for the population of Thies (Senegal) by copy-n-paste the code and replacing the population name across the whole code or you can implement a loop bash script to run it automatically one after another, so no need to copy-n-paste all the code. More in the 'Parallelization' section.
# There are more filtering options to clean the SNPs dataset. For example:



