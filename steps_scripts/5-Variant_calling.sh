#!/usr/bin/bash



# call variants
bcftools mpileup -f $reference_DIR/reference_parsed.fasta $alignments_DIR/sample1.aln.dedup.raln_indels.fm.bam | bcftools call -mv -Oz -o $variants_DIR/sample1.aln.dedup.raln_indels.fm.vcf.gz;
wait;
sleep 2;


# extract only SNP variants:
bcftools view -v snps  $variants_DIR/sample1.aln.dedup.raln_indels.fm.vcf.gz  -Oz -o $variants_DIR/sample1.aln.dedup.raln_indels.fm.snps.vcf.gz
wait;
sleep 2;


# filter SNPs
bcftools view -q 0.01:minor  $variants_DIR/sample1.aln.dedup.raln_indels.fm.snps.vcf.gz | bcftools filter -i '%QUAL>20' -Oz -o $variants_DIR/sample1.aln.dedup.raln_indels.fm.snps.fltrd.vcf.gz
wait;
sleep 2;


# get only biallelic SNPs
bcftools -m2 -M2  -Oz -o $variants_DIR/sample1.aln.dedup.raln_indels.snps.fm.fltrd.biallelic.vcf.gz  $variants_DIR/sample1.aln.dedup.raln_indels.snps.fltrd.vcf.gz
wait;
sleep 2;

# count final number of SNPs
bcftools index -t $variants_DIR/sample1.aln.dedup.raln_indels.snps.fltrd.vcf.gz
wait;
sleep 2;

# count SNPs:
bcftools index -n $variants_DIR/sample1.aln.dedup.raln_indels.snps.fltrd.vcf.gz
