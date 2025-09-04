## Step 5: variant calling

This script will call both type of variants (SNPs and INDELs) for a two samples that were collected in a specific population. In our 4 examples, 2 pairs belong to two different populations, Thies (Senegal)  and Tapachula (Mexico).

Having two samples Although these are two examples of two populations, the number of samples per population should be expected to be higher, so we can have a better calls quality  (identification) of each variant. See recommendations in our book chapter on this matter.

### Variant calling: files and commands


First, we have to create two TXT files with the names of the BAM files of each population separately. The files need to be created in the directory `variants_DIR`. So, let's move to the this directory.


```bash
cd ~  # our $HOME
workdir="$HOME/book_chapter_SNPs";  # in case you are using another terminal
variants_DIR=$workdir/variant_calling;
cd variants_DIR;
```

Next, create these file to make the calls.


```txt
# Create text files and add all SRA samples IDs. Here, only the ID is added without the hashtag symbol (it is used here as commment to show the ids). You can use any text editor you want, here we use "vi". A brief list of other suggestions: VS code, Codium, Sublime, Emacs, Nedit. 

username@server_name:~$ vi mexico.tapachula.txt
# SRR11196649
# SRR11196650

username@server_name:~$ vi senegal.thies.txt
# SRR11006725
# SRR11006726
```

Now, we perform the variant calling procedure.

```bash
# NOTE: this procedure has to be implemented inside the conda environment "variant_calling":
conda deactivate; 
conda acivate variant_calling;

# Inside the environment, set variables for infiles:
# reference_DIR : location of the reference genome (fasta) and associated index files.
# variants_DIR  : location where the file with the variants (variant calling file, VCF) will be saved.

# working directory
workdir="$HOME/book_chapter_SNPs";

reference_genome=$workdir/reference_genome/Aedes-aegypti-LVP_AGWG_CHROMOSOMES.AaegL5_2.parsed;
variants_DIR=$workdir/variant_calling;



#   /// MAIN STEPS ///

# [1] call variants
bcftools mpileup -f ${reference_genome}.fasta  -b $variants_DIR/mexico.tapachula.txt | bcftools call -mv -Oz -o $variants_DIR/mexico.tapachula.aln.dedup.raln_indels.fm.vcf.gz;
wait;
sleep 2;

# [2] Suggestion (recommended): Remove SNPs with a close proximity to indes of 10 bp:
bcftools filter --threads 3 --SnpGap 10 -Oz -o $variants_DIR/mexico.tapachula.closeprox_rmvd.vcf.gz  $variants_DIR/mexico.tapachula.aln.dedup.raln_indels.fm.vcf.gz
wait;


# [3] Extract only SNP variants:
bcftools view -v snps  $variants_DIR/mexico.tapachula.closeprox_rmvd.vcf.gz  -Oz -o $variants_DIR/mexico.tapachula.aln.dedup.raln_indels.fm.snps.vcf.gz
wait;
sleep 2;

# similar way: bcftools view --types snps 


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
# This procedure has to be repeated for the population of Thies (Senegal) by copy-n-paste the code and replacing the population name across the whole code or you can implement a loop bash script to run it automatically one after another, so no need to copy-n-paste all the code.
```

### MORE OPTIONS FOR HARD-FILTERING SNPS in VCF FORMAT

Filtering SNPs is not easy step due the diverse set of technical and data source instances, such as the variant caller program used that will defined the name of the options used for filtering (e.g. GATK, Freebayes, Platypus) or whether the *Reads data* to be analyzed comes from *low*, *medium* or *high* coverage data, respectively. Thus, define what set of criteria will be used for this procedure is key to further analyses and should be deeply discussed; sometimes maybe perform different cutoffs on the same VCF file can also performed, them compare them to identify the differences. A highly recomended discussion with worked filtering examples is provided by [Jon Puritz in the DDOCENT website](https://ddocent.com//filtering/), in which such protocols are applied specifically for uses RAD (Restriction-site Associated DNA), ddRAD, ezRAD, SE-RAD, and PE-RAD sequencing data. Also, a diverse set of useful options to filter SNPs with BCFTOOLS is posted below:


```text
a) DDOCENT:
    https://ddocent.com//filtering/ 

b) BCFtools cheat sheet from Ernesto Lowy:
    https://gist.github.com/elowy01/93922762e131d7abd3c7e8e166a74a0b

c) BCFTOOLS official website:
    https://samtools.github.io/bcftools/bcftools.html
    https://samtools.github.io/bcftools/howtos/filtering.html
    
d)  BCFTOOLS official website with some additional instructions to plot filtered data using pipelines with BCFTOOLS + PERL + GNUPLOT:
    https://www.htslib.org/workflow/filter.html


```


