
## UPDATE for VectorBase database

This file has the same information that in the other short command lines from the book [chapter](https://github.com/naborlozada/BookChapter_SNP_identification_tutorial/blob/main/steps_scripts/2-Reference_genome_index_files.sh), 
the only change is the updated information associated to the VectorBase database: to have access to any information/resource of this database you have to create a **free** account through a registration process. 
They share following message:

*VEuPathDB is evolving under a new organizational structure. In order to use VEuPathDB resources, you will now need to log into your free account. This helps us collect accurate user metrics to guide future development.*

Also, thy clarify:

*IMPORTANT: If you already registered in another site
(AmoebaDB, CryptoDB, FungiDB, GiardiaDB, MicrosporidiaDB, PiroplasmaDB, PlasmoDB, SchistoDB, ToxoDB, TrichDB, TriTrypDB, VectorBase or VEuPathDB)
you do NOT need to register again.*

After registration you will be able to see all resources and decide what version of the reference genome you want to download. The genome version of this example can be downloaded with no registration (so far). 

```bash 
# ------------------------------------------------------------------------------------------
# Make index files for the reference genome:
# ------------------------------------------------------------------------------------------

# variable:
# $reference     : full path location of the reference genome in fasta format. 
# $reference_DIR : location of the reference genome (fasta) and associated index files.
# $java22        : location of the binary java file (version 22).
# $picard        : location of Picard program.  


# Working directory:
cd $HOME/snps_identification/reference_genome

# get reference genome
wget https://vectorbase.org/common/downloads/Legacy%20VectorBase%20Files/Aedes-aegypti/Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.fa.gz
# Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.fa.gz      100%[=====================================================================================================================>] 389,25M  25,5MB/s    in 32s     


wget https://vectorbase.org/common/downloads/Legacy%20VectorBase%20Files/Aedes-aegypti/Aedes-aegypti-LVP_AGWG_BASEFEATURES_AaegL5.2.gff3.gz
# Aedes-aegypti-LVP_AGWG_BASEFEATURES_AaegL5.2.gff3.gz 100%[=====================================================================================================================>]   5,12M  2,98MB/s    in 1,7s    

# parse fasta headers
zcat Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.fa.gz | awk '{print $1}' > Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.parsed.fasta

# Verify CHRM IDS: get list of chromosomes and contigs IDs from FASTA and GFF files to compare IDs in both files: 
grep '>' Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.parsed.fasta | sed -E 's/^>//' > aaegL5.v5_2.chrms_contigs_ids.fasta_list.txt
zcat Aedes-aegypti-LVP_AGWG_BASEFEATURES_AaegL5.2.gff3.gz | grep -v '^\#' | awk '{print $1}' | sort -u > aaegL5.v5_2.chrms_contigs_ids.gff_list.txt
cat aaegL5.v5_2.chrms_contigs_ids.fasta_list.txt  aaegL5.v5_2.chrms_contigs_ids.gff_list.txt | sort | uniq -c | sort -n | awk '{if($1==1){print "single_ID\t"$0} else {print "matched_ID\t"$0} }' | cut -f 1| sort | uniq -c | sort -n
# 2310 matched_ID, no single chrm ids.

# directory
reference_genome=/home/username/book_variant_calling/steps/2_reference_genome/Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.parsed.fasta

# Proceed to make each index separately.
# index bwa
bwa index ${reference_genome.fasta;
# index samtools
samtools faidx ${reference_genome}.fasta;
# index picard
$java22 -jar $picard CreateSequenceDictionary.jar  REFERENCE=${reference_genome}.fasta  OUTPUT=${reference_genome}.dict;

