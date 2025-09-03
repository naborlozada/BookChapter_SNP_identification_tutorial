
## Step 2: Download & Parse Reference Genome Files

### UPDATE for VectorBase database

This file has the same information that in the other short command lines instructions from our book [chapter](https://link.springer.com/protocol/10.1007/978-1-0716-4583-3_2), however, the only change is the updated information associated to the VectorBase (VB) database: to have access to any information/resource of this database you have to create a **free** account through a registration process. This is a major update that this database has been introduced in the past month after major changes.
They share the following message:

*VEuPathDB is evolving under a new organizational structure. In order to use VEuPathDB resources, you will now need to log into your free account. This helps us collect accurate user metrics to guide future development.*

Also, they clarify:

*IMPORTANT: If you already registered in another site
(AmoebaDB, CryptoDB, FungiDB, GiardiaDB, MicrosporidiaDB, PiroplasmaDB, PlasmoDB, SchistoDB, ToxoDB, TrichDB, TriTrypDB, VectorBase or VEuPathDB)
you do NOT need to register again.*

After registration you will be able to see all resources and can decide what version of the reference genome you want to download (or do any other type of analyses you wish to do). The single individual genomes as reads sequences used as examples for our book chapter can be downloaded with no registration in VB (so far), and any other sequences were created from specific regions.

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
cd $wordir/snps_identification/reference_genome

# download reference genome: FASTA
wget https://vectorbase.org/common/downloads/Legacy%20VectorBase%20Files/Aedes-aegypti/Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.fa.gz
# Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.fa.gz      100%[=====================================================================================================================>] 389,25M  25,5MB/s    in 32s     

# download reference genome: GFF
wget https://vectorbase.org/common/downloads/Legacy%20VectorBase%20Files/Aedes-aegypti/Aedes-aegypti-LVP_AGWG_BASEFEATURES_AaegL5.2.gff3.gz
# Aedes-aegypti-LVP_AGWG_BASEFEATURES_AaegL5.2.gff3.gz 100%[=====================================================================================================================>]   5,12M  2,98MB/s    in 1,7s    

# Parse fasta headers
zcat Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.fa.gz | awk '{print $1}' > Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.parsed.fasta

# // Verify chromosome/contigs IDS //
# Get lists of chromosomes and contigs IDs from FASTA and GFF files. Next, compare both IDs lists, and find whether all IDs are present in both files ('matched_ID') or there are IDs that are uniquely present either in the fasta or gff files annotation ('single_ID').

# Get chromosome/contigs IDs from fasta file
grep '>' Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.parsed.fasta | sed -E 's/^>//' > aaegL5.v5_2.chrms_contigs_ids.fasta_list.txt
# Get chromosome/contigs IDs from gff file
zcat Aedes-aegypti-LVP_AGWG_BASEFEATURES_AaegL5.2.gff3.gz | grep -v '^\#' | awk '{print $1}' | sort -u > aaegL5.v5_2.chrms_contigs_ids.gff_list.txt

# Compare both lists
cat aaegL5.v5_2.chrms_contigs_ids.fasta_list.txt  aaegL5.v5_2.chrms_contigs_ids.gff_list.txt | sort | uniq -c | sort -n |  awk '{if($1==1){print "NOT_ALL_IDs_matched\t"$0} else {print "ALL_IDs_matched\t"$0} }' | cut -f 1 | sort | uniq -c | sort -n
# 2310 ALL_IDs_matched, no single chrm ids.

# Reference genome path
# Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.parsed.fasta
reference_genome=/home/username/book_variant_calling/steps/2_reference_genome/Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.parsed


# ** Proceed to make each index separately **
# index bwa
bwa index ${reference_genome}.fasta;
# index samtools
samtools faidx ${reference_genome}.fasta;
# index picard
$java22 -jar $picard CreateSequenceDictionary.jar  REFERENCE=${reference_genome}.fasta  OUTPUT=${reference_genome}.dict;


