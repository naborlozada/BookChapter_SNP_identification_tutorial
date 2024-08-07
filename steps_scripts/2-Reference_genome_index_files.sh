
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

wget https://vectorbase.org/common/downloads/release-52/AaegyptiLVP_AGWG/fasta/data/VectorBase-52_AaegyptiLVP_AGWG_Genome.fasta

# download the genome annotation of the reference genome
wget https://vectorbase.org/common/downloads/release-52/AaegyptiLVP_AGWG/gff/data/VectorBase-52_AaegyptiLVP_AGWG.gff



# First, process reference genome (fasta file):
cat $reference | awk '{print $1}' > reference_parsed.fasta


# Proceed to make each index separately.
# index bwa
bwa index $reference_DIR/reference_parsed.fasta;
# index samtools
samtools faidx $reference_DIR/reference_parsed.fasta;
# index picard
$java22 -jar $picard CreateSequenceDictionary.jar  REFERENCE=$reference_DIR/reference_parsed.fasta  OUTPUT=$reference_DIR/reference_parsed.dict;


