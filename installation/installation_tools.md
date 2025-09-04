### Installation of programs and packages

First, we define the structure of our working directory, then we proceed with the installation of all programs and packages.

1.CREATE `WORKING DIRECTORY` DIRECTORY STRUCTURE
```bash
# current location: $HOME // ~
cd $HOME
mkdir book_chapter_SNPs
cd book_chapter_SNPs;
mkdir programs snps_identification && cd snps_identification && mkdir reference_genome raw_reads mapped_unmapped alignments variant_calling;
cd $HOME;

# * WORKDIR STRUCTURE * 
tree book_chapter_SNPs/
# book_chapter_SNPs/
# ├── programs
# └── snps_identification
#     ├── alignments
#     ├── mapped_unmapped
#     ├── raw_reads
#     ├── reference_genome
#     └── variant_calling

# working directory
workdir="$HOME/book_chapter_SNPs";
```


2. INSTALLATION of PROGRAMS with MINICONDA

```bash
# move to programs dir:
cd $workdir/programs/

# Download and install miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh;

bash Miniconda3-latest-Linux-x86_64.sh;

# Dinished installation, then, create a conda environment with a name “variant_calling”.
# First, type each command one-by-one.
conda create -n variant_calling;
conda deactivate;
conda activate variant_calling;
conda config --add channels defaults;
conda config --add channels bioconda;
conda config --add channels conda-forge;
conda config --set channel_priority strict;

# Next, Install main programs via conda:
conda install -c bioconda gzip;
conda install -c bioconda bwa;
conda install -c bioconda fastqc;
conda install -c bioconda samtools;
conda install -c bioconda bcftools;
conda install -c bioconda vcftools;
conda install bioconda::bbmap;
conda install bioconda::parallel;
conda install bioconda::fastqc;
conda install bioconda::multiqc;

```

3. INSTALLATION of PROGRAMS via `wget` from SOURCE WEBSITE DEVELOPERS.

   3.1. Installation of **JAVA**.

```bash
# INSTALL JAVA version 22 for PICARD:
# correction in the link from the book chapter.
# Incorrect: https://download.oracle.com/java/22/latest/jdk-22_linux-x64_bin.tar.gz;
# Correct (`archive` instead of `latest`): 
wget https://download.oracle.com/java/22/archive/jdk-22_linux-x64_bin.tar.gz

tar -xzvf jdk-22_linux-x64_bin.tar.gz;
mv jdk-22 java22;

# location of executable files
cd java22/bin/;

# define location of picard:
java22="$workdir/programs/java22/bin/java";


# INSTALL version 8 for GATK:
wget https://javadl.oracle.com/webapps/download/AutoDL?BundleId=249840_43d62d619be4e416215729597d70b8ac
mv AutoDL?BundleId=249840_43d62d619be4e416215729597d70b8ac  java8.tar.gz
tar -zxvf java8.tar.gz;
mv jre1.8.0_411 java8

# test java:
cd java8/bin/
./java -jar -version
# java version "1.8.0_411"
# Java(TM) SE Runtime Environment (build 1.8.0_411-b09)
# Java HotSpot(TM) 64-Bit Server VM (build 25.411-b09, mixed mode)

# location of executable files
cd java8/bin/;

# define location of Java version 8:
java8="$workdir/book_chapter_SNPs/programs/java8/bin/java";

```

&emsp; &nbsp; 3.2. Installation of **PICARD**.

```bash
cd $workdir/programs/;

wget https://github.com/broadinstitute/picard/releases/download/2.27.0/picard.jar

# testing picard
$java22 -jar picard.jar

# define location of picard:
picard="$workdir/programs/picard.jar";
```

&emsp; &nbsp; 3.3. Installation of **GATK**.

```bash
cd $workdir/programs/;

# TDownload GATK v.3.8
wget https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2;

# decompress
tar -xvjf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2;
mv GenomeAnalysisTK-3.8-1-0-gf15c1c3ef gatk_v381;
cd gatk_v381/;

# testing GATK:
 $java8  -jar GenomeAnalysisTK.jar --help

# define location of GATK:
gatk_v381="$workdir/programs/gatk_v381/GenomeAnalysisTK.jar";
```

&emsp; &nbsp; 3.4. Installation of **TRIMMOMATIC**.

```bash
cd $workdir/programs

wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip;
unzip Trimmomatic-0.39.zip;
cd Trimmomatic-0.39;

# test trimmomatic
$java22 -jar trimmomatic-0.39.jar
```

&emsp; &nbsp; 3.5. Installation of **SRA-TOOLS**.

```bash
# create a separated conda environment from the named "variant_calling"
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda create -n sratools sra-tools
conda deactivate
conda activate sratools
conda install bioconda::sra-tools
```

Next, download infile datasets: References genome, reads from A. aegypti populations.
