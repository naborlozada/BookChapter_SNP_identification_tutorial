
# download & make fastQ files
# -------------------------------------------------------------------------------------------

# Download process can be performed by downloading 1-by-1 SRA file (typing each time the same command or in a loop reading a list), or in parallel. Although there are some apps that claim to do this process in parallel, 
# we recommend to do it  1-by-1.

prefetch SRR11006725
prefetch SRR11006726
prefetch SRR11196650
prefetch SRR11196649



# from SRA to paired-ended reads in Fastq format
fasterq-dump –progress --split 2 --verbose --threads 10 SRR11006725.sra
fasterq-dump –progress --split 2 --verbose --threads 10 SRR11006726.sra
fasterq-dump –progress --split 2 --verbose --threads 10 SRR11196650.sra
fasterq-dump –progress --split 2 --verbose --threads 10 SRR11196649.sra

# compress all files using gzip program. It can be done in GNU parallel.

ls *.fastq | parallel 'gzip {}'


# Quality control using FASTQC:

# create a directory for quality of fastq files and move inside, then, the call all FASTQ files with full path and do a LOOP of jobs to perform a quality test. 

mkdir quality_control;
cd quality_control;

FASTQ_QUALITY=$(pwd);

for $FILE in *.fastq.gz; do fastqc $FILE -o  $FASTQ_QUALITY/ ; done

# where, FILE represents a single FASTQ file, and FASTQ_QUALITY the current working directory.





