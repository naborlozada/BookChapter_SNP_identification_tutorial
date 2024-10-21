
## Step 1: get reads from SRA database using SRA IDs.

```bash
# /// get 4 single samples whole genomes ///

# WORKING DIRECTORY (pwd): /home/username/book_variant_calling/steps/1_get_sequences

# load 'sratools' conda's environment:
[me@server current_directory]$ conda activate sratools
# small change in the prompt command line:
(sratools) [username@server current_directory]$

# Download process can be performed by downloading 1-by-1 SRA file (typing each time the same command or in a loop reading a list), or in parallel. 
# Although there are some apps that claim to do this process in parallel, it might be possible that such processes can get some errors while the FASTQ file is downloaded and/or the final file has some corrupted information.
# Therefore, we recommend to do it 1-by-1. 
(sratools) [username@server current_directory]$ prefetch SRR11006725
(sratools) [username@server current_directory]$ prefetch SRR11006726
(sratools) [username@server current_directory]$ prefetch SRR11196650
(sratools) [username@server current_directory]$ prefetch SRR11196649


# It might be possible to get the reads in parallel using one fo the two following options: 1) GNU parallel function or 2) the PBS jobs system implemented throught the SLURM cluster's system.
# I did not try these options for this particular task, but although it might be posssible (with any of these functions), I think it might be risky since the multiple calls network request from a single IP address point can be suspicious for the recipient SRA server. However, you are free to try it. 
# We recommend to do it 1-by-1 using a loop, to do so just create a script to make the job: 

(sratools) [username@server current_directory]$ vi get_SRA_sequences.through_loops.sh

while IFS="" read -r SRAID || [ -n "$SRAID" ]
do
  echo "prefetch $SRAID";    # this is to know on what sample is currently reading/working the `bash` script. "time" just to know how long does it takes that process to be completed.
  time prefetch "$SRAID";
done < SRA_IDs_to_download.txt 

# Where: `SRA_IDs_to_download.txt` hash all SRA IDs, which in this case are four.
# NOTE: if you run this script using a linux screen system (virtual inside session), you have to load again the sratools' conda's environment. If you run it with the `nohup` command, there is no need to do reload anything.
```

### ERROR messages

Unfortunately, sometimes sequences cannot be downloaded due different errors, if one of them comes from the `sratools` tool package, like in the example below, it will be better to install the `sratools` separately.

```bash
# Error message example:

(sratools) [username@server current_directory]$ prefetch SRR11006725
# 2024-10-18T22:02:20 prefetch.2.8.0 sys: connection failed while opening file within cryptographic module - mbedtls_ssl_handshake returned -9984 ( X509 - Certificate verification failed, e.g. CRL, CA or signature check failed )
# 2024-10-18T22:02:20 prefetch.2.8.0 sys: mbedtls_ssl_get_verify_result returned 0x8 (  !! The certificate is not correctly signed by the trusted CA  )
# 2024-10-18T22:02:20 prefetch.2.8.0 err: path not found while resolving tree within virtual file system module - 'SRR11006725' cannot be found.


# /// Installation of 'sratools' ///

# Best option is to use the binary files: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit
# Once installed, use the full path directory where the `sratoolkit` was saved, like in the example path below:
/home/username/programs/bioinformatics/sratoolkit.3.0.5/bin/prefetch

# or, save it in a variable (for short command lines declaration) and used in the command lines or integrate this in the bash script from above:
sratoolkit_prefetch=/home/username/programs/bioinformatics/sratoolkit.3.0.5/bin/prefetch

# run it:
(sratools) [username@server current_directory]$ time $sratoolkit_prefetch SRR11006726 2>> tee get_fastq.SRR11006726.stderr.log &> get_fastq.SRR11006726.stderr2.log
# real    8m6,274s
# user    1m2,301s
# sys     0m14,558s


# SRA to FASTQ reads format:

# save the function `fasterq-dump` as a variable and use it as above, but in parallel (careful with the number of threads):
sratoolkit_fasterq_dump=/home/chavez/programs/bioinformatics/sratoolkit.3.0.5/bin/fasterq-dump

(base) [username@server current_directory]$ time $sratoolkit_fasterq_dump  --threads 10 --verbose  SRR11006726 2>> get_fastq.sra2fastq.SRR11006726.stderr.log &> get_fastq.sra2fastq.SRR11006726.stderr.log
# real    7m4,818s
# user    13m55,084s
# sys     1m58,881s

# Done, reads were transform to FASTQ. Now, just compress this big file with the gzip command. This can be don in parallel once all 4 genomes have been downloaded:
ls SRR*.fastq | parallel 'gzip {}' &

# Done! Now go to the step 2: download and index reference genome 

```

