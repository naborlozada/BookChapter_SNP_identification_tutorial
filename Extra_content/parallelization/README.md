## Parallelization 


### Introduction

Here, we share XXX bash scripts that can be used in a server (i.e. work stations) and/or HPC computer system. Briefly, these sytems have a high I/O speed, big space data storage, 
several number of cores/threads, a considerable number of nodes, and large amount of RAM memory in each node.
In other words, a super computer system (server/cluster) is designed to perform, at high speeds using multiple processors (CPUs) working in parallel: 

  &emsp; &nbsp; (a) a large amount of computational jobs, \
  &emsp; &nbsp; (b) complex computations (physics, math, chemistry, biology, among others), \
  &emsp; &nbsp; (c) database storage, running web applications, and online jobs supported by the infrastructure of the cluster.
     
**Importantly, this part of the tutorial not only requires that all users have a basic knowledge on how to use a HPC computer system (i.e. threads, memory, storage, speed), but also to be 
fully aware that these systems are not a personal computer, they are valuable and unique shared computational resources with other users. Therefore, please be cautious on the use 
of the number of threads, RAM memory usage and storage space for each single task across a pipeline.**

Also, before using a Server or HPC, contact your IT administrator of the research institute, unversity or web-cloud to know more about these systems. 
It is very likely that they have an available manual for users in a pdf file or website (it must be!). Also, conversations with other previous/current persons familiar with these systems 
are always very helpful.

### Scripts

#### — Servers system —

Server system bash scripts. The bash scripts are a compilation of the steps shown in this tutorial, which are divided en 4 major steps shown below. These scripts work basically as a 
loop function, reading each pair of FASTQ files and/or BAM files, and/or VCF files to apply all corresponding bioinformatics tools. Please, change all variables for properly for their 
corresponding files, paths and programs. They most shows the exact location of these files, so they can be processed properly.

1. Quality Control: FASTQC and SAMtools.

```bash
# run script using a nohup function or a screen linux session.
# nohup
nohup time bash 01_quality_control.all_samples.sh 2>> 01_quality_control.all_samples.stderr.log &> 01_quality_control.all_samples.stderr.log &

# screen
# open a new screen session and run it:
time bash 01_quality_control.all_samples.sh 2>&1 | tee 01_quality_control.all_samples.stderr.log
# exit of the screen session properly to left the script running in the background.
```
   
2. Processing samples: trimming and alignment of Raw reads.

```bash
# run script using a nohup function or a screen linux session.
# nohup
nohup time bash 02_trimming_n_alignment.server.sh 2>> 02_trimming_n_alignment.server.stderr.log &> 02_trimming_n_alignment.server.stderr.log &

# screen
# open a new screen session and run it:
time bash 02_trimming_n_alignment.server.sh 2>&1 | tee 02_trimming_n_alignment.server.stderr.log
# exit of the screen session properly to left the script running in the background.
```
   
3. Re-alignment of indels in dedupped alignments.

```bash
# nohup
nohup time bash 03_re_alignment.server.sh 2>> 03_re_alignment.server.stderr.log &> 03_re_alignment.server.stderr.log &

# screen
# open a new screen session and run it:
time bash 03_re_alignment.server.sh 2>&1 | tee 03_re_alignment.server.stderr.log
# exit of the screen session properly to left the script running in the background.
```
   
4. Variant calling and hard-filtering with BCFtools.

```bash
# nohup
nohup time bash 04_variant_calling.server.sh 2>> 04_variant_calling.server.stderr.log &> 04_variant_calling.server.stderr.log &

# screen
# open a new screen session and run it:
time bash 04_variant_calling.server.sh 2>&1 | tee 04_variant_calling.server.stderr.log
# exit of the screen session properly to left the script running in the background.
```

#### — Cluster system —

Cluster system bash scripts. The bash scripts are a compilation of the steps shown in this tutorial divided en 4 major steps shown below (as for Servers). These scripts work basically using
the SLUM system, and have all required parameters that make it functional. Please, change all variables for properly for their corresponding files, paths and programs. They most shows the 
exact location of these files, so they can be processed properly.

1. Quality Control: FASTQC and SAMtools.

```bash
# running script using sbatch function

sbatch 01_quality_control.cluster.sh
```
   
2. Processing samples: trimming and alignment of Raw reads.

```bash
# run script using a nohup function or a screen linux session.
sbatch 02_trimming_n_alignment.cluster.sh
```
   
3. Re-alignment of indels in dedupped alignments.

```bash
sbatch bash 03_trimming_n_alignment.cluster.sh
```
   
4. Variant calling and hard-filtering with BCFtools.

```bash
sbatch 04_trimming_n_alignment.cluster.sh
```

