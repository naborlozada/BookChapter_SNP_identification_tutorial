# Tutorial for SNPs identification in genomic sequences in insects.

_Repository with scripts and steps_

## From: Book Chapter: Identification of Single Nucleotide Polymorphism from insect genomic data

**Authors:** _Alejandro Nabor Lozada-Chavez and Mariangela Bonizzoni_. University of Pavia (Pavia, Italy).

**Citation:** 
* Lozada-Chávez AN, Bonizzoni M. _Identification of Single Nucleotide Polymorphism from Insect Genomic Data._ Methods Mol Biol. 2025;2935:29-49. doi: 10.1007/978-1-0716-4583-3_2. PMID: 40828273.
* Lozada-Chávez AN and Bonizzoni M. 2025. _Identification of Single Nucleotide Polymorphism from Insect Genomic Data (Chapter 2)_. In: _Insect Genomics: Methods and Protocols_. Methods in Molecular Biology, Vol. 2935. Bonizzoni M. & Ometto L. (eds). Human Press (Springer Publisher). doi: https://doi.org/10.1007/978-1-0716-4583-3_2._

----
### Description.

The series of scripts and command lines presented in this repository are fully described in [our book chapter](https://doi.org/10.1007/978-1-0716-4583-3_2) in the book **"Insect Genomics: Methods and Protocols"**. The guideline of this tutorial is divided in three major steps based on the GATK protocol: 1) preprocessing data and alignments, 2) refinement of alignments, and 3) variant calling and hard-filtering of Single Nucleotide Polymorphisms (SNPs).

I. Installation of tools: FASTQC, BWA, vcftools, bcftools, samtools, picard, GATK. 

II. The main content of this tutorial describes the following topics:
 * Quality control of the raw reads
 * Mapping of short reads to a reference genome
 * Realignment of insertion and deletions (indels)
 * Variant calling and hard-filtering of SNPs candidates

III. Extra content: recalibration of aligned reads and parallelization, in which the information of these topics are located in the directories "BQSR" and "Parallelization".

IV. Final comments.

### Directories:

1. **Installation**: Command lines from the book chapter.
2. **Scripts_tutorial**:\
   2.1. Single_scripts: single sample command lines by each step.
   2.2. Scripts_descriptions: description of each single step. 
3. **BQSR**: recalibration of aligned reads.\
    3.1 Steps and description.
4. **Parallelization**: SLURM & GNU parallel as scripts followed by descriptions of each step.

Working directory structure:

```text
/home/
  └── username/ 
    └── snps_identification/ 
      ├── alignments/ 
      ├── mapped_n_unmapped/ 
      ├── raw_reads/ 
      ├── reference_genome/
      └── variant_calling/
```


