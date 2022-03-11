***********************
What does this do?
***********************

This program takes an input of a reads file in fasta format coupled with a query file (again, in fasta format), constructs contigs from the reads file, and then outputs the largest constructed contig which contains the initial query as "Alleles.fasta", the alignment of sequence reads (from the reads file) to assembled contigs as "Alleles.aln", which is a tab delimited file, columns purpose/names are explained in the "Alleles.aln" section. 

What is included in this repository??
===================================
1. **Python program for Sequence assembly:** 
These components include a FASTA file reader function , a query comparison function (from a query string generated from a FASTA file) [not yet complete], and the rest of the components are currently stubbed out as they are still being made. 

2. **Default Files:** 
READS.fasta and QUERY.fasta are the files that are default for the program to use as opposed to user inputs of other fasta query/reads files, and are designed to be held in the same directory such that without user input they become the default inputs of the included scripts. 

************************
Getting Started
************************

This program requires Python version 3 or later.

Dependencies
-------------
Included in this repository are two files for use, a READS.fasta file and a QUERY.fasta file. These are included to be used as default inputs, but users are prompted to input the path and filename of their FASTA files of interest. 
This workflow also requires the following prerequesites installed (thus far):
pandas 
numpy


Jupyter Notebook
==========================================
The Jupyter notebook provided is the primary component of this build, which currently has two working functions: file reading, and query comparison, the latter of which is debugged. 

Alleles.aln
-------------
ALLELES.ALN is a tab-delimited file describing the alignment of sequence reads to the assembled contig(s) in the ALLELES.FASTA file. 
The columns in ALLELES.ALN are as follows: 

SSEQID name of sequencing read (from READS.fastq.gz)

QSEQID name of contig matched (from ALLELES.fasta)

SSTART starting coordinate in sequencing read sseqid that matches qseq

SEND ending coordinate in sequencing read  sseqid that matches qseq

QSTART starting coordinate in contig that matches sseq

QEND ending coordinate in contig that matches sseq