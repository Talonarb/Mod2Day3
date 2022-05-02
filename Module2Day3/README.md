## What does this program do? 

## Short version
Takes an input of two fasta files, one as a query, the other as reads, and assembles a contig that is searched using the query sequence, outputting a tsv file for the
contig coverage from the read file, as well as a fasta file of the largest assembled contig.

## Long version 
This program takes an input of a reads file in fasta format coupled with a query file (again, in fasta format), 
constructs contigs from the reads file, and then outputs the largest constructed contig which contains the initial query as "Alleles.fasta", 
the alignment of sequence reads (from the reads file) to assembled contigs as "Alleles.aln", which is a tab delimited file, columns purpose/names are explained in the "Alleles.aln" 
section. 
The process by which these tasks are completed are as follows: 

files are read in based on user or default inputs (in the case of no user response)

the fasta reads file is used to create a DeBruijn graph

 the graph is traversed via an implementation of hierholzer's algorithm
 
 after traversal the path to of the Eulerian circuit is used to construct a contig 
 
 the contig is searched for query matches using Boyer-Moore (via https://github.com/BenLangmead/ads1-notebooks/blob/master/2.01_BoyerMoore.ipynb , with alterations)
 
 the contig coverage is searched by using the fasta reads to match, with the same Boyer-Moore algorithm, and matches are recorded in a tsv file titled "ALLELES.aln"



#### Command Line Usage
```commandline
$ python .\main.py 
```
The user will be prompted for various inputs, including: directory of the reads/query files, names of the reads/query files, 
as well as kmer size and starting node for the Hierholzer algorithm. Default values are provided for all, 
however for Kmer-size the current run is using kmers of size 10. This has yet to complete, we are on hour 18. Results pending. 

## Input
1. READS.fasta
   - fasta file, lines equating to sequence ID and sequence (the sequence ID being identifiable via the ">" before it)

2. QUERY.fasta
   - fasta file, lines equating to sequence ID and sequence (the sequence ID being identifiable via the ">" before it)

## Output
Alleles.fasta
-------------
Alleles.fasta is a fasta-format file containing the largest constructed contig 

Alleles.aln
-------------
ALLELES.ALN is a tab-delimited file describing the alignment of sequence reads to the assembled contig(s) in the ALLELES.FASTA file. 
The columns in ALLELES.ALN are as follows: 

SSEQID name of sequencing read (from READS.fasta)

QSEQID name of contig matched (from ALLELES.fasta)

SSTART starting coordinate in sequencing read sseqid that matches qseq

SEND ending coordinate in sequencing read  sseqid that matches qseq

QSTART starting coordinate in contig that matches sseq

QEND ending coordinate in contig that matches sseq
