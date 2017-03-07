# Codes and functions for the multiscale analysis of the Escherichia coli 3C data 
This page presents and explains the different codes and functions developped for the analysis of 3C data of the organism *Escherichia coli* presented in the article **Cooperation of a condensin complex and nucleoid-associated proteins for the multiscale conformation of a bacterial chromosome** by Virginia S. Lioy, Axel Cournac, Martial Marbouty, Stéphane Duigou, Julien Mozziconacci, Olivier Espéli, Frédéric Boccard, Romain Koszul.
The codes presented here allow to exactly reproduce the various plots present in the main manuscript and supplementary data. 


For queries or help getting these running, you can send email or open an issue at the github repository.

### Table of contents

* [Dependencies](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#dependencies)
* [Raw data extraction and alignment](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#raw-data-extraction-and-alignment)
* [Building of the contacts map](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#building-of-the-contacts-map)
* [Normalization of the data](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#normalization-of-the-data)
* [Computation of genomic distance law](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#computation-of-genomic-distance-law)
* [Computation of correlation matrices](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#computation-of-correlation-matrices)
* [Directional Index tool to detect TADs](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#directional-index-tool-to-detect-tads)
* [Decomposition into eigen vectors](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#decomposition-into-eigen-vectors)
* [Use of sparse formalism](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#use-of-sparse-formalism)
* [Miscellaneous](https://github.com/axelcournac/3C_tutorial/blob/master/README.md#miscellaneous)

### Dependencies

Scripts and codes can be run on OS X and other Unix-based systems. It basically requires to have Python installed on your machine which is commonly installed on Unix-based systems. 
For windows, you can have a look to https://www.python.org/downloads/windows/. Then, a few python modules are necessary for diverses operations on arrays and vizualisation. 

#### Python (>=2.7)
* Numpy
* Matplotlib (>=1.0)
* Scipy
* Biopython

#### External programs

* `Bowtie2 ` / [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* 'Pymol' / [Pymol](https://www.pymol.org/)
* `Shrec3D` / [Shrec3D](https://sites.google.com/site/julienmozziconacci/)


## Raw data extraction and alignment
#### Data extraction
Data can be dowloaded on Short Read Archive server at the following address **http://www.ncbi.nlm.nih.gov/sra**.

A SRA executable called fastq-dump from SRA can be used to extract and split both mates of a library (to use it, you can go with your terminal to the directory containg the executables files by using the bash command cd).Then the program can be used like this:  /fastq-dump library_identification --split-3 -O /path_to_a_directory

```bash
./fastq-dump SRR639031 --split-3 -O /run/media/axel/RSG3/IMR90_data/
```

#### Alignment

For the alignement step, we will use the sofware Bowtie2 and an iterative procedure like the one of [hiclib] (https://bitbucket.org/mirnylab/hiclib). 

Here, some lines that can be used to do this latter task:

```bash
#  Keeping only the columns of the sam file that contain necessary information:
awk '{print $1,$3,$4,$2,$5;}' p1.sam > p1.sam.0
awk '{print $1,$3,$4,$2,$5;}' p2.sam > p2.sam.0

# Sort according to the read identification to have both mates in the same order
# if sort does not have -V option try -d
sort -V -k1 p1.sam.0 > p1.sam.0.sorted
sort -V -k1 p2.sam.0 > p2.sam.0.sorted

# Pairing of both mates in a single file
paste p1.sam.0.sorted p2.sam.0.sorted > p1_p2_merged

# Removal of intermediar files
rm p1.sam.0.sorted
rm p2.sam.0.sorted

# Filtering of paires of reads that both have a Mapping Quality above 30
awk '{if($1 eq $6 && $5>= 30 && $10 >= 30) print $2,$3,$4,$7,$8,$9}'  p1_p2_merged  > output_alignment_idpt.dat

# Removal of intermediar file
rm p1_p2_merged
```





