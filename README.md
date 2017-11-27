#ProteomeGenerator

##Introduction

ProteomeGenerator is an open, modular, and scalable framework for reference guided and de novo proteogenomic database generation written in the Snakemake workflow management system. The workflow consists of a base file that defines the project details and input samples, and pgm, a modular set of rules sourced as an include.


##Requirements

ProteomeGenerator depends on multiple free and open source tools. The following software and all dependencies must be installed:

* [Snakemake](https://snakemake.readthedocs.io/en/stable/)
* [STAR](https://github.com/alexdobin/STAR)
* [Picard](http://broadinstitute.github.io/picard/)
* [bedtools](http://bedtools.readthedocs.io/en/latest/)
* [samtools](http://samtools.sourceforge.net)
* [StringTie](https://ccb.jhu.edu/software/stringtie/)
* [TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki)
* [gffread](https://github.com/gpertea/gffread)
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)

##Configuration

##Running ProteomeGenerator