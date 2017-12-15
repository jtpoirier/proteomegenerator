# ProteomeGenerator

## Introduction

ProteomeGenerator is an open, modular, and scalable framework for reference guided and de novo proteogenomic database generation written in the Snakemake workflow management system. The workflow consists of a base file that defines the project details and input samples, and pgm, a modular set of rules sourced as an include.


## Requirements

ProteomeGenerator depends on multiple free and open source tools. The following software and all dependencies must be installed:

* [Snakemake](https://snakemake.readthedocs.io/en/stable/)
* [STAR](https://github.com/alexdobin/STAR)
* [Picard](http://broadinstitute.github.io/picard/)
* [samtools](http://samtools.sourceforge.net)
* [StringTie](https://ccb.jhu.edu/software/stringtie/)
* [TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki)
* [gffread](https://github.com/gpertea/gffread)
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)

## Configuration

ProteomeGenerator is configured by setting a series of variables to the explicit location of each of the required tools above in order to avoid conflicting versions that may be available in the path.

## Running ProteomeGenerator

ProteomeGenerator can be run locally or in a variety of cluster environments. The below example demonstrates how to run ProteomeGenerator on an LSF cluster head node from within screen in case the connection to the head node is lost.

```bash
screen -S pg
snakemake --snakefile Snakefile-K0562 --cluster "bsub -J {params.J} -n {params.n} -R {params.R} -W 4:00 -o {params.o} -eo {params.eo}" --jn {rulename}.{jobid}.sj -j 50 -k --latency-wait 60 --ri
```