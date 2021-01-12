# Metoid
Pipeline for metagenomic organism identification

## Installation:
1) Install [`nextflow`](https://nf-co.re/usage/installation) (version >= 20.04.0)
2) Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/)
3) Clone repository
4) docker build -t metoid_ks .
5) Set default run configurations in `nextflow.config`

## Running the pipeline:
The typical command for running the pipeline is as follows:

 ```bash
 nextflow run main.nf '--reads=data/*.bam' --bam -profile docker
  ```
This will launch the pipeline with single-end reads as input and the docker configuration profile. See below for more information about profiles.

 ```bash
 nextflow run main.nf '--reads=data/*{1,2}.fastq' --pairedEnd -profile docker
  ```
  This will launch the pipeline with paired-end reads as input.
  
   ```bash
 nextflow run main.nf '--reads=data/*.fastq.gz' --longReads -profile docker
  ```
This will launch the pipeline with long-read data as input.

## Main arguments:

### `-profile`

* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)

### `-build`

* Build a custom kraken database


