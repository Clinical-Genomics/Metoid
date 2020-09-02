# Metoid
Pipeline for metagenomic organism identification

Installation:
- Clone repository
- docker build -t metoid_ks .
- Set default run configurations in nextflow.config

Usage:
- nextflow run main.nf '--reads=data/*.bam' --bam -profile docker
- nextflow run main.nf '--reads=data/*.bam' --bam --build -profile docker (if the kraken database will be rebuilt and updated)

Requires:
NextFlow 20.04: Allows the output of process env variables
