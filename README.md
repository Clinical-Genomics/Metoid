# Metoid
Pipeline for metagenomic organism identification

For the docker image:
- Create environment.yml and Dockerfile
- docker build -t metoid .
- docker run -ti metoid
- sudo nextflow run main.nf --reads="data/*.bam" --bam -profile docker

NextFlow 20.04: Allows the output of process env variables
