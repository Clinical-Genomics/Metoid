# Metoid
Pipeline for metagenomic organism identification

For the docker image:
- Create environment.yml and Dockerfile
- docker build -t metoid_ks .
- docker run -ti metoid_ks
- nextflow run main.nf '--reads=data/*.bam' --bam --kraken2_build-with-docker metoid_ks (if the kraken database will be rebuilt and updated)




NextFlow 20.04: Allows the output of process env variables
