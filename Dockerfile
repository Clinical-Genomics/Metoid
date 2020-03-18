FROM ubuntu

MAINTAINER Sofia Stamouli <sofia.stamouli@folkhalsomyndigheten.se>


RUN apt-get update --fix-missing && \
  apt-get install -q -y python wget unzip samtools


RUN apt-get update && apt-get install -y pigz

RUN apt-get update && apt-get install -y software-properties-common


