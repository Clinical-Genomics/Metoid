FROM continuumio/miniconda3:4.8.2

COPY environment.yml /
RUN conda env create -f /environment.yml
RUN echo "source activate sofia" > ~/.bashrc
ENV PATH /opt/conda/envs/sofia/bin:$PATH
