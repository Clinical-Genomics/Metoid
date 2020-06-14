FROM continuumio/miniconda3:4.8.2

RUN apt-get update --fix-missing && apt install -yq make gcc g++ gfortran git

COPY environment.yml /

RUN conda env create -f environment.yml && conda clean --all 
RUN echo "source activate metoid_ks" > ~/.bashrc
ENV PATH /opt/conda/envs/metoid_ks/bin:$PATH

RUN apt-get install -y curl grep sed dpkg && \
    TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb && \
    apt-get clean

ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD [ "/bin/bash" ]

