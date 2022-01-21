FROM continuumio/miniconda3

RUN mkdir /proj /tools

COPY packages.txt /tools/packages.txt

RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

RUN conda update -n base conda && \
    conda create -n var_call --file /tools/packages.txt 

RUN apt-get update && \
    apt-get -y install build-essential libz-dev liblzma-dev libbz2-dev libcurl4-gnutls-dev

RUN cd /tools && \
    wget https://github.com/samtools/bcftools/releases/download/1.14/bcftools-1.14.tar.bz2 && \
    tar -xf *.tar.bz2 && \
    cd bcftools-1.14 && \
    make

ENV BCF=/tools/bcftools-1.14

WORKDIR /proj

SHELL ["conda", "run", "-n", "var_call", "/bin/bash", "-c"]

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "var_call"]
CMD ["snakemake", "-c 4", "-p"]
