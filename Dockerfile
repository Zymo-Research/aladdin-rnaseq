FROM nfcore/base:2.1

COPY environment.yml /

RUN conda env create -f /environment.yml && conda clean -a

COPY assets/multiqc_plugins /opt/multiqc_plugins

SHELL ["/bin/bash", "--login", "-c"]

RUN conda activate rnaseq && cd /opt/multiqc_plugins && python setup.py install

ENV PATH /opt/conda/envs/rnaseq/bin:$PATH
