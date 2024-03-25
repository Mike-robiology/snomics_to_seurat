FROM continuumio/miniconda3
COPY conf/conda_env.yml .
RUN mamba create -f conda_env.yml
RUN echo "conda activate snomics_to_seurat" > ~/.bashrc
ENV PATH /opt/conda/envs/snomics_to_seurat/bin:$PATH