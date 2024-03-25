# get base
FROM condaforge/mambaforge:24.1.2-0

# install env
COPY conf/conda_env.yml .
RUN mamba env create -f conda_env.yml

# activate env
RUN echo "conda activate snomics_to_seurat" > ~/.bashrc
ENV PATH /opt/conda/envs/snomics_to_seurat/bin:$PATH