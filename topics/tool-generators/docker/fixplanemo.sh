#!/bin/bash
planemo/con/bin/conda init  bash
. ~/.bashrc \
/planemo/con/bin/conda activate base \
conda install -y -c bioconda -c conda-forge configparser galaxyxml
